#include "TROOT.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMath.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TArrayD.h"
#include "TArrayI.h"
#include "TNamed.h"

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <ctime>

#include "ConfigValue.C"

/*
	Changelog:
	2.0: removed reset on histograms for image (spedup processing by ~100 times)
	2.1: added noise 
	3.0: Fixed reading so it didn't skip single-deposit events (previously it would skip over those)
	Changed pixel size to correspond to physical dimensions
	3.1: made some changes to allow damic1k simulation to work (histograms are now vectors instead of arrays, dynamic number of CCDs, etc.)
	3.2: fixed bug with infotree reading X/Y/Z position from the *next* entry prior to diffusion
	3.3: Added unique event ID and count number of CCDs hit per unique IDs
	3.4: Moved TRandom3 to a global generator called r3 (looked like it was taking lots of time to initizalize the RNG everytime findnext was called).
	3.5: Added Joao Process Categorization code
	3.6: Fixed bug of skipping last event
	3.7: Added "surrounding" pixels and npixels for LL processing
	3.8: Fixed bug with generation of unique ID not actually being unique
	3.9: Minor bugfix, tree for production inside CCD should have been ProductCCD not InsideP
	4.0: Added noise to LL surround
	4.1: Added 1x100 capability 
	4.2: Added DAMIC42 processing capability (mostly number of pixels)
	5.0: Added support for new DamicSimu output format
	5.1: Modified the way we generate unique_event_id slightly (now random for lowest 9 digits, simulation eventID for digits 9+)
	5.2: Fixed bug with reading X/Y/Z coordinates for decay from next entry, again, hopefully this time for good (*knock on wood*). Note that damic1k simulation is *not* fixed by this
	5.3: Set the number of pixels to the correct DAMIC100 value
	5.4: Added some additional simulation variables into the finfo tree
	5.5: Added some more simulation parameters into the finfo tree (using TNamed's)
	5.6: Made added an option (writeExtra) to disable some variables that massively inflate the output file size
*/

Double_t ReadSimVersion=5.5;

//Physical CCD simulation
Int_t num_ccds = 8;
Double_t zD = 675.; //thickness of the CCD in um
Int_t nx = 4116; // Fixed # of bins
Int_t ny = 4148;
const Int_t window = 250;
Double_t sigma_value=0.0; // noise in e-
const Double_t seed = 4.0;
const double hwid=0.0;
// Conversion Parameters
const double etoev = 0.00377; // keV/e-

Double_t px = 15.; //Default pixel size in um - overwritten by simulation
Double_t py = 15.;



// Image Histogram
std::vector<TH2D*> him;
std::vector<TH2D*> zdepth;
std::vector<TH2D*> Etrue;

// Vector of Values
vector<Double_t>* Vx;
vector<Double_t>* Vy;
vector<Double_t>* Ve;
vector<Double_t>* VE;
vector<Double_t>* Vz;
vector<Double_t>* Vn;
vector<Int_t>* Vp1;
vector<Int_t>* Vp2;
vector<Double_t>*VpxPart;
vector<Double_t>*VpyPart;
vector<Int_t>*VpIns;

// Search seed vector for FindCluster
vector<vector<pair<Int_t,Int_t> > >* S = new vector<vector<pair<Int_t,Int_t> > >();

// Final Histogram
std::vector<TH1F*> hFinal;
const Double_t minbin = 0.;
const Double_t maxbin = 70.;
const Int_t numbins = 100;

// Original decay locations
Double_t Xi, Yi, Zi;
Long64_t unique_event_id; //Unique event ID accross all simulations
Int_t num_ccds_hit;

//Helper function to translate CCD number (1-8) to extID(1,2,3,4,6,11,12)
Int_t CCDNumToExtID(Int_t ccdnum) {
	//DAMIC100 configuration
	if (num_ccds==8) {
		switch(ccdnum) {
		case 1:
			return 12;
			break;
		case 2:
			return 11;
			break;
		case 3:
			return 6;
			break;
		case 4:
			return 5;
			break;
		case 5:
			return 4;
			break;
		case 6:
			return 3;
		case 7:
			return 2;
			break;
		case 8:
			return 1;
			break;
		default:
			return -1;
			break;
		}
		return -1;
	}
	//Toy configuration (damic1k, mostly)
	else {
		return ccdnum;
	}
}

TString GetSimSVNVersion(TString fname) {
	/*
		Attemps to find the sim SVN version from a file assuming standard simulation directory format
	 */
	TString dirname=fname.Remove(fname.Last('/'), fname.Length());
	std::ifstream file(dirname+"/../../../SVN_VERSION");
	if (! file.good()) {
		std::cout << "SIM_VERSION file not found" << std::endl;
		return TString("NULL");
	}
	else {
		TString version;
		file >> version;
		std::cout << "Found simulation version:" << version <<std::endl;
		return version;
	}
}

Double_t TransformToCCD(Double_t X, Double_t pixX){
  Double_t X1 = X;
   X1 = X1 * 10000;
   X1 = (X1/pixX) +1;
   X1 = floor(X1);
   return X1;
}

// Diffuses Electrons
void DiffuseDep(vector<Double_t>* xloc, vector<Double_t>* yloc, vector<Double_t>* zloc, vector<Double_t>* ene, vector<Int_t>* CCD, Double_t F, Double_t A, Double_t b){

	Int_t nevents = xloc->size();
	Int_t CCDId;
	for(int i=0; i<nevents; i++){

		Double_t z = ((*zloc)[i]) * 10000; // cm -> um
		Double_t x = ((*xloc)[i]) * 10000;
		Double_t y = ((*yloc)[i]) * 10000;

		Double_t var = -A*TMath::Log(1. - b*z);
		if (var > 0) {
			var=TMath::Sqrt(var);
		} else {
			var=0;
		}

		Double_t x_shift=0;
		Double_t y_shift=0;

		CCDId=(*CCD)[i];
		// Calculate # of observed electrons (w/ Fano)
		Int_t eobs = TMath::Nint(r3.Gaus((*ene)[i]/etoev,TMath::Sqrt(F*(*ene)[i]/etoev)));
		Double_t x_int=(x/px)+1;
		Double_t y_int=(y/py)+1;

		/* Etrue[CCDId]->Fill(x_int,y_int,eobs); */
		/* zdepth[CCDId]->Fill(x_int,y_int,z*eobs/10000.); */
		TH2D* current_him=him[CCDId];
		for(Int_t j=0; j<eobs; j++){

			if(var>0){
				y_shift=r3.Gaus(0,var);
				x_shift=r3.Gaus(0,var);
			} else {
				y_shift=0;
				x_shift=0;
			}

			Double_t bx = ((x+x_shift)/px) + 1;
			Double_t by = ((y+y_shift)/py) + 1;

			Int_t curr_content = current_him->GetBinContent(bx,by);
			Etrue[CCDId]->SetBinContent(bx, by, Etrue[CCDId]->GetBinContent(bx,by)+1);
			zdepth[CCDId]->SetBinContent(bx, by, zdepth[CCDId]->GetBinContent(bx,by)+(z/10000.));
			/* Etrue[CCDId]->Fill(bx,by,1.); */
			/* zdepth[CCDId]->Fill(bx,by,z/10000); */
			current_him->SetBinContent( bx, by, curr_content+1);

			if (curr_content < 1)
				(*S)[CCDId].push_back(make_pair(bx,by));
		}
	}

}

void findnext(Int_t x, Int_t y, Int_t CCD){

	Double_t bin_content = him[CCD]->GetBinContent(x,y);
	Double_t E_bin_content = Etrue[CCD]->GetBinContent(x,y);
	Double_t z_bin_content = zdepth[CCD]->GetBinContent(x,y);

	if ((x < 1) || (x > nx) || (y < 1) || (y > ny) || (bin_content<0.01)){
		Etrue[CCD]->SetBinContent(x,y,0.);
		zdepth[CCD]->SetBinContent(x,y,0.);

		return;
	}

	bin_content += r3.Gaus(0,sigma_value);

	if (bin_content > seed*sigma_value){

		Vx->push_back((Double_t)x);
		Vy->push_back((Double_t)y);
		Ve->push_back(bin_content*etoev);
		VE->push_back(E_bin_content*etoev);
		Vz->push_back((Double_t)(z_bin_content/E_bin_content));
		Vn->push_back(bin_content);
	}

	him[CCD]->SetBinContent(x,y,0.);
	Etrue[CCD]->SetBinContent(x,y,0.);
	zdepth[CCD]->SetBinContent(x,y,0.);

	findnext(x-1,y,CCD);
	findnext(x-1,y-1,CCD);
	findnext(x-1,y+1,CCD);
	findnext(x,y-1,CCD);
	findnext(x,y+1,CCD);
	findnext(x+1,y-1,CCD);
	findnext(x+1,y,CCD);
	findnext(x+1,y+1,CCD);

}

void FindClusters(TTree *tree, TArrayD* pixel_x, TArrayD* pixel_y, TArrayD* pixel_val, TArrayD* pixel_sigma, TArrayD* pixel_ped, Int_t &extid, Int_t &npixels, Double_t &Xdep, Double_t &Ydep, Double_t &Zdep, TArrayD* pixel_sime, TArrayD* pixel_simz, TArrayD* pixel_simn, TArrayI* processCluster1, TArrayI* processCluster2, TArrayI* ProduceInside, TArrayD* xpart, TArrayD*ypart, Bool_t make_extra_pix, Bool_t one_by_onehundred, Int_t neighbor, bool writeExtra){
	Int_t binx = 0;
	Int_t biny = 0;

	Double_t wk_pixel = 0;
	Double_t dep_ene = 0;
	num_ccds_hit=0;
	for (int CCD = 1; CCD<=num_ccds; CCD++) {
		Int_t size_seed=(*S)[CCD].size();
		if ( size_seed > 0) {
			num_ccds_hit+=1;
		}
	}
	for (int CCD = 1; CCD<=num_ccds; CCD++){

		Int_t size_seed = (*S)[CCD].size();
		extid = CCDNumToExtID(CCD);
		Xdep = Xi;
		Ydep = Yi;
		Zdep = Zi;
		while (size_seed > 0){

			binx = ((*S)[CCD])[size_seed-1].first;
			biny = ((*S)[CCD])[size_seed-1].second;

			wk_pixel = him[CCD]->GetBinContent(binx,biny);
			if (wk_pixel > 0.01) { //it's electrons so quantized...
				findnext(binx,biny,CCD);
				dep_ene = accumulate(Ve->begin(), Ve->end(), 0.0); //removed *etoev
				hFinal[CCD]->Fill(dep_ene);
				hFinal[0]->Fill(dep_ene);
				if (Vx->size() > 0) { //Only add to our tree if we've got >0 pixels
					// Add surrounding pixels for LL
					if (make_extra_pix){
						vector<pair<Int_t,Int_t> >* P = new vector<pair<Int_t,Int_t> >();
						npixels = Vx->size();
						for (Int_t i = 0; i<npixels; i++){
							P->push_back(make_pair(((Int_t)(Vx->at(i))),((Int_t)(Vy->at(i)))));
						}
						for (Int_t i=0; i<npixels; i++){
							Int_t currx = ((Int_t)(Vx->at(i)));
							Int_t curry = ((Int_t)(Vy->at(i)));
							for (Int_t a=-neighbor; a<=neighbor; a++){
								if ((curry+a < 1) || (curry+a >ny))
									continue;
								if (one_by_onehundred && a!=0)
									continue;
								for (Int_t b=-neighbor; b<=neighbor; b++){
									if (a==0 && b==0)
										continue;
									if ((currx+b < 1) || (currx+b > nx))
										continue;
									pair<Int_t,Int_t> p = make_pair(currx+b,curry+a);
									if(find(P->begin(), P->end(), p) == P->end()) {
										P->push_back(p);
										Vx->push_back(currx+b);
										Vy->push_back(curry+a);
										if (sigma_value > 0.)
											Ve->push_back(r3.Gaus(0,sigma_value)*etoev);
										else
											Ve->push_back(0.);	
									}
								}
							}
						}
					}
					
					pixel_x->Set(Vx->size(), &(*Vx)[0]);
					pixel_y->Set(Vy->size(), &(*Vy)[0]);
					pixel_val->Set(Ve->size(), &(*Ve)[0]);
					pixel_sime->Set(VE->size(),&(*VE)[0]);
					pixel_simz->Set(Vz->size(),&(*Vz)[0]);
					pixel_simn->Set(Vn->size(),&(*Vn)[0]);

					if (writeExtra) {
						xpart->Set(VpxPart->size(), &(*VpxPart)[0]);
						ypart->Set(VpyPart->size(), &(*VpyPart)[0]);
						processCluster1->Set(Vp1->size(), &(*Vp1)[0]);
						processCluster2->Set(Vp2->size(), &(*Vp2)[0]);
						ProduceInside->Set(VpIns->size(), &(*VpIns)[0]);
					}
					//Setup out ped/sigma w/ empty values
					pixel_sigma->Set(Vx->size());
					pixel_ped->Set(Vx->size());
					for (Int_t pix=0; pix < Vx->size(); pix++) {
						pixel_sigma->AddAt(sigma_value, pix);
						pixel_ped->AddAt(0.0, pix);
					}
					tree->Fill();
				}
				Vx->clear();
				Vy->clear();
				Ve->clear();
				VE->clear();
				Vz->clear();
				Vn->clear();
			}
			size_seed-=1;

		}
		//Note: these are VERY SLOW
		/* him[CCD]->Reset(); */
		/* zdepth[CCD]->Reset(); */
		/* Etrue[CCD]->Reset(); */
		(*S)[CCD].clear();
	}
}

void ReadSimFast(TString cfgfile="default.cfg", TString InputFile="CopperTopPlate44PV60a27z.root", TString OutputFile="Slow.root", TString version="NULL", bool drawHistograms=false,bool writeExtra=true ){

	//Initialize vectors

	vector<Double_t>* VxRaw = new vector<Double_t>();
	vector<Double_t>* VyRaw = new vector<Double_t>();
	vector<Double_t>* VzRaw = new vector<Double_t>();
	vector<Double_t>* VeRaw = new vector<Double_t>();
	vector<Int_t>* VcRaw = new vector<Int_t>();

	vector<Int_t>* Vp1Raw = new vector<Int_t>();
	vector<Int_t>* Vp2Raw = new vector<Int_t>();
	vector<Double_t>* VpxPartRaw = new vector<Double_t>();
	vector<Double_t>* VpyPartRaw = new vector<Double_t>();
	vector<Int_t>* VpInsRaw = new vector<Int_t>();


	Vx = new vector<Double_t>();
	Vy = new vector<Double_t>();
	Ve = new vector<Double_t>();
	VE = new vector<Double_t>();
	Vz = new vector<Double_t>();
	Vn = new vector<Double_t>();

	Vp1 = new vector<Int_t>();
	Vp2 = new vector<Int_t>();
	VpxPart = new vector<Double_t>();
	VpyPart = new vector<Double_t>();
	VpIns = new vector<Int_t>();

	TTree *tree = new TTree("clusters_tree","clusters_tree");

	TArrayD *pixel_x = new TArrayD();
	TArrayD *pixel_y = new TArrayD();
	TArrayD *pixel_sime = new TArrayD();
	TArrayD *pixel_simz = new TArrayD();
	TArrayD *pixel_simn = new TArrayD();
	TArrayD *pixel_val = new TArrayD();
	TArrayD *pixel_sigma = new TArrayD();
	TArrayD *pixel_ped = new TArrayD();

	TArrayI *inside = new TArrayI();
	TArrayD *xpart = new TArrayD();
	TArrayD *ypart = new TArrayD();
	TArrayI* processCluster1 =  new TArrayI();
	TArrayI* processCluster2 = new TArrayI();


	Int_t extid;
	Double_t Xdep, Ydep, Zdep, Xi_last, Yi_last, Zi_last;
	Int_t npixels;
	Bool_t damic1k=false;
	TParameter<double>* HWID=new TParameter<double>("HWID", hwid);

	tree->Branch("pixel_x", &pixel_x);
	tree->Branch("pixel_y", &pixel_y);
	tree->Branch("pixel_val", &pixel_val);
	tree->Branch("pixel_sime", &pixel_sime);
	tree->Branch("pixel_simz", &pixel_simz);
	tree->Branch("pixel_simn", &pixel_simn);
	tree->Branch("pixel_sigma", &pixel_ped);
	tree->Branch("pixel_ped", &pixel_ped);
	tree->Branch("Xdep", &Xdep);
	tree->Branch("Ydep", &Ydep);
	tree->Branch("Zdep", &Zdep);
	
	if (writeExtra) {
		tree->Branch("Xpart", &xpart);
		tree->Branch("Ypart", &ypart);
		tree->Branch("Process1",&processCluster1);
		tree->Branch("Process2",&processCluster2);
		tree->Branch("ProduceInside", &inside);
	}
	tree->Branch("HWID", HWID);
	tree->Branch("EXTID", &extid);
	tree->Branch("unique_event_id", &unique_event_id);
	tree->Branch("num_ccds_hit", &num_ccds_hit);
	tree->Branch("npixels", &npixels);
	// Load Input Files
	TChain* OutputTree = new TChain("OutPut");
	TFile *f = new TFile(OutputFile,"RECREATE");
	f->cd();
	OutputTree->Add(InputFile);
	const Int_t entries = OutputTree->GetEntries();
	cout << entries << endl;

	TChain* InfoTree = new TChain("InfoAna");
	//	InfoTree->Add(InputFile);
	TChain *CCDInfo = new TChain("CCDInfo");
	CCDInfo->Add(InputFile);


	Int_t eventNum, CCDNum, ccdConfig;
	Double_t sizeX, sizeY, sizeZ; // CCD Sizes
	Double_t XCoord, YCoord, ZCoord, EnergyDeposit;
	Int_t EventID;
	Int_t InsideP;
	Int_t processNum1, processNum2, partID;

	CCDInfo->SetBranchAddress("sizeX", &sizeX);
	CCDInfo->SetBranchAddress("sizeY", &sizeY);
	CCDInfo->SetBranchAddress("Thick", &sizeZ);
	CCDInfo->SetBranchAddress("CCDConfiguration", &ccdConfig);
	CCDInfo->GetEntry(0);
	Double_t sim_version=-1;
	if (CCDInfo->GetLeaf("SimVersion")) {
		sim_version=CCDInfo->GetLeaf("SimVersion")->GetValue();
	}

	//Hardcoded for the damic1k simulation, bad Ryan
	if (sim_version>999 || ccdConfig > 8) {
		std::cout << "Using 1k processing code..." << std::endl;
		damic1k=true;
		nx=6000;
		ny=6000;
		num_ccds=ccdConfig;
	} else if ((sim_version > 99 && sim_version < 999)) {
		std::cout << "Using Damic42 processing code" << std::endl;
		nx=4000;
		ny=2000;
		num_ccds=ccdConfig;
	}

	//Config
	LoadConfigFile(cfgfile);

	// Fano Factor
	Double_t F = ConfigValue("sim_fano_factor");
	
	// LL parameters
	Bool_t make_extra_pix = ConfigValue("make_extra_pix");
	Bool_t one_by_onehundred = ConfigValue("sim_onebyhundred");
	Int_t neighbor = ConfigValue("make_extra_pix");

	// Diffusion Parameters
	Double_t A = ConfigValue("sim_diff_A");
	Double_t b = ConfigValue("sim_diff_b");

	zD=sizeZ*1000; //Gets the thickness of our CCD
	if (one_by_onehundred)
		ny = ny/100;
	px = 10000.*sizeX/nx; // Set new pixel size
	py = 10000.*sizeY/ny;
	
	// Sigma
	sigma_value = ConfigValue("sim_noise");

	//Check for some legacy file formats
	//TODO: maybe alter to output versions?
	Double_t output_version;
	if (OutputTree->GetLeaf("ProductCCD")) {
		output_version=2.0;
	} else {
		cout << "Warning, no Inside Flag " << endl;
	}
	if (OutputTree->GetLeaf("ParticleID")) {
		output_version=3.0;

	}

	cout << "Using DamicSimu output_version " << output_version << endl;

	OutputTree->SetBranchAddress("XCoord",&XCoord);
	OutputTree->SetBranchAddress("YCoord",&YCoord);
	OutputTree->SetBranchAddress("ZCoord",&ZCoord);
	OutputTree->SetBranchAddress("EnergyDeposit",&EnergyDeposit);
	OutputTree->SetBranchAddress("CCDNum", &CCDNum);
	OutputTree->SetBranchAddress("ProcessNum1",&processNum1);
	OutputTree->SetBranchAddress("ProcessNum2",&processNum2);

	if(output_version >= 2.0){
		OutputTree->SetBranchAddress("ProductCCD", &InsideP);
	}

	if (output_version >= 3.0) {
		OutputTree->SetBranchAddress("EventID",&eventNum);
		OutputTree->SetBranchAddress("ParticleID",&partID);
		OutputTree->SetBranchAddress("PrimaryX",&Xi_last);
		OutputTree->SetBranchAddress("PrimaryY",&Yi_last);
		OutputTree->SetBranchAddress("PrimaryZ",&Zi_last);
	} else {
		OutputTree->SetBranchAddress("IDPrim",&eventNum);
		OutputTree->SetBranchAddress("PartIDElec",&partID);

		InfoTree->Add(InputFile);
		InfoTree->SetBranchAddress("X",&Xi_last);
		InfoTree->SetBranchAddress("Y",&Yi_last);
		InfoTree->SetBranchAddress("Z",&Zi_last);
	}

	//Chain Processing
	OutputTree->GetEntry(0);
	Int_t eventNumPrev = eventNum;
	Int_t CCDNumPrev = CCDNum;
	Int_t IDNumPrev = partID;

	clock_t begin = clock();

	for (int i = 0; i<=num_ccds; i++){
		hFinal.push_back(new TH1F(Form("hFinal_%d",i),Form("CCD: %d",i),numbins,minbin,maxbin));
		him.push_back(new TH2D(Form("him_%d",i),Form("CCD: %d",i),nx,0,nx,ny,0,ny));
		Etrue.push_back(new TH2D(Form("Etrue_%d",i),Form("CCD: %d",i),nx,0,nx,ny,0,ny));
		zdepth.push_back(new TH2D(Form("zdepth_%d",i),Form("CCD: %d",i),nx,0,nx,ny,0,ny));
		S->push_back(vector<pair<Int_t,Int_t> >());
	}

	Int_t tempProcess1 = -1;
	Int_t tempProcess2 = -1;
	Int_t IDpartPrev = -1;
	Int_t numberParts = 0;
	//Loop over our energydeposit events
	Int_t n = 0;
	Bool_t pushedback = false;
	while (n <= entries){
		OutputTree->GetEntry(n);
		if ((eventNum == eventNumPrev || !pushedback) && n<entries) {
			Xi=Xi_last;
			Yi=Yi_last;
			Zi=Zi_last;
			VxRaw->push_back(XCoord);
			VyRaw->push_back(YCoord);
			VzRaw->push_back(ZCoord);
			VeRaw->push_back(EnergyDeposit);
			VcRaw->push_back(CCDNum);

			if (writeExtra) {
				if(IDpartPrev != partID){
					Double_t Xco = TransformToCCD(XCoord,px);
					Double_t Yco = TransformToCCD(YCoord,py);
					VpxPartRaw->push_back(Xco);
					VpyPartRaw->push_back(Yco);
					if(output_version >= 2.0){
						VpInsRaw->push_back(InsideP);
					}
					else{
						VpInsRaw->push_back(4);
					}
					Vp1Raw->push_back(processNum1);
					Vp2Raw->push_back(processNum2);
					IDpartPrev = partID;
					numberParts ++;
				}
			}
			n++;
			pushedback = true;
		}

		else {
			unique_event_id=r3.Integer(1000000000)+eventNumPrev*1000000000;
			DiffuseDep(VxRaw, VyRaw, VzRaw, VeRaw, VcRaw, F, A, b);
			
			if (writeExtra) {
				Vp1 = Vp1Raw;
				Vp2 = Vp2Raw;
				VpIns = VpInsRaw;
				VpxPart = VpxPartRaw;
				VpyPart = VpyPartRaw;
			}
			FindClusters(tree, pixel_x, pixel_y, pixel_val, pixel_sigma, pixel_ped, extid, npixels, Xdep, Ydep, Zdep, pixel_sime, pixel_simz, pixel_simn, processCluster1, processCluster2, inside, xpart, ypart, make_extra_pix, one_by_onehundred, neighbor, writeExtra);

			VxRaw->clear();
			VyRaw->clear();
			VzRaw->clear();
			VeRaw->clear();
			VcRaw->clear();

			if (writeExtra) {
				VpInsRaw->clear();
				VpxPartRaw->clear();
				VpyPartRaw->clear();
				Vp1Raw->clear();
				Vp2Raw->clear();
			}

			if (damic1k) {
				//Note: I'm not sure this actually works
				//TODO: rewrite so it works properly, this is a placeholder
				InfoTree->GetEntry(eventNumPrev);
			}
			else {
				if (output_version < 3.0) {
					InfoTree->GetEntry(n-1);
					Xi=Xi_last;
					Yi=Yi_last;
					Zi=Zi_last;
				}
			}
			eventNumPrev = eventNum;
			CCDNumPrev = CCDNum;
			pushedback = false;
			IDpartPrev = -1;
			if (n>=entries) n++;
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Time: " << elapsed_secs << endl;


	TTree *finfo = new TTree("finfo", "finfo");

	CCDInfo->GetEntry(0);
	Int_t num_event=-1;
	finfo->Branch("SIM_VERSION", &version);
	finfo->Branch("ReadSimVersion", &ReadSimVersion);
	if (CCDInfo->GetLeaf("NumEvent")) {
		num_event=(Int_t) CCDInfo->GetLeaf("NumEvent")->GetValue();
	}
	finfo->Branch("NumEvent",&num_event);
	finfo->Branch("sim_npix_x", &nx);
	finfo->Branch("sim_npix_y", &ny);
	finfo->Branch("sim_ccd_x", &sizeX);
	finfo->Branch("sim_ccd_y", &sizeY);
	finfo->Branch("sim_ccd_z", &sizeZ);
	
	//To store the strings of various outputs
	TList *vars=new TList();
	vars->Add(new TNamed(TString("SIMFPATH"), InputFile));
	TString sim_svn_version = GetSimSVNVersion(InputFile);
	vars->Add(new TNamed(TString("SIMSVNVERSION"), sim_svn_version));
	finfo->Branch(vars);

	finfo->Fill();
	finfo->Write();

	if (drawHistograms) {
		TCanvas* c = new TCanvas();
		c->cd();
		c->Draw();
		hFinal[0]->SetLineWidth(3);
		hFinal[0]->Draw();
		hFinal[0]->GetXaxis()->SetTitle("Energy [keV]");
		hFinal[0]->GetYaxis()->SetTitle("Clustered Events");
		hFinal[0]->SetName("Total");
		hFinal[0]->SetTitle("Total");
		hFinal[0]->Write();
		for (int i = 1; i<=num_ccds; i++){
			hFinal[i]->SetLineColor(i);
			hFinal[i]->SetName(Form("CCD%d",i));
			hFinal[i]->SetTitle(Form("CCD %d",i));
			hFinal[i]->Draw("same");
			hFinal[i]->Write();
		}
		c->BuildLegend();
		hFinal[0]->SetTitle("All CCDs");
		c->Write();
	}
	tree->Write();
	f->Close();
}

int main(int, char **) {
	ReadSimFast();
	return 0;
}

