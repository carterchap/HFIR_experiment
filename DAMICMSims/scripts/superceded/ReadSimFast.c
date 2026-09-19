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

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <ctime>

#include "ConfigValue.C"

//Physical CCD simulation
const Int_t num_ccds = 8;
const Double_t px = 15; //pixel size in um
const Double_t py = 15;
const Double_t zD = 675; //thickness of the CCD in um
const Int_t nx = 4000;
const Int_t ny = 4000;
const Int_t window = 250;
const Double_t sigma_value = 3.0; // noise in e-
const double hwid=0.0;
// Conversion Parameters
const double etoev = 0.00377; // keV/e-


// Image Histogram
TH2D* him[num_ccds+1];

// Vector of Values
vector<Double_t>* Vx;
vector<Double_t>* Vy;
vector<Double_t>* Ve;

// Search seed vector for FindCluster
vector<vector<pair<Int_t,Int_t> > >* S = new vector<vector<pair<Int_t,Int_t> > >();

// Final Histogram
TH1F* hFinal[num_ccds+1];
const Double_t minbin = 0.;
const Double_t maxbin = 70.;
const Int_t numbins = 100;

//Helper function to translate CCD number (1-8) to extID(1,2,3,4,6,11,12)
Int_t CCDNumToExtID(Int_t ccdnum) {
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
// Diffuses Electrons
void DiffuseDep(vector<Double_t>* xloc, vector<Double_t>* yloc, vector<Double_t>* zloc, vector<Double_t>* ene, vector<Int_t>* CCD, Double_t F, Double_t A, Double_t b){
    
    TRandom3 r(0);

    Int_t nevents = xloc->size();
    
    for(int i=0; i<nevents; i++){

        Double_t z = ((*zloc)[i])*10000; // cm -> um
        Double_t x = ((*xloc)[i])*10000;
        Double_t y = ((*yloc)[i])*10000;
        
        Double_t var = -A*TMath::Log(1. - b*z);
        TVector3 v(0,0,0);
        
        // Calculate # of observed electrons (w/ Fano)
        Int_t eobs = TMath::Nint(r.Gaus((*ene)[i]/etoev,TMath::Sqrt(F*(*ene)[i]/etoev))+0.5);
        
        for(Int_t j=0; j<eobs; j++){
            
            if(var>0){
                
                v.SetX(r.Gaus(0,TMath::Sqrt(var)));
                v.SetY(r.Gaus(0,TMath::Sqrt(var)));
            }
            
            Int_t bx = (Int_t)((x+v.X())/px + 1);
            Int_t by = (Int_t)((y+v.Y())/py + 1);
            
            Int_t curr_content = him[(*CCD)[i]]->GetBinContent(bx,by);
            
            him[(*CCD)[i]]->SetBinContent( bx, by, curr_content+1);
            
            if (curr_content < 1)
                (*S)[(*CCD)[i]].push_back(make_pair(bx,by));
        }
    }
    
}

void findnext(Int_t x, Int_t y, Int_t CCD){
    
    Double_t bin_content = him[CCD]->GetBinContent(x,y);
    
    if ((x < 1) || (x > nx) || (y < 1) || (y > ny) || (bin_content<0.01)){
        return;
    }
    
    Vx->push_back((Double_t)x);
    Vy->push_back((Double_t)y);
    Ve->push_back(bin_content);

    him[CCD]->SetBinContent(x,y,0);
    
    findnext(x-1,y,CCD);
    findnext(x-1,y-1,CCD);
    findnext(x-1,y+1,CCD);
    findnext(x,y-1,CCD);
    findnext(x,y+1,CCD);
    findnext(x+1,y-1,CCD);
    findnext(x+1,y,CCD);
    findnext(x+1,y+1,CCD);
    
}

void FindClusters(TTree *tree, TArrayD* pixel_x, TArrayD* pixel_y, TArrayD* pixel_val, TArrayD* pixel_sigma, TArrayD* pixel_ped, Int_t &extid){
    Int_t binx = 0;
    Int_t biny = 0;
    
    Double_t wk_pixel = 0;
    Double_t dep_ene = 0;

    for (int CCD = 1; CCD<=num_ccds; CCD++){
			
        Int_t size_seed = (*S)[CCD].size();
				extid = CCDNumToExtID(CCD);
        while (size_seed > 0){
        
            binx = ((*S)[CCD])[size_seed-1].first;
            biny = ((*S)[CCD])[size_seed-1].second;

            wk_pixel = him[CCD]->GetBinContent(binx,biny);
            if (wk_pixel > 0.01) { //it's electrons so quantized...
                findnext(binx,biny,CCD);
                dep_ene = accumulate(Ve->begin(), Ve->end(), 0.0)*etoev;
                hFinal[CCD]->Fill(dep_ene);
                hFinal[0]->Fill(dep_ene);
								if (Vx->size() > 0) { //Only add to our tree if we've got >0 pixels
									pixel_x->Set(Vx->size(), &(*Vx)[0]);
									pixel_y->Set(Vy->size(), &(*Vy)[0]);
									pixel_val->Set(Ve->size(), &(*Ve)[0]);
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
            }
            size_seed-=1;
        
        }
    
        him[CCD]->Reset();
        (*S)[CCD].clear();
    }
}

void ReadSimFast(TString cfgfile="default.cfg", TString InputFile="sur_pb210_merged.root", TString OutputFile="Pb210.root", TString version="NULL", bool drawHistograms=false ){
    
    //Initialize vectors

    vector<Double_t>* VxRaw = new vector<Double_t>();
    vector<Double_t>* VyRaw = new vector<Double_t>();
    vector<Double_t>* VzRaw = new vector<Double_t>();
    vector<Double_t>* VeRaw = new vector<Double_t>();
    vector<Int_t>* VcRaw = new vector<Int_t>();


    Vx = new vector<Double_t>();
    Vy = new vector<Double_t>();
    Ve = new vector<Double_t>();

		TTree *tree = new TTree("clusters_tree","clusters_tree");
		
		TArrayD *pixel_x = new TArrayD();
		TArrayD *pixel_y = new TArrayD();
		TArrayD *pixel_val = new TArrayD();
		TArrayD *pixel_sigma = new TArrayD();
		TArrayD *pixel_ped = new TArrayD();
		Int_t extid;
		TParameter<double>* HWID=new TParameter<double>("HWID", hwid);

		tree->Branch("pixel_x", &pixel_x);
		tree->Branch("pixel_y", &pixel_y);
		tree->Branch("pixel_val", &pixel_val);
		tree->Branch("pixel_sigma", &pixel_ped);
		tree->Branch("pixel_ped", &pixel_ped);		
		tree->Branch("HWID", HWID);
		tree->Branch("EXTID", &extid);
		
    for (int i = 0; i<=num_ccds; i++){
        hFinal[i]= new TH1F("",Form("CCD: %d",i),numbins,minbin,maxbin);
        him[i] = new TH2D("",Form("CCD: %d",i),nx,0,nx,ny,0,ny);
        S->push_back(vector<pair<Int_t,Int_t> >());
    }
    
    //Config
    LoadConfigFile(cfgfile);
    
    // Fano Factor
    Double_t F = ConfigValue("sim_fano_factor");
    
    // Diffusion Parameters
    Double_t A = ConfigValue("sim_diff_A");
    Double_t b = ConfigValue("sim_diff_b");
    
    // Load Input Files
    TChain* OutputTree = new TChain("OutPut");
    OutputTree->Add(InputFile);
    const Int_t entries = OutputTree->GetEntries();
    cout << entries << endl;
    
    //Chain Processing
    OutputTree->GetEntry(0);
    Int_t eventNumPrev = OutputTree->GetLeaf("IDPrim")->GetValue();
    Int_t CCDNumPrev = OutputTree->GetLeaf("CCDNum")->GetValue();

    Int_t eventNum, CCDNum;
    Double_t XCoord, YCoord, ZCoord, EnergyDeposit;

    OutputTree->SetBranchAddress("XCoord",&XCoord);
    OutputTree->SetBranchAddress("YCoord",&YCoord);
    OutputTree->SetBranchAddress("ZCoord",&ZCoord);
    OutputTree->SetBranchAddress("EnergyDeposit",&EnergyDeposit);
    OutputTree->SetBranchAddress("IDPrim",&eventNum);
    OutputTree->SetBranchAddress("CCDNum", &CCDNum);
    
    clock_t begin = clock();
    
    for (int n = 0; n < entries; n++){
        
        OutputTree->GetEntry(n);
        
        if (eventNum == eventNumPrev){
            VxRaw->push_back(XCoord);
            VyRaw->push_back(YCoord);
            VzRaw->push_back(ZCoord);
            VeRaw->push_back(EnergyDeposit);
            VcRaw->push_back(CCDNum);
        }
        else {
					DiffuseDep(VxRaw, VyRaw, VzRaw, VeRaw, VcRaw, F, A, b);
					VxRaw->clear();
					VyRaw->clear();
					VzRaw->clear();
					VeRaw->clear();
					VcRaw->clear();
					FindClusters(tree, pixel_x, pixel_y, pixel_val, pixel_sigma, pixel_ped, extid);
					eventNumPrev = eventNum;
					CCDNumPrev = CCDNum;
        }
        
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "Time: " << elapsed_secs << endl;
    
    TFile *f = new TFile(OutputFile,"RECREATE");
		TTree *finfo = new TTree("finfo", "finfo");
		finfo->Branch("SIM_VERSION", &version);
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
