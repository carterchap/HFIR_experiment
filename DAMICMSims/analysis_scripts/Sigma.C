#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TList.h"
#include "TMath.h"
#include "TLatex.h"
#include "TRegexp.h"
#include "TLegend.h"
#include "TLeaf.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "DAMICPartTable.h"

enum ReadMode_t {ALL, //All files independently
								 VOLUME, //by simulation volumed
								 ISOTOPE, //by radiactive isotope
								 PART,  //by physical part (e.g. copper box, copper modules)
								 ISO_SOURCE, //radiogenic or primeval
								 ISO_CHAIN, //Source chain of isotope
								 PART_CHAIN};


DAMICPartTable *partTable= new DAMICPartTable();
TString dir_prefix="../";
TString data_dir="/data/simulation_output/shorts/";

TString readModeToString(ReadMode_t mode, Bool_t plotData=false) {
	TString name;
	if (mode==ALL) {
		 name=TString("all");
	} else if (mode==VOLUME) {
		name=TString("volume");
	} else if (mode==ISOTOPE ) {
		name=TString("isotope");
	} else if (mode==PART) {
		name=TString("part");
	} else if (mode==ISO_SOURCE) {
		name=TString("iso-source");
	}	else if (mode==ISO_CHAIN) {
		name=TString("iso-chain");
	} else if (mode==PART_CHAIN) {
		name=TString("part-chain");
	}
	else {
		name=TString("error");
	}
	if (plotData) {
		name+=TString("_data");
	}
	return name;
}

Int_t getLine(TString read_volume, TString read_isotope, vector<TH1F*> &spectra, ReadMode_t mode, Int_t mask=0) {
	TString key_string;
	if (mode==ALL) return spectra.size();
	if (mode==VOLUME) key_string=read_volume;
	if (mode==ISOTOPE) key_string=read_isotope;
	if (mode==PART) key_string=partTable->GetPartName(read_volume);
	if (mode==ISO_SOURCE) key_string=partTable->GetIsotopeSource(read_isotope);
	if (mode==ISO_CHAIN) key_string=partTable->GetIsotopeChain(read_isotope);
	if (mode==PART_CHAIN) {
		key_string=partTable->GetPartName(read_volume);
		key_string+="_"+partTable->GetIsotopeChain(read_isotope);
	}
	for (int i=0; i < spectra.size(); i++) {
		TString title=spectra[i]->GetTitle();
		if (title.Contains(key_string) && (mask==0 || title.Contains(partTable->GetMaskName(mask)))) return i;
	}
	return spectra.size();
}

Int_t getNumEvents(TString read_volume, TString read_isotope, Double_t &read_sim_version) {
	TChain *finfo_chain=new TChain("finfo");
	finfo_chain->Add(dir_prefix+read_volume+"/*"+read_isotope+"*");
	Int_t NumEvent=-1;
	for (Int_t i=0; i < finfo_chain->GetEntries(); i++) {
		finfo_chain->GetEntry(i);
		if (finfo_chain->GetLeaf("NumEvent")) {
			NumEvent+=finfo_chain->GetLeaf("NumEvent")->GetValue();
		};
		if (finfo_chain->GetLeaf("ReadSimVersion")) {
			read_sim_version=finfo_chain->GetLeaf("ReadSimVersion")->GetValue();
		} else {
			read_sim_version=-1.;
		}
	}
	return NumEvent;
}

void Sigma(TString actFile="base_activity.in", ReadMode_t readMode=ISOTOPE, Bool_t plotRest=true, Bool_t plotData=false, Bool_t act_from_in=false, Int_t nbins=500, Double_t minX=0.5, Double_t maxX=1000.5, Double_t etokev=0.00377, Double_t CCDmass = 5.8e-3, Int_t EXTNUM=1) {
    
	Int_t min2d = 0;
	Int_t max2d = 4000;
	Int_t nbins2d = 100;
	
	Int_t nbins1d = 100;
	
	Int_t nbinsSig = 60;
	Int_t minSig = 0;
	Int_t maxSig = 3;
	
	const Int_t NCCDs=8;
	//For the plot legend
	const Int_t max_legend_size=20;
	Double_t leg_font=.02;
    
	Int_t tenkevbin = TMath::Nint((10.-minX)/((maxX-minX)/nbins));
    
	//	cout << tenkevbin << endl;
    
	TString saveFolder = "../sim_results_pdf/";

	TChain *data = new TChain("clusters_tree");

	data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoON_30000s-IntW200-OS_1x1_run3/short/*.root");
	data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoON_30000s-IntW200-OS_1x1_run4/short/*.root");
	data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoON_30000s-IntW200-OS_1x1_run5/short/*.root");
	data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoON_30000s-IntW200-OS_1x1_run6/short/*.root");

	data->SetAlias("run3","RUNID>=239 && RUNID<=516");
	data->SetAlias("run4","RUNID>=65 && RUNID<=152");
	data->SetAlias("run5","(RUNID>=224 && RUNID<239) ||(RUNID>=642 && RUNID<=663)");
	data->SetAlias("run6","RUNID>=793");
	data->SetAlias("rnhi","(RUNID>=245 && RUNID<=250) || (RUNID>=1090 && RUNID<=1094)");
	data->SetAlias("rnlo","!rnhi");
	data->SetAlias("polyon","RUNID>=269 || RUNID<239");
	data->SetAlias("polyof","RUNID<269 && RUNID>=239");
	data->SetAlias("culo","RUNID>272 || RUNID<239");
	data->SetAlias("cuhi","RUNID<=272 && RUNID>=239");
	data->SetAlias("eraw","charge_total.fVal*1.2e-3");
	data->SetAlias("eneS","charge_total.fVal*0.001*( 1.09*(EXTID==1) + 1.05*(EXTID==2) + 1.07*(EXTID==3) + 0.97*(EXTID==4) + 1.10*(EXTID==6) + 0.92*(EXTID==11) + 1.11*(EXTID==12) )");
	data->SetAlias("eneC","charge_total.fVal*0.001*( 1.10*(EXTID==1) + 1.03*(EXTID==2) + 1.11*(EXTID==3) + 1.00*(EXTID==4) + 1.22*(EXTID==6) + 0.99*(EXTID==11) + 1.11*(EXTID==12) )");
	data->SetAlias("ene","(run5==1)*eneC + (run5==0)*eneS");
	data->SetAlias("x","charge_mean_x.fVal");
	data->SetAlias("y","charge_mean_y.fVal");
	data->SetAlias("maxx","charge_max_x.fVal");
	data->SetAlias("maxy","charge_max_y.fVal");
	data->SetAlias("rmsx","charge_rms_x.fVal");
	data->SetAlias("rmsy","charge_rms_y.fVal");
	data->SetAlias("rmsxy","sqrt((rmsx*rmsx+rmsy*rmsy)/2)");
	data->SetAlias("mask","mask_edge.fVal>0");
	data->SetAlias("nsat","raw_nsat.fVal");
	data->SetAlias("valid","!mask && !(EXTID==1 && maxx==6831 && maxy==3525) && !(EXTID==12 && maxx==4441 && maxy==1706) && !(EXTID==12 && maxx==4441 && maxy==1707) && !(EXTID==6 && maxx==8152 && maxy==3399)");
	data->SetAlias("np","size_npixels.fVal");
	data->SetAlias("sz","size_slength.fVal*size_slength.fVal/(TMath::Sqrt(1.+size_aspect_ratio.fVal*size_aspect_ratio.fVal)*TMath::Sqrt(1.+1./(size_aspect_ratio.fVal*size_aspect_ratio.fVal)))");
	data->SetAlias("alpha","(np/sz > 0.45 && ene > 900) || (ene>550 && np/sz>0.75) || ene>5000");
	data->SetAlias("plasma","rmsx/rmsy>0.75 &&  alpha");
	data->SetAlias("bloom","rmsx/rmsy<=0.75 && alpha");
	data->SetAlias("face","16 - ( (EXTID<10)*(alpha*2*EXTID-bloom) + (EXTID>10)*(alpha*2*(EXTID-4)-bloom) )");
	
	TH1F *htemp;
	data->Draw("EXPTIME.fVal/24/3600 >> htemp","EXTID==1","GOFF");
	htemp=(TH1F*)gROOT->FindObject("htemp");
	Double_t expot = htemp->GetMean() * htemp->Integral();
	Int_t nbins_data=nbins/10;
	Double_t data_scale=nbins_data/(maxX-minX)/expot;
	vector<TChain*> chains;
	vector<TH1F*> spectra;
	vector<TH1F*> spectra_indiv;
    
	TString read_volume, read_isotope;
    
	Double_t read_activity;
	vector<Double_t> activity;
    
	Double_t read_mass, read_ndecays;
    
	Int_t colors[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,44,45,46}; 
    
	vector<Double_t> scales;
    
	ifstream file("activity/"+actFile);
    
	actFile.ReplaceAll("_activity.in","");
	
	Int_t num_lines=0;
	Int_t num_files=0;

	TH1F* temp = new TH1F("temp","temp",nbins,minX,maxX);
	//TODO: rewrite this to read lines a bit cleaner
	while (file >> read_volume >> read_isotope >> read_activity >> read_mass >> read_ndecays) {
		Double_t read_sim_version=-1;
		Int_t mask=0;
		//		cout << read_volume + "_" + read_isotope << endl;
		chains.push_back(new TChain("clusters_tree", read_volume + "_" + read_isotope));
		Int_t num_events=getNumEvents(read_volume, read_isotope, read_sim_version);
		if (num_events > 0) {
			read_ndecays=num_events;
		} else {
			std::cout << "Num events not in file, using .in file ndecays..." << std::endl;
		}
		if (!act_from_in) {
			read_activity=partTable->GetActivity(read_volume, read_isotope);
			read_mass=partTable->GetMass(read_volume);
		}
		std::cout << read_volume << " " << read_isotope << " " << read_activity << " " << read_mass <<std::endl;
		scales.push_back((nbins/(maxX-minX))*read_activity*read_mass/read_ndecays);
		if (read_volume.Contains('%')) {
			mask=TString(TString(read_volume("%[0-9]+")).Strip(TString::kLeading,'%')).Atoi();
			read_volume=TString(read_volume(".*%")).Strip(TString::kTrailing,'%');
		}
		int nf = chains[num_files]->Add(dir_prefix+read_volume+"/*"+read_isotope+"*");
		if (nf < 1)
			cout << "Warning: Files not found for specified isotope" << endl;
		Int_t line=getLine(read_volume, read_isotope, spectra, readMode, mask); 
		if (line >= num_lines) {
			line=num_lines;
			TString spec_title;
			//TODO: move this to it's own function / the getline function
			if (readMode==ALL) spec_title=read_volume+"_"+read_isotope;
			if (readMode==VOLUME) spec_title=read_volume;
			if (readMode==ISOTOPE) spec_title=read_isotope;
			if (readMode==PART) spec_title=partTable->GetPartName(read_volume);
			if (readMode==ISO_SOURCE) spec_title=partTable->GetIsotopeSource(read_isotope);
			if (readMode==ISO_CHAIN) spec_title=partTable->GetIsotopeChain(read_isotope);
			if (readMode==PART_CHAIN) {
				spec_title=partTable->GetPartName(read_volume);
				spec_title+="_"+partTable->GetIsotopeChain(read_isotope);
			}
			if (mask!=0) spec_title+=" "+partTable->GetMaskName(mask);
			spectra.push_back(new TH1F("", spec_title,nbins, minX,maxX));
			spectra_indiv.push_back(new TH1F("", spec_title, nbins, minX,maxX));
			num_lines++;
		}
		TLeaf *maskleaf=chains[num_files]->GetLeaf("mask");
		TString draw_cmd;
		if (read_sim_version > 0) {
			draw_cmd="charge_total.fVal >> temp";
		} else {
			draw_cmd=Form("charge_total.fVal*%f >> temp",etokev);
		}
		if (maskleaf==NULL) {
			chains[num_files]->Draw(draw_cmd,Form("EXTID==%d",EXTNUM),"GOFF");
			spectra_indiv[line]->Add(temp,scales[num_files]);
			temp->Reset();
			chains[num_files]->Draw(draw_cmd,Form("EXTID!=%d",EXTNUM),"GOFF");			
			spectra[line]->Add(temp, scales[num_files]*1./7.);
			temp->Reset();
		} else {			
			chains[num_files]->Draw(draw_cmd,Form("EXTID==%d && mask==%d",EXTNUM, mask),"GOFF");
			spectra_indiv[line]->Add(temp,scales[num_files]);
			temp->Reset();
			chains[num_files]->Draw(draw_cmd,Form("EXTID!=%d && mask==%d",EXTNUM, mask),"GOFF");			
			spectra[line]->Add(temp, scales[num_files]*1./7.);
			temp->Reset();
		};
		num_files++;
	}

	TH1F* all_spectra_indiv = new TH1F("all_spectra_indiv",actFile+Form(", CCD# %d",EXTNUM),nbins,minX,maxX);
	TH1F* all_spectra =       new TH1F("all_spectra",actFile+", Other CCDs",nbins,minX,maxX);

	for (int i=0; i<num_lines; i++) {

		spectra_indiv[i]->GetXaxis()->SetTitle("Energy [keV]");
		spectra_indiv[i]->GetYaxis()->SetTitle("Count [/day/keV]");
		all_spectra_indiv->Add(spectra_indiv[i],1);
		spectra_indiv[i]->SetTitle((TString)spectra_indiv[i]->GetTitle()+Form(", %f dru",spectra_indiv[i]->Integral(1,tenkevbin)/tenkevbin/CCDmass));

		spectra[i]->GetXaxis()->SetTitle("Energy [keV]");
		spectra[i]->GetYaxis()->SetTitle("Count [/day/keV/CCD]");
		spectra[i]->SetTitle((TString)spectra[i]->GetTitle()+Form(", %f dru",spectra[i]->Integral(1,tenkevbin)/tenkevbin/CCDmass));
		all_spectra->Add(spectra[i],1);
	}
	gStyle->SetOptStat(0);
    
    
	TH2F* temp2d = new TH2F("temp2d","temp2d",nbins2d,min2d,max2d,nbins2d,min2d,max2d);
	TH1F* temp1d = new TH1F("temp1d","temp1d",nbins1d,min2d,max2d);
	TH1F* temp1dx = new TH1F("temp1dx","temp1dx",nbins1d,4275,8275);
	TH1F* tempSig = new TH1F("tempSig","tempSig",nbinsSig,minSig,maxSig);
    

	TCanvas *c6 = new TCanvas();
	c6->cd();
	c6->Draw();
	TH1F *sY[num_files];
	TH1F *all_sY = new TH1F("all_sY",actFile + ", Sigma Y, Other CCDs",nbinsSig, minSig, maxSig);
	TH1F *DATA_sY = new TH1F("DATA_sY","DATA Sigma Y",nbinsSig, minSig, maxSig);
	for (int i=0; i<num_files; i++) {
		sY[i] = new TH1F("","",nbinsSig, minSig, maxSig);
		chains[i]->Draw("charge_rms_y.fVal >> tempSig",Form("EXTID!=%d && charge_total.fVal<20 && charge_total.fVal>0.004",EXTNUM),"GOFF");
		sY[i]->Add(tempSig,scales[i]);
		tempSig->Reset();
		all_sY->Add(sY[i],1);
	}
	data->Draw("rmsy >> tempSig",Form("EXTID!=%d && ene<20 && ene>0.004 && valid",EXTNUM),"GOFF");
	DATA_sY->Add(tempSig,1);
	tempSig->Reset();
	all_sY->Scale(1./all_sY->Integral());
	DATA_sY->Scale(1./DATA_sY->Integral());
	all_sY->GetXaxis()->SetTitle("Sigma Y");
	all_sY->GetYaxis()->SetTitle("Normalized Count");
	all_sY->SetLineWidth(3);
	all_sY->Draw();
	DATA_sY->SetLineColor(2);
	DATA_sY->SetLineWidth(3);
	DATA_sY->Draw("same");
	TLegend *l6 = c6->BuildLegend();
	l6->SetBorderSize(0);
	all_sY->SetTitle(actFile + ", Other CCDs, Sigma Y");
	//c6->SaveAs(saveFolder+actFile+"_OtherCCD_sigmaY.svg");

    
};
