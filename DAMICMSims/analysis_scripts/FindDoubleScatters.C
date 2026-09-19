#include "TChain.h"
#include "TROOT.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "common.C"
#include <vector>
#include <iostream>

//TODO: rewrite this to avoid single-electron events from diffusion
void FindDoubleScatters(TString fname, bool show_spectrum=true, bool open_viewer=false, Int_t nbins=500, Double_t xmin=0.1, Double_t xmax=1000, Double_t CCDmass=20e-3, Double_t activity=1., Double_t act_mass=1., Int_t ndecays=1, Int_t num_ccds=1) {
	TChain *chain=new TChain("clusters_tree");
	chain->Add(fname);
	if (!chain->GetBranch("unique_event_id")) {
		std::cout << "No unique event id branch found, rerun processing code" << std::endl;
	}
	Long64_t id;
	Int_t tenkevbin=TMath::Nint((10.-xmin)/((xmax-xmin)/nbins));
	TParameter<double>* charge_total=NULL;
	std::vector<Long64_t> ids;
	std::vector<Double_t> interactions;
	std::vector<Double_t> energy_deposit;
	std::vector<std::vector<Int_t> > event_numbers;
	chain->SetBranchAddress("unique_event_id", &id);
	chain->SetBranchAddress("charge_total", &charge_total);
	TH1D *num_interactions=new TH1D("num_interactions", "num_interactions",10,0,10);
	TH1D *spectrum=new TH1D("spectrum", "spectrum", nbins,xmin,xmax);
	Int_t events=0;
	Int_t nentries=chain->GetEntries();
	bool new_id=true;
	Double_t charge;
	std::cout <<"starting loop..." << std::endl;
	for (Int_t i=0; i < nentries; i++ ) {
		new_id=true;
		chain->GetEntry(i);
		charge=charge_total->GetVal();
		for (Int_t j=0; j < ids.size();j++ ) {
			if (ids[j]==id) {
				interactions[j]+=1;
				energy_deposit[j]+=charge;
				event_numbers[j].push_back(i);
				new_id=false;
			}
		}
		if (new_id) {
			ids.push_back(id);
			interactions.push_back(1);
			energy_deposit.push_back(charge);
			event_numbers.push_back(std::vector<Int_t>());
		}
	}
	for (Int_t i=0; i < interactions.size(); i++ ) {
		num_interactions->Fill(interactions[i]);
	}
	TCanvas *c1 = new TCanvas("c1", "c1");
	num_interactions->Draw();
	TCanvas *c2 = new TCanvas("c2", "c2");
	c2->cd();
	spectrum->FillN(energy_deposit.size(), &energy_deposit[0],0);
	Double_t scale=1.*(nbins/(xmax-xmin))*activity*act_mass/(ndecays*num_ccds);
	spectrum->Scale(scale);
	spectrum->Draw();
	spectrum->GetXaxis()->SetTitle("Energy [kev]");
	spectrum->GetYaxis()->SetTitle("Count");
	std::cout << "dru: " << spectrum->Integral(1,tenkevbin)/tenkevbin/CCDmass << std::endl;
	c2->SetLogy();
	TH1D* htemp=new TH1D("htemp", "htemp",nbins,xmin,xmax);
	chain->Draw("charge_total.fVal >> htemp","","GOFF");
	std::cout << "Unsummed DRU: " << htemp->Integral(1,tenkevbin)/tenkevbin/CCDmass*scale << std::endl;
	
							
	std::cout << "Fraction of events with >1 interaction:" << 1-num_interactions->Integral(0,2)*1./num_interactions->Integral() << std::endl;

	if (open_viewer) {
		std::getchar();
		TArrayD* pixel_x=NULL;
		TArrayD* pixel_y=NULL;
		TArrayD* pixel_val=NULL;
		chain->SetBranchAddress("pixel_x",&pixel_x);
		chain->SetBranchAddress("pixel_y",&pixel_y);
		chain->SetBranchAddress("pixel_val",&pixel_val);
		TString last_input="";
		Int_t id_counter=0;
		Int_t clust_counter=0;
		Int_t max_counter=ids.size();
		Int_t max_clust_counter;
		TH2F* hview=NULL;
		c1->cd();
		while (last_input!="q" && id_counter<max_counter && id_counter >=0) {
			Long64_t ID=ids[id_counter];
			max_clust_counter=event_numbers[id_counter].size();
			chain->GetEntry(event_numbers[id_counter][clust_counter]);
			delete hview;
			hview=ArraysToTH2F(pixel_x, pixel_y, pixel_val);
			hview->Draw("COLZ");
			c1->Modified();
			c1->Update();
			std::cout << "Viewing id:  " << ID << ", cluster:" << clust_counter << std::endl;
			std::cin >> last_input;
			if (last_input=="p"){ 
				if (clust_counter >0) clust_counter--;
				else {
					id_counter--;
					clust_counter=event_numbers[id_counter].size()-1;
				}
			} else if (last_input=="n") {
				if (clust_counter < max_clust_counter-1) clust_counter++;
				else {
					clust_counter=0;
					id_counter++;
				}
			}
		}
	}
}
