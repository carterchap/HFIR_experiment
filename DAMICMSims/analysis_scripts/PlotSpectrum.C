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
#include "TError.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "DAMICPartTable.h"

enum ReadMode_t {ALL, //All files independently
								 VOLUME, //by simulation volumed
								 ISOTOPE, //by radiactive isotope
								 PART,  //by physical part (e.g. copper box, copper modules)
								 ISO_SOURCE, //radiogenic or primeval
								 ISO_CHAIN, //Source chain of isotope
								 PART_CHAIN, 
								 MEASUREMENT,
								 PART_MEASUREMENT,
								 PROPOSAL//Plot for the proposal 
};


DAMICPartTable *partTable= new DAMICPartTable();
TString dir_prefix="../";
TString data_dir="/data/simulation_output/shorts/1x100/";

TChain *data = new TChain("clusters_tree");
Bool_t loaded = false;

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
	} else if (mode==MEASUREMENT) {
		name=TString("measurement");
	} else if (mode==PART_MEASUREMENT) {
		name=TString("part_measurement");
	} else if (mode==PROPOSAL) {
		name=TString("proposal");
	}
	else {
		name=TString("error");
	}
	if (plotData) {
		name+=TString("_data");
	}
	return name;
}	

Int_t getLine(TString read_volume, TString read_isotope, vector<TH1D*> &spectra, ReadMode_t mode, Int_t mask=0) {
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
	if (mode==MEASUREMENT) key_string=partTable->GetIsotopeMeasurement(read_volume, read_isotope);
	if (mode==PART_MEASUREMENT) {
		key_string=partTable->GetPartName(read_volume);
		key_string+="_"+partTable->GetIsotopeMeasurement(read_volume, read_isotope);
	}
	if (mode==PROPOSAL) key_string=partTable->GetProposalCategory(read_volume, read_isotope);
	for (unsigned int i=0; i < spectra.size(); i++) {
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
	delete finfo_chain;
	return NumEvent;
}

bool valid_ccd(Int_t ccd) {
	if (ccd==1 || ccd==2 || ccd==3 ||ccd==4||ccd==6||ccd==11||ccd==12) {
		return true;
	} else {
		return false;
	}
}

void loadChain(Int_t bin){

	if(loaded) return;
	loaded = true;

	if(bin==100){ //1x100
	
		data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoOFF_100000s-IntW800_OS_1x100_run1/short/*.root");
		data->Add("/data/damic/snolab/processed/DAMIC_SNOLAB_RUN1/Jan2017/cryoOFF_100000s-IntW800_OS_1x100_run2/short/*.root");
	
		data->SetAlias("eraw","charge_total.fVal*2.6e-4");
		data->SetAlias("ene","eraw*( 1.039*(EXTID==1) + 0.962*(EXTID==2) + 0.999*(EXTID==3) + 0.986*(EXTID==4) + 1.034*(EXTID==6) + 0.982*(EXTID==11) + 0.993*(EXTID==12) )");
		data->SetAlias("enell","LL_eneh.fVal*2.6e-4*(  1.039*(EXTID==1) + 0.962*(EXTID==2) + 0.999*(EXTID==3) + 0.986*(EXTID==4) + 1.034*(EXTID==6) + 0.982*(EXTID==11) + 0.993*(EXTID==12) )");
		data->SetAlias("sigma","LL_sigma.fVal");
		data->SetAlias("x","charge_mean_x.fVal");
		data->SetAlias("y","charge_mean_y.fVal");
		data->SetAlias("maxx","charge_max_x.fVal");
		data->SetAlias("maxy","charge_max_y.fVal");
		data->SetAlias("rmsx","charge_rms_x.fVal");
		data->SetAlias("rmsy","charge_rms_y.fVal");
		data->SetAlias("rmsxy","sqrt((rmsx*rmsx+rmsy*rmsy)/2)");
		data->SetAlias("mask","mask_edge.fVal>0 || y<2.5");
		data->SetAlias("nsat","raw_nsat.fVal");
		data->SetAlias("ll","LL_ll.fVal");
		data->SetAlias("valid","!mask && nsat==0");
		data->SetAlias("validll","valid && LL_efit.fVal!=-1000 && ll<100 && ene>0.5 && ene<14.5");
		data->SetAlias("np","size_npixels.fVal");
		data->SetAlias("sz","size_slength.fVal*size_slength.fVal/(TMath::Sqrt(1.+size_aspect_ratio.fVal*size_aspect_ratio.fVal)*TMath::Sqrt(1.+1./(size_aspect_ratio.fVal*size_aspect_ratio.fVal)))");

	}
	else{ //default to 1x1
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
	} 

}

void PlotSpectrum(TString actFile="flex_activity.in", ReadMode_t readMode=ISOTOPE, Bool_t plotRest=true, Bool_t plotData=true, Bool_t act_from_in=false, Int_t binning = 1, Int_t nbins=500, Double_t minX=0.5, Double_t maxX=1000.5, Double_t etokev=0.00377, Double_t CCDmass = 5.8e-3, Int_t EXTNUM=1) {
    
	Int_t min2d = 0;
	Int_t max2d = 4000;
	Int_t nbins2d = 100;
	
	Int_t nbins1d = 100;
	
	Int_t nbinsSig = 20;
	Int_t minSig = 0;
	Int_t maxSig = 5;
	
	Double_t sigXYene = 20;
	
	const Int_t NCCDs=8;
	//For the plot legend
	const Int_t max_legend_size=100;
	Double_t leg_font=.02;
    
	Int_t tenkevbin = TMath::Nint((10.-minX)/((maxX-minX)/nbins));
    
	//	cout << tenkevbin << endl;
    
	TString saveFolder = "../sim_results_pdf/";
	
	loadChain(binning);
	
	TH1D *htemp;
	data->Draw("EXPTIME.fVal/24/3600 >> htemp","EXTID==1","GOFF");
	htemp=(TH1D*)gROOT->FindObject("htemp");
	Double_t expot = htemp->GetMean() * htemp->Integral();
	Int_t nbins_data=nbins;
	if (maxX > 200) nbins_data/=10;
	Double_t data_scale=nbins_data/(maxX-minX)/expot;
	vector<TChain*> chains;
	vector<TH1D*> spectra;
	vector<TH1D*> spectra_indiv;
	vector<TString> draw_cmds;
    
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

	TH1D* temp = new TH1D("temp","temp",nbins,minX,maxX);
	TString	draw_cmd="charge_total.fVal >> temp";
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
		if (!act_from_in && (read_isotope!="210a2z" && read_isotope!="210a3z")) {
		  read_activity=partTable->GetActivity(read_volume, read_isotope);
			Double_t temp_mass=partTable->GetMass(read_volume);
			if (temp_mass >0) {
				read_mass=temp_mass;
			}
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
			if (readMode==PART_MEASUREMENT) {
				spec_title=partTable->GetPartName(read_volume);
				spec_title+="_"+partTable->GetIsotopeMeasurement(read_volume, read_isotope);
			}
			if (readMode==MEASUREMENT) spec_title=partTable->GetIsotopeMeasurement(read_volume, read_isotope);
			if (readMode==PROPOSAL) spec_title=partTable->GetProposalCategory(read_volume, read_isotope);
			if (mask!=0) spec_title+=" "+partTable->GetMaskName(mask);
			spectra.push_back(new TH1D("", spec_title,nbins, minX,maxX));
			spectra_indiv.push_back(new TH1D("", spec_title, nbins, minX,maxX));
			num_lines++;
		}
		TLeaf *maskleaf=chains[num_files]->GetLeaf("mask");
 		if (read_sim_version > 0) {
			draw_cmd="charge_total.fVal >> temp";
		} else {
			draw_cmd=Form("charge_total.fVal*%f >> temp",etokev);
		}
		draw_cmds.push_back(draw_cmd);
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

	TH1D* all_spectra_indiv = new TH1D("all_spectra_indiv",actFile+Form(", CCD# %d",EXTNUM),nbins,minX,maxX);
	//TH1D* all_spectra_indiv = new TH1D("all_spectra_indiv","PNNL CCD",nbins,minX,maxX);
	TH1D* all_spectra =       new TH1D("all_spectra",actFile+", Other CCDs",nbins,minX,maxX);

	for (int i=0; i<num_lines; i++) {

		spectra_indiv[i]->GetXaxis()->SetTitle("Energy [keV]");
		spectra_indiv[i]->GetXaxis()->CenterTitle();
		spectra_indiv[i]->GetYaxis()->SetTitle("Count [/day/keV]");
		spectra_indiv[i]->GetYaxis()->CenterTitle();
		all_spectra_indiv->Add(spectra_indiv[i],1);
		spectra_indiv[i]->SetTitle((TString)spectra_indiv[i]->GetTitle()+Form(", %f dru",spectra_indiv[i]->Integral(1,tenkevbin)/tenkevbin/CCDmass));

		spectra[i]->GetXaxis()->SetTitle("Energy [keV]");
		spectra[i]->GetYaxis()->SetTitle("Count [/day/keV/CCD]");
		spectra[i]->SetTitle((TString)spectra[i]->GetTitle()+Form(", %f dru",spectra[i]->Integral(1,tenkevbin)/tenkevbin/CCDmass));
		all_spectra->Add(spectra[i],1);
	}
	gStyle->SetOptStat(0);
    
    
	TCanvas *c1 = new TCanvas();
	c1->cd();

	TH1D* data_ene=new TH1D("data_ene", "data",nbins_data,minX,maxX);
	if (plotData) {
		data->Draw("ene >> data_ene", Form("EXTID!=%d",EXTNUM),"GOFF");
		//		data_ene=(TH1D*) gROOT->FindObject("data_temp");
		data_ene->Sumw2();
		//		data_ene->Scale(all_spectra->GetMaximum()/data_ene->GetMaximum());
		data_ene->Scale(data_scale*1./7);
		//		data_ene->SetLineStyle(7);
		data_ene->SetLineWidth(2);
		data_ene->Draw("E1");
	}

	if (plotData) {
		all_spectra->Draw("hist SAME");
		data_ene->GetXaxis()->SetTitle("Energy [keV]");
		data_ene->GetYaxis()->SetTitle("Count [/day/keV/CCD]");
		data_ene->GetXaxis()->CenterTitle();
		data_ene->GetYaxis()->CenterTitle();
		
	} else {
		all_spectra->Draw("hist");
		all_spectra->GetXaxis()->SetTitle("Energy [keV]");
		all_spectra->GetYaxis()->SetTitle("Count [/day/keV/CCD]");

	}
	cout << "Total :" << all_spectra->Integral() << endl;
	for (Int_t i=0; i < num_lines; i++ ) {
		spectra[i]->Draw("HIST SAME");
		spectra[i]->SetLineWidth(2);
		spectra[i]->SetLineColor(colors[i%27]);
	}

	Double_t y1=.88-(.38/12*num_lines);
	Double_t x1=.6;

	c1->SetLogy();
	if (num_lines < max_legend_size) {
		TLegend *l1 = c1->BuildLegend(x1,y1,.88,.88);
		l1->SetTextSize(leg_font);
		l1->SetBorderSize(0);
	} else {
		cout << "Hiding legend..." << endl;
	}
	TLatex *t1 = new TLatex(0.25,0.8,Form("dru (.5-10 keV): %f",all_spectra->Integral(1,tenkevbin)/tenkevbin/CCDmass));  
	t1->SetNDC(kTRUE);
	t1->SetTextSize(0.035);
	t1->Draw();
	c1->SaveAs(saveFolder+actFile+"_"+readModeToString(readMode,plotData)+"_OtherCCDs_spectrum.svg");
    
	TCanvas *c2 = new TCanvas();
	c2->cd();
	TH1D* data_ene_indiv=new TH1D("data_ene_indiv", "data", nbins_data, minX, maxX);
	//TH1D* data_ene_indiv=new TH1D("data_ene_indiv", "Data", nbins_data, minX, maxX);
	if (plotData) {
		data->Draw("ene >> data_ene_indiv", Form("EXTID==%d",EXTNUM),"GOFF");
		data_ene_indiv->Sumw2();
		data_ene_indiv->Scale(data_scale);
		//		data_ene_indiv->SetLineStyle(7);
		data_ene_indiv->SetLineWidth(2);
		data_ene_indiv->Draw("E1");
	}
	if (plotData) {
		all_spectra_indiv->Draw("hist SAME");
		data_ene_indiv->GetXaxis()->SetTitle("Energy [keV]");
		data_ene_indiv->GetYaxis()->SetTitle("Count [/day/keV/CCD]");
		data_ene_indiv->GetXaxis()->CenterTitle();
		data_ene_indiv->GetYaxis()->CenterTitle();
	} else {
		all_spectra_indiv->Draw("hist");
		all_spectra_indiv->GetXaxis()->SetTitle("Energy [keV]");
		all_spectra_indiv->GetYaxis()->SetTitle("Count [/day/keV/CCD]");	
	}

	for (Int_t i=0; i < num_lines; i++ ) {
		spectra_indiv[i]->Draw("HIST SAME");
		spectra_indiv[i]->SetLineWidth(2);
		spectra_indiv[i]->SetLineColor(colors[i%27]);        
		if (TString(spectra_indiv[i]->GetTitle()).Contains("Upper")) {
			std::cout << "Setting line to dashed" << std::endl;
			spectra_indiv[i]->SetLineStyle(9);
		}
	}

	TLatex *t2 = new TLatex(0.25,0.8,Form("dru (.5-10 keV): %f",all_spectra_indiv->Integral(1,tenkevbin)/tenkevbin/CCDmass));  
	c2->SetLogy();
	if (num_lines < max_legend_size) {
		TLegend *l2 = c2->BuildLegend(x1,y1,.88,.88);
		l2->SetTextSize(leg_font);
		l2->SetBorderSize(0);
	}
	t2->SetNDC(kTRUE);
	t2->SetTextSize(0.035);
	t2->Draw();
	c2->SaveAs(saveFolder+actFile+"_"+readModeToString(readMode,plotData)+"_PNNL_spectrum.svg");

    
	TH2F* temp2d = new TH2F("temp2d","temp2d",nbins2d,min2d,max2d,nbins2d,min2d,max2d);
	TH1D* temp1d = new TH1D("temp1d","temp1d",nbins1d,min2d,max2d);
	TH1D* temp1dx = new TH1D("temp1dx","temp1dx",nbins1d,4275,8275);
	TH1D* tempSig = new TH1D("tempSig","tempSig",nbinsSig,minSig,maxSig);

	TCanvas *c_byext=new TCanvas();
	c_byext->cd();
	c_byext->Draw();
	Double_t dru_by_ext[12];
	Double_t dru_errors[12];
	TH1F *dru_extension=new TH1F("dru_extension", "DRU By Extension", 12,1,12);
	dru_extension->Sumw2();
	gErrorIgnoreLevel = kFatal;
	for (Int_t i=0; i < 12; i++) {
		dru_by_ext[i]=0;
		dru_errors[i]=0;
	};
	for (Int_t i=0; i <num_files;i++) {
		for (Int_t k=0; k < 12; k++) { 
			if (valid_ccd(k+1)) {
				chains[i]->Draw(draw_cmds[i],Form("EXTID==%d", k+1),"GOFF");
				dru_by_ext[k]+=temp->Integral(1,tenkevbin)/tenkevbin/CCDmass*scales[i];
				//TODO: double checked this math
				dru_errors[k]+=TMath::Sqrt(temp->Integral(1,tenkevbin))/tenkevbin/CCDmass*scales[i];
			}
		}
	}

	for (int j=0; j<12;j++) {
		if (valid_ccd(j+1)) {
			std::cout << "Ext " << j+1 << " DRU: " << std::setprecision(4) << dru_by_ext[j] << " +/- " << std::setprecision(2) << dru_errors[j] << std::endl;
			dru_extension->SetBinContent(j+1, dru_by_ext[j]);
			dru_extension->SetBinError(j+1, dru_errors[j]);
		}
	};
	
	dru_extension->Draw("E1");
if (plotRest){

    
	TCanvas *c3 = new TCanvas();
	c3->cd();
	c3->Draw();
	TH2F *XY_indiv[num_files];
	TH2F *all_XY_indiv = new TH2F("all_XY_indiv",Form(actFile + ", CCD# %d",EXTNUM),nbins2d,min2d,max2d,nbins2d,min2d,max2d);
	for (int i=0; i<num_files; i++) {
		XY_indiv[i] = new TH2F("","",nbins2d,min2d,max2d,nbins2d,min2d,max2d);
		chains[i]->Draw("charge_mean_y.fVal:charge_mean_x.fVal >> temp2d",Form("EXTID==%d",EXTNUM),"GOFF");
		XY_indiv[i]->Add(temp2d,scales[i]);
		temp2d->Reset();
		all_XY_indiv->Add(XY_indiv[i],1);
	}
	all_XY_indiv->GetXaxis()->SetTitle("x");
	all_XY_indiv->GetYaxis()->SetTitle("y");
	all_XY_indiv->Draw("COLZ");
	c3->SaveAs(saveFolder+actFile+"_PNNL_XY.svg");
    
	TCanvas *c4 = new TCanvas();
	c4->cd();
	c4->Draw();
	TH2F *XY[num_files];
	TH2F *all_XY = new TH2F("all_XY",actFile + ", Other CCDs",nbins2d,min2d,max2d,nbins2d,min2d,max2d);
	for (int i=0; i<num_files; i++) {
		XY[i] = new TH2F("","",nbins2d,min2d,max2d,nbins2d,min2d,max2d);
		chains[i]->Draw("charge_mean_y.fVal:charge_mean_x.fVal >> temp2d",Form("EXTID!=%d",EXTNUM),"GOFF");
		XY[i]->Add(temp2d,scales[i]);
		temp2d->Reset();
		all_XY->Add(XY[i],1);
	}
	all_XY->GetXaxis()->SetTitle("X [pixels]");
	all_XY->GetYaxis()->SetTitle("Y [pixels]");
	all_XY->Draw("COLZ");
	c4->SaveAs(saveFolder+actFile+"_OtherCCDs_XY.svg");
	
	TCanvas *c5 = new TCanvas();
	c5->cd();
	c5->Draw();
	TH1D *dX[num_files];
	TH1D *dY[num_files];
	TH1D *all_dX = new TH1D("all_dX",actFile + ", X, Other CCDs",nbins1d, min2d, max2d);
	TH1D *all_dY = new TH1D("all_dY",actFile + ", Y, Other CCDs",nbins1d, min2d, max2d);
	TH1D *DATA_dX = new TH1D("DATA_dX","DATA X",nbins1d, min2d, max2d);
	TH1D *DATA_dY = new TH1D("DATA_dY","DATA Y",nbins1d, min2d, max2d);
	for (int i=0; i<num_files; i++) {
		dX[i] = new TH1D("","",nbins1d, min2d, max2d);
		chains[i]->Draw("charge_mean_x.fVal >> temp1d",Form("EXTID!=%d",EXTNUM),"GOFF");
		dX[i]->Add(temp1d,scales[i]);
		temp1d->Reset();
		all_dX->Add(dX[i],1);
		dY[i] = new TH1D("","",nbins1d, min2d, max2d);
		chains[i]->Draw("charge_mean_y.fVal >> temp1d",Form("EXTID!=%d",EXTNUM),"GOFF");
		dY[i]->Add(temp1d,scales[i]);
		temp1d->Reset();
		all_dY->Add(dY[i],1);
	}

	data->Draw("charge_mean_x.fVal >> temp1dx",Form("EXTID!=%d",EXTNUM),"GOFF");
	for (int i=1; i<=nbins1d; i++){
		DATA_dX->SetBinContent(i,temp1dx->GetBinContent(i));
	}
	temp1dx->Reset();
	data->Draw("charge_mean_y.fVal >> temp1d",Form("EXTID!=%d",EXTNUM),"GOFF");
	DATA_dY->Add(temp1d,1);
	temp1d->Reset();
	all_dX->Scale(1./all_dX->Integral());
	all_dY->Scale(1./all_dY->Integral());
	DATA_dX->Scale(1./DATA_dX->Integral());
	DATA_dY->Scale(1./DATA_dY->Integral());
	all_dX->GetXaxis()->SetTitle("X,Y");
	all_dX->GetYaxis()->SetTitle("Normalized Count");
	all_dX->SetLineWidth(2);
	all_dX->Draw();

	all_dX->GetYaxis()->SetRangeUser(0.,0.03);
	all_dY->SetLineColor(2);
	all_dY->SetLineWidth(2);
	all_dY->Draw("same");
	DATA_dY->SetLineColor(2);
	DATA_dX->SetLineWidth(3);
	DATA_dY->SetLineWidth(3);
	DATA_dX->SetLineStyle(2);
	DATA_dY->SetLineStyle(2);
	DATA_dX->Draw("same");
	DATA_dY->Draw("same");
	TLegend *l5 = c5->BuildLegend();
	l5->SetBorderSize(0);
	all_dX->SetTitle(actFile + ", Other CCDs, X&Y Dist.");
	c5->SaveAs(saveFolder+actFile+"_OtherCCD_1dXY.svg");
	
	TCanvas *c6 = new TCanvas();
	c6->cd();
	c6->Draw();
	TH1D *sX[num_files];
	TH1D *sY[num_files];
	TH1D *all_sX = new TH1D("all_sX",actFile + ", Sigma X, Other CCDs",nbinsSig, minSig, maxSig);
	TH1D *all_sY = new TH1D("all_sY",actFile + ", Sigma Y, Other CCDs",nbinsSig, minSig, maxSig);
	TH1D *DATA_sX = new TH1D("DATA_sX","DATA Sigma X",nbinsSig, minSig, maxSig);
	TH1D *DATA_sY = new TH1D("DATA_sY","DATA Sigma Y",nbinsSig, minSig, maxSig);
	for (int i=0; i<num_files; i++) {
		sX[i] = new TH1D("","",nbinsSig, minSig, maxSig);
		chains[i]->Draw("charge_rms_x.fVal >> tempSig",Form("EXTID!=%d && charge_total.fVal>0.004 && charge_total.fVal<%f" ,EXTNUM, sigXYene),"GOFF");
		sX[i]->Add(tempSig,scales[i]);
		tempSig->Reset();
		all_sX->Add(sX[i],1);
		sY[i] = new TH1D("","",nbinsSig, minSig, maxSig);
		chains[i]->Draw("charge_rms_y.fVal >> tempSig",Form("EXTID!=%d && charge_total.fVal>0.004 && charge_total.fVal<%f",EXTNUM, sigXYene),"GOFF");
		sY[i]->Add(tempSig,scales[i]);
		tempSig->Reset();
		all_sY->Add(sY[i],1);
	}
	data->Draw("charge_rms_x.fVal >> tempSig",Form("EXTID!=%d && ene<%f && valid",EXTNUM, sigXYene),"GOFF");
	DATA_sX->Add(tempSig,1);
	tempSig->Reset();
	data->Draw("charge_rms_y.fVal >> tempSig",Form("EXTID!=%d && ene<%f && valid",EXTNUM, sigXYene),"GOFF");
	DATA_sY->Add(tempSig,1);
	tempSig->Reset();
	all_sX->Scale(1./all_sX->Integral());
	all_sY->Scale(1./all_sY->Integral());
	DATA_sX->Scale(1./DATA_sX->Integral());
	DATA_sY->Scale(1./DATA_sY->Integral());
	all_sX->GetXaxis()->SetTitle("Sigma X,Y");
	all_sX->GetYaxis()->SetTitle("Normalized Count");
	all_sX->SetLineWidth(2);
	all_sX->Draw();
	all_sY->SetLineColor(2);
	all_sX->SetLineWidth(2);
	all_sY->Draw("same");
	DATA_sY->SetLineColor(2);
	DATA_sX->SetLineWidth(3);
	DATA_sY->SetLineWidth(3);
	DATA_sX->SetLineStyle(2);
	DATA_sY->SetLineStyle(2);
	DATA_sX->Draw("same");
	DATA_sY->Draw("same");
	TLegend *l6 = c6->BuildLegend();
	l6->SetBorderSize(0);
	all_sX->SetTitle(Form(actFile + ", Sigma XY with ene < %f, Other CCD", sigXYene));
	c6->SaveAs(saveFolder+actFile+"_OtherCCD_sigmaXY.svg");
	
	TCanvas *c7 = new TCanvas();
	c7->cd();
	c7->Draw();
	TH1D *PNNL_dX = new TH1D("PNNL_dX",actFile + ", X, PNNL",nbins1d, min2d, max2d);
	TH1D *PNNL_dY = new TH1D("PNNL_dY",actFile + ", Y, PNNL",nbins1d, min2d, max2d);
	TH1D *PDATA_dX = new TH1D("PDATA_dX","DATA X, PNNL",nbins1d, min2d, max2d);
	TH1D *PDATA_dY = new TH1D("PDATA_dY","DATA Y, PNNL",nbins1d, min2d, max2d);
	for (int i=0; i<num_files; i++) {
		chains[i]->Draw("charge_mean_x.fVal >> temp1d",Form("EXTID==%d",EXTNUM),"GOFF");
		PNNL_dX->Add(temp1d,scales[i]);
		temp1d->Reset();
		chains[i]->Draw("charge_mean_y.fVal >> temp1d",Form("EXTID==%d",EXTNUM),"GOFF");
		PNNL_dY->Add(temp1d,scales[i]);
		temp1d->Reset();
	}
	data->Draw("charge_mean_x.fVal >> temp1dx",Form("EXTID==%d",EXTNUM),"GOFF");
	for (int i=1; i<=nbins1d; i++){
		PDATA_dX->SetBinContent(i,temp1dx->GetBinContent(i));
	}
	temp1dx->Reset();
	data->Draw("charge_mean_y.fVal >> temp1d",Form("EXTID==%d",EXTNUM),"GOFF");
	PDATA_dY->Add(temp1d,1);
	temp1d->Reset();
	PNNL_dX->Scale(1./PNNL_dX->Integral());
	PNNL_dY->Scale(1./PNNL_dY->Integral());
	PDATA_dX->Scale(1./PDATA_dX->Integral());
	PDATA_dY->Scale(1./PDATA_dY->Integral());
	PNNL_dY->SetLineColor(2);
	PNNL_dX->SetLineWidth(2);
	PNNL_dY->SetLineWidth(2);
	PNNL_dX->Draw();
	PNNL_dX->GetYaxis()->SetRangeUser(0.,0.03);
	PNNL_dY->Draw("same");
	PDATA_dY->SetLineColor(2);
	PDATA_dX->SetLineWidth(3);
	PDATA_dY->SetLineWidth(3);
	PDATA_dX->SetLineStyle(2);
	PDATA_dY->SetLineStyle(2);
	PDATA_dX->Draw("same");
	PDATA_dY->Draw("same");
	TLegend *l7 = c7->BuildLegend();
	l7->SetBorderSize(0);
	all_dX->SetTitle(actFile + ", X&Y Dist.");
	c7->SaveAs(saveFolder+actFile+"_PNNL_1dXY.svg");
	
	TCanvas *c8 = new TCanvas();
	c8->cd();
	c8->Draw();
	TH1D *PNNL_sX = new TH1D("PNNL_sX",actFile + ", Sigma X, PNNL",nbinsSig, minSig, maxSig);
	TH1D *PNNL_sY = new TH1D("PNNL_sY",actFile + ", Sigma Y, PNNL",nbinsSig, minSig, maxSig);
	TH1D *PDATA_sX = new TH1D("PDATA_sX","DATA, Sigma X, PNNL",nbinsSig, minSig, maxSig);
	TH1D *PDATA_sY = new TH1D("PDATA_sY","DATA, Sigma Y, PNNL",nbinsSig, minSig, maxSig);
	for (int i=0; i<num_files; i++) {
		chains[i]->Draw("charge_rms_x.fVal >> tempSig",Form("EXTID==%d && charge_total.fVal>0.004 && charge_total.fVal<%f" ,EXTNUM, sigXYene),"GOFF");
		PNNL_sX->Add(tempSig,scales[i]);
		tempSig->Reset();
		chains[i]->Draw("charge_rms_y.fVal >> tempSig",Form("EXTID==%d && charge_total.fVal>0.004 && charge_total.fVal<%f" ,EXTNUM, sigXYene),"GOFF");
		PNNL_sY->Add(tempSig,scales[i]);
		tempSig->Reset();	
	}
	data->Draw("charge_rms_x.fVal >> tempSig",Form("EXTID!=%d && ene<%f && valid",EXTNUM, sigXYene),"GOFF");
	PDATA_sX->Add(tempSig,1);
	tempSig->Reset();
	data->Draw("charge_rms_y.fVal >> tempSig",Form("EXTID!=%d && ene<%f && valid",EXTNUM, sigXYene),"GOFF");
	PDATA_sY->Add(tempSig,1);
	tempSig->Reset();
	PNNL_sX->Scale(1./PNNL_sX->Integral());
	PNNL_sY->Scale(1./PNNL_sY->Integral());
	PDATA_sX->Scale(1./PDATA_sX->Integral());
	PDATA_sY->Scale(1./PDATA_sY->Integral());
	PNNL_sY->SetLineColor(2);
	PNNL_sX->SetLineWidth(2);
	PNNL_sY->SetLineWidth(2);
	PNNL_sX->Draw();
	PNNL_sY->Draw("same");
	PDATA_sY->SetLineColor(2);
	PDATA_sX->SetLineWidth(3);
	PDATA_sY->SetLineWidth(3);
	PDATA_sX->SetLineStyle(2);
	PDATA_sY->SetLineStyle(2);
	PDATA_sX->Draw("same");
	PDATA_sY->Draw("same");
	TLegend *l8 = c8->BuildLegend();
	l8->SetBorderSize(0);
	PNNL_sX->SetTitle(Form(actFile + ", Sigma XY with ene < %f, PNNL", sigXYene));
	c8->SaveAs(saveFolder+actFile+"_PNNL_sigmaXY.svg");
   
}
    
};
