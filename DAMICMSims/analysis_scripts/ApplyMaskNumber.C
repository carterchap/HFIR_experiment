#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <iostream>
#include <TCanvas.h>
#include <TList.h>
#include <vector>
#include <stdio.h>
#include <TList.h>
#include <TGraph.h>
#include <TMath.h>
#include <TROOT.h>
#include <TEllipse.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TChain.h>
#include <iostream>
#include <fstream>

const char *ext=".root";

Bool_t InsideEllipse(Double_t x, Double_t y, Double_t X, Double_t Y, Double_t R){
    
    Double_t checkval;
    
    checkval = ((x-X)*(x-X) + (y-Y)*(y-Y))/R/R;
    
    if (checkval <= 1.)
        return 1;
    
    return 0;
}

void ApplyMaskNumber(TString maskFile="masks/base_screws.mask", TString inputDirectory="/data/simulation_output/CopperBasePlatePV/", Double_t defaultRad = 1.75){

    Double_t x,y,z; //Read from maskFile
    Int_t maskN;
    
    vector<Double_t> Xmask, Ymask, Zmask; //Stores from maskFile
    vector<Int_t> maskNum;
    
    Double_t Xdep, Ydep, Zdep; //Read from inputFile
    Bool_t check_inside;
    
    Int_t mask; //variable for the new branch

    // Load Mask File
    Int_t num_lines = 0;
    ifstream file(maskFile);
    while (file >> x >> y >> z >> maskN){
        num_lines+=1;
        Xmask.push_back(x);
        Ymask.push_back(y);
        Zmask.push_back(z);
        maskNum.push_back(maskN);
        //cout << "Screw: " << Xmask[num_lines-1] << ", " << Ymask[num_lines-1] << ", " << Zmask[num_lines-1] << endl;
    }
    
    // Ellipses for plotting
    TEllipse *el[num_lines];
    for (Int_t i = 0; i<num_lines; i++){
        el[i] = new TEllipse(Xmask[i],Ymask[i],defaultRad);
    }
    
    TSystemDirectory dir(inputDirectory, inputDirectory);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
              
                // Load Input Files
                TFile *f = new TFile(inputDirectory+fname, "update");
                TTree *ClustersTree = (TTree*)f->Get("clusters_tree");
                TBranch *newBranch = ClustersTree->Branch("mask", &mask,"mask/I");
                const Int_t entries = ClustersTree->GetEntries();
                
                ClustersTree->SetBranchAddress("Xdep",&Xdep);
                ClustersTree->SetBranchAddress("Ydep",&Ydep);
                ClustersTree->SetBranchAddress("Zdep",&Zdep);
                
                for (Int_t n = 0; n<entries; n++){
                    ClustersTree->GetEntry(n);
                    mask = 0;
                    for (Int_t i = 0; i<num_lines; i++){
                        check_inside = InsideEllipse(Xdep, Ydep, Xmask[i],Ymask[i],defaultRad);
                        if (check_inside){
                            //cout << "Inside: " << Xdep << ", " << Ydep << endl;
                            mask = maskNum[i];
                        }
                    }
                    newBranch->Fill();
                }
                
                ClustersTree->Write("",TObject::kOverwrite);
                
                cout << fname << " :" << entries << endl;
                
                /*TCanvas *c1 = new TCanvas();
                c1->cd();
                ClustersTree->Draw("Ydep:Xdep","","GOFF");
                TGraph *gr = new TGraph(ClustersTree->GetSelectedRows(), ClustersTree->GetV2(), ClustersTree->GetV1());
                gr->Draw("AP*");
                gr->GetXaxis()->SetRangeUser(-60,60);
                gr->GetYaxis()->SetRangeUser(-60,60);
                for (Int_t i = 0; i<num_lines; i++){
                    el[i]->SetFillColorAlpha(0,0);
                    el[i]->Draw("same");
                }*/
                
            }
        }
    }
    //gSystem->FreeDirectory(inputDirectory);
    
    

    
}
