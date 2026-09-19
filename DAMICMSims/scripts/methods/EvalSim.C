#include "TH2F.h"
#include "TList.h"
#include "TMath.h"

void EvalSim(TH2F* cluster, TString prefix, TList* vals, Bool_t vid){

    Double_t total = 0;
    Double_t nfill = 0;
    Int_t nx = cluster->GetNbinsX();
    Int_t ny = cluster->GetNbinsY();
    
    for(Int_t i=1; i<=nx; i++)
        for(Int_t j=1; j<=ny; j++)
        {
        
        Double_t cont = cluster->GetBinContent(i,j);        
        if(cont>0){
            
            total+=cont;
            nfill++;
        }
    }
    
    TString vname = prefix;
    vname += "_total";
    AddVariable(vals, vname, total);
    
    if(nfill>0) total /= nfill;
    
    vname = prefix;
    vname += "_mean";
    AddVariable(vals, vname, total);
    
    if(vid){
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        
        cluster->Draw("COLZ");
    }
}
