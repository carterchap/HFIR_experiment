#include "TH2F.h"
#include "TList.h"
#include "TMath.h"

void EvalRaw(TH2F* cluster, TH2F* pedestal, TList* vals, Bool_t vid){

    TH2F* raw = (TH2F*) cluster->Clone();
    raw->Add(pedestal);
    
    Int_t nsat = 0;
    
    Int_t nbins_x = pedestal->GetNbinsX();
    Int_t nbins_y = pedestal->GetNbinsY();
    double max = -1e9;
    double min = 1e9;
    
    for(Int_t i=1; i<=nbins_x; i++)
        for(Int_t j=1; j<=nbins_y; j++){
            
            Double_t bc = raw->GetBinContent(i,j);
            if(bc>65534) nsat++; //hard-coded from now but need to find way to access number of bits from here
            
            if(bc>max) max=bc;
            if(bc<min) min=bc;
        }
    
    AddVariable(vals, "raw_nsat", nsat);
    
    if(vid) {
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        
        raw->SetAxisRange(min,max,"Z");
        raw->Draw("COLZ");
    }
    
    else
        raw->Delete();

}
