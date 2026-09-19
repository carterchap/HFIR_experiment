#include "TH2F.h"
#include "TList.h"
#include "TMath.h"

void EvalSigma(TH2F* sigma, TList* vals, Bool_t vid){

    Double_t maximum = 0;
    Double_t minimum = 1E9;
    Double_t mean = 0;
    Double_t rms = 0;
    
    Double_t total = 0;
    Double_t sum = 0;
    Double_t sum2 = 0;
    
    Int_t nbins_x = sigma->GetNbinsX();
    Int_t nbins_y = sigma->GetNbinsY();
    
    for(Int_t i=1; i<=nbins_x; i++)
        for(Int_t j=1; j<=nbins_y; j++){
            
            Double_t bc = sigma->GetBinContent(i,j);
            if(bc>0){
                
                sum += bc;
                sum2 += bc*bc;
                total++;
                if(bc>maximum) maximum = bc;
                if(bc<minimum) minimum = bc;
            }
        }
    
    mean = sum/total;
    rms = TMath::Sqrt(sum2/total - mean*mean);
    
    AddVariable(vals, "sigma_maximum", maximum);
    AddVariable(vals, "sigma_minimum", minimum);
    AddVariable(vals, "sigma_mean", mean);
    AddVariable(vals, "sigma_rms", rms);
    
    if(vid) {
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        
        sigma->SetAxisRange(minimum,maximum,"Z");
        sigma->Draw("COLZ");
    }

}
