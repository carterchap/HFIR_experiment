#include "TH2F.h"
#include "TList.h"
#include "TMath.h"

void EvalPedestal(TH2F* pedestal, TList* vals, Bool_t vid){

    Double_t maximum = 0;
    Double_t minimum = 1E9;
    Double_t mean = 0;
    Double_t rms = 0;
    
    Double_t total = 0;
    Double_t sum = 0;
    Double_t sum2 = 0;
    
    Int_t nbins_x = pedestal->GetNbinsX();
    Int_t nbins_y = pedestal->GetNbinsY();
    
    for(Int_t i=1; i<=nbins_x; i++)
        for(Int_t j=1; j<=nbins_y; j++){
            
            Double_t bc = pedestal->GetBinContent(i,j);
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
    
    AddVariable(vals, "pedestal_maximum", maximum);
    AddVariable(vals, "pedestal_minimum", minimum);
    AddVariable(vals, "pedestal_mean", mean);
    AddVariable(vals, "pedestal_rms", rms);
    
    if(vid) {
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        pedestal->SetAxisRange(minimum,maximum,"Z");
        pedestal->Draw("COLZ");
    }

}
