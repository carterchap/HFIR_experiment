#include "TList.h"
#include "TParameter.h"
#include "TArray.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <algorithm>

Int_t min_entries_fit = 100;

void AddVariable(TList* vals, TString vname, Double_t vval){
    
    if(vals->Contains(vname)){
        
        TParameter<double>* vpar = (TParameter<double>*) vals->FindObject(vname);
        vpar->SetVal(vval);
    }
    
    else{
        
        TParameter<double>* vpar = new TParameter<double>(vname,vval);
        vals->Add(vpar);
    }
}

TH2F* ArraysToTH2F(TArrayD* x, TArrayD* y, TArrayD* val, Int_t npixels, Bool_t LL){
    
    Int_t min_x = 9999999;
    Int_t max_x = 0;
    Int_t min_y = 9999999;
    Int_t max_y = 0;
    
    Int_t array_lim = 0;
    if ((npixels == 0) || LL)
        array_lim = x->GetSize();
    else
        array_lim = npixels;
    
    //for(Int_t i=0; i<x->GetSize(); i++){
    // Added by KR
    for (Int_t i=0; i<array_lim; i++){
    
        if(x->At(i)>max_x) max_x = (Int_t) x->At(i);
        if(x->At(i)<min_x) min_x = (Int_t) x->At(i);
        if(y->At(i)>max_y) max_y = (Int_t) y->At(i);
        if(y->At(i)<min_y) min_y = (Int_t) y->At(i);
    }
    
    
    TH2F* h2 = new TH2F("h2","h2",max_x-min_x+1,min_x-1,max_x,max_y-min_y+1,min_y-1,max_y);
   
    
    // Modified by KR
    for(Int_t i=0; i<array_lim; i++){
        h2->SetBinContent(((Int_t) x->At(i))-min_x+1, ((Int_t) y->At(i))-min_y+1, val->At(i));
    }
    return h2;
}

void SetMinEntriesFit(Int_t val){
    
    min_entries_fit = val;
    std::cout << "Min entries for Gaussian fit set to " << min_entries_fit << std::endl;
}

//Histogram contents must be in counts
void FitToGaus(TH1F* h, Double_t* mean, Double_t* mean_unc, Double_t* sigma, Double_t* sigma_unc, Double_t* chi2_ndf){
    
    //Gaussian function to fit number distribution
    Int_t nbinsx = h->GetNbinsX();
    
    h->SetAxisRange(h->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nbinsx),"X");
    const Int_t hentries = h->Integral();
    Double_t* values = new Double_t[hentries];
    
    //set starting parameters to reasonable values from histogram
    Double_t mean_guess, sigma_guess;
    Double_t maximum = 0;
    Int_t k = 0;
    
    for(Int_t i=1; i<=nbinsx; i++){
        
        Double_t bc = h->GetBinContent(i);
        if(bc>maximum) maximum = bc;
        
        for(Int_t j=0; j<bc; j++){
            
            values[k] = h->GetBinCenter(i);
            k++;
        }
    }
       
    std::sort(values,values+hentries);
    
    if(hentries%2==1)
        mean_guess = values[(hentries-1)/2];
    else
        mean_guess = (values[hentries/2] + values[hentries/2-1])/2.;
    
    //Calculate MAD and use to guess RMS
    for(Int_t i=0; i<hentries; i++)
        values[i] = TMath::Abs(values[i]-mean_guess);
    
    std::sort(values,values+hentries);
    
    if(hentries%2==1)
        sigma_guess = values[(hentries-1)/2];
    else
        sigma_guess = (values[hentries/2] + values[hentries/2-1])/2.;
    
    sigma_guess *= 1.4826;
    
    delete[] values;
    
    /*
    //Old way to calculate guess for RMS
    Int_t sbin = 0;
    Double_t mean_bin = h->FindBin(mean_guess);
    Double_t running_total_up = h->GetBinContent(mean_bin);
    Double_t running_total_down = h->GetBinContent(mean_bin);
    
    while(running_total_up<0.34*hentries+0.5 && running_total_down<0.34*hentries+0.5){
        sbin++;
        
        if(sbin>=mean_bin && sbin+mean_bin>nbinsx) break;
        if(mean_bin-sbin>=1) running_total_down += h->GetBinContent(mean_bin-sbin);
        if(mean_bin+sbin<=nbinsx) running_total_up += h->GetBinContent(mean_bin+sbin);
    }

    sigma_guess = (mean_guess - h->GetBinCenter(mean_bin-sbin));
    */
    
    //std::cout << maximum << " " << mean_guess << " " << sigma_guess << std::endl;

    Double_t min_range = mean_guess - 2*sigma_guess;
    Double_t max_range = mean_guess + 2*sigma_guess;
    
    //When there are too few entries in the histogram do not do the fit
    //below is robust estimate of mean and sigma
    if(hentries<min_entries_fit){
        
        h->SetAxisRange(min_range,max_range,"X");
        *mean = h->GetMean();
        *mean_unc = h->GetRMS()/TMath::Sqrt(h->Integral());
        *sigma = sigma_guess;
        *sigma_unc = -1;
        *chi2_ndf = -1;
        return;
    }
    
    TF1* g = new TF1("g","gaus(0)",h->GetBinCenter(1), h->GetBinCenter(nbinsx));
    g->SetLineColor(2);
    
    g->SetParameter(0, maximum);
    g->SetParameter(1, mean_guess);
    g->SetParameter(2, sigma_guess);
    
    //Set parameter limits
    g->SetParLimits(0, 0, 3*maximum);
    g->SetParLimits(1, mean_guess-2*sigma_guess, mean_guess+2*sigma_guess);
    g->SetParLimits(2, 0, 3*sigma_guess);
    
    //perform the fit
    h->Fit("g","QL","",min_range,max_range);
    
    *mean = g->GetParameter(1);
    *mean_unc = g->GetParError(1);
    *sigma = g->GetParameter(2);
    *sigma_unc = g->GetParError(2);
    *chi2_ndf = g->GetChisquare()/g->GetNDF();
    if(*chi2_ndf==TMath::Infinity()) *chi2_ndf = -1;
    g->Delete();
    
    //std::cout << *mean << " " << *sigma << std::endl;
}
