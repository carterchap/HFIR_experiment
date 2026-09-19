#include "TH2F.h"
#include "TList.h"
#include "TMarker.h"

void EvalCharge(TH2F* cluster, TList* vals, Bool_t vid, Bool_t is_g4_sim=false){

    AddVariable(vals, "center_x", cluster->GetXaxis()->GetBinCenter(cluster->GetNbinsX()/2.+1));
    AddVariable(vals, "center_y", cluster->GetYaxis()->GetBinCenter(cluster->GetNbinsY()/2.+1));
    
    AddVariable(vals, "charge_total", cluster->Integral());
    AddVariable(vals, "charge_mean_x", cluster->GetMean(1));
    AddVariable(vals, "charge_mean_y", cluster->GetMean(2));
    AddVariable(vals, "charge_rms_x", cluster->GetRMS(1));
    AddVariable(vals, "charge_rms_y", cluster->GetRMS(2));
    AddVariable(vals, "charge_skew_x", cluster->GetSkewness(1));
    AddVariable(vals, "charge_skew_y", cluster->GetSkewness(2));
    
    Double_t qmax = cluster->GetMaximum();
    Double_t qmin = cluster->GetMinimum();
    
    AddVariable(vals, "charge_maximum", qmax);
    
    Int_t binx, biny, binz;
    cluster->GetBinXYZ(cluster->GetMaximumBin(),binx,biny,binz);
    
    AddVariable(vals, "charge_max_x",
    cluster->GetXaxis()->GetBinCenter(binx)+0.5);
    AddVariable(vals, "charge_max_y",
    cluster->GetYaxis()->GetBinCenter(biny)+0.5);
    
    if (!is_g4_sim){
    
    	Int_t imax = qmax+0.5;
    	Int_t imin = qmin-0.5;
    
    	TH1F* charge_dist = new TH1F("charge_dist","charge_dist", imax-imin, qmin, qmax);
    	Int_t nbins_x = cluster->GetNbinsX();
   		Int_t nbins_y = cluster->GetNbinsY();
    
    	for(Int_t i=1; i<=nbins_x; i++)
    	    for(Int_t j=1; j<=nbins_y; j++){
            
     	       Double_t bc = cluster->GetBinContent(i,j);
     	       if(bc>0)  charge_dist->Fill(bc);
     	   }
    	
    	AddVariable(vals, "charge_mean", charge_dist->GetMean());
   	 	AddVariable(vals, "charge_rms", charge_dist->GetRMS());
   	 	AddVariable(vals, "charge_skew", charge_dist->GetSkewness());
    
   	 	charge_dist->Delete();
    }

    if(vid){
        
        TMarker* charge_mean = new TMarker(cluster->GetMean(1),cluster->GetMean(2),23);
        TMarker* charge_max = new TMarker(cluster->GetXaxis()->GetBinCenter(binx)+0.5,cluster->GetYaxis()->GetBinCenter(biny)+0.5,22);
        charge_mean->SetMarkerSize(2);
        charge_max->SetMarkerSize(2);
        
        charge_mean->Draw("same");
        charge_max->Draw("same");
    }
}

