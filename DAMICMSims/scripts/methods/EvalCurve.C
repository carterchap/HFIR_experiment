#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TF1.h"

void EvalCurve(TH2F* cluster, TList* vals, Bool_t vid)
{
	Int_t nbinsx = cluster->GetNbinsX();
    Int_t nbinsy = cluster->GetNbinsY();
	Double_t correlationfactor = cluster->GetCorrelationFactor(1,2);
	TF1 *f1 = new TF1("f1","pol1(0)",0,10000);
	
    //fit to a line
    Double_t slope = cluster->GetCovariance()/cluster->GetRMS(1)/cluster->GetRMS(1);
    Double_t intercept = cluster->GetMean(2) - slope*cluster->GetMean(1);
    
    //set these parameters in f1
    f1->SetParameters(intercept,slope);
    
	Double_t totrms = 0;
    Double_t totbc = 0;
    
	for(Int_t i=1;i<=nbinsx;i++)
		for(Int_t j=1;j<=nbinsy;j++)
			if(cluster->GetBinContent(i,j)>0) 
			{	
				Double_t bc = cluster->GetBinContent(i,j);
                Double_t x = cluster->GetXaxis()->GetBinCenter(i);
				Double_t y = cluster->GetYaxis()->GetBinCenter(j);
                
                Double_t l = y - slope*x - intercept;
                Double_t s = l/TMath::Sqrt(1.+slope*slope);
                
                totbc += bc;
                totrms += bc*s*s;
            }
        
	
    Double_t rms = TMath::Sqrt(totrms/totbc);
    
	AddVariable(vals, "curve_track_rms",rms);
    AddVariable(vals, "curve_correlation_factor",correlationfactor);
    
    if(vid) f1->Draw("same");
    else f1->Delete();
}
