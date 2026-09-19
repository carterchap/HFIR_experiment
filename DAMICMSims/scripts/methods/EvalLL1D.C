#include "TH2F.h"
#include "TArrayD.h"
#include "TTree.h"
#include "TFile.h"
#include "TPrincipal.h"
#include "TBrowser.h"
#include "TF1.h"
#include "TF2.h"
#include "TList.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TGraph.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TParameter.h"
#include "TBox.h"

#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <vector>

TH1F* just1Dcluster;

//to fit for sigma
double LogLikelihood1D(const double *xx){
    
    const Double_t x = xx[0];
    const Double_t s = xx[1];
    const Double_t q = xx[2];
    const Double_t o = xx[3];
    
    Double_t bot = TMath::Sqrt2()*s;
    
    Double_t ll = 0;
    
    for(Int_t i=1; i<=nx; i++){
        
        Double_t obs = just1Dcluster->GetBinContent(i) - o;
        
        if(obs!=0){
            
            //expected number of electrons
            Double_t ne = q*0.5*c*(TMath::Erf((i+bxlow+0.5-x)/bot) - TMath::Erf((i+bxlow-0.5-x)/bot));
            
            if(ne==0 && has_noise) {
                
                ll -= (obs)*(obs)/(2.*sadc*sadc) + TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
            }
            
            else if (ne>0 && has_noise){
                
                Double_t kp = maxk(obs,ne);
                
                //round to integer and find maximum contribution to ll
                Int_t ko = kp + 0.5;
                if(ko<0) ko = 0;
                Double_t llmax = ko*TMath::Log(ne) - ne - TMath::LnGamma(ko+1) - (obs-ko/c)*(obs-ko/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                
                //now look around
                Int_t k = ko;
                Double_t llk = llmax;
                Double_t lli = llmax;
                
                //to the left
                while(llmax-llk<4){
                    
                    k--;
                    if(k<0) break; //smallest value of k==0
                    llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1) - (obs-k/c)*(obs-k/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                    lli += TMath::Log(1. + TMath::Exp(llk-lli));
                }
                
                //to the right
                k=ko;
                llk = llmax;
                while(llmax-llk<4){
                    
                    k++;
                    llk = k*TMath::Log(ne) - ne - TMath::LnGamma(k+1) - (obs-k/c)*(obs-k/c)/(2.*sadc*sadc) - TMath::Log(TMath::Sqrt(2.*TMath::Pi()*sadc*sadc));
                    lli += TMath::Log(1. + TMath::Exp(llk-lli));
                }
                
                //now increase the log likelihood
                ll += lli;
            }
            
            else{
                
                Int_t k = (Int_t)(obs*c + 0.5);
                ll += k*TMath::Log(ne) - ne - TMath::LnGamma(k+1);
            }

        }
    }
    
    return -1.*ll;
    
}

void Fit1D(TH2F* cluster, TH2F* sigma, Double_t meanx, Int_t LLwindow, Bool_t zero_base, Bool_t is_g4_sim, Double_t ene_lim, Double_t sim_noise){
    
    nx = 2*LLwindow + 1;
    
    Int_t center_bin_x = cluster->GetXaxis()->FindBin(meanx);
    Int_t bin_x_max = TMath::Min(center_bin_x+LLwindow, cluster->GetNbinsX());
    
    Int_t bx_start = center_bin_x - LLwindow;
    bxlow = cluster->GetXaxis()->GetBinLowEdge(bx_start);
    
    just1Dcluster = new TH1F("H1D","",nx,bxlow,bxlow+nx);
    
    Double_t qtot = 0;
    Double_t ibc = 0;
    
    Double_t qmax = 0;
    Int_t imax = 0;
    
    Int_t numbins = 0;
    for (Int_t i=bx_start; i <= bin_x_max; i++){
        
        //cout << cluster->GetBinContent(i,1) << endl;
        ibc = cluster->GetBinContent(i,1) - base;
        if (has_noise && is_g4_sim)
        	sadc += (sim_noise*ev_per_e/1000.);
        else if (has_noise)
            sadc += sigma->GetBinContent(i,1);
        else
            sadc += 0;
        qtot += ibc;
        just1Dcluster->SetBinContent(i-bx_start+1, ibc+base);
        
        if(ibc>qmax){
            qmax = ibc;
            imax = i;
        }
        numbins++;
    }
    
    sadc /= numbins;
    
    //cout << sadc << "; " << qtot << "; " << base << endl;
    
    //To perform the fit
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex");
    
    // set tolerance , etc...
    min->SetMaxFunctionCalls(10000); // for Minuit/Minuit2
    min->SetTolerance(0.01);
    min->SetPrintLevel(0);
    
    // create funciton wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&LogLikelihood1D,4);
    min->SetErrorDef(0.5);
    min->SetFunction(f);
    
    // Set the free variables to be minimized!
    min->SetLimitedVariable(0,"x",meanx,0.01,meanx-2, meanx+2);
    min->SetLimitedVariable(1,"s",0.5, 0.01, 0.01, 3);
    
   	if (is_g4_sim) 
   		min->SetLimitedVariable(2,"q",qtot,1E-5,0.,1.2*ene_lim);
   	else 
   		min->SetLimitedVariable(2,"q",qtot,1E-3,0.,3.*qtot);
    
    if(zero_base) 
    	min->SetFixedVariable(3,"o",0);
    else 
    	min->SetLimitedVariable(3,"o",base,1E-3,base-sadc,base+sadc);
    
    // do the minimization
    Bool_t success = min->Minimize();
    
    const double *xs = min->X();
    const double *xe = min->Errors();
    Double_t ll = min->MinValue();
    
    qfit = xs[2];
    meanxout = xs[0]-0.5;
    LLout = TMath::IsNaN(ll) ? -1001 : ll;
    sigout = xs[1];
    sigoutunc = xe[1];
    base = xs[3];
    
    //Calculate energy variable inspired from eneh in pointFitScan
    Double_t q1org = 0;
    Double_t qfrac1org = 0;
    Double_t maxdorg = 2*smax; //may not be best for all smax
    if (maxdorg<0.51) maxdorg = 0.51;
    Double_t bot = TMath::Sqrt2()*sigout;
    
    for(Int_t i=1; i<=nx; i++){
        
        Double_t xlo = just1Dcluster->GetXaxis()->GetBinLowEdge(i);
        Double_t xhi = just1Dcluster->GetXaxis()->GetBinUpEdge(i);
        Double_t bcx = (xhi+xlo)/2.;
        Double_t frac = 0.5*(TMath::Erf((xhi-meanx)/bot)-TMath::Erf((xlo-meanx)/bot));
        Double_t bcADU = just1Dcluster->GetBinContent(i);
        
        if(qfit*frac > 4*sadc || TMath::Abs(meanxout-bcx) < maxdorg){

            q1org+=bcADU;
            qfrac1org+=frac;
        }
    }
    
    qh = q1org/qfrac1org;
    
    delete min;
}

void EvalLL1D(TH2F* cluster_large, TH2F* cluster, TH2F* sigma, TList* vals, Bool_t vid, Int_t LLwindow, Double_t LL_conv_factor, Double_t ene_lim, Double_t ene_min, Bool_t zero_base, Bool_t has_noise_cfg, Bool_t is_g4_sim_cfg, Double_t npixels, Double_t full_size, Double_t sim_noise){
    
    c = LL_conv_factor * 1000 / ev_per_e;
    
    has_noise = has_noise_cfg; // Set global (I know, I know...)
    
    Double_t small_charge_mean_x = cluster->GetMean(1); // Compute params. for small cluster
    Double_t small_charge_mean_y = cluster->GetMean(2);
    Double_t small_charge_total = cluster->Integral();
    
    if (!zero_base){
    	base = (cluster_large->Integral() - cluster->Integral())/(full_size-npixels);
    	//cout << base << "; " << cluster_large->Integral() << "; " << cluster->Integral() << "; " << small_charge_total*LL_conv_factor << endl;
    }
    else
    	base = 0;
    
    meanyout = small_charge_mean_y;
    q3 = -1000;
    q1 = -1000;
    qh = -1000;
    
    if ((small_charge_total*LL_conv_factor > ene_lim) || (small_charge_total*LL_conv_factor < ene_min) || cluster_large->GetNbinsY()>1){
        meanxout = small_charge_mean_x;
        LLout = -1000;
        sigout = -1000;
        sigoutunc = -1000;
        qfit = -1000;
        vid = false;
    }
    else{
        Fit1D(cluster_large, sigma, small_charge_mean_x, LLwindow, zero_base, is_g4_sim_cfg, ene_lim, sim_noise);
    }
    
    AddVariable(vals, "LL_meanx", meanxout);
    AddVariable(vals, "LL_meany", meanyout);
    AddVariable(vals, "LL_ll", LLout);
    AddVariable(vals, "LL_sigma", sigout);
    AddVariable(vals, "LL_sigmaunc", sigoutunc);
    AddVariable(vals, "LL_efit", qfit);
    AddVariable(vals, "LL_ene9", q3);
    AddVariable(vals, "LL_ene1", q1);
    AddVariable(vals, "LL_eneh", qh);
    
    if (vid){
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        
        TF1* fll = new TF1("fll","[0]*TMath::Gaus(x,[1],[2],1)+[3]",just1Dcluster->GetXaxis()->GetBinLowEdge(1),just1Dcluster->GetXaxis()->GetBinUpEdge(just1Dcluster->GetNbinsX()));
        fll->SetParameter(0,c*qfit);
        fll->SetParameter(1,meanxout);
        fll->SetParameter(2,sigout);
        fll->SetParameter(3,c*base);
        fll->SetLineStyle(2);
        fll->SetNpx(1000);
        just1Dcluster->Scale(c);
        TH1F* draw1Dfill = (TH1F*) just1Dcluster->Clone("D1D");
        draw1Dfill->SetLineWidth(0);
        draw1Dfill->SetFillStyle(3001);
        draw1Dfill->SetFillColor(2);
        
        double max = 0;
        
        for(int i=1; i<=just1Dcluster->GetNbinsX(); i++){
            
            if(just1Dcluster->GetBinContent(i)>max) max = just1Dcluster->GetBinContent(i);
            just1Dcluster->SetBinError(i,TMath::Sqrt(just1Dcluster->GetBinContent(i)+4)); //assume 2e- of noise
            draw1Dfill->SetBinContent(i,fll->Integral(draw1Dfill->GetXaxis()->GetBinLowEdge(i),draw1Dfill->GetXaxis()->GetBinUpEdge(i)));
        }
        
        max = max + 3*TMath::Sqrt(max);
        draw1Dfill->SetAxisRange(-10,max,"Y");
        draw1Dfill->Draw();
        just1Dcluster->Draw("same");
        fll->Draw("same");
    }

}

