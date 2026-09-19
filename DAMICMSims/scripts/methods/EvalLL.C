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

Double_t ev_per_e = 3.77;
Double_t smax = 1.3; // Taken from pointFitScan

TH2F* justLLcluster;
Int_t nx, ny;
Double_t bxlow, bylow;
Double_t c; // conv factor to e-
Double_t base, sadc;
Double_t q1, q3, qh, qfit, sigout, sigoutunc, meanxout, meanyout, LLout;
Bool_t has_noise = false;

Double_t maxk(Double_t obs, Double_t ne){
	//find k that offers maximum contribution
	Double_t kp = c*c*sadc*sadc*(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc);
	Int_t iters = 0;
	Int_t iiters = 0;
	while( TMath::Abs(kp/c/c/sadc/sadc + TMath::Log(kp) - TMath::Log(ne/2.) - TMath::EulerGamma() - obs/c/sadc/sadc) > 0.1 ){
        
		kp /= (kp/c/c/sadc/sadc + TMath::Log(kp))/(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc);
		iters++;
		if(iters>6){
			kp = c*c*sadc*sadc*(TMath::Log(ne/2.) + TMath::EulerGamma() + obs/c/sadc/sadc) - c*c*sadc*sadc*TMath::Log(kp);
			iters=0;
			iiters++;
		}
		//worst case
		if(iiters>2){
			kp=0;
			break;
		}
	}
	return kp;
}

//to fit for sigma
double LogLikelihood(const double *xx){
    
	const Double_t x = xx[0];
	const Double_t y = xx[1];
	const Double_t s = xx[2];
	const Double_t q = xx[3];
	const Double_t o = xx[4];
    
	Double_t bot = TMath::Sqrt2()*s;
    
	Double_t ll = 0;
    
	for(Int_t i=1; i<=nx; i++){
		for(Int_t j=1; j<=ny; j++){
            
			Double_t obs = justLLcluster->GetBinContent(i,j) - o;
            
			if(obs!=0){
                
				//expected number of electrons
				Double_t ne = q*0.25*c*(TMath::Erf((i+bxlow+0.5-x)/bot) - TMath::Erf((i+bxlow-0.5-x)/bot))*(TMath::Erf((j+bylow+0.5-y)/bot) - TMath::Erf((j+bylow-0.5-y)/bot));
                
                
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
	}
    
	return -1.*ll;
    
}

void Fit(TH2F* cluster, TH2F* sigma, Double_t meanx, Double_t meany, Int_t LLwindow, Bool_t zero_base, Bool_t is_g4_sim, Double_t ene_lim, Double_t sim_noise){
    
	nx = 2*LLwindow + 1;
	ny = nx;
    
	Int_t center_bin_x = cluster->GetXaxis()->FindBin(meanx);
	Int_t center_bin_y = cluster->GetYaxis()->FindBin(meany);
    
	Int_t bin_x_max = TMath::Min(center_bin_x+LLwindow, cluster->GetNbinsX());
	Int_t bin_y_max = TMath::Min(center_bin_y+LLwindow, cluster->GetNbinsY());
    
	//cout << center_bin_x << " " << center_bin_y << endl;
    
	Int_t bx_start = center_bin_x - LLwindow;
	bxlow = cluster->GetXaxis()->GetBinLowEdge(bx_start);
	Int_t by_start = center_bin_y - LLwindow;
	bylow = cluster->GetYaxis()->GetBinLowEdge(by_start);
    
	//cout << bxlow << " " << bylow << endl;
    
	justLLcluster = new TH2F("","",nx,0,nx,ny,0,ny);
    
	Double_t qtot = 0;
	Double_t ibc = 0;
    
	Double_t qmax = 0;
	Int_t imax = 0;
	Int_t jmax = 0;
    
	Int_t numbins = 0;
	for (Int_t i=bx_start; i <= bin_x_max; i++){
		for (Int_t j=by_start; j <= bin_y_max; j++){
			ibc = cluster->GetBinContent(i,j) - base;
			if (has_noise && is_g4_sim)
				sadc += (sim_noise*ev_per_e/1000.); 
			else if (has_noise)
				sadc += sigma->GetBinContent(i,j);
			else
				sadc += 0;
			qtot += ibc;
			justLLcluster->SetBinContent(i-bx_start+1, j-by_start+1,ibc+base);
            
			if(ibc>qmax){
				qmax = ibc;
				imax = i;
				jmax = j;
			}
			numbins++;
		}
	}
    
	sadc /= numbins;
    
	//cout << sadc << "; " << qtot << endl;
    
	//To perform the fit
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex");
    
	// set tolerance , etc...
	min->SetMaxFunctionCalls(10000); // for Minuit/Minuit2
	min->SetTolerance(0.01);
	min->SetPrintLevel(0);
    
	// create funciton wrapper for minmizer
	// a IMultiGenFunction type
	ROOT::Math::Functor f(&LogLikelihood,5);
	min->SetErrorDef(0.5);
	min->SetFunction(f);
    
	// Set the free variables to be minimized!
	min->SetLimitedVariable(0,"x",meanx,0.01,meanx-1, meanx+1);
	min->SetLimitedVariable(1,"y",meany,0.01,meany-1, meany+1);
	min->SetLimitedVariable(2,"s",0.5, 1E-2, 1E-2, 5.);
   
	if (is_g4_sim) 
		min->SetLimitedVariable(3,"q",qtot, 1E-5, 0., 1.2*ene_lim);
	else 
		min->SetLimitedVariable(3,"q",qtot,1E-3,0.,3.*qtot);
    
	if(zero_base) 
		min->SetFixedVariable(4,"o",0);
	else 
		min->SetLimitedVariable(4,"o",base,1E-3,base-sadc,base+sadc);
    
	// do the minimization
	Bool_t success = min->Minimize();
    
	const double *xs = min->X();
	const double *xe = min->Errors();
	Double_t ll = min->MinValue();
    
	qfit = xs[3];
	meanxout = xs[0]-0.5;
	meanyout = xs[1]-0.5;
	LLout = TMath::IsNaN(ll) ? -1001 : ll;
	sigout = xs[2];
	sigoutunc = xe[2];
    
	Int_t dummy;
	justLLcluster->GetBinXYZ(cluster->FindBin(meanxout,meanyout),imax,jmax,dummy);
    
	//Calculate q3, i.e. the charge from the 3x3 pixel square centered about best-fit meanx pixel
	q3 = 0;
	Double_t qfrac3 = 0;
	q1 = 0;
	Double_t qfrac1 = 0;
	Double_t maxd = 2*smax; //may not be best for all smax
	if(maxd<0.71) maxd = 0.71;
	//variable to contain corrected charge from highest signal-to-background pixels
	qh = 0;
	Double_t qfrach = 0;
    
	Double_t bot = TMath::Sqrt2()*sigout;
    
	for(Int_t i=1; i<=nx; i++)
		for(Int_t j=1; j<=ny; j++){
            
			Double_t xlo = justLLcluster->GetXaxis()->GetBinLowEdge(i)+bxlow+0.5;
			Double_t xhi = justLLcluster->GetXaxis()->GetBinUpEdge(i)+bxlow+0.5;
			Double_t ylo = justLLcluster->GetYaxis()->GetBinLowEdge(j)+bylow+0.5;
			Double_t yhi = justLLcluster->GetYaxis()->GetBinUpEdge(j)+bylow+0.5;
            
			Double_t bcx = (xhi+xlo)/2.;
			Double_t bcy = (yhi+ylo)/2.;
            
			Double_t frac = 0.25*(TMath::Erf((xhi-meanxout)/bot)-TMath::Erf((xlo-meanxout)/bot))*(TMath::Erf((yhi-meanyout)/bot)-TMath::Erf((ylo-meanyout)/bot));
			Double_t bc = justLLcluster->GetBinContent(i,j) - xs[4];
            
			if((TMath::Abs((imax-bx_start+1)-i)<=1 || TMath::Abs((jmax-by_start+1)-j)<=1) && !TMath::IsNaN(bc)){
                
				q3+=bc;
				qfrac3+=frac;
			}
            
			if(TMath::Sqrt((meanxout-bcx)*(meanxout-bcx) + (meanyout-bcy)*(meanyout-bcy)) < maxd){
                
				q1+=bc;
				qfrac1+=frac;
			}
            
			if((xs[3]*frac > 4*sadc) || (TMath::Sqrt((meanxout-bcx)*(meanxout-bcx) + (meanyout-bcy)*(meanyout-bcy)) < maxd)){
                
				qh+=bc;
				qfrach+=frac;
			}
		}
    
	//correct for only the charge expected in 3x3 square to obtain total
	q3/=qfrac3;
	q1/=qfrac1;
	qh/=qfrach;
    
	//code was written to return mean x and y according to bin number being bin center.
	//Here shift back to use bin up edge as bin number, as general convention in our code.
	//cout << success << "; LL: " << ll << "; X:" << meanx << "; Y:" << meany << "; S:" << xs[2] << "+-" << xe[2] << "; q:" << xs[3]*c*ev_per_e/1000. << endl;
	//cout << "q3: " << q3*c*ev_per_e/1000. << ", q1: " << q1*c*ev_per_e/1000. << ", qh: " << qh*c*ev_per_e/1000. << endl;
    
	delete min;
    
}

void EvalLL(TH2F* cluster_large, TH2F* cluster, TH2F* sigma, TList* vals, Bool_t vid, Int_t LLwindow, Double_t LL_conv_factor, Double_t ene_lim, Double_t ene_min, Bool_t zero_base, Bool_t has_noise_cfg, Bool_t is_g4_sim_cfg, Double_t npixels, Double_t full_size, Double_t sim_noise){
    
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
    	
    
	if ((small_charge_total*LL_conv_factor > ene_lim) || (small_charge_total*LL_conv_factor < ene_min)){
		meanxout = small_charge_mean_x;
		meanyout = small_charge_mean_y;
		LLout = -1000;
		sigout = -1000;
		sigoutunc = -1000;
		qfit = -1000;
		q3 = -1000;
		q1 = -1000;
		qh = -1000;
	}
	else{
		Fit(cluster_large, sigma, small_charge_mean_x, small_charge_mean_y, LLwindow, zero_base, is_g4_sim_cfg, ene_lim, sim_noise);
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
		TCanvas *cview3=NULL;
		if(!gROOT->GetListOfCanvases()->FindObject("cview3"))  {
			cview3 = new TCanvas("cview3","cview3",1400,500);
			cview3->Divide(2,1);
		}				
		else
			cview3= (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3");
        
		TH1F* just1DclusterX = (TH1F*) cluster_large->ProjectionX();
		TH1F* just1DclusterY = (TH1F*) cluster_large->ProjectionY();
		TF1* fllX = new TF1("fllX","[0]*TMath::Gaus(x,[1],[2],1)+[3]",just1DclusterX->GetXaxis()->GetBinLowEdge(1),just1DclusterX->GetXaxis()->GetBinUpEdge(just1DclusterX->GetNbinsX()));
		TF1* fllY = new TF1("fllY","[0]*TMath::Gaus(x,[1],[2],1)+[3]",just1DclusterY->GetXaxis()->GetBinLowEdge(1),just1DclusterY->GetXaxis()->GetBinUpEdge(just1DclusterY->GetNbinsX()));

		fllX->SetParameter(0,c*qfit);
		fllX->SetParameter(1,meanxout);
		fllX->SetParameter(2,sigout);
		fllX->SetParameter(3,c*base);
		fllX->SetLineStyle(2);
		fllX->SetNpx(1000);
		just1DclusterX->Scale(c);

		fllY->SetParameter(0,c*qfit);
		fllY->SetParameter(1,meanyout);
		fllY->SetParameter(2,sigout);
		fllY->SetParameter(3,c*base);
		fllY->SetLineStyle(2);
		fllY->SetNpx(1000);
		just1DclusterY->Scale(c);

		TH1F* draw1DfillX = (TH1F*) just1DclusterX->Clone("D1DX");
		draw1DfillX->SetLineWidth(0);
		draw1DfillX->SetFillStyle(3001);
		draw1DfillX->SetFillColor(2);
		draw1DfillX->SetTitle("X-axis Projection");

		TH1F* draw1DfillY = (TH1F*) just1DclusterY->Clone("D1DY");
		draw1DfillY->SetLineWidth(0);
		draw1DfillY->SetFillStyle(3001);
		draw1DfillY->SetFillColor(2);
		draw1DfillY->SetTitle("Y-axis Projection");
        

		double maxX = 0;
		double maxY = 0;
		for(int i=1; i<=just1DclusterX->GetNbinsX(); i++){            
			if(just1DclusterX->GetBinContent(i)>maxX) maxX = just1DclusterX->GetBinContent(i);
			just1DclusterX->SetBinError(i,TMath::Sqrt(just1DclusterX->GetBinContent(i)+4)); //assume 2e- of noise
			draw1DfillX->SetBinContent(i,fllX->Integral(draw1DfillX->GetXaxis()->GetBinLowEdge(i),draw1DfillX->GetXaxis()->GetBinUpEdge(i)));
		}
		for(int i=1; i<=just1DclusterY->GetNbinsX(); i++){            
			if(just1DclusterY->GetBinContent(i)>maxY) maxY = just1DclusterY->GetBinContent(i);
			just1DclusterY->SetBinError(i,TMath::Sqrt(just1DclusterY->GetBinContent(i)+4)); //assume 2e- of noise
			draw1DfillY->SetBinContent(i,fllY->Integral(draw1DfillY->GetXaxis()->GetBinLowEdge(i),draw1DfillY->GetXaxis()->GetBinUpEdge(i)));
		}
		
		cview3->cd(1);
		maxX = maxX + 3*TMath::Sqrt(maxX);
		draw1DfillX->SetAxisRange(-10,maxX,"Y");
		draw1DfillX->Draw();
		just1DclusterX->Draw("same");
		fllX->Draw("same");

		cview3->cd(2);
		maxY = maxY + 3*TMath::Sqrt(maxY);
		draw1DfillY->SetAxisRange(-10,maxY,"Y");
		draw1DfillY->Draw();
		just1DclusterY->Draw("same");
		fllY->Draw("same");

	}

}

