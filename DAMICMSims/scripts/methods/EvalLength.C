#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include <vector>
#include "TF1.h"

TH2F *mark(TH2F * hist_easy,Int_t binx,Int_t biny,Int_t npixels)
{

	vector<int>* kx = new vector<int>[npixels+1];
	vector<int>* ky = new vector<int>[npixels+1];
	Int_t i=0;
	kx[0].push_back(binx);
	ky[0].push_back(biny);
	hist_easy->SetBinContent(binx,biny,i);
	do	
    {
	for(size_t j =0;j<kx[i].size();j++)
	{
		binx=kx[i][j];
		biny=ky[i][j];
 		if(hist_easy->GetBinContent(binx+1,biny)==-1)  {kx[i+1].push_back(binx+1); ky[i+1].push_back(biny);  hist_easy->SetBinContent(binx+1,biny,i+1);     }
 		if(hist_easy->GetBinContent(binx-1,biny)==-1)  {kx[i+1].push_back(binx-1); ky[i+1].push_back(biny);  hist_easy->SetBinContent(binx-1,biny,i+1);     }
 		if(hist_easy->GetBinContent(binx,biny-1)==-1)  {kx[i+1].push_back(binx);   ky[i+1].push_back(biny-1);hist_easy->SetBinContent(binx,biny-1,i+1);     }
 		if(hist_easy->GetBinContent(binx,biny+1)==-1)  {kx[i+1].push_back(binx);   ky[i+1].push_back(biny+1);hist_easy->SetBinContent(binx,biny+1,i+1);     }
		
    }
	i++;
    }while(kx[i].size()>0);
    
    delete [] kx;
    delete [] ky;

	return hist_easy;
	
}
	

void EvalLength(TH2F* cluster, TList* vals, Bool_t vid)
{

	Int_t nbinsx = cluster->GetNbinsX();
    Int_t nbinsy = cluster->GetNbinsY();
	
	//Create two 2D hist with binvalue =-1 whever the track is 
	TH2F *hist_easy = NULL;
	TH2F *hist_easy1 = NULL;
	hist_easy=(TH2F *)cluster->Clone();
	Int_t npixels=0;
	for(Int_t i=1;i<=nbinsx;i++)
	{	
		for(Int_t j=1;j<=nbinsy;j++)
		{
			if(cluster->GetBinContent(i,j)>0) {hist_easy->SetBinContent(i,j,-1);npixels++;}
			}
	}
	
	hist_easy1=(TH2F *)hist_easy->Clone();
	
	
	
	//Begin finding the end of the track
	Int_t track_end_x,track_end_y,track_start_x,track_start_y;
	Int_t binx,biny,binz;
	Int_t track_length=0;
	cluster->GetBinXYZ(cluster->GetMaximumBin(),binx,biny,binz);
	if(npixels==1){track_end_x=binx+cluster->GetXaxis()->GetBinCenter(1)-0.5;
		track_start_x=binx+cluster->GetXaxis()->GetBinCenter(1)+0.5;
		track_end_y=biny+cluster->GetYaxis()->GetBinCenter(1)-0.5;
		track_start_y=biny+cluster->GetYaxis()->GetBinCenter(1)+0.5;
		track_length=1;}
	
	else{
	
	
	hist_easy=mark(hist_easy,binx,biny,npixels);
	hist_easy->GetBinXYZ(hist_easy->GetMaximumBin(),binx,biny,binz);
	track_start_x=binx+cluster->GetXaxis()->GetBinCenter(1)-0.5;
	track_start_y=biny+cluster->GetYaxis()->GetBinCenter(1)-0.5;
	//Finding the start of the track;
	hist_easy1=mark(hist_easy1,binx,biny,npixels);
	
	hist_easy1->GetBinXYZ(hist_easy1->GetMaximumBin(),binx,biny,binz);
	track_end_x=binx+cluster->GetXaxis()->GetBinCenter(1)-0.5;
	track_end_y=biny+cluster->GetYaxis()->GetBinCenter(1)-0.5;
	track_length=hist_easy1->GetMaximum()+1;}
	
    AddVariable(vals, "track_startx",track_start_x);
    AddVariable(vals, "track_starty",track_start_y);
    AddVariable(vals, "track_endx",track_end_x);
	AddVariable(vals, "track_endy",track_end_y);
	AddVariable(vals, "track_length",track_length);
	
	if(vid)
        {
            if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
                new TCanvas("cview3","cview3",700,500);
            else
                ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
                
            hist_easy1->Draw("COLZ");
          TMarker *tm = new TMarker(track_start_x,track_start_y,23);
            tm->SetMarkerSize(2);
            tm->Draw("same");
            TMarker *tm1 = new TMarker(track_end_x,track_end_y,29);
          tm1->SetMarkerSize(2);
            tm1->Draw("same");
	}
    
    else{

        hist_easy->Delete();
        hist_easy1->Delete();
    }
}
