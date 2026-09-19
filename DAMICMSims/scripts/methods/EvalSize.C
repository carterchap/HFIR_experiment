#include "TH2F.h"
#include "TList.h"
#include "TMath.h"
#include "TArray.h"
#include "TGraph.h"
#include "TVector3.h"
#include <iostream>

Bool_t contains(TArrayI* ax, TArrayI* ay, Int_t vx, Int_t vy){
    
    Int_t size = ax->GetSize();
    for(Int_t i=0; i<size; i++)
        if(ax->At(i)==vx && ay->At(i)==vy) return true;
    return false;
}

Double_t DoLength(TH2F* cluster, TArrayI* x, TArrayI* y, Int_t* flag, Double_t* end_points){
    
    //OPTIONS
    Int_t edge_dist = ConfigValue("edge_dist");
    Double_t covered_frac = ConfigValue("covered_frac");
    //======================================//
    
    Int_t nbinsx = cluster->GetNbinsX();
    Int_t nbinsy = cluster->GetNbinsY();
    
    //to handle one pixel clusters
    if(nbinsx==1 || nbinsy==1){
        
        *flag=0;
        return nbinsx*nbinsy;
    }
    
    Int_t ci = x->GetSize();
    Int_t cx = 0;
    Int_t cy = 0;
    Int_t maxb = 0;
    Int_t mayb = 0;
    
    //Find starting point
    Int_t nbins_edge = 999999;
    Int_t nn=0;
    Int_t tx=0;
    Int_t ty=0;
    Bool_t up = false;
    Double_t smax = 0;
    
    for(Int_t i=0; i<=nbinsx; i++){
        if(cluster->GetBinContent(i,1)>smax){
            maxb=i;
            mayb=1;
            smax=cluster->GetBinContent(i,1);
        }
        
        if(cluster->GetBinContent(i,2)>0) nn++;
        if(cluster->GetBinContent(i,1)>0){
            nn++;
            if(cluster->GetBinContent(i,2)>0 && cluster->GetBinContent(i-1,2)*cluster->GetBinContent(i,2)*cluster->GetBinContent(i+1,2)==0){
                tx=i;
                ty=1;
                up=true;
            }
        }
    }
        
    if(nn<nbins_edge&&nn>0&&up){
        nbins_edge=nn;
        cx=tx;
        cy=ty;
    }
    
    nn=0; up=false;
    for(Int_t i=0; i<=nbinsx; i++){
        if(cluster->GetBinContent(i,nbinsy)>smax){
            maxb=i;
            mayb=nbinsy;
            smax=cluster->GetBinContent(i,nbinsy);
        }
        
        if(cluster->GetBinContent(i,nbinsy-1)>0) nn++;
        if(cluster->GetBinContent(i,nbinsy)>0){
            nn++;
            if(cluster->GetBinContent(i,nbinsy-1)>0 && cluster->GetBinContent(i-1,nbinsy-1)*cluster->GetBinContent(i,nbinsy-1)*cluster->GetBinContent(i+1,nbinsy-1)==0){
                tx=i;
                ty=nbinsy;
                up=true;
            }
        }
    }
    
    if(nn<nbins_edge&&nn>0&&up){
        nbins_edge=nn;
        cx=tx;
        cy=ty;
    }
    
    nn=0; up=false;
    for(Int_t j=0; j<=nbinsy; j++){
        if(cluster->GetBinContent(1,j)>smax){
            maxb=1;
            mayb=j;
            smax=cluster->GetBinContent(1,j);
        }
        
        if(cluster->GetBinContent(2,j)>0) nn++;
        if(cluster->GetBinContent(1,j)>0){
            nn++;
            if(cluster->GetBinContent(2,j)>0 && cluster->GetBinContent(2,j-1)*cluster->GetBinContent(2,j)*cluster->GetBinContent(2,j+1)==0){
                tx=1;
                ty=j;
                up=true;
            }
        }
    }
        
    if(nn<nbins_edge&&nn>0&&up){
        nbins_edge=nn;
        cx=tx;
        cy=ty;
    }
    
    nn=0; up=false;
    for(Int_t j=0; j<=nbinsy; j++){
        if(cluster->GetBinContent(nbinsx,j)>smax){
            maxb=nbinsx;
            mayb=j;
            smax=cluster->GetBinContent(nbinsx,j);
        }
        
        if(cluster->GetBinContent(nbinsx-1,j)>0) nn++;
        if(cluster->GetBinContent(nbinsx,j)>0){
            nn++;
            if(cluster->GetBinContent(nbinsx-1,j)>0 && cluster->GetBinContent(nbinsx-1,j-1)*cluster->GetBinContent(nbinsx-1,j)*cluster->GetBinContent(nbinsx-1,j+1)==0){
                tx=nbinsx;
                ty=j;
                up=true;
            }
        }
    }
        
    if(nn<nbins_edge&&nn>0&&up){
        nbins_edge=nn;
        cx=tx;
        cy=ty;
    }
    
    if(cx==0 && cy==0){

        cx=maxb;
        cy=mayb;
    }
    
    //at this point cx and cy have the starting point, which is either the last pixel along the boundary of the canvas that has the smallest number of "edge" pixels or (if no suitabl edge pixels found) the maximum point along the boundary.
	end_points[0]=cx;
	end_points[1]=cy;
	//cout<<"cx = "<<cx<<", cy = "<<cy<<endl;
	
    
    Int_t nextx = 0; Int_t nexty = 0;
    Int_t minx = 999999; Int_t miny = 999999; Int_t maxx = 0; Int_t maxy = 0;
    
    //Add starting point to array
    if(ci==0){
        
        ci++;
        x->Set(ci);
        y->Set(ci);
        x->AddAt(cx,ci-1);
        y->AddAt(cy,ci-1);
        if(cx>maxx) maxx=cx;
        if(cy>maxy) maxy=cy;
        if(cx<minx) minx=cx;
        if(cy<miny) miny=cy;
    }
    
    Double_t length = 0;
    TVector3 cdir(0,0,0);
    TVector3 pdir(0,0,0);
    TVector3 sdir(0,0,0);
    TVector3 ipoint(cx,cy,0);
    
    Int_t vec_check = 4;
    
    //make initial direction toward center
    cdir.SetX(nbinsx/2-cx);
    cdir.SetY(nbinsy/2-cy);
    pdir=cdir;
    sdir=pdir;
    
    while(true){
        
        Double_t max = 0;
        for(Int_t i=cx-1; i<=cx+1; i++)
            for(Int_t j=cy-1; j<=cy+1; j++){
                
                TVector3 v(i-cx,j-cy,0);
                TVector3 iv(i-ipoint.X(),j-ipoint.Y(),0);
                
                if(i<=0 || j<=0 || i>nbinsx || j>nbinsy ||
                   contains(x,y,i,j)) continue;
                
                if(
                   ((iv.Dot(v)<0&&vec_check>0)   && !(v.Dot(cdir)>0 && v.Dot(pdir)>0 && v.Dot(sdir)>0)) ||
                   (cdir.Dot(v)<0&&vec_check>1) ||
                   (pdir.Dot(v)<0&&vec_check>2) ||
                   (sdir.Dot(v)<0&&vec_check>3)
                   )  //this is to prevent a U turn
                    continue;
                
                Double_t bc = cluster->GetBinContent(i,j);
                Int_t bn = 1;
                if(bc>0){
                    Int_t ii=i+v.X(); Int_t jj=j+v.Y();
                    while(ii>0 && jj>0 && ii<=nbinsx && jj<=nbinsy){
                        Double_t bcc = cluster->GetBinContent(ii,jj);
                        if(bcc==0 || contains(x,y,ii,jj)) break;
                        bc+=bcc;
                        bn++;
                        ii+=v.X(); jj+=v.Y();
                    }
                }
                
                if(bc>max){
                    max = bc;
                    nextx = i;
                    nexty = j;
                }
                
            }
        
        if(max==0){
            
            if(((maxx>nbinsx-edge_dist && maxy>nbinsy-edge_dist && minx<=edge_dist && miny<=edge_dist)
               && (cluster->Integral(minx,maxx,miny,maxy)/cluster->Integral(1,nbinsx,1,nbinsy) > covered_frac)
               ) || vec_check==0 )
                break;
            else{ //somehow stuck in middle of TH2, try to get back on track by losening vector constraints
                vec_check--;
                continue;
            }
        }
        
        vec_check=4;
        sdir=pdir;
        pdir=cdir;
        cdir.SetX(nextx-cx);
        cdir.SetY(nexty-cy);
        length += cdir.Mag();
        
        cx=nextx;
        cy=nexty;
            
        ci++;
        x->Set(ci);
        y->Set(ci);
        x->AddAt(cx,ci-1);
        y->AddAt(cy,ci-1);
        if(cx>maxx) maxx=cx;
        if(cy>maxy) maxy=cy;
        if(cx<minx) minx=cx;
        if(cy<miny) miny=cy;
    }

	end_points[2]=cx;
	end_points[3]=cy;
    
    //Let's do the flagging
    if((minx<=edge_dist && miny<=edge_dist && maxx>nbinsx-edge_dist && maxy>nbinsy-edge_dist) && (cluster->Integral(minx,maxx,miny,maxy)/cluster->Integral(1,nbinsx,1,nbinsy) > covered_frac))
        *flag=0;
    else
        *flag=1;
    
    return length;
    
}

void EvalSize(TH2F* cluster, TList* vals, Bool_t vid){

    Int_t npixels = 0;
    Int_t nbinsx = cluster->GetNbinsX();
    Int_t nbinsy = cluster->GetNbinsY();
	//cout<<"start x =" <<cluster->GetBinCenter(20)<<endl;

    
    for(Int_t i=1; i<=nbinsx; i++)
        for(Int_t j=1; j<=nbinsy; j++)
            if(cluster->GetBinContent(i,j)>0) npixels++;
    
    AddVariable(vals, "size_npixels",npixels);
    AddVariable(vals, "size_slength",TMath::Sqrt(nbinsx*nbinsx + nbinsy*nbinsy));
    AddVariable(vals, "size_aspect_ratio", (Double_t) nbinsy/nbinsx);
    
    TArrayI* x = new TArrayI();
    TArrayI* y = new TArrayI();

    Int_t flag;
    Double_t end_points[4] = {0,0,0,0};
    Double_t clength = DoLength(cluster,x,y,&flag,end_points);
    AddVariable(vals, "size_clength", clength);
    AddVariable(vals, "size_cflag", flag);
    if(cluster->GetEntries()>2){
	end_points[0]+=cluster->GetXaxis()->GetBinCenter(1)-0.5;
	end_points[1]+=cluster->GetYaxis()->GetBinCenter(1)-0.5;
	end_points[2]+=cluster->GetXaxis()->GetBinCenter(1)-0.5;
	end_points[3]+=cluster->GetYaxis()->GetBinCenter(1)-0.5;}
	else {end_points[0]=cluster->GetXaxis()->GetBinCenter(1)-0.5;
		 end_points[1]=cluster->GetYaxis()->GetBinCenter(1)-0.5;
		end_points[2]=cluster->GetXaxis()->GetBinCenter(1)+0.5;
		end_points[3]=cluster->GetYaxis()->GetBinCenter(1)+0.5;}
	
	AddVariable(vals, "track_startx",end_points[0]);
	AddVariable(vals, "track_starty",end_points[1]);
	AddVariable(vals, "track_endx",end_points[2]);
	AddVariable(vals, "track_endy",end_points[3]);

    //Make histogram of along track
    TH1F* ddx = NULL;
    
    if(x->GetSize()>0){
        
        ddx = new TH1F("ddx","ddx",x->GetSize()-1,0,x->GetSize()-1);
        TVector3 prev_pos(0,0,0);
        Double_t prev_cont = 0;
        for(Int_t k=0; k<x->GetSize(); k++){
            
            Int_t idx = x->At(k);
            Int_t jdx = y->At(k);
            
            Double_t cont = cluster->GetBinContent(idx,jdx);
            TVector3 dir(idx-prev_pos.X(),jdx-prev_pos.Y(),0);
            
            if(k>0) ddx->SetBinContent(k,(cont+prev_cont)/(2.*dir.Mag()));
            
            prev_cont=cont;
            prev_pos.SetX(idx);
            prev_pos.SetY(jdx);
        }
        
        //Add variables associated with ddx
        AddVariable(vals, "size_ddx_mean", ddx->Integral()/ddx->GetNbinsX());
        AddVariable(vals, "size_ddx_max", ddx->GetMaximum());
        AddVariable(vals, "size_ddx_max_x", cluster->GetXaxis()->GetBinCenter(x->At(ddx->GetMaximumBin()))+0.5);
        AddVariable(vals, "size_ddx_max_y", cluster->GetYaxis()->GetBinCenter(y->At(ddx->GetMaximumBin()))+0.5);
        AddVariable(vals, "size_ddx_min", ddx->GetMinimum());
        AddVariable(vals, "size_ddx_min_x", cluster->GetXaxis()->GetBinCenter(x->At(ddx->GetMinimumBin()))+0.5);
        AddVariable(vals, "size_ddx_min_y", cluster->GetYaxis()->GetBinCenter(y->At(ddx->GetMinimumBin()))+0.5);
    }

    else{
        
        AddVariable(vals, "size_ddx_mean", -1);
        AddVariable(vals, "size_ddx_max", -1);
        AddVariable(vals, "size_ddx_max_x", -1);
        AddVariable(vals, "size_ddx_max_y", -1);
        AddVariable(vals, "size_ddx_min", -1);
        AddVariable(vals, "size_ddx_min_x", -1);
        AddVariable(vals, "size_ddx_min_y", -1);
    }
    
    //==
    //For visual represenation with viewer
    if(vid && x->GetSize()>0){
        for(Int_t i=0; i<x->GetSize(); i++){
            
            x->AddAt(x->At(i)+cluster->GetXaxis()->GetXmin(),i);
            y->AddAt(y->At(i)+cluster->GetYaxis()->GetXmin(),i);
        }
        
        TGraph* g = new TGraph(x->GetSize(),x->GetArray(),y->GetArray());
        g->Print();
        g->SetLineWidth(2);
        g->Draw("L*");
        
        if(!gROOT->GetListOfCanvases()->FindObject("cview3"))
            new TCanvas("cview3","cview3",700,500);
        else
            ( (TCanvas*) gROOT->GetListOfCanvases()->FindObject("cview3"))->cd();
        
        ddx->Draw();
    }
    //==
    
    if(!vid && x->GetSize()>0) ddx->Delete();
    delete x;
    delete y;
}
