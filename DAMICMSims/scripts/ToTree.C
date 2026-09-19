#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "ConfigValue.C"
#include "methods/Methods.C"

void ToTree(TString fname, TString cfgname=""){
    
    LoadConfigFile(cfgname);
    
    TFile* f = new TFile(fname, "UPDATE");
    
    if(f->IsZombie()){
        std::cout << "Could not open " << fname << std::endl;
        return;
    }
    
    TTree* t = NULL;
    t = (TTree*) f->Get("clusters_tree");
    if(t==NULL){
        std::cout << "Could not open clusters_tree" << std::endl;
        std::cout << "Run FindClusters first!" << std::endl;
        f->Close();
        return;
    }
    
    //to know if simulated file
    TParameter<double>* hwid = NULL;
    t->SetBranchAddress("HWID", &hwid);
    t->GetEntry(0);
    bool is_sim = false;
    if(hwid!=NULL && hwid->GetVal()<0) is_sim = true;
    
    //arrays to store object
    TArrayD* pixel_x = NULL;
    TArrayD* pixel_y = NULL;
    TArrayD* pixel_val = NULL;
    TArrayD* pixel_sigma = NULL;
    TArrayD* pixel_ped = NULL;
    
    TArrayD* pixel_sime = NULL;
    TArrayD* pixel_simz = NULL;
    //TArrayD* pixel_simx = NULL;
    //TArrayD* pixel_simy = NULL;
    TArrayD* pixel_simn = NULL;
    TArrayD* xpart = NULL;
    TArrayD* ypart = NULL;
    TArrayI* processCluster1 = NULL;
    TArrayI* processCluster2 = NULL;
    TArrayI* inside = NULL;

    // Added by KR
    // to know if it has extra pixels
    Bool_t has_extra_pix = false;
    Bool_t doLL = false;
    Bool_t LLcorr = false;
    Bool_t one_by_onehundred = false;
    Bool_t zero_base = false;
    Bool_t has_noise = false;
    Bool_t is_g4_sim = false;
    Int_t npixels = 0;
    Int_t LLwindow = 5;
    Double_t ene_lim = 0;
    Double_t ene_min = 0;
    Double_t sim_noise = 0;
    Double_t LL_conv_factor = 0.000625;
    if (t->GetListOfLeaves()->Contains("npixels")){
        has_extra_pix = true;
        t->SetBranchAddress("npixels", &npixels);
        doLL = ConfigValue("doLL");
        LL_conv_factor = ConfigValue("LL_conv_factor");
        LLwindow = ConfigValue("LLwindow");
        ene_lim = ConfigValue("ene_lim");
		ene_min = ConfigValue("ene_min");
        has_noise = ConfigValue("has_noise");
        one_by_onehundred = ConfigValue("one_by_onehundred");
        zero_base = ConfigValue("zero_base");
        is_g4_sim = ConfigValue("is_g4_sim");
    	sim_noise = ConfigValue("sim_noise");
    }
    
    t->SetBranchAddress("pixel_x", &pixel_x);
    t->SetBranchAddress("pixel_y", &pixel_y);
    t->SetBranchAddress("pixel_val", &pixel_val);
    t->SetBranchAddress("pixel_sigma", &pixel_sigma);
    t->SetBranchAddress("pixel_ped", &pixel_ped);
    
    //for simulated
    if(is_sim || is_g4_sim) {
        
        t->SetBranchAddress("pixel_sime", &pixel_sime);
        t->SetBranchAddress("pixel_simz", &pixel_simz);
        //t->SetBranchAddress("pixel_simx", &pixel_simx);
        //t->SetBranchAddress("pixel_simy", &pixel_simy);
        t->SetBranchAddress("pixel_simn", &pixel_simn);
        if (is_g4_sim){
        	t->SetBranchAddress("Xpart",&xpart);
        	t->SetBranchAddress("Ypart", &ypart);
        	t->SetBranchAddress("Process1", &processCluster1);
        	t->SetBranchAddress("Process2", &processCluster2);
        	t->SetBranchAddress("ProduceInside", &inside);
        }
    }
    
    //List to store values for cluster
    TList* vals = new TList();
    //TTree t owrite
    TTree* t2w = NULL;
    
    Int_t nentries = t->GetEntries();
    
    for(Int_t i=0; i<nentries; i++){
        
        t->GetEntry(i);
        TH2F* cluster = ArraysToTH2F(pixel_x, pixel_y, pixel_val, npixels, false); // Modified by KR
        cluster->SetNameTitle("cluster","cluster");
        TH2F* sigma = ArraysToTH2F(pixel_x, pixel_y, pixel_sigma, npixels, false); // Modified by KR
        sigma->SetNameTitle("sigma","sigma");
        TH2F* pedestal = ArraysToTH2F(pixel_x, pixel_y, pixel_ped, npixels, false); // Modified by KR
        pedestal->SetNameTitle("pedestal","pedestal");
        TH2F* cluster_large = ArraysToTH2F(pixel_x, pixel_y, pixel_val, npixels, doLL); // Added by KR
       	cluster_large->SetNameTitle("cluster_large","cluster_large");
        TH2F* sigma_large = ArraysToTH2F(pixel_x, pixel_y, pixel_sigma, npixels, doLL); // Added by KR
        sigma_large->SetNameTitle("sigma_large","sigma_large");

        Methods(cluster, sigma, pedestal, vals, 0, is_g4_sim);
        MethodsLL(cluster_large, cluster, sigma_large, vals, 0, LLwindow, LL_conv_factor, ene_lim, ene_min, one_by_onehundred, zero_base, has_noise, is_g4_sim, npixels, pixel_x->GetSize(), sim_noise);
        
        // for debug
        /*if (i == 515){
        
        	cout << npixels << endl;
        
            TCanvas *c1 = new TCanvas();
            c1->Draw();
            c1->cd();
            TH2F *hc = (TH2F*)cluster->Clone();
            hc->SetDirectory(0);
            hc->Draw("COLZ");
            
            TCanvas *c2 = new TCanvas();
            c2->Draw();
            c2->cd();
            TH2F *hcl = (TH2F*)cluster_large->Clone();
            hcl->SetDirectory(0);
            hcl->Draw("COLZ");
            
            //
          
        } */
        
        if(is_sim || is_g4_sim){
            TH2F* sime = ArraysToTH2F(pixel_x, pixel_y, pixel_sime, npixels, false); // Modified by KR
            sime->SetNameTitle("sime","sime");
            TH2F* simz = ArraysToTH2F(pixel_x, pixel_y, pixel_simz, npixels, false); // Modified by KR
            simz->SetNameTitle("simz","simz");
            TH2F* simn = ArraysToTH2F(pixel_x, pixel_y, pixel_simn, npixels, false); // Modified by KR
            simn->SetNameTitle("simn","simn");
            
            //TH2F* simx = ArraysToTH2F(pixel_x, pixel_y, pixel_simx);
            //simx->SetNameTitle("simx","simx");
            //TH2F* simy = ArraysToTH2F(pixel_x, pixel_y, pixel_simy);
            //simy->SetNameTitle("simy","simy");
            
            MethodsSim(sime, simz, simn, vals, 0);
            
            if (is_g4_sim){
              MethodsSimG4(vals, pixel_x, pixel_y, xpart, ypart, processCluster1, processCluster2, inside);
            }
            
            sime->Delete();
            simz->Delete();
            simn->Delete();
            //simx->Delete();
            //simy->Delete();
        }
        
        //If first time add Branches to tree
        if(i==0){
         
            TObjArray* barray0 = t->GetListOfBranches();
            //Do not clone branches from a previous run of ToTree.C
            for(Int_t j=0; j<vals->GetSize(); j++){
                
                if(barray0->Contains(vals->At(j)->GetName()))
                    ( (TBranch*) barray0->FindObject(vals->At(j)->GetName()) )->SetStatus(0);
            }
            
            t2w = t->CloneTree();
            //Add branches to new tree
            t2w->Branch(vals);
        }
        
        //Fill new branches
        for(Int_t j=0; j<vals->GetSize(); j++){
            
            TBranch* b = t2w->GetBranch(vals->At(j)->GetName());
            b->Fill();
        }
        
        cluster->Delete();
        cluster_large->Delete();
        sigma->Delete();
        sigma_large->Delete();
        pedestal->Delete();
    }
    
    t->Delete();
    t2w->SetNameTitle("clusters_tree","clusters_tree");
    t2w->Write("",TObject::kWriteDelete);
    f->Close();
    
    std::cout << "ToTree: Done!" << std::endl;
}
