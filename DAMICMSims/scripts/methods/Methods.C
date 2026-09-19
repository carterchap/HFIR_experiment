#include "TH2F.h"
#include "TList.h"

#include "../methods/common.C"
#include "../methods/EvalCharge.C"
#include "../methods/EvalSize.C"
#include "../methods/EvalCurve.C"
#include "../methods/EvalSigma.C"
#include "../methods/EvalPedestal.C"
#include "../methods/EvalRaw.C"
#include "../methods/EvalSim.C"
#include "../methods/EvalLength.C"
#include "../methods/EvalLL.C"
#include "../methods/EvalLL1D.C"
#include "../methods/EvalSimG4.C"

void Methods(TH2F* cluster, TH2F* sigma, TH2F* pedestal, TList* vals, Int_t vid, Bool_t is_g4_sim = false){

    //IN ORDER
    EvalCharge(cluster, vals, vid==1, is_g4_sim);
    EvalSize(cluster, vals, vid==2);
    EvalCurve(cluster,vals, vid==3);
    EvalSigma(sigma, vals, vid==4);
    EvalPedestal(pedestal, vals, vid==5);
    EvalRaw(cluster, pedestal, vals, vid==6);
    EvalLength(cluster,vals,vid==7);
}

void MethodsSim(TH2F* sime, TH2F* simz, TH2F* simn, TList* vals, Int_t vid){

    EvalSim(sime, "sime", vals, vid==10);
    EvalSim(simz, "simz", vals, vid==13);
    EvalSim(simn, "simn", vals, vid==11);
}

// Added by KR
void MethodsLL(TH2F* cluster_large, TH2F* cluster, TH2F* sigma, TList* vals, Int_t vid, Int_t LLwindow, Double_t LL_conv_factor, Double_t ene_lim, Double_t ene_min, Bool_t one_by_onehundred, Bool_t zero_base, Bool_t has_noise, Bool_t is_g4_sim, Double_t npixels, Double_t full_size, Double_t sim_noise){
       
    if(one_by_onehundred)
        EvalLL1D(cluster_large, cluster, sigma, vals, vid==8, LLwindow, LL_conv_factor, ene_lim, ene_min, zero_base, has_noise, is_g4_sim, npixels, full_size, sim_noise);
    
    else
        EvalLL(cluster_large, cluster, sigma, vals, vid==8, LLwindow, LL_conv_factor, ene_lim, ene_min, zero_base, has_noise, is_g4_sim, npixels, full_size, sim_noise);
    
}

// Added by Joao
void MethodsSimG4(TList* vals, TArrayD* Xas, TArrayD* Yas, TArrayD* Xpart, TArrayD* Ypart, TArrayI* p1, TArrayI* p2, TArrayI* in ){
  //clock_t begin = clock();
  vector<Int_t>* Number = EvalSimG4Num(Xas,Yas,Xpart,Ypart);
  //clock_t begin2 = clock();
  //cout << begin2-begin << endl;
  EvalSimG4Prop(vals, p1, p2, in, Number);
}
