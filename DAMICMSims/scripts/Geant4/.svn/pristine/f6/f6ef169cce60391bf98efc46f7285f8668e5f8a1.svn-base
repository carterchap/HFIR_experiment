#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include <vector>
#include <fstream>
#include <iostream>

struct Decay {

    Int_t idx = 0; //index where we are in tree
    Double_t tau = 0; //lifetime parameter for decay

    TTree* tree = NULL; //tree containing individual decays
    vector<Double_t> time; //time for individual decays
};

struct OTree{

    Int_t ndep = 0;
    Double_t time = 0;
    TTree* depo = NULL;
    TTree* info = NULL;
};

//exposure in seconds
//test.in list of files with details
void ReadGeant(Double_t exposure, Double_t pixel_size_x=0.0015, Double_t pixel_size_y=0.0015, TString input="test.in", TString opath="out/simage"){

    //General CCD params
    Double_t thickness;
    Int_t nbinsx, nbinsy;

    //Load input trees
    TString ifi, typ; Double_t tim;
    ifstream opt(input);

    vector< vector<Decay> > inputv; //vector containing input decays
    vector<Decay> vt;

    while(opt >> ifi >> typ >> tim){

        if(typ=="p" && vt.size()>0){

            inputv.push_back(vt);
            vt.clear();
        }

        TFile* irf = new TFile(ifi,"READ");

        if(vt.size()==0){
            //This is the first tree
            //Load CCD parameters assume the same for all files
            //They better be!

            Double_t Thick, sizeX, sizeY;
            TTree* ft = (TTree*) irf->Get("CCDInfo");
            ft->SetBranchAddress("Thick",&Thick);
            ft->SetBranchAddress("sizeX",&sizeX);
            ft->SetBranchAddress("sizeY",&sizeY);
            ft->GetEntry(0);

            thickness = Thick;
            nbinsx = sizeX/pixel_size_x-1e-3; nbinsx++;
            nbinsy = sizeY/pixel_size_y-1e-3; nbinsy++;
        }

        TTree* t = (TTree*) irf->Get("OutPut");
        Decay d; d.tree = t; d.idx = 0; d.tau = tim;
        vt.push_back(d);
    }

    if(vt.size()>0){

        inputv.push_back(vt);
        vt.clear();
    }
    opt.close();

    std::cout << "Read input files." << std::endl;

    //Now set the times for the decays
    TRandom3 r(0);

    for(size_t m=0; m<inputv.size(); m++)
        for(size_t n=0; n<inputv.at(m).size(); n++){

            //Set times for parent
            if(n==0){

                Double_t pt = 0; //previous time
                Int_t pp = 0; //previous primary

                Int_t cp; inputv.at(m).at(n).tree->SetBranchAddress("IDPrim",&cp);

                for(int k=0; k<inputv.at(m).at(n).tree->GetEntries(); k++){
                    inputv.at(m).at(n).tree->GetEntry(k);

                    Int_t pi = cp - pp; //primary step
                    Double_t ti = pi==0 ? 0 : inputv.at(m).at(n).tau * pi; //time step

                    //if(pi!=0) std::cout << ti + pt << std::endl;

                    inputv.at(m).at(n).time.push_back( ti + pt );
                    pp = cp;
                    pt = ti + pt;
                }
            }

            //Set times for daughters
            else{

            }
        }

    //Now input vector has the decay chains
    double timeout = false;
    Int_t rid = 1;

    //For this test do two runs only
    while(!timeout){

        //File name prefit for this run
        TString oname = opath;
        if(rid<1000) oname+="0";
        if(rid<100 ) oname+="0";
        if(rid<10  ) oname+="0";
        oname+=rid; oname+="_";

        //Create files for this run
        const int nccd = 12;

        //Go to top directory
        gROOT->cd();

        //Vector to contain output trees - one per CCD per exposure
        vector<OTree> vccd;
        Double_t dep_x, dep_y, dep_z, dep_e, dep_t;

        for(int k=0; k<nccd; k++){

            OTree ot;
            //Create sinfo tree
            Int_t simid = k+1;

            //Here are the binnig options
            //Assume the same for all CCDs
            ot.info = new TTree("sinfo","sinfo");
            ot.info->Branch("simid", &simid, "simid/I");
            ot.info->Branch("nbins_x", &nbinsx, "nbins_x/I");
            ot.info->Branch("nbins_y", &nbinsy, "nbins_y/I");
            ot.info->Branch("bw_x", &pixel_size_x, "bw_x/D");
            ot.info->Branch("bw_y", &pixel_size_y, "bw_y/D");
            ot.info->Branch("bw_z", &thickness, "bw_z/D");
            ot.info->Fill();

            //Now deposits tree
            ot.depo = new TTree("deposits","deposits");
            ot.depo->Branch("dep_x",&dep_x,"dep_x/D");
            ot.depo->Branch("dep_y",&dep_y,"dep_y/D");
            ot.depo->Branch("dep_z",&dep_z,"dep_z/D");
            ot.depo->Branch("dep_e",&dep_e,"dep_e/D");
            ot.depo->Branch("dep_t",&dep_t,"dep_t/D");
            vccd.push_back(ot);
        }

        std::cout << "Created vectors with output trees." << std::endl;

        //Time for this run to end
        Double_t tend = rid*exposure;

        //Now loop over the decays
        for(size_t m=0; m<inputv.size(); m++)
            for(size_t n=0; n<inputv.at(m).size(); n++){

                if(n==0){
                    //Deal with parent
                    Int_t CCDNum; inputv.at(m).at(n).tree->SetBranchAddress("CCDNum",&CCDNum);
                    Double_t EnergyDeposit; inputv.at(m).at(n).tree->SetBranchAddress("EnergyDeposit",&EnergyDeposit);
                    Double_t XCoord; inputv.at(m).at(n).tree->SetBranchAddress("XCoord",&XCoord);
                    Double_t YCoord; inputv.at(m).at(n).tree->SetBranchAddress("YCoord",&YCoord);
                    Double_t ZCoord; inputv.at(m).at(n).tree->SetBranchAddress("ZCoord",&ZCoord);

                    size_t k1 = inputv.at(m).at(n).idx;
                    Double_t time = inputv.at(m).at(n).time.at(k1);


                    while( time < tend ){

                        inputv.at(m).at(n).tree->GetEntry(k1);

                        //Now store depostis in output file
                        dep_e = EnergyDeposit;
                        dep_x = XCoord;
                        dep_y = YCoord;
                        dep_z = ZCoord;
                        dep_t = time;

                        if(CCDNum>0 && CCDNum<=nccd){

                            vccd.at(CCDNum-1).depo->Fill();
                            vccd.at(CCDNum-1).ndep++;
                        }

                        else{

                            std::cout << "Invalid CCD number! – skipping" << std::endl;
                        }

                        k1++;
                        if(k1>=inputv.at(m).at(n).time.size()){ timeout=true; break; }
                        time = inputv.at(m).at(n).time.at(k1);
                    }

                    inputv.at(m).at(n).idx = k1;
                }

                else{
                    //Here we will have to deal with daughters
                }
            }

        std::cout << "Writing output files." << std::endl;

        //Write files for this run
        for(int k=0; k<nccd; k++)
            if(vccd.at(k).ndep>0){

                TString fname = oname; fname += k+1; fname += ".pds";
                TFile* fout = new TFile(fname,"RECREATE");
                vccd.at(k).info->Write();
                vccd.at(k).depo->Write();
                fout->Close();
                gROOT->cd();
            }

        rid++;
    }

}
