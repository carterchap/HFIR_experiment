#include <TROOT.h>
#include <TObject.h>
#include <TStyle.h>
#include <TF1.h>
#include <TAxis.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TMath.h>
#include <THnSparse.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include <string>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include "TRandom3.h"
#include <algorithm>

//#include "ReadMCNP.h"

using namespace std;

//some general values
const double electron_to_eV = 3.77; // //Mean energy to ionize e in Si: keV/e-
THnSparse *img_ediff = NULL;
THnSparse *img_edep = NULL;
THnSparse *img_zdep = NULL;
THnSparse *img_npart = NULL;
THnSparse *img_adcu = NULL;
THnSparse *flag_noise = NULL;
THnSparse *img_adcu_1x100 = NULL;
THnSparse *img_ediff_1x100 = NULL;
THnSparse *flag_noise_1x100 = NULL;

THnSparse *img_ediff_cont = NULL;
THnSparse *img_zdiff_cont = NULL;
THnSparse *img_adcu_cont = NULL;
THnSparse *img_zadcu_cont = NULL;
THnSparse *flag_noise_cont = NULL;
THnSparse *img_adcu_cont_1x100 = NULL;
THnSparse *img_ediff_cont_1x100 = NULL;
THnSparse *flag_noise_cont_1x100 = NULL;

TFile *fout;

//int ccd_size_x = 4150;
//int ccd_size_y = 4150;

TTree *clustersTrue;
TTree *clustersRec;

TTree *clustersRec100; 

//vars clustertree
int ClusterID = 0;
int ccd_size_x, ccd_size_y, Nccds;
int PrimaryPart, NTruePix, NRecPix, clCCDid, clEventID, NSatPixels;
double TrueEnergy, TruePosX, TruePosY, TrueDX, TrueDY, TrueMeanX, TrueMeanY, TrueRMSX, TrueRMSY, TrueQmax, TrueQmaxX, TrueQmaxY, TrueMeanDepth, TrueRMSDepth,
  TrueXmin, TrueXmax, TrueYmin, TrueYmax; 
double RecEnergy, RecPosX, RecPosY, RecDX, RecDY, RecMeanX, RecMeanY, RecRMSX, RecRMSY, RecQmax, RecQmaxX, RecQmaxY, RecMeanDepth, RecRMSDepth,
  RecXmin, RecXmax, RecYmin, RecYmax;

vector<int> TruePixels_x, TruePixels_y; vector<double> TruePixels_z; vector<double> TruePixels_E;
vector<int> RecPixels_x, RecPixels_y; vector<double> RecPixels_z; vector<double> RecPixels_E;


void Usage(string codename="SimulateImages") {
  cout << endl << endl
       << "*******************************************************"  << endl << endl
    << " simulation of the charge diffusion on top of a image " << endl
    //    << " Usage: ./" << codename << " -s <sim_opt> -f <filesim> -i <imagefile> -n <Nsims> -d <dirfiles>" << endl
    << " Usage  ./" << codename << " <configfile>  [default config.txt]" << endl << endl
    << " * SIMMODE: geant4, uni(form) " << endl
    << " * SIMFILE: geant4_filename (mandatory if geant4 option is chosen), neglected if uniform is chosen " << endl
       << "    - number of clusters : used in case of Uniform mode, neglected for geant4 case" << endl << endl
    << " * IMGMODE: img, new (for existing image or new image) " <<endl 
    << "    - if img : inputdir name, filename are mandatory. Filename can be * if all files in the dir should be runed " <<endl
    << "    - if new : a flat image is built for given pedestal and background model. 0 are used otherwise. " << endl
    << " * DIFFUSION parameters [optional, default values : a = *** and b = **** , offset = *** " << endl
    << " * DIFFMODEL: simple [default], only one available " << endl
    << " * Output : FORMAT : fits, root (root as default) " << endl
    << " * OUTFILE can be optionally provided. An automatic filename as: <filein_out.ext> will be used. " << endl
    << " * OUTDIR directory can be optionally provided. The files will be written in the current folder. " << endl
    << " * DEBUG (0, 1, 2) : 0 = only clusters Tree stored (defaul); 1 = clusters tree,  true_edep (in eV) and img after diff+noise+RO (in ADCu) store; 2 = images stored at different stages (edep, diff only, diff + RO)" << endl
    << " * BINFLAG (0, 1, 2) : 0 = only 1x1; 1 = 1x100; 2 = 1x1 and 1x100" << endl
    << " * SEED : used to initialize the random generator [default : timeStamp] "  << endl <<  endl
    << " * THRESHOLD : pixel-signal threshold for cluster search (in ADCu)  "  << endl <<  endl
    << "* M. Settimo - Subatech " << endl << " Aug. 2018 - v0.1 (internal) " << endl << endl;
   
  exit(0);
}



//global var for definition 
//-------------------------------
string simmode, imgmode, simfile, simdir, imgdir,  imgfile, outdir, outfile, format, diffmodel;
int seed, nsim,  sizex, sizey, saturation_adcu;
double depth, diff_A, diff_b, diff_offset, pedestal_mean, pedestal_sigma, e_factor; 
bool diff_a_flag, diff_b_flag, diff_offset_flag;
int debug, binflag, bin_size;
double fano_factor = 0.16;
int nreads = 1 ; //nreads for noise studies

int n_ampli = 1;
double darkcurrent = 0.001; //in e-/pix/day
double readtime = 0.00001; //in s
double Thr_adcu = -1.00;

void CheckConsistency() {

  //  cout << endl << endl << " Checking configuration values " << endl;
  if (simmode=="" || imgmode=="" || outfile=="")  {
    cout << " Error. SIMMODE or IMGMODE or OUTFILE are missed. Please check the configuration file. Exiting "<<endl;
    Usage();
  }

  if (simmode!="geant4" && simmode!="uni"){
    cout << " SIMMODE " << simmode << " unknown. Exiting "  << endl;
     exit(0);
  } 

   if (imgmode!="new" && imgmode!="ext"){
     cout << "IMGMODE: " << imgmode <<  " unknown. Exiting "  << endl;
     exit(0);
  } 
    
  if (simmode=="geant4" && (simfile=="" || simdir=="") ){
    cout << " Error: SIMFILE and SIMDIR mandatory with geant4 option " << endl;
    Usage();
  }

  if (simmode=="uni" && nsim<=0){
    cout << " Error: NCLUSTERS > 0 needed with SIMMODE = uni" << endl;
    Usage();
  }

  if (imgmode=="new" && pedestal_mean<0){
    cout << " Error: PEDESTAL >= 0 needed with IMGMODE = new" << endl;
    Usage();
  }

  if (imgmode=="new" && pedestal_sigma<0){
    cout << " Error: NOISE PIXELS (SIGMA) in e- >= 0 needed with IMGMODE = new" << endl;
    Usage();
  }

  if (binflag<0 || binflag>2) {
    cout << " Error: BINning storing flag: POSSIBLE VALUES 0, 1, 2 " << endl;
    Usage();
  }

if (debug <0 || debug>2) {
    cout << " Error: DEBUG storing options: POSSIBLE VALUES 0, 1, 2 " << endl;
    Usage();
  }

   if (imgmode=="new" && saturation_adcu<0){
    cout << " Error: SATURATION >= 0 needed with IMGMODE = new" << endl;
    Usage();
  }
 
  if (imgmode=="ext" && (imgfile=="" || imgdir=="") ){
    cout << " Error: IMGFILE and IMGDIR mandatory with IMGMODE = ext " << endl;
    Usage();
  }

  if (outfile=="") {
     cout << " Error: OUTFILE mandatory "  << endl;
     Usage();
  }

   if (darkcurrent<0) {
     cout << " Error: DARKCURRENT mandatory "  << endl;
     Usage();
  }

    if (Thr_adcu<0) {
     cout << " WARNING: Thr_adcu for cluster search not defined. Setting to default value (= 5*sigma_noise [ADCu])  "  << endl;
     Thr_adcu = 5*pedestal_sigma*electron_to_eV/(e_factor*1000); //5sigma
  }
   if (n_ampli==3 || n_ampli>4 || n_ampli <=0) {
     cout << " Error: AMPLI mandatory. Possible values = 1, 2, 4 "  << endl;
     Usage();
  }
 
  if (format=="" || !(format=="fits" || format=="root")) {
    format = "root";
    cout << " Warning : OUTFORMAT unknown. Set to default value " << format << endl;
  } 
 
  if (outdir=="") {
    cout << " Warning: OUTDIR not provided, Storing in the local directory "  << endl;
  }
  
  if (!diff_a_flag || !diff_b_flag || !diff_offset_flag) {
    cout << " Default diffusion parameters used : a = " << diff_A <<  "; b = " << diff_b << " diff_off " << diff_offset <<  endl;
  }
}

//----------------------------------

void ReadConfigFile(string cfgname) {
  //read and parse the config file 
  ifstream file(cfgname.c_str());

  if (!file.is_open()) {
   cout << " File " << cfgname << " does not exist or is corrupted. " << endl;
   Usage();
  }
  const char comment = '#';
  cout << endl<< endl;
  for (std::string line; std::getline(file, line); )
    {
       std::istringstream iss(line);
       std::string id, eq, val, val2, val3;     

       if (!(iss >> id) ) 
	 continue;      
       else {
	 if (id[0] == comment) {
	 //this is a comment 
	 continue;
	 }
	  else if (!(iss >> val) ) {
	    continue;
	  }
	 if (id == "SIZE"){ 
	   iss >> val2 >> val3;
	   sizex = stoi(val);
	   sizey = stoi(val2);
	   depth = stof(val3);
	 }
	 //sim:
	 else if (id == "SIMMODE" || id == "simmode") {
	   simmode = val;
	   cout << "simmode: " << simmode << endl;
	 }
	 else if (id == "SIMFILE" || id == "simfile") {
  	   simfile = val;
	   cout << "simfile: " << simfile << endl;
	 }
	 else if (id == "SIMDIR" || id == "simdir") {
  	   simdir = val;
	   cout << "simdir: " << simdir << endl;
	 }
	 else if (id == "NCLUSTERS" ) {
  	   nsim = stoi(val);
	   cout << "clusters: " << nsim << endl;
	 }
	 //input image
	 else if (id == "IMGMODE" || id == "imgmode") {
	   imgmode = val;
	   cout << "imgmode: " << imgmode << endl;
	 }
	 else if (id == "IMGFILE" || id == "imgfile") {
	   imgfile = val;
	   cout << "imgfile: " << imgfile << endl;
	 }
	 else if (id == "IMGDIR" || id == "imgdir") {
  	   imgdir = val;
	   cout << "imgdir: " << imgdir << endl;
	 }
	 else if (id == "PEDESTAL") { 
	    pedestal_mean = stof(val);
	    cout << "pedestal_mean = " << pedestal_mean << endl;
	 }
	 else if (id == "DEBUG") { //values : 0, 1, 2 
	    debug = stoi(val);
	    cout << "debug = " << debug << endl;
	 }
	 else if (id == "BINFLAG") { //values : 0, 1, 2 
	    binflag = stoi(val);
	    cout << "binflag = " << binflag << endl;
	 }
	 else if (id == "SIGMA") { 
	    pedestal_sigma = stof(val);
	    cout << "pedestal_sigma = " << pedestal_sigma << endl; //in e-
	 }
	 else if (id == "SATURATION") { 
	   saturation_adcu = stoi(val);
	    cout << "saturation_adcu = " << saturation_adcu << endl; //in adcu
	 }
	  else if (id == "DARKCURRENT") { 
	   darkcurrent = stof(val);
	    cout << "dark current = " << darkcurrent << endl; //in e-
	 }
	  else if (id == "AMPLI") { 
	   n_ampli = stoi(val);
	    cout << "n_ampli = " << n_ampli << endl; 
	 }
	  else if (id == "THR_ADCU") {
	    Thr_adcu = stof(val);
	    cout << "Threshold for cluster search = " << Thr_adcu << "[ADCu] "<< endl;
	 }
	  else if (id == "READTIME") { 
	   readtime = stof(val);
	   cout << "readtime = " << readtime << " s " << endl; //in s
	 }
	 //outfiles
         else if (id == "OUTFORMAT" || id == "outformat") {
	   format = val;
	   cout << "out format: " << format << endl;
	 }
	  else if (id == "OUTFILE" || id == "outfile") {
	   outfile = val;
	   cout << "outfile: " << outfile << endl;
	 }
	 else if (id == "OUTDIR" || id == "outdir") {
  	   outdir = val;
	   cout << "outdir: " << outdir << endl;
	 }
	 //diffusion
	 else if (id == "DIFF_MODEL" ) {
	   //       diffmodel = val;
	   // cout << "diff model: " << diffmodel << endl;
	   cout << "Warning : only diffusion_model =  simple is currently implemented. skipping assignement"  << endl;
	 }
	  else if (id == "DIFF_A") {
	    diff_A = stof(val);
	    diff_a_flag = true;
	    cout << "diff a = " << diff_A << endl;
	 }
	  else if (id == "DIFF_B") {
	    diff_b = stof(val);
	    diff_b_flag = true;
	    cout << "diff b = " << diff_b << endl;
	 }
	 else if (id == "DIFF_OFFSET") {
	    diff_offset = stof(val);
	    diff_offset_flag = true;
	    cout << "diff offset = " << diff_offset << endl;
	 }
	 else if (id == "SEED" || id == "seed")  {
	   seed = stoi(val);
	   cout << "seed " << seed << endl;
	 }
	 else if (id == "e_FACTOR" || id == "e_factor")  {
	   e_factor = stof(val);
	   cout << "e_factor " << e_factor << endl;
	 }
	 else if (id == "NPIXELX" || id == "PIXELSIZEX" ) {
  	   ccd_size_x = stoi(val);
	   cout << "CCD x pixel size: " << ccd_size_x << endl;
	 }
	 else if (id == "NPIXELY" || id == "PIXELSIZEY" ) {
  	   ccd_size_y = stoi(val);
	   cout << "CCD y pixel size: " << ccd_size_y << endl;
	 }
	 else if (id == "Nccds" || id == "NCCD" || id == "NCCDS") {
	 	 Nccds = stoi(val);
	 	 cout << "Number of CCDs: " << Nccds << endl;
	 }

       //    else if (id[1] == "PATH" ) {
       // }
       // else if (id[1] == "FILENAME" ) {
       // }
       
       } 
    }//end for

  CheckConsistency();
  return;
}


int CalculateElectrons(double Edep_eV) {
  return Edep_eV/3.77;
}




//=====================================
void DiffuseElectrons( double x_p, double y_p, double z_p, double ene, int CCDId, vector<double> &binxdiff,  vector<double> &binydiff)
{
  TRandom3 *r3 = new TRandom3();
  binxdiff.clear();
  binydiff.clear();
  double var = 0;
  Double_t pix_sizex = 15.; //Default pixel size in um - overwritten by simulation
  Double_t pix_sizey = 15.;
  //MS : add option  if(sim_flip_CCD==false)  ???  //false means the depth is calculated from the top ?? 

  if(z_p > diff_offset)
  {
  	  var = -diff_A*TMath::Log(1. - diff_b*(z_p-diff_offset));
  }
  
  if(var > 0) 
  {
  	  var=TMath::Sqrt(var);
  }
  else
  {
  	  var=0;
  }
  
  double x_shift=0;
  double y_shift=0;
 
  // Calculate # of observed electrons (w/ Fano)
  // TMath:Nint  Round to nearest integer. Rounds half integers to the nearest even integer
  int eobs = TMath::Nint(r3->Gaus(ene/electron_to_eV, TMath::Sqrt(fano_factor*(ene/electron_to_eV))));
  
  for(int j=0; j<eobs; j++)
  {
  	  if(var>0)
  	  {
  	  	  y_shift=r3->Gaus(0,var);
  	  	  x_shift=r3->Gaus(0,var);
  	  }
  	  else 
  	  {
  	  	  y_shift=0;
  	  	  x_shift=0;
  	  }
  	  
  	  binxdiff.push_back((int)((x_p + x_shift)/pix_sizex )); //MS //org: ((x_p + x_shift)/pix_sizex) + 1;
  	  binydiff.push_back((int)((y_p + y_shift)/pix_sizey ));
  	  
  	  //if (var>0) cout << " var " << var << "  x_p " <<  x_p << " " << x_shift <<" " << " pix_sizex " 
  	  //<< pix_sizex <<  " binxdiff.push_back(" <<  (int)((x_p + x_shift)/pix_sizex ) << " ) " << endl;

  	}//end j<eobs
  delete r3;
}//end function




//-------------------------------

void ReadImage(){
    //read img size from file header
  //    sizex = .. 
  }

//-------------------------------

void ConvertToFits(){
  
  //add header, convert to fits
    //  fstream fout;
  //fout.open(filenameout.c_str(), std::fstream::out);

}

void WriteOut() {
  //create file name + write out in root and convert to fits if needed 
}

void GenerateEdep() {
  //generate Edep uniformly or with another generator
}


void FindClusters(int nEventID, int primPDG) 
{
	
  	//read img_adcu_cont and img_edep for Rec and True clusters respectively
  	Int_t* x = new Int_t[3];//x,y,ccdid
  	int clusterid = 0;
  	int last_cid_found = 0; //just a counter
  	bool foundNeigh = false;
  	unsigned long touch = 0;
  	int ncrowns = 2; //numb or rings which we look as touching pixels

  	std::vector<unsigned long long>::iterator it;

   	//========================
  	//------- true image ------
  	//=========================
  	//true infos on the clusters 
  	vector<unsigned long long> vpos; //max value 6560006000 = 6e9 //long long uses 64b
  	vector<int> vcid; //clusters id
  	vector<double> clusterEne; //cumulative energy of the found clusters, summing up the energy.
  	vector<int> clusterNpix; //counter of the npix belonging the a given clusters
  	vector<int> clusterCCDid; //counter of the nCCD 
  	vector<int> clusterMinX; //Xmin
  	vector<int> clusterMaxX; //Xmax 
  	vector<int> clusterMinY; //Ymin
   	vector<int> clusterMaxY; //Ymax
   	vector<double> clusterTmpX; //temporary vector to store the X coord
   	vector<double> clusterTmpY; //temporary vector to store the X coord
   	vector<double> clusterMeanX; //weighted average X 
   	vector<double> clusterMeanY; //weighted average X 
   	vector<double> clusterMeanZ; //weighted average X 
   	vector<double> clusterRMSX;  //RMS X 
   	vector<double> clusterRMSY; //RMS Y 
   	vector<double> clusterRMSZ; //RMS Z 
   	vector<double> clusterQmaxX; //X-coord of Qmax 
   	vector<double> clusterQmaxY; //Y-coord of Qmax
   	vector<double> clusterQmax;  //max charge (if more pixels are saturate it is the first one found)
   	//temp vectors used to construct the clusters pixels vectors 
   	vector<int> tmp_cluster_id;
   	vector<int> tmp_pixels_x;
   	vector<int> tmp_pixels_y;
   	vector<double> tmp_pixels_z;
   	vector<double> tmp_pixels_E;
   	int index = 0;
   	
   	last_cid_found=clusterid=index=0;
   	
   	
   	for (Long64_t l = 0; l < img_edep->GetNbins(); ++l) 
   	{
   		Double_t ni = img_npart->GetBinContent(l, x);    
   		Double_t zi = img_zdep->GetBinContent(l, x);
   		Double_t vi = img_edep->GetBinContent(l, x);
   
    	foundNeigh = false;
    	if (vi>0) 
    	{
    		tmp_pixels_x.push_back(x[0]);
      		tmp_pixels_y.push_back(x[1]);
      		//      cout << " zi " << zi << " " << " ni "  << ni << " " << zi/ni << endl;
      		zi/=ni;
      		tmp_pixels_z.push_back(zi);
      
      		tmp_pixels_E.push_back(vi);
      
      		//look for existing (already found) neighboors
      		if (vpos.size()>0) 
      		{
      			for (int m = -ncrowns; m<ncrowns+1; ++m) 
      			{ 
      				for (int n = -ncrowns; n<ncrowns+1; ++n) 
      				{
      					if (m==0 && n==0) continue;
      					touch = (x[0]+m) + (x[1]+n)*10000 + x[2]*1e8;
      					it = std::find(vpos.begin(), vpos.end(), touch);
      					if (it != vpos.end()) 
      					{
      						index = std::distance(vpos.begin(), it);
      						//by pix
      						foundNeigh = true;
      						clusterid = vcid.at(index);
      						tmp_cluster_id.push_back(clusterid);  //add a new entry for each pix 
      						//by clusters
      						clusterEne.at(clusterid-1)+= vi/1000.; //in keV
      						clusterNpix.at(clusterid-1)++;
      						if (vi/1000. > clusterQmax.at(clusterid-1))
      						{ 
      							clusterQmax.at(clusterid-1) = vi/1000.;
      							clusterQmaxX.at(clusterid-1) = x[0];
      							clusterQmaxY.at(clusterid-1) = x[1];
      						}
      						
      						if (x[0]<clusterMinX.at(clusterid-1)) 
							{
								clusterMinX.at(clusterid-1) = x[0];
							}
	    					
	    					if (x[0]>clusterMaxX.at(clusterid-1)) 
							{
								clusterMaxX.at(clusterid-1) = x[0];
							}

							if (x[1]<clusterMinY.at(clusterid-1))
							{
								clusterMinY.at(clusterid-1) = x[1];
							}

							if (x[1]>clusterMaxY.at(clusterid-1))
							{
								clusterMaxY.at(clusterid-1) = x[1];
							}
	    
	    					clusterTmpX.at(clusterid-1) += (x[0]*vi/1000.);
	    					clusterTmpY.at(clusterid-1) += (x[1]*vi/1000.);
	    					clusterMeanX.at(clusterid-1) += (x[0]);
	    					clusterMeanY.at(clusterid-1) += (x[1]);
	    					clusterMeanZ.at(clusterid-1) += zi; //from img_zdep
	    					clusterRMSX.at(clusterid-1) += (x[0]*x[0]);
	    					clusterRMSY.at(clusterid-1) += (x[1]*x[1]);
	    					clusterRMSZ.at(clusterid-1) += (zi*zi);
	    					break;
	    				}
	    			}
	    			if (foundNeigh)
					{
						break;
					}
				}//end loop over x,y crowns;
			}
			//no neighboors found, creating new cluster 
			if (foundNeigh == false)
			{
				last_cid_found++;
				clusterid = last_cid_found;
				tmp_cluster_id.push_back(clusterid);  //add a new entry for each pix 
				clusterEne.push_back(vi/1000.);
				clusterNpix.push_back(1);
				clusterCCDid.push_back(x[2]); //dont need further  updates
				clusterMinX.push_back(x[0]);
				clusterMaxX.push_back(x[0]);	  
				clusterMinY.push_back(x[1]);
				clusterMaxY.push_back(x[1]);	  
				clusterTmpX.push_back(x[0]*vi/1000.);
				clusterTmpY.push_back(x[1]*vi/1000.);
				clusterMeanX.push_back(x[0]);
				clusterMeanY.push_back(x[1]);
				clusterMeanZ.push_back(zi);
				clusterRMSX.push_back(x[0]*x[0]);
				clusterRMSY.push_back(x[1]*x[1]);
				clusterRMSZ.push_back(zi*zi);
				clusterQmax.push_back(vi/1000.);
				clusterQmaxX.push_back(x[0]);
				clusterQmaxY.push_back(x[1]);
			}
			
			//store this pixel info
			vpos.push_back(x[2]*1e8 + x[1]*10000 + x[0] );//pixel centers (x,y) 
			vcid.push_back(clusterid); //clusterid the pixel belongs to
		} //end if vi>0    
		
	}//end loop over imp_edep
	
	//calculate generic informations for each cluster and fill tree infos
	for (int k=0; k<clusterEne.size(); ++k) 
	{
		TruePixels_x.clear();
		TruePixels_y.clear();
		TruePixels_z.clear();
		TruePixels_E.clear();
		
		clEventID = nEventID;
		clCCDid = clusterCCDid.at(k);
		ClusterID = k+1;
		PrimaryPart = primPDG ;
		//   Time = 0;
		TrueEnergy = clusterEne.at(k);
		TruePosX = clusterTmpX.at(k)/clusterEne.at(k);
		TruePosY = clusterTmpY.at(k)/clusterEne.at(k);
		TrueDX = clusterMaxX.at(k)-clusterMinX.at(k);
		TrueDY = clusterMaxY.at(k)-clusterMinY.at(k);
		NTruePix = clusterNpix.at(k);

    	TrueMeanX = clusterMeanX.at(k)/NTruePix;
    	TrueMeanY = clusterMeanY.at(k)/NTruePix;
    	TrueMeanDepth = clusterMeanZ.at(k)/NTruePix;
    	
    	TrueRMSX = TMath::Sqrt(clusterRMSX.at(k)/NTruePix);
    	TrueRMSY = TMath::Sqrt(clusterRMSY.at(k)/NTruePix);
    	TrueRMSDepth = TMath::Sqrt(clusterRMSZ.at(k)/NTruePix);
    	
    	TrueXmin = clusterMinX.at(k);
    	TrueXmax = clusterMaxX.at(k);
    	TrueYmin = clusterMinY.at(k);
    	TrueYmax = clusterMaxY.at(k);
    	TrueQmax = clusterQmax.at(k);
    	TrueQmaxX = clusterQmaxX.at(k);
    	TrueQmaxY = clusterQmaxY.at(k);    

    	//fill pixels_vectors
    	int npix = tmp_cluster_id.size();  
    	for (int jj = 0; jj< npix; ++jj) 
    	{
    		if (tmp_cluster_id.at(jj) == ClusterID) 
    		{
    			TruePixels_x.push_back(tmp_pixels_x.at(jj));
    			TruePixels_y.push_back(tmp_pixels_y.at(jj));
    			TruePixels_z.push_back(tmp_pixels_z.at(jj));
    			TruePixels_E.push_back(tmp_pixels_E.at(jj));
    		}
    	}//end jj
    	clustersTrue->Fill();
    }
    
    //========================
    //========================
  	//------- Rec image ------
  	//=========================

  	//rec infos on the clusters 
  	vector<unsigned long long> vposrec; //max value 6560006000 = 6e9 //long long uses 64b
  	vector<int> vcidrec; //clusters id
  	vector<int> clusterRecNpix; //counter of the npix belonging the a given clusters
  	vector<int> clusterRecNSatpix; //counter of the saturated pix belonging the a given clusters
  	vector<int> clusterRecCCDid;
  	vector<double> clusterRecEne; //cumulative energy of the found clusters, summing up the energy.
    vector<double> clusterRecMinX; 
    vector<double> clusterRecMaxX; 
    vector<double> clusterRecMinY; 
    vector<double> clusterRecMaxY; 
    vector<double> clusterRecTmpX; 
    vector<double> clusterRecTmpY;
    vector<double> clusterRecMeanX; 
    vector<double> clusterRecMeanY;
    vector<double> clusterRecMeanZ;
    vector<double> clusterRecRMSX; 
    vector<double> clusterRecRMSY;
    vector<double> clusterRecRMSZ;
    vector<double> clusterRecQmaxX; //counter of the npix belonging the a given clusters
    vector<double> clusterRecQmaxY; //counter of the npix belonging the a given clusters
    vector<double> clusterRecQmax; //counter of the npix belonging the a given clusters
    //temp vectors used to construct the clusters pixels vectors 
    vector<int> tmp_reccluster_id;
    vector<int> tmp_recpixels_x;
    vector<int> tmp_recpixels_y;
    vector<double> tmp_recpixels_z;
    vector<double> tmp_recpixels_E;

  	last_cid_found=clusterid=0;

  	for (Long64_t l = 0; l < img_adcu_cont->GetNbins(); ++l) 
  	{
  		Double_t zi = img_zadcu_cont->GetBinContent(l, x);
  		Double_t vi = img_adcu_cont->GetBinContent(l, x);

  		foundNeigh = false;
  		//bins in the diffuse and cont histo are not  ordered in x,y so the "step-by-step" 
  		//online method cannot work. A merge is applied to correct for this issue.
  		if (vi>Thr_adcu) 
  		{
  			tmp_recpixels_x.push_back(x[0]);
  			tmp_recpixels_y.push_back(x[1]);
  			tmp_recpixels_z.push_back(zi);
  			tmp_recpixels_E.push_back(vi); //in adcu
  			
  			if(vposrec.size()>0)
  			{
  				//look for existing (already found) neighboors
  				for (int m = -ncrowns; m<ncrowns+1; ++m) 
  				{ 
  					for (int n = -ncrowns; n<ncrowns+1; ++n) 
  					{
  						if (m==0 && n==0)
						{
							continue;
						}

						touch = (x[0]+m) + (x[1]+n)*10000 + x[2]*1e8;
						
						it = std::find(vposrec.begin(), vposrec.end(), touch);
						
						if (it != vposrec.end()) 
						{
							index = std::distance(vposrec.begin(), it);
							
							foundNeigh = true;
							clusterid = vcidrec.at(index);
							tmp_reccluster_id.push_back(clusterid);  //add a new entry for each pix 
							
							//by clusters
							clusterRecEne.at(clusterid-1)+= vi*e_factor; //in keV //old vi/1000.
							clusterRecNpix.at(clusterid-1)++;
							
							if (vi>=saturation_adcu)
							{
								clusterRecNSatpix.at(clusterid-1)++;
							}
							
							if (vi*e_factor > clusterRecQmax.at(clusterid-1))
							{ 
								clusterRecQmax.at(clusterid-1) = vi*e_factor;
								clusterRecQmaxX.at(clusterid-1) = x[0];
								clusterRecQmaxY.at(clusterid-1) = x[1];
							}
							
							if (x[0]<clusterRecMinX.at(clusterid-1))
							{
								clusterRecMinX.at(clusterid-1) = x[0];
							}
							
							if (x[0]>clusterRecMaxX.at(clusterid-1))
							{
								clusterRecMaxX.at(clusterid-1) = x[0];
							}

							if (x[1]<clusterRecMinY.at(clusterid-1))
							{
								clusterRecMinY.at(clusterid-1) = x[1];
							}
							
							if (x[1]>clusterRecMaxY.at(clusterid-1))
							{
								clusterRecMaxY.at(clusterid-1) = x[1];
							}

							clusterRecTmpX.at(clusterid-1) += (x[0]*vi*e_factor);
							clusterRecTmpY.at(clusterid-1) += (x[1]*vi*e_factor);
							clusterRecMeanX.at(clusterid-1) += (x[0]);
							clusterRecMeanY.at(clusterid-1) += (x[1]);
	    					
	    					if (zi!=-9999)
							{
								clusterRecMeanZ.at(clusterid-1) += zi;
							}
	    
						    clusterRecRMSX.at(clusterid-1) += (x[0]*x[0]);
						    clusterRecRMSY.at(clusterid-1) += (x[1]*x[1]);
						    if (zi!=-9999)
							{
								clusterRecRMSZ.at(clusterid-1) += (zi*zi);
							}
							break;
						}
					}
					
					if (foundNeigh)
					{
						break;
					}
				
				}//end loop over x,y crowns;
			} //end if size   
			//no neighboors found, creating new cluster

			if (foundNeigh == false)
			{
				last_cid_found++;
				clusterid = last_cid_found;
				tmp_reccluster_id.push_back(clusterid);  
				clusterRecEne.push_back(vi*e_factor); //keV
				clusterRecNpix.push_back(1);
				clusterRecNSatpix.push_back(0);
				clusterRecCCDid.push_back(floor(x[2])); //dont need any update 
				clusterRecMinX.push_back(x[0]);
				clusterRecMaxX.push_back(x[0]);
				clusterRecMinY.push_back(x[1]);
				clusterRecMaxY.push_back(x[1]);
				clusterRecTmpX.push_back(x[0]*vi*e_factor);
				clusterRecTmpY.push_back(x[1]*vi*e_factor);
				clusterRecMeanX.push_back(x[0]);
				clusterRecMeanY.push_back(x[1]);
				
				if (zi!=-9999)
				{
					clusterRecMeanZ.push_back(zi); 
				}
				else
				{
					clusterRecMeanZ.push_back(0); 
				}

	  			clusterRecRMSX.push_back(x[0]*x[0]);
				clusterRecRMSY.push_back(x[1]*x[1]);
	  			if (zi!=-9999)
				{
					clusterRecRMSZ.push_back(zi*zi);
				}
				else
				{
					clusterRecRMSZ.push_back(0); 
				}
				clusterRecQmax.push_back(vi*e_factor);
				clusterRecQmaxX.push_back(x[0]);
	  			clusterRecQmaxY.push_back(x[1]);
			}
			//store this pixel info

			vposrec.push_back(x[2]*1e8 + x[1]*10000 + x[0] );//pixel centers (x,y) 
			vcidrec.push_back(clusterid); //clusterid the pixel belongs to
		
		} //end if vi>0    
	}//end loop over img_adcu
	
	//merging close-by-events?
	//to be applied only to clusters having PosX,PosY touching each other
	//calculate generic informations for each cluster and fill tree infos
	
	vector<unsigned long long> mergeMin;
	vector<unsigned long long> mergeMax;
	unsigned long long touch1 = 0;
	unsigned long long touch2 = 0;
	vector<int> merged;
	int size = clusterRecEne.size();
	merged.resize(size);
	bool HasMerged=0;
	int newClusterID = 0;
	int ClusterID1 = 0;
	int ClusterID2 = 0;

  	for (int k=0; k<size; ++k) 
  	{
  		if (merged.at(k) == 1)
		{
			continue; //if already merged no need to do anything
		}
    	HasMerged=0;
  
    	clEventID = nEventID;
    	PrimaryPart = primPDG ;
    	clCCDid = clusterRecCCDid.at(k);
    	if (clusterRecEne.at(k)>saturation_adcu)   {
	  NSatPixels++;	}
    	RecEnergy = clusterRecEne.at(k); //in keV
    	RecPosX = clusterRecTmpX.at(k)/clusterRecEne.at(k);
    	RecPosY = clusterRecTmpY.at(k)/clusterRecEne.at(k);
    	RecDX = clusterRecMaxX.at(k)-clusterRecMinX.at(k);
    	RecDY = clusterRecMaxY.at(k)-clusterRecMinY.at(k);
    	NRecPix = clusterRecNpix.at(k);
    	ClusterID1 = k+1;

    	if (k == size-1) 
    	{ //for last event need the vector of pix to be filled here
    		//fill pixels_vectors
    		int npix1 = tmp_reccluster_id.size();
    		//   int ncl = 0;	  
    		for (int jj = 0; jj< npix1; ++jj) 
    		{
    			if (tmp_reccluster_id.at(jj) == ClusterID1) 
    			{ //add to RecVector the pixels of cluster k and j	       
    				RecPixels_x.push_back(tmp_recpixels_x.at(jj));
    				RecPixels_y.push_back(tmp_recpixels_y.at(jj));
    				RecPixels_z.push_back(tmp_recpixels_z.at(jj));
    				RecPixels_E.push_back(tmp_recpixels_E.at(jj));	      	        
    				// ncl++;
    			}
    		}//end jj
    	} 
    	else 
    	{
    		for (int j=k+1; j<size; ++j)
    		{
    			if (merged.at(j) == 1) continue; //needed to further check if already merged 
    			if (clusterRecCCDid.at(j)!=clCCDid) continue;
    			double RecPosX2 = clusterRecTmpX.at(j)/clusterRecEne.at(j);
    			double RecPosY2 = clusterRecTmpY.at(j)/clusterRecEne.at(j);
    			double RecDX2 = clusterRecMaxX.at(j)-clusterRecMinX.at(j);
    			double RecDY2 = clusterRecMaxY.at(j)-clusterRecMinY.at(j);
    			ClusterID2 = -1;
    			if (TMath::Abs(RecPosX - RecPosX2) < RecDX+RecDX2+1 && TMath::Abs(RecPosY - RecPosY2) < RecDY+RecDY2+1) 
    			{
    				merged.at(j) = 1;
    				HasMerged=1;
    				ClusterID2 = j+1;

    				clusterRecTmpX.at(k)+=clusterRecTmpX.at(j);
    				clusterRecTmpY.at(k)+=clusterRecTmpY.at(j);
    				clusterRecMeanX.at(k)+=clusterRecMeanX.at(j);
    				clusterRecMeanY.at(k)+=clusterRecMeanY.at(j);
    				clusterRecMeanZ.at(k)+=clusterRecMeanZ.at(j);
    				
    				clusterRecRMSX.at(k)+=clusterRecRMSX.at(j);
    				clusterRecRMSY.at(k)+=clusterRecRMSY.at(j);
    				clusterRecRMSZ.at(k)+=clusterRecRMSZ.at(j);
    				
    				clusterRecEne.at(k)+=clusterRecEne.at(j);
    				clusterRecNpix.at(k)+=clusterRecNpix.at(j);
    				clusterRecNSatpix.at(k)+=clusterRecNSatpix.at(j);
    				
    				if (clusterRecQmax.at(j)>clusterRecQmax.at(k)) 
    				{
    					clusterRecQmax.at(k)=clusterRecQmax.at(j);
    					clusterRecQmaxX.at(k)=clusterRecQmaxX.at(j);
    					clusterRecQmaxY.at(k)=clusterRecQmaxY.at(j);
    				}

    				if (clusterRecMaxX.at(j)>clusterRecMaxX.at(k))
					{
						clusterRecMaxX.at(k)=clusterRecMaxX.at(j);
					}
				    if (clusterRecMaxY.at(j)>clusterRecMaxY.at(k))
					{
						clusterRecMaxY.at(k)=clusterRecMaxY.at(j);
					}
					if (clusterRecMinX.at(j)<clusterRecMinX.at(k))
					{
						clusterRecMinX.at(k)=clusterRecMinX.at(j);
					}
					if (clusterRecMinY.at(j)<clusterRecMinY.at(k))
					{ 
						clusterRecMinY.at(k)=clusterRecMinY.at(j);
					}
				
				}//end if
				
				//fill pixels_vectors

				int npix1 = tmp_reccluster_id.size();
				
				for (int jj = 0; jj< npix1; ++jj) 
				{
					if (tmp_reccluster_id.at(jj) == ClusterID1 || tmp_reccluster_id.at(jj) == ClusterID2) 
					{ //add to RecVector the pixels of cluster k and j	       
						RecPixels_x.push_back(tmp_recpixels_x.at(jj));
						RecPixels_y.push_back(tmp_recpixels_y.at(jj));
						RecPixels_z.push_back(tmp_recpixels_z.at(jj));
						RecPixels_E.push_back(tmp_recpixels_E.at(jj));	      	        
						// ncl++;
					}
				}//end jj	   
			}//end for  j
		}
		
		if (HasMerged==1) 
		{ 
			k--; continue;
		} //to allow for more than 2 pieces
		newClusterID++;
		ClusterID = newClusterID;
		RecEnergy = clusterRecEne.at(k); //in keV
		RecPosX = clusterRecTmpX.at(k)/clusterRecEne.at(k);
		RecPosY = clusterRecTmpY.at(k)/clusterRecEne.at(k);
		RecDX = clusterRecMaxX.at(k)-clusterRecMinX.at(k);
		RecDY = clusterRecMaxY.at(k)-clusterRecMinY.at(k);
		NRecPix = clusterRecNpix.at(k);
		NSatPixels = clusterRecNSatpix.at(k);
		
		RecXmin = clusterRecMinX.at(k);
		RecXmax = clusterRecMaxX.at(k);
		RecYmin = clusterRecMinY.at(k);
		RecYmax = clusterRecMaxY.at(k);
		
		RecMeanX = clusterRecMeanX.at(k)/NRecPix;
		RecMeanY = clusterRecMeanY.at(k)/NRecPix;
		RecMeanDepth = clusterRecMeanZ.at(k)/NRecPix;
		
		RecRMSX = TMath::Sqrt(clusterRecRMSX.at(k)/NRecPix);
		RecRMSY = TMath::Sqrt(clusterRecRMSY.at(k)/NRecPix);
		RecRMSDepth = TMath::Sqrt(clusterRecRMSZ.at(k)/NRecPix);
		
		RecQmax = clusterRecQmax.at(k);
		RecQmaxX = clusterRecQmaxX.at(k);
		RecQmaxY = clusterRecQmaxY.at(k);
		
		clustersRec->Fill();
		//clear needs to be after filling - not at beginning of the loop
		RecPixels_x.clear();
		RecPixels_y.clear();
		RecPixels_z.clear();
		RecPixels_E.clear();
	
	}//end k

  	// =====================================
  	// 1x100
  	// ====================================
  	if (binflag > 0) 
  	{  
  		last_cid_found=clusterid=0;
  		tmp_recpixels_x.clear();
  		tmp_recpixels_y.clear();
  		tmp_recpixels_z.clear();
  		tmp_recpixels_E.clear();
  		tmp_reccluster_id.clear();
  		clusterRecEne.clear();
  		clusterRecNpix.clear();
  		
    	vposrec.clear();
    	vcidrec.clear();
    	clusterRecNSatpix.clear();
    	clusterRecQmax.clear();
    	clusterRecCCDid.clear();
    	clusterRecQmaxY.clear();
    	clusterRecQmaxX.clear();
    	clusterRecMinX.clear();
    	clusterRecMinY.clear();
    	clusterRecMaxX.clear();
    	clusterRecMaxY.clear();
    	clusterRecTmpX.clear();
    	clusterRecTmpY.clear();
    	clusterRecMeanX.clear();
    	clusterRecMeanY.clear();
    	clusterRecMeanZ.clear();
    	clusterRecRMSX.clear();
    	clusterRecRMSY.clear();
    	clusterRecRMSZ.clear();
    	for (Long64_t l = 0; l < img_adcu_cont_1x100->GetNbins(); ++l) 
    	{
    		Double_t zi = 0; //img_zadcu_cont_1x100->GetBinContent(l, x);
    		Double_t vi = img_adcu_cont_1x100->GetBinContent(l, x);
    		
    		foundNeigh = false;
    		//bins in the diffuse and cont histo are not  ordered in x,y 
    		//so the "step-by-step" online method cannot work. A merge is applied 
    		//to correct for this issue.
    		if (vi>Thr_adcu) 
    		{
    			tmp_recpixels_x.push_back(x[0]);
    			tmp_recpixels_y.push_back(x[1]);
    			tmp_recpixels_z.push_back(zi);
    			tmp_recpixels_E.push_back(vi*e_factor);

				if (vposrec.size()>0) 
				{
					//look for existing (already found) neighboors
					for (int m = -ncrowns; m<ncrowns+1; ++m) 
					{ 
						for (int n = -ncrowns; n<ncrowns+1; ++n) 
						{
							if (m==0 && n==0) continue;

							touch = (x[0]+m) + (x[1]+n)*10000 + x[2]*1e8;
							it = std::find(vposrec.begin(), vposrec.end(), touch);

							if (it != vposrec.end()) 
							{
								index = std::distance(vposrec.begin(), it);
								foundNeigh = true;
								clusterid = vcidrec.at(index);
								tmp_reccluster_id.push_back(clusterid);  //add a new entry for each pix 
								
								//by clusters
								clusterRecEne.at(clusterid-1)+= vi*e_factor; //keV
								clusterRecNpix.at(clusterid-1)++;
								
								if (vi>=saturation_adcu)
								{
									clusterRecNSatpix.at(clusterid-1)++;
								}
								
								if (vi*e_factor > clusterRecQmax.at(clusterid-1))
								{ 
									clusterRecQmax.at(clusterid-1) = vi*e_factor;
									clusterRecQmaxX.at(clusterid-1) = x[0];
									clusterRecQmaxY.at(clusterid-1) = x[1];
								}
								
								if (x[0]<clusterRecMinX.at(clusterid-1))
								{
									clusterRecMinX.at(clusterid-1) = x[0];
								}
								if (x[0]>clusterRecMaxX.at(clusterid-1)) 
								{
									clusterRecMaxX.at(clusterid-1) = x[0];
								}
								if (x[1]<clusterRecMinY.at(clusterid-1))
								{
									clusterRecMinY.at(clusterid-1) = x[1];
								}
								if (x[1]>clusterRecMaxY.at(clusterid-1))
								{
									clusterRecMaxY.at(clusterid-1) = x[1];
								}

								clusterRecTmpX.at(clusterid-1) += (x[0]*vi*e_factor); 
								clusterRecTmpY.at(clusterid-1) += (x[1]*vi*e_factor);
		
								clusterRecMeanX.at(clusterid-1) += (x[0]);
								clusterRecMeanY.at(clusterid-1) += (x[1]);
								if (zi!=-9999)
								{
							          	clusterRecMeanZ.at(clusterid-1) += zi;
								}
		
								clusterRecRMSX.at(clusterid-1) += (x[0]*x[0]);
								clusterRecRMSY.at(clusterid-1) += (x[1]*x[1]);
								if (zi!=-9999)
								{
									clusterRecRMSZ.at(clusterid-1) += (zi*zi);
								}
								break;
							}
						}
						
						if (foundNeigh) break;
					}//end loop over x,y crowns;
				} //end if size   
				//no neighboors found, creating new cluster

				if (foundNeigh == false)
				{
					last_cid_found++;
					clusterid = last_cid_found;
					tmp_reccluster_id.push_back(clusterid);  
					clusterRecEne.push_back(vi*e_factor);
					clusterRecNpix.push_back(1);
					clusterRecNSatpix.push_back(0);
					clusterRecCCDid.push_back(floor(x[2])); //dont need any update 
					clusterRecMinX.push_back(x[0]);
					clusterRecMaxX.push_back(x[0]);
					clusterRecMinY.push_back(x[1]);
					clusterRecMaxY.push_back(x[1]);
					clusterRecTmpX.push_back(x[0]*vi*e_factor);
					clusterRecTmpY.push_back(x[1]*vi*e_factor);

					clusterRecMeanX.push_back(x[0]);
					clusterRecMeanY.push_back(x[1]);
					if (zi!=-9999)
					{
						clusterRecMeanZ.push_back(zi); 
					}
					else
					{
						clusterRecMeanZ.push_back(0); 
					}

					clusterRecRMSX.push_back(x[0]*x[0]);
					clusterRecRMSY.push_back(x[1]*x[1]);
					if (zi!=-9999)
					{
						clusterRecRMSZ.push_back(zi*zi);
					}
					else
					{
						clusterRecRMSZ.push_back(0); 
					}
					
					clusterRecQmax.push_back(vi*e_factor);
					clusterRecQmaxX.push_back(x[0]);
					clusterRecQmaxY.push_back(x[1]);
				}
				
				//store this pixel info
				vposrec.push_back(x[2]*1e8 + x[1]*10000 + x[0] );//pixel centers (x,y) 
				vcidrec.push_back(clusterid); //clusterid the pixel belongs to
			
			} //end if vi>0    
		}//end loop over img_adcu
		
		//merging close-by-events?
		//to be applied only to clusters having PosX,PosY touching each other
		//calculate generic informations for each cluster and fill tree infos
		
		vector<unsigned long long> mergeMin;
		vector<unsigned long long> mergeMax;
		unsigned long long touch1 = 0;
		unsigned long long touch2 = 0;
		vector<int> merged;
		int size = clusterRecEne.size();
		merged.resize(size);
		bool HasMerged=0;
		int newClusterID = 0;
		int ClusterID1 = 0;
		int ClusterID2 = 0;
		
		for (int k=0; k<size; ++k) 
		{
			if (merged.at(k) == 1) continue; //if already merged no need to do anything
			HasMerged=0;
			clEventID = nEventID;
			PrimaryPart = primPDG ;
			clCCDid = clusterRecCCDid.at(k);
			if (clusterRecEne.at(k)>saturation_adcu) NSatPixels++;
			RecEnergy = clusterRecEne.at(k); //in eV
			RecPosX = clusterRecTmpX.at(k)/clusterRecEne.at(k);
			RecPosY = clusterRecTmpY.at(k)/clusterRecEne.at(k);
			RecDX = clusterRecMaxX.at(k)-clusterRecMinX.at(k);
			RecDY = clusterRecMaxY.at(k)-clusterRecMinY.at(k);
			NRecPix = clusterRecNpix.at(k);
			ClusterID1 = k+1;
			
			if (k == size-1) 
			{ //for last event need the vector of pix to be filled here
				//fill pixels_vectors
				int npix1 = tmp_reccluster_id.size();
				//   int ncl = 0;	  
				for (int jj = 0; jj< npix1; ++jj) 
				{
					if (tmp_reccluster_id.at(jj) == ClusterID1) 
					{ //add to RecVector the pixels of cluster k and j	       
						RecPixels_x.push_back(tmp_recpixels_x.at(jj));
						RecPixels_y.push_back(tmp_recpixels_y.at(jj));
						RecPixels_z.push_back(tmp_recpixels_z.at(jj));
						RecPixels_E.push_back(tmp_recpixels_E.at(jj));	      	        
						// ncl++;
					}
				}//end jj
			} 
			else 
			{
				for (int j=k+1; j<size; ++j)
				{
					if (merged.at(j) == 1) continue; //needed to further check if already merged 
					if (clusterRecCCDid.at(j)!=clCCDid) continue;
					
					double RecPosX2 = clusterRecTmpX.at(j)/clusterRecEne.at(j);
					double RecPosY2 = clusterRecTmpY.at(j)/clusterRecEne.at(j);
					double RecDX2 = clusterRecMaxX.at(j)-clusterRecMinX.at(j);
					double RecDY2 = clusterRecMaxY.at(j)-clusterRecMinY.at(j);
					ClusterID2 = -1;
					if (TMath::Abs(RecPosX - RecPosX2) < RecDX+RecDX2+1 && TMath::Abs(RecPosY - RecPosY2) < RecDY+RecDY2+1) 
					{
						merged.at(j) = 1;
						HasMerged=1;
						ClusterID2 = j+1;
						
						clusterRecTmpX.at(k)+=clusterRecTmpX.at(j);
						clusterRecTmpY.at(k)+=clusterRecTmpY.at(j);
						clusterRecMeanX.at(k)+=clusterRecMeanX.at(j);
						clusterRecMeanY.at(k)+=clusterRecMeanY.at(j);
						clusterRecMeanZ.at(k)+=clusterRecMeanZ.at(j);
						
						clusterRecRMSX.at(k)+=clusterRecRMSX.at(j);
						clusterRecRMSY.at(k)+=clusterRecRMSY.at(j);
						clusterRecRMSZ.at(k)+=clusterRecRMSZ.at(j);
						
						clusterRecEne.at(k)+=clusterRecEne.at(j);
						clusterRecNpix.at(k)+=clusterRecNpix.at(j);
						clusterRecNSatpix.at(k)+=clusterRecNSatpix.at(j);
						
						if (clusterRecQmax.at(j)>clusterRecQmax.at(k)) 
						{
							clusterRecQmax.at(k)=clusterRecQmax.at(j);
							clusterRecQmaxX.at(k)=clusterRecQmaxX.at(j);
							clusterRecQmaxY.at(k)=clusterRecQmaxY.at(j);
						}

						if (clusterRecMaxX.at(j)>clusterRecMaxX.at(k)) clusterRecMaxX.at(k)=clusterRecMaxX.at(j);
						if (clusterRecMaxY.at(j)>clusterRecMaxY.at(k)) clusterRecMaxY.at(k)=clusterRecMaxY.at(j);
						if (clusterRecMinX.at(j)<clusterRecMinX.at(k)) clusterRecMinX.at(k)=clusterRecMinX.at(j);
						if (clusterRecMinY.at(j)<clusterRecMinY.at(k)) clusterRecMinY.at(k)=clusterRecMinY.at(j);
					
					}//end if
					
					//fill pixels_vectors
					int npix1 = tmp_reccluster_id.size();
					for (int jj = 0; jj< npix1; ++jj) 
					{
						if (tmp_reccluster_id.at(jj) == ClusterID1 || tmp_reccluster_id.at(jj) == ClusterID2) 
						{ //add to RecVector the pixels of cluster k and j	       
							RecPixels_x.push_back(tmp_recpixels_x.at(jj));
							RecPixels_y.push_back(tmp_recpixels_y.at(jj));
							RecPixels_z.push_back(tmp_recpixels_z.at(jj));
							RecPixels_E.push_back(tmp_recpixels_E.at(jj));	      	        
							// ncl++;
						}
					}//end jj	   
				}//end for  j
			}
			
			if (HasMerged==1) 
			{ 
				k--; continue;
			} //to allow for more than 2 pieces
			
			newClusterID++;
			ClusterID = newClusterID;
			RecEnergy = clusterRecEne.at(k); //in eV
			RecPosX = clusterRecTmpX.at(k)/clusterRecEne.at(k);
			RecPosY = clusterRecTmpY.at(k)/clusterRecEne.at(k);
			RecDX = clusterRecMaxX.at(k)-clusterRecMinX.at(k);
			RecDY = clusterRecMaxY.at(k)-clusterRecMinY.at(k);
			NRecPix = clusterRecNpix.at(k);
			NSatPixels = clusterRecNSatpix.at(k);
			
			RecXmin = clusterRecMinX.at(k);
			RecXmax = clusterRecMaxX.at(k);
			RecYmin = clusterRecMinY.at(k);
			RecYmax = clusterRecMaxY.at(k);
			
			RecMeanX = clusterRecMeanX.at(k)/NRecPix;
			RecMeanY = clusterRecMeanY.at(k)/NRecPix;
			RecMeanDepth = clusterRecMeanZ.at(k)/NRecPix;
			
			RecRMSX = TMath::Sqrt(clusterRecRMSX.at(k)/NRecPix);
			RecRMSY = TMath::Sqrt(clusterRecRMSY.at(k)/NRecPix);
			RecRMSDepth = TMath::Sqrt(clusterRecRMSZ.at(k)/NRecPix);
			
			RecQmax = clusterRecQmax.at(k);
			RecQmaxX = clusterRecQmaxX.at(k);
			RecQmaxY = clusterRecQmaxY.at(k);
			
			clustersRec100->Fill();
			//clear needs to be after filling - not at beginning of the loop
			RecPixels_x.clear();
			RecPixels_y.clear();
			RecPixels_z.clear();
			RecPixels_E.clear();
		
		}//end k
	}//end if binflag 
} //end FindClusters



void ReadGeant4Input() 
{
	string str = simdir + "/" + simfile;
	TFile *fin = new TFile(str.c_str());

	TRandom3 *rand = new TRandom3();
	TTree *RunInfo = (TTree*)fin->Get("RunInfo"); 
	TTree *CCDOut = (TTree*)fin->Get("CCDOut");
	//TTree *CCDOut = (TTree*)fin->Get("TrackOut");
	TTree *EventOut = (TTree*)fin->Get("EventOut");
	
	char *title = new char[100];
	char *grname = new char[100];

  	int NEvts, NCCDs;// Seed, CCDVersion, CCDFrameVersion,
  	//  CableVersion, VesselVersion, CryoVersion, ShieldingVersion;

  	// Read Brunch RunInfor to get the number of simulated events (and Nccds!)
  	RunInfo->SetBranchAddress("NEvts", &NEvts);
  	RunInfo->SetBranchAddress("NCCDs", &NCCDs);
  	//  RunInfo->SetBranchAddress("Seed", &Seed);
  	/*  RunInfo->SetBranchAddress("CCDVersion", &CCDVersion);
  	RunInfo->SetBranchAddress("CCDFrameVersion", &CCDFrameVersion);
  	RunInfo->SetBranchAddress("CableVersion", &CableVersion);
  	RunInfo->SetBranchAddress("VesselVersion", &VesselVersion);
  	RunInfo->SetBranchAddress("CryoVersion", &CryoVersion);
	RunInfo->SetBranchAddress("ShieldingVersion", &ShieldingVersion);
  	*/
  	
  	//--- Read the branch CCDOut where the energy depositions are stored
	//org version 
	//	int EventID, pdg, trackid, parentid, CCDid;
  	//double posx, posy, posz, edep_tot, time;
	
	//vector-version
	vector<double> *posx = NULL;
	vector<double> *posy = NULL;
	vector<double> *posz = NULL;
	vector<double> *edep_tot = NULL;
	vector<double> *time = NULL;
	int  EventID;  
	vector<int> *pdg = nullptr;
	vector<int> *trackid=nullptr;
	vector<int> *parentid=nullptr;
	vector<int> *CCDid=nullptr;
  
  
  	CCDOut->SetBranchAddress("EventID", &EventID);
  	CCDOut->SetBranchAddress("CCDid", (&CCDid));
  	CCDOut->SetBranchAddress("posx", (&posx));
  	CCDOut->SetBranchAddress("posy", (&posy));
  	CCDOut->SetBranchAddress("posz", (&posz));
  	CCDOut->SetBranchAddress("Edep", (&edep_tot));
  	CCDOut->SetBranchAddress("time", (&time));  

	// Get the simulated EventId and its PDG number
	int EvtPDG= 0; int EvtID = 0;
   	EventOut->SetBranchAddress("EventID", &EvtID);
   	EventOut->SetBranchAddress("pdg", &EvtPDG);
  
  	// AMPLIFICATION ??
	int Npix_line = ccd_size_x;
  	if(n_ampli == 2 || n_ampli == 4)
	{
		Npix_line = ccd_size_x/2.;
	}
	double readtime_line = readtime * Npix_line;  
 
  	//----- output
  	int fDim = 3; //x,y,ccdid;
  	int fBins = ccd_size_x; //size of ccd in pixels
  	Int_t* bins = new Int_t[fDim];
  	Double_t *xmin = new Double_t[fDim];
  	Double_t *xmax = new Double_t[fDim];

  	Int_t* bins1x100 = new Int_t[fDim];
  	Double_t *xmin1x100 = new Double_t[fDim];
  	Double_t *xmax1x100 = new Double_t[fDim];
 

  	cout << " " << endl;
  	cout << "Using " << Nccds << " CCDs with size "<< ccd_size_x << "x" << ccd_size_x << " pixels" << endl;
  
  	for(Int_t d = 0; d < fDim; ++d) 
  	{
    	if(d < 2)
    	{
    		bins[d] = fBins;
    		xmin[d] = 0;
    		xmax[d] = fBins;
    	}
    	else 
    	{
    		bins[d] = Nccds;
    		xmin[d] = 1;
    		xmax[d] = Nccds+1;
    	}
    }
  
    bins1x100[0] = fBins;
    xmin1x100[0] = 0;
    xmax1x100[0] = fBins;
    
    bins1x100[1] = (int)floor(fBins/100);
    xmin1x100[1] = 0;
    xmax1x100[1] = (int)floor(fBins/100);
					  
    bins1x100[2] = Nccds;
    xmin1x100[2] = 1;
    xmax1x100[2] = Nccds+1;
    sprintf(title, " ; coord x [pix]; coord y [pix];"); 
    
	// pre-step images to build reconstructed cluesters: clustersRec, clustersRec1x100
  	img_edep = new THnSparseF("img_edep", title, fDim, bins, xmin, xmax);
 	img_zdep = new THnSparseF("img_zdep", title, fDim, bins, xmin, xmax);
  	img_npart = new THnSparseF("img_npart", title, fDim, bins, xmin, xmax);
  
  	img_ediff = new THnSparseF("img_ediff", title, fDim, bins, xmin, xmax);
  	img_adcu = new THnSparseF("img_adcu", title, fDim, bins, xmin, xmax);
  	flag_noise = new THnSparseF("flag_noise", title, fDim, bins, xmin, xmax);

  	img_ediff_1x100 = new THnSparseF("img_ediff_1x100", title, fDim, bins1x100, xmin1x100, xmax1x100);
  	img_adcu_1x100 = new THnSparseF("img_adcu_1x100", title , fDim, bins1x100, xmin1x100, xmax1x100);
  	flag_noise_1x100 = new THnSparseF("flag_noise_1x100", title, fDim, bins1x100, xmin1x100, xmax1x100);


  	//continuous readout
  	img_ediff_cont = new THnSparseF("img_ediff_cont", title, fDim, bins, xmin, xmax);
  	img_zdiff_cont = new THnSparseF("img_zdiff_cont", title, fDim, bins, xmin, xmax);
  	img_adcu_cont = new THnSparseF("img_adcu_cont", title, fDim, bins, xmin, xmax);
  	img_zadcu_cont = new THnSparseF("img_zadcu_cont", title, fDim, bins, xmin, xmax);
  	flag_noise_cont = new THnSparseF("flag_noise_cont", title, fDim, bins, xmin, xmax);

  	img_ediff_cont_1x100 = new THnSparseF("img_ediff_cont_1x100", title, fDim, bins1x100, xmin1x100, xmax1x100);
  	img_adcu_cont_1x100 = new THnSparseF("img_adcu_cont_1x100", title , fDim, bins1x100, xmin1x100, xmax1x100);
  	flag_noise_cont_1x100 = new THnSparseF("flag_noise_cont_1x100", title, fDim, bins1x100, xmin1x100, xmax1x100);
  
  	Int_t* x = new Int_t[fDim];
  	int pix_size = 15;//um
  	double bw_x = img_edep->GetAxis(0)->GetBinWidth(1)* pix_size;
  	double bw_y = img_edep->GetAxis(1)->GetBinWidth(1)* pix_size; 
  	Int_t nentries = CCDOut->GetEntries();

  	char *filename = new char[500];
  	int EventID_last = 0;
	//        nentries = 10000;
   	Double_t *xloop = new Double_t[3];
   	double exposure_day = 8./24.;
	int perc_pritout = (int)(nentries/20);
  	for(int i=0; i<nentries; i++)
  	{
	   if(i% perc_pritout == 0)
		{
			cout << " Nentries :: " << i << "/" << nentries << endl;
		}
 
	posx->clear();	posy->clear(); posz->clear();  edep_tot->clear(); time->clear();

    	CCDOut->GetEntry(i);
	if (edep_tot->size()==0) continue;
	//	if (i%100==0)cout << "pippo " << i << " EvtiD : " << EventID << " size " << edep_tot->size() << endl;		
    	/*if(i == 0)
		{
			EventID_last = EventID;
		}*/
    	
    	if(EventID!=EventID_last) 
    	{
    		fout->cd();
    		//------------
       		// writedown the img_ediff previously filled
      		for(Long64_t l = 0; l < img_ediff->GetNbins(); ++l) 
      		{	
      			Double_t vi = img_ediff->GetBinContent(l, x);
      			// convert eV to electrons
      			vi = vi/electron_to_eV;
      			
      			xloop[0] = x[0];
      			xloop[1] = x[1];
      			xloop[2] = x[2];
	
			double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma)); //in e-
			double dc = rand->Poisson(darkcurrent*exposure_day); //in e-/pix/expo
				
			img_adcu->Fill(xloop, min(saturation_adcu, TMath::Nint((vi + noise + dc)*electron_to_eV/(e_factor*1000))) );
			flag_noise->Fill(xloop, 1);
					
			//cout << saturation_adcu << " " << vi << " " << noise << " " << dc << " " << electron_to_eV << " " << e_factor*1000 << endl;
			//cout << TMath::Nint((vi + noise + dc)*electron_to_eV/(e_factor*1000)) << endl;

			// fill the img_ediff_1x100 based on 1x1 previously filled
			if (binflag > 0) 
				{ 
					xloop[0] = x[0];
					xloop[1] = TMath::Nint(x[1]/100);
	
					noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
	  				dc = rand->Poisson(darkcurrent*exposure_day*100);
	  				img_ediff_1x100->Fill(xloop, vi+noise+dc);  //fill 1x100
	  				flag_noise_1x100->Fill(xloop, 1);
				}
	
				//add noise to 1x1 around the clusters 
				for(int bix=-15; bix <= 15; bix++) 
				{
					for(int biy=-15; biy <= 15; biy++) 
					{
						xloop[0] = x[0] + bix;
						xloop[1] = x[1] + biy;
						//xloop[2] = x[2];
						if(xloop[0] > xmax[0] || xloop[0] < xmin[0] || xloop[1]>xmax[1] || xloop[1]<xmin[1]) 
						{
							continue;
						}
						// Long64_t  index = img_ediff->GetBin(xloop, "false");
						Long64_t  index = flag_noise->GetBin(xloop);
						double  flag_pix = flag_noise->GetBinContent(index);
						
						if(flag_pix<=0) 
						{ //0 charges, RO noise not yet added
							double noise = rand->Gaus(pedestal_mean, TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
							dc = rand->Poisson(darkcurrent*exposure_day);
							img_adcu->Fill(xloop, min(saturation_adcu, TMath::Nint((noise + dc)*electron_to_eV/(e_factor*1000))));
							flag_noise->Fill(xloop, 1);
						}
					}
				}//end x, y loop
			}//end img_ediff entries
			
			//------------
			// add noise to 1x100
			if (binflag > 0) 
			{ 
				for (Long64_t l = 0; l < img_ediff_1x100->GetNbins(); ++l) 
				{
					Double_t vi = img_ediff_1x100->GetBinContent(l, x);
	      			// convert eV to electrons
    	  			vi = vi/electron_to_eV;

					xloop[0] = x[0];
					xloop[1] = x[1];
					xloop[2] = x[2];
					
					double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
					double dc = rand->Poisson(darkcurrent*exposure_day*100);
					img_adcu_1x100->Fill(xloop, min(saturation_adcu, TMath::Nint((vi + noise + dc)*electron_to_eV/(e_factor*1000))));
					flag_noise_1x100->Fill(xloop, 1);	

					for (int bix=-15; bix <= 15; bix++) 
					{
						for (int biy=-1; biy <= 1; biy++) 
						{
							xloop[0] = x[0] + bix;
							xloop[1] = x[1] + biy;
							//xloop[2] = x[2];
							if (xloop[0] > xmax1x100[0] || xloop[0] < xmin1x100[0] || xloop[1]>xmax1x100[1] || xloop[1]<xmin1x100[1]) 
							{
								continue;
							}
							Long64_t  index = flag_noise_1x100->GetBin(xloop, "false");
							double  flag_pix = flag_noise_1x100->GetBinContent(index);
							if (flag_pix<=0) 
							{ 
								//noise not yet added --> Should the DC include diffusion? 
								double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
								double dc = rand->Poisson(darkcurrent*exposure_day*100);
								img_adcu_1x100->Fill(xloop, min(saturation_adcu, TMath::Nint((noise + dc)*electron_to_eV/(e_factor*1000))));
								flag_noise_1x100->Fill(xloop, 1);
							}
						}
					}
				}
			} //end binflag

      		//=========================================
      		//CONTINUOUS READOUT -- 1x1 diffusion
      		//=========================================
      		for (Long64_t l = 0; l < img_ediff_cont->GetNbins(); ++l) 
      		{
      			Double_t vi = img_ediff_cont->GetBinContent(l, x);
      			// vi in units of electrons
      			vi = vi/electron_to_eV;
      			Double_t zi = img_zdiff_cont->GetBinContent(l, x);
      			//normalisation to the N. of electrons in the bin since img_zdiff
      			//img_zdiff is filled once per electron;
      			zi/=vi;
      			xloop[0] = x[0];
      			xloop[1] = x[1];
      			xloop[2] = x[2];

				double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma)); //in e-
				double dc = rand->Poisson(darkcurrent*exposure_day); //in e-/pix/expo
				img_adcu_cont->Fill(xloop, min(saturation_adcu, TMath::Nint((vi + noise + dc)*electron_to_eV/(e_factor*1000))) );
				img_zadcu_cont->Fill(xloop, zi); //fill with zi normalized 
				flag_noise_cont->Fill(xloop, 1);
				// fill the img_ediff_1x100 based on 1x1 previously filled
				xloop[0] = x[0];
				xloop[1] = TMath::Nint(x[1]/100);
				//	xloop[2] = x[2];
	
				noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
				dc = rand->Poisson(darkcurrent*exposure_day*100);
				img_ediff_cont_1x100->Fill(xloop, (vi+noise+dc)*electron_to_eV);  //fill 1x100 in electrons in eV
				flag_noise_cont_1x100->Fill(xloop, 1);
				
				//add noise to 1x1 around the clusters 
				for (int bix=-15; bix <= 15; bix++) 
				{
					for (int biy=-15; biy <= 15; biy++) 
					{
						xloop[0] = x[0] + bix;
						xloop[1] = x[1] + biy;
						// xloop[2] = x[2];
						if (xloop[0] > xmax[0] || xloop[0] < xmin[0] || xloop[1]>xmax[1] || xloop[1]<xmin[1])
						{
							continue;
						}
						Long64_t  index = flag_noise_cont->GetBin(xloop);
						double  flag_pix = flag_noise_cont->GetBinContent(index);

	    				if (flag_pix<=0) 
	    				{ //0 charges, RO noise not yet added
	    					double noise = rand->Gaus(pedestal_mean, TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
	    					dc = rand->Poisson(darkcurrent*exposure_day);
	    					img_adcu_cont->Fill(xloop, min(saturation_adcu, TMath::Nint((noise + dc)*electron_to_eV/(e_factor*1000))));
	    					img_zadcu_cont->Fill(xloop, -9999); //for the bins with noise only set z to a faka value easy to discard. 
	    					
	    					flag_noise_cont->Fill(xloop, 1);
	    				}
	    			}
	    		}//end x, y loop
	    	}//end img_ediff_cont entries
	    	
	    	//------------
	    	//// add noise to 1x100
	    	if (binflag>0) 
	    	{
	    		for (Long64_t l = 0; l < img_ediff_cont_1x100->GetNbins(); ++l) 
	    		{
	    			Double_t vi = img_ediff_cont_1x100->GetBinContent(l, x);
				vi=vi/electron_to_eV;
	    			xloop[0] = x[0];
	    			xloop[1] = x[1];
	    			xloop[2] = x[2];

					double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
					double dc = rand->Poisson(darkcurrent*exposure_day*100);
					img_adcu_cont_1x100->Fill(xloop, min(saturation_adcu, TMath::Nint((vi + noise + dc)*electron_to_eV/(e_factor*1000))));
					flag_noise_cont_1x100->Fill(xloop, 1);	
	  
	  				for (int bix=-15; bix <= 15; bix++) 
	  				{
	  					for (int biy=-1; biy <= 1; biy++) 
	  					{
	  						xloop[0] = x[0] + bix;
	  						xloop[1] = x[1] + biy;
	  						// xloop[2] = x[2];
	  						if (xloop[0] > xmax1x100[0] || xloop[0] < xmin1x100[0] || xloop[1]>xmax1x100[1] || xloop[1]<xmin1x100[1]) 
							{
								continue;
							}
							Long64_t  index = flag_noise_cont_1x100->GetBin(xloop, "false");
							double  flag_pix = flag_noise_cont_1x100->GetBinContent(index);
							if (flag_pix<=0)
							{ //noise not yet added --> Should the DC include diffusion? 
								double noise = rand->Gaus(pedestal_mean,TMath::Sqrt(nreads*pedestal_sigma*pedestal_sigma));
								double dc = rand->Poisson(darkcurrent*exposure_day*100);
								img_adcu_cont_1x100->Fill(xloop, min(saturation_adcu, TMath::Nint((noise + dc)*electron_to_eV/(e_factor*1000))));
								flag_noise_cont_1x100->Fill(xloop, 1);
							}
						}
					}
				}
			}//end binflag
			
			sprintf(title, "Evt: %d ; coord x [pix]; coord y [pix];", (int)EventID_last);
			
			img_edep->SetTitle(title);
			img_zdep->SetTitle(title);
			img_ediff->SetTitle(title);
			img_adcu->SetTitle(title);
			img_zadcu_cont->SetTitle(title);
			img_ediff_cont->SetTitle(title);
			img_adcu_cont->SetTitle(title);

	      	if (binflag > 0)   
	      	{
	      		img_ediff_1x100->SetTitle(title);
	      		img_adcu_1x100->SetTitle(title);
	      		img_ediff_cont_1x100->SetTitle(title);
	      		img_adcu_cont_1x100->SetTitle(title);
	      	}

			//*************************
			//Find Clusters
			//************************
			if (img_edep->GetNbins() > 0) 
			{
				FindClusters(EventID_last, EvtPDG); 
			}

			//*********************************
			//Write out and reset for next evt
			//********************************
			//1x1
			if ((binflag == 0 || binflag == 2)) 
			{ 
				sprintf(grname, "img_edep_%d", (int)EventID_last);
				img_edep->SetName(grname);
				if(debug > 0)
				{
					img_edep->Write();
				}

	  			sprintf(grname, "img_zdep_%d", (int)EventID_last);
	  			img_zdep->SetName(grname);
	  			if(debug > 1 )
				{
					img_zdep->Write();
				}
	  
	  			sprintf(grname, "img_ediff_%d", (int)EventID_last);
	  			img_ediff->SetName(grname);
	  			if (debug > 1)
				{
					img_ediff->Write();
				}
	  
	  
	  			sprintf(grname, "img_adcu_%d", (int)EventID_last);
	  			img_adcu->SetName(grname);
	  			if (debug > 1)
	  			{
	  				img_adcu->Write();
				}
	  
	  			//continuous readout
	  			sprintf(grname, "img_adcu_cont_%d", (int)EventID_last);
	  			img_adcu_cont->SetName(grname);
	  			if (debug > 0)
				{
					img_adcu_cont->Write();
				}
	  
	  			sprintf(grname, "img_ediff_cont_%d", (int)EventID_last);
	  			img_ediff_cont->SetName(grname);
	  			if (debug > 1)
				{
					img_ediff_cont->Write();
				}
			}
			
			//1x100
			if (binflag > 0) 
			{ 
				sprintf(grname, "img_ediff_1x100_%d", (int)EventID_last);
				img_ediff_1x100->SetName(grname);
				if (debug > 0)
				{
					img_ediff_1x100->Write();
				}
      
				sprintf(grname, "img_adcu_1x100_%d", (int)EventID_last);
				img_adcu_1x100->SetName(grname);     
				if (debug > 1)
				{
					img_adcu_1x100->Write();
				}
      
				//continuous readout
				sprintf(grname, "img_ediff_cont_1x100_%d", (int)EventID_last);
				img_ediff_cont_1x100->SetName(grname);
				if (debug > 0)
				{
					img_ediff_cont_1x100->Write();
				}
      
				sprintf(grname, "img_adcu_cont_1x100_%d", (int)EventID_last);
				img_adcu_cont_1x100->SetName(grname);     
				if (debug > 0)
				{
					img_adcu_cont_1x100->Write();
				}
			}
			
			img_edep->Reset();
			img_zdep->Reset();
			img_npart->Reset();
			img_ediff->Reset();

		    img_adcu->Reset();
		    flag_noise->Reset();
		    img_ediff_1x100->Reset();
		    img_adcu_1x100->Reset();
		    flag_noise_1x100->Reset();
		    
		    img_ediff_cont->Reset();
		    img_zdiff_cont->Reset();
		    img_adcu_cont->Reset();
		    img_zadcu_cont->Reset();
		    
		    flag_noise_cont->Reset();
		    img_ediff_cont_1x100->Reset();
		    img_adcu_cont_1x100->Reset();
		    flag_noise_cont_1x100->Reset();
		    
		    EventID_last = EventID;
	} //end if EventID != EventID_last
	
	for (int jj= 0; jj<edep_tot->size(); jj++)  {
   		// Go to the pixel bin given a position (posx,posy) in mm
   		// and assuming a pixel size equal to pix_size microns
    	int binx = TMath::Nint(posx->at(jj)*1000. / pix_size);
    	int biny = TMath::Nint(posy->at(jj)*1000. / pix_size);
	
	  
	  //check pre-existing content 
	  double fValue[3] = {(double)binx, (double)biny, (double)CCDid->at(jj)};
	  vector<double> binxdiff; // vector filled by diff. 
	  vector<double> binydiff; // vector filled by diff. 
	  
	  // position in units of microns

	  DiffuseElectrons(posx->at(jj)*1000., posy->at(jj)*1000., posz->at(jj)*1000, edep_tot->at(jj), CCDid->at(jj), binxdiff, binydiff); //pass coord. in um
	  //    cout << " pos z : " << posz << endl ;
	  img_edep->Fill(fValue, edep_tot->at(jj)); //fill 1 electron per time
	  img_zdep->Fill(fValue, posz->at(jj)); //fill 1 electron per time
	  img_npart->Fill(fValue, 1); //fill 1 electron per time
	
    	//noise/ped/
    	int nelec = binxdiff.size();
    	double time0 = 0;

   	//y-coordinate shift due to time, in time equivalent wrt 0.
    	//int ytime1x100 = floor((time - time0)/readtime_line/100.); //y-coordinate shift due to
    	//time equivalent wrt 0.
    	int ytime = 0.0;
	
    	if(readtime_line>0.0)
		{
		  int ytime = floor((time->at(jj) - time0)/readtime_line); 
		}

    	for (int k=0; k<nelec; ++k) 
    	{
    		fValue[0] = binxdiff.at(k);
    		fValue[1] = binydiff.at(k);
			
	       	//fill 1 electron per time converted in eV
      		img_ediff->Fill(fValue, electron_to_eV);  
      		// img_ediff->Fill(fValue, 1.0);
      		
      		//continuous readout here --> 
    
      		// for (int k=0; k<binxdiff.size(); ++k)
      		// fValue[0] = binxdiff.at(k);
      		//
      		fValue[1] = binydiff.at(k) + ytime; 
     		
     		//fill 1 electron per time converted in eV
      		// XXX img_ediff->Fill(fValue, electron_to_eV);
      		img_ediff_cont->Fill(fValue, electron_to_eV);
      		//img_ediff_cont->Fill(fValue, 1.0 );
      		
      		//fill one per electron. Note that it needs to be
      		//normalized for the numb of electrons in that bin. see later
      		img_zdiff_cont->Fill(fValue, posz->at(jj));
      	}//end loop binx_diff;
      }//end jj
    } //end nentries
	delete rand;
}


//-------------------------------

int main(int argc, char **argv) 
{
	TRandom3 *rand = new TRandom3();
	seed = 0;
	diff_a_flag = diff_b_flag = false;
	format = "root"; //default value
	diffmodel = "simple" ;
	diff_A = 0.14; //default values.. to be checked 
	diff_b = 0.001;   //default values.. to be checked
	diff_offset = 0;
	string filein = "";
	
	if(argc>2)
	{
		Usage();
	}
	else if(argc==1)
	{
		filein = "config.txt";
	}
	else
	{
		filein = argv[1];
	}
	
	//read config file 
	ReadConfigFile(filein);
	if (seed > 0)
	{
		rand->SetSeed(seed);
	}
  
	Double_t adc_unit = electron_to_eV/e_factor;
	clEventID=0;

  	//clusters trees - 1x1 (default) 
  	clustersTrue = new TTree("clustersTrue","clusters tree from Edep image");
  	clustersTrue->Branch("EventID",&clEventID);
  	clustersTrue->Branch("CCDID",&clCCDid);
  	clustersTrue->Branch("ClusterID",&ClusterID);
  	clustersTrue->Branch("PrimaryPart",&PrimaryPart);
  	clustersTrue->Branch("MeanDepth",&TrueMeanDepth);
  	clustersTrue->Branch("RMSDepth",&TrueRMSDepth);
  	clustersTrue->Branch("Thr_adcu",&Thr_adcu);
  	clustersTrue->Branch("Energy",&TrueEnergy);
  	clustersTrue->Branch("PosX",&TruePosX);
  	clustersTrue->Branch("PosY",&TruePosY);
  	clustersTrue->Branch("MeanX",&TrueMeanX);
  	clustersTrue->Branch("MeanY",&TrueMeanY);
  	clustersTrue->Branch("RMSX",&TrueRMSX);
  	clustersTrue->Branch("RMSY",&TrueRMSY);
  	clustersTrue->Branch("Xmin",&TrueXmin);
  	clustersTrue->Branch("YminY",&TrueYmin);
  	clustersTrue->Branch("Xmax",&TrueXmax);
  	clustersTrue->Branch("Ymax",&TrueYmax);
  	clustersTrue->Branch("Qmax",&TrueQmax);
  	clustersTrue->Branch("QmaxX",&TrueQmaxX);
  	clustersTrue->Branch("QmaxY",&TrueQmaxY);
  	clustersTrue->Branch("DX",&TrueDX);
  	clustersTrue->Branch("DY",&TrueDY);
  	clustersTrue->Branch("NPix",&NTruePix);
  	clustersTrue->Branch("Pixels_x",&TruePixels_x);
  	clustersTrue->Branch("Pixels_y",&TruePixels_y);
  	clustersTrue->Branch("Pixels_z",&TruePixels_z);
  	clustersTrue->Branch("Pixels_E",&TruePixels_E);

  	//add info on the particle initiating the cluster? if it's a gamma or e- or something else?
  	////how to define in case of more particles depositing energy in the same pixel?
  	
  	clustersRec = new TTree("clustersRec","clusters tree from image (diff,DC, noise and continuous RO) in ADCu");
  	clustersRec->Branch("EventID",&clEventID);
  	clustersRec->Branch("CCDID",&clCCDid);
  	clustersRec->Branch("ClusterID",&ClusterID);
  	clustersRec->Branch("PrimaryPart",&PrimaryPart);
  	// clustersRec->Branch("Time",&time);
  	clustersRec->Branch("Thr_adcu",&Thr_adcu);
  	clustersRec->Branch("Energy",&RecEnergy);
  	clustersRec->Branch("MeanDepth",&RecMeanDepth);
  	clustersRec->Branch("RMSDepth",&RecRMSDepth);
  	clustersRec->Branch("PosX",&RecPosX);
  	clustersRec->Branch("PosY",&RecPosY);
  	clustersRec->Branch("MeanX",&RecMeanX);
  	clustersRec->Branch("MeanY",&RecMeanY);
  	clustersRec->Branch("RMSY",&RecRMSY);
  	clustersRec->Branch("RMSX",&RecRMSX);  
  	clustersRec->Branch("Xmin",&RecXmin);
  	clustersRec->Branch("Ymin",&RecYmin);
  	clustersRec->Branch("Xmax",&RecXmax);
  	clustersRec->Branch("Ymax",&RecYmax);
  	clustersRec->Branch("Qmax",&RecQmax);
  	clustersRec->Branch("QmaxX",&RecQmaxX);
  	clustersRec->Branch("QmaxY",&RecQmaxY);
  	clustersRec->Branch("DX",&RecDX);
  	clustersRec->Branch("DY",&RecDY);
  	clustersRec->Branch("NPix",&NRecPix);
  	clustersRec->Branch("NSatPixels",&NSatPixels);
  	clustersRec->Branch("Pixels_x",&RecPixels_x);
  	clustersRec->Branch("Pixels_y",&RecPixels_y);
  	//  clustersRec->Branch("Pixels_z",&RecPixels_z); //does not make too much sense as the z per pixels as it cannot be evaluated for noisy pix. 
  	clustersRec->Branch("Pixels_E",&RecPixels_E);

  
  	//=================================
  	//clusters trees - 1x100
  	//==================================
    clustersRec100 = new TTree("clustersRec100","clusters tree from image (diff,DC, noise and continuous RO) in ADCu, 1x100");
    clustersRec100->Branch("EventID",&clEventID);
    clustersRec100->Branch("CCDID",&clCCDid);
    clustersRec100->Branch("ClusterID",&ClusterID);
    clustersRec100->Branch("PrimaryPart",&PrimaryPart);
    clustersRec100->Branch("Thr_adcu",&Thr_adcu);
    clustersRec100->Branch("Energy",&RecEnergy);
    clustersRec100->Branch("MeanDepth",&RecMeanDepth);
    clustersRec100->Branch("RMSDepth",&RecRMSDepth);
    clustersRec100->Branch("PosX",&RecPosX);
    clustersRec100->Branch("PosY",&RecPosY);
    clustersRec100->Branch("MeanX",&RecMeanX);
    clustersRec100->Branch("MeanY",&RecMeanY);
    clustersRec100->Branch("RMSY",&RecRMSY);
    clustersRec100->Branch("RMSX",&RecRMSX);  
    clustersRec100->Branch("Xmin",&RecXmin);
    clustersRec100->Branch("Ymin",&RecYmin);
    clustersRec100->Branch("Xmax",&RecXmax);
    clustersRec100->Branch("Ymax",&RecYmax);
    clustersRec100->Branch("Qmax",&RecQmax);
    clustersRec100->Branch("QmaxX",&RecQmaxX);
    clustersRec100->Branch("QmaxY",&RecQmaxY);
    clustersRec100->Branch("DX",&RecDX);
    clustersRec100->Branch("DY",&RecDY);
    clustersRec100->Branch("NPix",&NRecPix);
    clustersRec100->Branch("NSatPixels",&NSatPixels);
    clustersRec100->Branch("Pixels_x",&RecPixels_x);
    clustersRec100->Branch("Pixels_y",&RecPixels_y);
    clustersRec100->Branch("Pixels_E",&RecPixels_E);

    if(imgmode == "ext") 
    {
    	ReadImage();
    }
	
	string outname = outfile + simfile;
  	fout = new TFile(outname.c_str(),"recreate");
	
	if(simmode == "uni")
	{
		GenerateEdep(); 
	}
	else if(simmode == "geant4")
	{
		ReadGeant4Input();
	}

  	fout->cd();
  	clustersTrue->Write();
  	clustersRec->Write();
  	
  	if(binflag > 0)
  	{ 
  		clustersRec100->Write();
	}
  	fout->Close();
	delete rand;
  	return 1;
}//end main
//-------------------------------
