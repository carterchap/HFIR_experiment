#include <fstream>
#include <iostream>
#include "TMap.h"
#include "TObjString.h"
#include "TTree.h"
#include "TList.h"
#include "TParameter.h"
#include "TNamed.h"

TMap* config = new TMap();
TList* volumesUsed = new TList();



void readFileConfig(TString fFile){

  string rawline;
  ifstream in(fFile);
  Bool_t involumes = false;

  while(getline(in,rawline)){

    TString line = rawline;
    if (line.BeginsWith("#Volumes")) involumes = true;
    if(line.BeginsWith("#")) continue;

    TString cname = "";
    TString cval = "";
    Double_t cvalD = 0;


    if(line.Length()==0) continue;
    while(line.First(" ")==0) line.Remove(0,1);
    cname=line(0,line.First(" "));
    line.Remove(0,cname.Length()+1);
    while(line.First(" ")==0) line.Remove(0,1);
    cval=line;

    if(!config->Contains(cname)){

      TObjString* oname = new TObjString(cname);
      TObjString* oval = new TObjString(cval);
      cvalD = cval.Atof();
      cout << involumes << endl;
      if (cvalD != 0 && involumes){
        volumesUsed->Add(oname);
      }
      config->Add(oname, oval);
    }

  }

  in.close();
}

Double_t ConfigValue(TString cname){

  if(config->Contains(cname)){

    TString s = ((TObjString*) config->GetValue(cname))->GetString();
    if(s=="true") return 1;
    else if(s=="false") return 0;
    else return s.Atof();
  }

  else{

    std::cout << cname << "is not a config parameter!" << std::endl;
    return 0;
  }
}



void createMacroContaminent(TString fName, TString fConfig ){
  TString base = "/damic/gun/";
  TString energyM = "energy/";
  TString directionM = "direction/";
  TString volumeM = "position/";

  ofstream macroFile(fName);
  readFileConfig("contaminant.cfg");

  if (macroFile){
    // Mandatory
    macroFile << "/control/verbose 0" << endl;
    macroFile << "/run/verbose 0" << endl;
    macroFile << "/run/initialize" << endl;
    macroFile << "/tracking/verbose 0" << endl;
    macroFile << "/event/verbose 0" << endl;

    // Particle
    TString ParticleCommand = base+ "particle ion";
    macroFile  << ParticleCommand << endl;

    Double_t ZConta = ConfigValue("Z");
    Double_t AConta = ConfigValue ("A");
    TString IonCmd = base + "ion ";

    macroFile << IonCmd << ZConta << " " << AConta << endl;

    // Energy
    TString energyCmd = base + energyM;
    energyCmd += "mono 0 eV";

    macroFile << energyCmd << endl;

    // Direction
    TString direction = base + directionM;
    TString directionOne = direction+"oned";
    TString directionX = direction+"onedX 1";
    TString directionY = direction+"onedY 0";
    TString directionZ = direction+"onedZ 0";

    macroFile << directionOne << endl;
    macroFile << directionX << endl;
    macroFile << directionY << endl;
    macroFile << directionZ << endl;

    // Volume
    cout << volumesUsed->GetSize() << endl;
    TString doVol = base + volumeM + "dovolume";
    macroFile << doVol << endl;
    for (int i =0; i < volumesUsed->GetSize(); i++){
      TString nameVol = volumesUsed->At(i)->GetName();
      Double_t activity = ConfigValue(nameVol);
      TString Add = "addvolume ";
      TString Command  = base + volumeM + Add + nameVol+ " ";
      macroFile << Command << activity << endl;
    }


    // Contaminant Number
    TString run = "/run/beamOn ";
    Double_t numberRun = ConfigValue("run_number");
    macroFile << run << numberRun << endl;

    macroFile.close();
  }

  else{
    cout << " Error opening macro file " << endl;
  }

}
