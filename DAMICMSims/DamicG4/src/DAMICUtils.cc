/**
 * @file DAMICUtils.cc
 * @author: Mariangela Settimo
 * @date 2018 DAMIC - Subatech  
 */

#include "DAMICUtils.hh"

#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


// Method to check if given volume exist in the geometry
G4bool FindVolume(G4String VolName)
{
  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = 0;
  G4String theRequiredVolumeName = VolName;
  PVStore      = G4PhysicalVolumeStore::GetInstance();
  G4int      i = 0;
  G4bool found = false;
  while (!found && i<G4int(PVStore->size())) {
    tempPV = (*PVStore)[i];
    found  = tempPV->GetName() == theRequiredVolumeName;
    if (!found)
      {i++;}
  }
  // found = true then the volume exists else it doesnt.
  if(found == true)
    {
      return true;
    }
  else
    {
      return false;
    }

}

// Method to get the process Id
G4int GetProcessID(G4String procname)
{
  if (procname == "Generator") return 0;
  //neutron/anti-neutron
  else if (procname == "nCapture") return 1;
  else if (procname == "nFission") return 2;
  else if (procname == "neutronInelastic") return 3;
  else if (procname == "anti_neutronInelastic") return 4;  

  //e+/e-
  else if (procname == "eIoni") return 10;
  else if (procname == "eBrem") return 11;
  else if (procname == "annihil") return 12;
  else if (procname == "electronNuclear") return 13;
  else if (procname == "positronNuclear") return 14;
  
  //phot
  else if (procname == "conv") return 20;
  else if (procname == "phot") return 21;
  else if (procname == "compt") return 22;

  //muons  
  else if (procname == "muIoni") return 30;
  else if (procname == "muPairProd") return 31;
  else if (procname == "muonNuclear") return 32;
  else if (procname == "muBrems") return 33;
  else if (procname == "muMinusCaptureAtRest") return 34;
  
  //hadrons 
  else if (procname == "hIoni") return 40;
  else if (procname == "ionIoni") return 41;
  else if (procname == "hadElastic") return 42;
  else if (procname == "hBertiniCaptureAtRest") return 43;
  else if (procname == "hFritiofCaptureAtRest") return 44;
  else if (procname == "hBertiniCaptureAtRest") return 45;  
  else if (procname == "hBrems") return 46;
  else if (procname == "hPairProd") return 47;
  
  //inelastic
  else if (procname == "protonInelastic") return 50;
  else if (procname == "anti_protonInelastic") return 51;
  
  else if (procname == "pi-Inelastic") return 52;
  else if (procname == "pi+Inelastic") return 53;

  else if (procname == "kaon0SInelastic") return 54;
  else if (procname == "kaon0LInelastic") return 55;

  else if (procname == "kaon+Inelastic") return 56;
  else if (procname == "kaon-Inelastic") return 57;

  else if (procname == "tInelastic") return 58; 
  else if (procname == "anti_tritonInelastic") return 59;
  
  else if (procname == "dInelastic") return 60;
  else if (procname == "anti_deuteronInelastic") return 61;
  
  else if (procname == "sigma+Inelastic") return 62;
  else if (procname == "anti_sigma+Inelastic") return 63;

  else if (procname == "sigma-Inelastic") return 64;
  else if (procname == "anti_sigma-Inelastic") return 65;

  else if (procname == "xi-Inelastic") return 66;
  else if (procname == "anti_xi-Inelastic") return 67;

  else if (procname == "xi0Inelastic") return 68;
  else if (procname == "anti_xi0Inelastic") return 69;

  else if (procname == "omega-Inelastic") return 70;
  else if (procname == "anti_omega-Inelastic") return 71;  

  else if (procname == "lambdaInelastic") return 72;
  else if (procname == "anti-lambdaInelastic") return 73;  

  else if (procname == "alphaInelastic") return 74;
  else if (procname == "anti_alphaInelastic") return 75;
  
  else if (procname == "He3Inelastic") return 76;
  else if (procname == "anti_He3Inelastic") return 77;

  else if (procname == "ionInelastic") return 78;

  //others
  else if (procname == "CoulombScat") return 80;
  else if (procname == "msc") return 81;
  
  //radioactive 
  else if (procname == "RadioactiveDecay") return 90;
  else if (procname == "Decay") return 91;

  else if (procname == "Transportation") return 100;
  else if (procname == "unknown") return 900;
  else return 1000;
                 
}

int GetVolumeID(G4String volName)  {
  
  if (volName== "WorldLV"  ) return 1 ; 
  else if (volName ==   "physLab"  ) return     2;
   /// Base 
   else if (volName ==  "LowerTubeAll")   return    10;
   else if (volName ==  "LowerTube")      return    15;
   else if (volName ==  "PrinTubLV")      return    20;
   else if (volName ==  "ColdFingerLV")   return    20;
   else if (volName ==  "TopPlate" )      return    21;
   else if (volName ==  "RearPlate")      return    22;
   else if (volName ==  "endCoverPlate" ) return    23;
   else if (volName ==  "SidePlate")      return    24;
   else if (volName ==  "InsulatorLV"   ) return    25;
   else if (volName ==  "endCoverBar"   ) return    26;
   else if (volName ==  "UpperPlateLV"  ) return    27;
   else if (volName ==  "SpacerFinalLV" ) return    28;
  //modules 
   else if (volName ==  "Module44LV"    )  return   40;
   else if (volName ==  "CopperBasePlate") return   41;
   else if (volName ==  "CopperMountBar")  return   42;
   else if (volName ==  "CopperTopPlate" ) return   43;
   else if (volName ==  "Sub44LV"  )       return   43;
   else if (volName ==  "CCDSensor"  )     return   44;
  //CCDsens
   else if (volName ==  "TopLayer"  )      return   45;
   else if (volName ==  "BotLayer"  )      return   46;
   else if (volName ==  "Sens" )           return   47;
   else if (volName ==  "SiliconSubstrate" ) return     48; 
   else if (volName ==  "Module44WSiLV"  )   return     49;
   else if (volName ==  "Sub44LVWSi")        return     50;
  //Flex
   else if (volName ==  "FlexCable")      return    60;
   else if (volName ==  "Flex2Cable"    ) return    61;
   else if (volName ==  "Flex3Cable"    ) return    62;
   else if (volName ==  "Flex4Cable"    ) return    63;
  //Lead
   else if (volName ==  "LeadBlockInBox") return    80;
   else if (volName ==  "BlockSupport"  ) return    81;
   else if (volName ==  "CopperOnLeadB" ) return    82;
   else if (volName ==  "SpacerPlateLV" ) return    83;
   else if (volName ==  "GoodLead" )      return    84;
   else if (volName ==  "SpacerPlateLV" ) return    85;
   else if (volName ==  "AncientLead2"  ) return    86;
   else if (volName ==  "BoxMountingPlateLV" ) return   87;
   else if (volName ==  "BoxTopLead"    )      return   88;
  /// Full Plane
   else if (volName ==  "BoxBottomLead" ) return    90;
   else if (volName ==  "InsideLead1"   ) return    91;
   else if (volName ==  "InsideLead2"   ) return    92;
   else if (volName ==  "ShieldingLead" ) return    93;
   ///ShieldingLeadPart
   else if (volName ==  "ShieldingLeadPart"  ) return   94;
   else if (volName ==  "LeadCastleFrame"    ) return   95;
   else if (volName ==  "PartLeadCastle")      return   96; 
  //Flanges
   else if (volName ==  "LowerEnd" )      return    100;
   else if (volName ==  "LowerFlangeAll") return    101;    
   else if (volName ==  "UpperFlange"   ) return    102;
   else if (volName ==  "Assembly1")      return    103;
   else if (volName ==  "Assembly3")      return    104;
  //Assembly 1  - same LV name... ?
   else if (volName ==  "RestraintSheet")   return  105;
   else if (volName ==  "RestraintBlocks" ) return  106;
   else if (volName ==  "CornerLead"    )   return  107; 
   else return -999;
}

int GetVolID(G4String volName)
{

    std::map<G4String, int> volNamesID = {
        {   "WorldLV"                   ,   1   }   ,
        {   "physLab"                   ,   2   }   ,
        /// Base 
        {   "LowerTubeAll"              ,   10  }   ,
        {   "LowerTube"                 ,   15  }   ,
        {   "PrinTubLV"                 ,   20  }   ,
    {   "ColdFingerLV"              ,   20  }   ,
    {   "TopPlate"                  ,   21  }   ,
    {   "RearPlate"                 ,   22  }   ,
    {   "endCoverPlate"             ,   23  }   ,
    {   "SidePlate"                 ,   24  }   ,
    {   "InsulatorLV"               ,   25  }   ,
    {   "endCoverBar"               ,   26  }   ,
    {   "UpperPlateLV"              ,   27  }   ,
    {   "SpacerFinalLV"             ,   28  }   ,
    //modules 
    {   "Module44LV"                ,   40  }   ,
    {   "CopperBasePlate"           ,   41  }   ,
    {   "CopperMountBar"            ,   42  }   ,
    {   "CopperTopPlate"             ,  43  }   ,
        {   "Sub44LV"                   ,   43  }   ,
    {   "CCDSensor"                   , 44  }   ,
    //CCDsens
    {   "TopLayer"                   ,  45  }   ,
    {   "BotLayer"                   ,  46  }   ,
    {   "Sens"                       ,  47  }   ,
    {   "SiliconSubstrate"           ,  48  }   ,   
    {   "Module44WSiLV"              ,  49  }   ,
    {   "Sub44LVWSi"                 ,  50  }   ,
        //Flex
    {   "FlexCable"                 ,   60  }   ,
    {   "Flex2Cable"                ,   61  }   ,
    {   "Flex3Cable"                ,   62  }   ,
    {   "Flex4Cable"                ,   63  }   ,
    //Lead
    {   "LeadBlockInBox"            ,   80  }   ,
    {   "BlockSupport"              ,   81  }   ,
    {   "CopperOnLeadB"             ,   82  }   ,
    {   "SpacerPlateLV"             ,   83  }   ,
    {   "GoodLead"                  ,   84  }   ,
    {   "SpacerPlateLV"             ,   85  }   ,
    {   "AncientLead2"              ,   86  }   ,
    {   "BoxMountingPlateLV"        ,   87  }   ,
    {   "BoxTopLead"                ,   88  }   ,
       /// Full Plane
        {   "BoxBottomLead"             ,   90  }   ,
        {   "InsideLead1"               ,   91  }   ,
        {   "InsideLead2"               ,   92  }   ,
        {   "ShieldingLead"             ,   93  }   ,
        ///ShieldingLeadPart
        {   "ShieldingLeadPart"         ,   94  }   ,
        {   "LeadCastleFrame"           ,   95  }   ,
        {   "PartLeadCastle"            ,   96  }   ,   
    //Flanges
    {   "LowerEnd"                  ,   100 }   ,
    {   "LowerFlangeAll"            ,   101 }   ,   
        {   "UpperFlange"               ,   102 }   ,
        {   "Assembly1"                 ,   103 }   ,
    {   "Assembly3"                 ,   104 }   ,
    //Assembly 1  - same LV name... ?
        {   "RestraintSheet"            ,   105 }   ,
        {   "RestraintBlocks"           ,   106 }   ,
        {   "CornerLead"                ,   107 }   ,
  
    };

    if( !volNamesID[volName] ) 
      return -999;
    else
      return volNamesID[volName];
}


