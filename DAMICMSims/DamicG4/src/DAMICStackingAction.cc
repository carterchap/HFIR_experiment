#include "DAMICStackingAction.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <vector>
#include "G4DecayTable.hh"
#include "G4RunManager.hh"

DAMICStackingAction::DAMICStackingAction():G4UserStackingAction(),
    _partName("ez"), 
    _partEnergyKin(0.), 
    _partEnergyTot(0.)
{
    _particlesKill = {};
    _limitNucleus = -1;
    _messenger = new DAMICStackingActionMessenger(this);
    _idpartPrimPartPDG = -1;
}

DAMICStackingAction::~DAMICStackingAction()
{ 
    if (_messenger != nullptr)
    {
        delete _messenger;
        _messenger = nullptr;
    }
}

void DAMICStackingAction::SetPartName(G4String name)
{
    _partName = name;
}

void DAMICStackingAction::SetPartEnergyKin(G4double nrj)
{
    _partEnergyKin = nrj;
}

void DAMICStackingAction::SetPartEnergyTot(G4double nrj)
{
    _partEnergyTot = nrj;
}

void DAMICStackingAction::SetLimitNucleus(G4int pdgenco)
{
    _limitNucleus = pdgenco;
}

void DAMICStackingAction::AddParticlesKill(G4String Particle)
{
    G4ParticleDefinition* Part = G4ParticleTable::GetParticleTable()->FindParticle(Particle);
    _particlesKill.push_back(Part->GetPDGEncoding());

}

G4bool DAMICStackingAction::DoKill(G4int PartPDG)
{
    G4bool Kill = false;
    G4int entries = _particlesKill.size();
    for (int i = 0; i< entries; i++)
    {
        if(PartPDG == _particlesKill[i])
        {
            Kill = true;
        }
    }
    return Kill;
}


G4bool DAMICStackingAction::DoKillNucleus(G4int PartPDG)
{
    if(_limitNucleus == -1)
    {
        return false;
    }
    if(PartPDG < 1000000000)
    {
        return false;
    }
    if(PartPDG == 1000020040)
    {
        return false;
    }
    
    //G4cout << PartPDG << G4endl;
    int Z = floor(PartPDG/10000)- 100000;
    int A = floor(PartPDG/10)-100000000 - Z*1000 ;
    int Zbase = floor(_limitNucleus/10000)- 100000;
    int Abase = floor(_limitNucleus/10)- 100000000 - Zbase*1000;
    int unity = PartPDG - 1000000000 - Z *10000 - A *10;
    //G4bool alphaC = A >= Abase+4;
    //G4bool firstC = A < Abase;
    //G4bool secondC = Z < Zbase-1;
    //G4bool thirdC = Z > Zbase +1;
    G4bool sameZ = Z ==Zbase;
    G4bool sameA = A == Abase;
    G4bool nulUnity = unity ==0 ;
    /*if(alphaC)
    {
        return false;
    }
    if(firstC)
    {
        return true;
    }
    else if(secondC || thirdC)
    {
        return true;
    }*/
    if(sameZ && sameA && nulUnity)
    {
        return true;
    }
    else
    {
        return false;
    }
}

G4bool DAMICStackingAction::DoOneDecay(G4int PartPDG, G4int IDpart)
{ // construction

    //G4cout << " -------" << PartPDG << "," << _idpartPrimPartPDG << G4endl;   
    if(PartPDG < 1000000000)
    {
        return false;
    }
    if(PartPDG == 1000020040)
    {
        return false;
    }
    if(IDpart == 0)
    {
        return false;
    }
    G4int diffID = abs(_idpartPrimPartPDG/10 - PartPDG/10);
    if(diffID == 0)
    {
        return false;
    }
    if((PartPDG % 10) == 0)
    {
        return true;
    }
    return false;
}

void DAMICStackingAction::ResetNewRun()
{
    _particlesKill = {};
    _limitNucleus = -1;
}

void DAMICStackingAction::UnsetDoOneDecay(G4bool val)
{
    // if val (/damicm/stack/fullchain is true: UnsetDoOneChain
    //  if not, only DoOneDecay
    _oneDecay = not val;
}

G4ClassificationOfNewTrack DAMICStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
    G4String Name = aTrack->GetDefinition()->GetParticleName();
    SetPartName(Name);
    //G4cout << "je rentre Stack" << G4endl;
    G4double energy = aTrack->GetTotalEnergy();
    G4double energykin = aTrack->GetKineticEnergy();
    SetPartEnergyKin(energykin);
    SetPartEnergyTot(energy);
    G4int IDpart = aTrack->GetParentID();
    G4int ProcessCreatorSub = -1;
    G4int ProcessCreatorType = -1;
    G4String ProcessName = "NULL";
    G4String nameVolume = "NULL";

    /* Test */
    G4int PDGEnco = aTrack->GetDefinition()->GetPDGEncoding();
    /*if(PDGEnco == 1000912340)
    {
        G4cout << "  decay at 234Pa !!!!!!!!! " << G4endl;
    }*/

    if(IDpart == 0)
    {
        _idpartPrimPartPDG = PDGEnco;
    }
    /*if(IDpart !=0 && aTrack->GetCreatorProcess()->GetProcessName() != "ionIoni" 
         && aTrack->GetCreatorProcess()->GetProcessName() != "eIoni")
    {
        G4cout << Name << G4endl;
        G4cout << energykin  << "   " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
    }*/
    //G4cout << Name << G4endl;
    //G4cout << aTrack->GetDefinition()->GetPDGEncoding() << G4endl;
    

    if(DoKill(PDGEnco))
    {
        return fKill;
    }
    if(DoKillNucleus(PDGEnco))
    {
        return fKill;
    }
    //if((Name !="e-" ) && (Name!="gamma"))
    //      G4cout << PDGEnco << G4endl;
    
    if(_oneDecay)
    {
        if(DoOneDecay(PDGEnco, IDpart))
        {
            return fKill;
        }
    }

    //G4cout << Name << G4endl;
    if(IDpart != 0)
    {
        ProcessCreatorSub = aTrack->GetCreatorProcess()->GetProcessSubType();
        ProcessCreatorType = aTrack->GetCreatorProcess()->GetProcessType();
        //ProcessName= aTrack->GetCreatorProcess()->GetProcessName();
    }
    else 
    {
        ProcessCreatorSub = -1;
        ProcessCreatorType = -1;
        nameVolume = "NULL";
    }

    // Kill electrons produced by ioni below 10 keV
    /*if(ProcessCreatorSub == 2 && ProcessCreatorType == 2 && energykin < 10*keV){
      return fKill;
    }*/

    // Check volume from particles not produced by ioni.
    if(IDpart != 0 && ProcessCreatorSub != 2 && ProcessCreatorType != 2 )
    {
        nameVolume = aTrack->GetVolume()->GetName();
    }
    // Kill Brem in CCD
    /*if( ProcessCreatorSub == 3  && energykin < 1*keV && nameVolume == "SensPartPV" && ProcessCreatorType == 2)
    {
      return fKill;
    }
    // Kill Brem not in CCD
    if( ProcessCreatorSub == 3  && energykin < 1*keV && nameVolume != "SensPartPV" && ProcessCreatorType == 2)
    {
      return fKill;
    }*/

    /*if(energykin<1*eV && Name=="e-"){ // Kill particles with less than 1 eV  kinetic energy
        return fKill;
    }*/

    /*if(Name =="gamma"){
      G4cout << energykin  << "   " << ProcessName << G4endl;
    }*/
    // if((Name !="e-" ) && (Name!="gamma"))
    //   G4cout << Name << G4endl;


    //G4cout << "je sors Stack" << G4endl;
    //G4RunManager::GetRunManager()->rndmSaveThisEvent();
    return fUrgent;
}
