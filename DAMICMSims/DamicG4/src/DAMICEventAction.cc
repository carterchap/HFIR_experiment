#include "DAMICEventAction.hh"
#include "DAMICRunAction.hh"
#include "DAMICAnalysisManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "DAMICCCDSD.hh"
#include "DAMICCCDHit.hh"
#include "G4SDManager.hh"

DAMICEventAction::DAMICEventAction(DAMICRunAction* runAction) :
    G4UserEventAction(),
    fRunAction(runAction)
{
    ID = 0;
    NamePrimary = "NULL";
    EnergyPrimary = 0;
    AtomicNumber = 0;
    AtomicMass = 0;
}

DAMICEventAction::~DAMICEventAction()
{}


void DAMICEventAction::BeginOfEventAction(const G4Event* Event)
{
  //  G4RunManager* runman = G4RunManager::GetInstance();

  //G4Run * run = runman->GetCurrentRun();
  //run->SetNumberOfEventToBeProcessed(_nentries);

    DAMICAnalysisManager* man = DAMICAnalysisManager::GetInstance();
    //clearing all vectors
    man->clear_vectors();

    G4int EventID = Event->GetEventID();
    SetID(EventID);

    SetNamePrimary(Event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName());
    SetAtomicNumber( Event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicNumber());
    SetAtomicMass(Event->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetAtomicMass());
    SetEnergyPrimary(Event->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy());

    if(EventID%10000 == 0)
    {
        G4cout << EventID << G4endl;
    }

    fRunAction->primary_volid = 0;
    fRunAction->triggerTime = 0;
    CCDCollID = -1;
    S_hits = -1;
    if ( CCDCollID == -1 )
    {
        G4SDManager *SDman = G4SDManager::GetSDMpointer();
        CCDCollID = SDman->GetCollectionID("CCDColl");
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DAMICEventAction::EndOfEventAction(const G4Event* Event)
{

   // Particle generated info
    G4int nvertex =  Event->GetNumberOfPrimaryVertex();
    DAMICAnalysisManager* man = DAMICAnalysisManager::GetInstance();
    G4int EventID = Event->GetEventID();

    for (int i=0; i<nvertex; ++ i) {
      G4ThreeVector momentum = Event->GetPrimaryVertex(i)->GetPrimary()->GetMomentum();
      G4double mass = Event->GetPrimaryVertex(i)->GetPrimary()->GetMass();
      G4double totMomentum = pow(momentum[0]*momentum[0] + momentum[1]*momentum[1] + momentum[2]*momentum[2], 0.5);
      G4double kinEnergy = pow(totMomentum*totMomentum + mass*mass, 0.5 ) - mass;
      G4ThreeVector finalMomentum(0,0,0);

      // XXX FIXME: Substitute this for a ZERO in term of a double
      // (std::fabs(totMomentum) < 1e-20??)
      if(totMomentum!=0)
	{
	  finalMomentum = G4ThreeVector(momentum[0]/totMomentum,momentum[1]/totMomentum,momentum[2]/totMomentum);
	}
         man->FillEventOut(
			Event->GetPrimaryVertex(i)->GetPrimary()->GetPDGcode(),
			Event->GetPrimaryVertex(i)->GetPrimary()->GetCharge(),
			fRunAction->primary_volid, kinEnergy/CLHEP::eV,
			Event->GetPrimaryVertex(i)->GetPosition(),
			finalMomentum,
			fRunAction->triggerTime);
    }


    //=============================
    // address hits collections
    DAMICCCDHitsCollection* SHC = nullptr;
    G4HCofThisEvent* HCE = Event->GetHCofThisEvent();
    if(HCE)
    {
      SHC = (DAMICCCDHitsCollection*)(HCE->GetHC(CCDCollID));
    }

    if(SHC)
    {
        S_hits = SHC->entries();
    }

    G4double energy = 0;
    G4double track = 0;
    G4double time = 0;
    G4double ltime = 0;
    G4double depth_in_ccd = 0;
    G4double posx = 0;
    G4double posy = 0;
    G4double gposx = 0;
    G4double gposy = 0;
    G4double gposz = 0;

    int klast=0;
    int ncount = 0;
    for (int i=0; i<S_hits; ++i)
    {
        G4int ccdID = (*SHC)[i]->GetCCDNum(); //missed [0-nmax_ccd]
        G4int pixID = (*SHC)[i]->GetPixID();  //to be fixed for damicm [0-36Mpix]
        G4int pID   = (*SHC)[i]->GetTrackID();
        G4int pdg   = (*SHC)[i]->GetPDGNumber();
        G4int mID   = (*SHC)[i]->GetMotherID();

        energy = time = ltime = track = ncount = depth_in_ccd = posx = posy = gposx = gposy = gposz = 0;
        for (int k=0; k<S_hits; ++k)
        {
            if ((*SHC)[k]->GetCCDNum() == ccdID && (*SHC)[k]->GetPixID() == pixID && (*SHC)[k]->GetTrackID() == pID)
            {
                energy += (*SHC)[k]->GetEnergyDeposit();
                
                time += (*SHC)[k]->GetTime();
                ltime += (*SHC)[k]->GetLocalTime();

                posx += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetPos().x();
                posy += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetPos().y();
                depth_in_ccd += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetPos().z();

                gposx += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetGlobPos().x();
                gposy += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetGlobPos().y();
                gposz += ((*SHC)[k]->GetEnergyDeposit())*(*SHC)[k]->GetGlobPos().z();

                klast = k;
                ncount++;
            }

        } //end for k

        if(ncount>0)
        {
            time /= (double)ncount;
            ltime /= (double)ncount;
            //time -= fRunAction->triggerTime;
        }

        if(energy>0)
        {
            depth_in_ccd /= energy;
            posx /= energy;
            posy /= energy;

            gposz /= energy;
            gposx /= energy;
            gposy /= energy;
        }
        else
        {
            depth_in_ccd = -1;
            posx = -1;
            posy = -1;

            gposz = -1;
            gposx = -1;
            gposy = -1;
        }
        i = klast;
        //test what happens if a particle deposit part of the energy in one ccd and part in another?
        // XXX test what happens in case of backscattering? is klast the correct value?

        if(energy>0)
        {
            man->FillCCDOut(pdg,pID,mID,ccdID,
                    posx/CLHEP::mm,posy/CLHEP::mm,depth_in_ccd/CLHEP::mm,
                    gposx/CLHEP::mm,gposy/CLHEP::mm,gposz/CLHEP::mm,
                    energy/CLHEP::eV,time/CLHEP::s,ltime/CLHEP::s);
        }
    }

    // Fill all trees
    man->FillTrees(EventID);
}



void DAMICEventAction::SetID(G4int id)
{
    ID = id;
}

void DAMICEventAction::SetNamePrimary(G4String name)
{
    NamePrimary = name;
}

void DAMICEventAction::SetEnergyPrimary(G4double nrj)
{
    EnergyPrimary = nrj;
}
void DAMICEventAction:: SetAtomicNumber ( G4int Z)
{
    AtomicNumber = Z;
}

void DAMICEventAction::SetAtomicMass(G4int A)
{
    AtomicMass = A;
}
