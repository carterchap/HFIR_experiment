#ifndef DAMICCCDSD_h
#define DAMICCCDSD_h 1

#include "DAMICCCDHit.hh"
#include "DAMICUtils.hh"

#include "G4VSensitiveDetector.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class DAMICRunAction;

class DAMICCCDSD : public G4VSensitiveDetector
{
    public:
        DAMICCCDSD(G4String name);
        virtual ~DAMICCCDSD();

        virtual void Initialize(G4HCofThisEvent*HCE);
        virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);

        G4int GetCCDNumber(G4double globZ);
    private:
        DAMICCCDHitsCollection* _hitsCollection;
        G4int _hitID;
        G4int _eventID;
        const G4double _pixSize=15e-3;

        void _setEventID(G4int evtID) {_eventID = evtID;};
        G4int _getEventID() {return _eventID;};


};


#endif
