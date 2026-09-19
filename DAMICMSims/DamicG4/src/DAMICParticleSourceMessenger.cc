#include <fstream>
#include <iomanip>
#include <iostream>

#include "DAMICParticleSourceMessenger.hh"
#include "DAMICParticleSource.hh"

#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"
#include "G4AnalysisUtilities.hh"
#include <cmath>

DAMICParticleSourceMessenger::DAMICParticleSourceMessenger(DAMICParticleSource* fPartSource):
    fParticleGun(fPartSource),
    fShootIon(false)
{
    particleTable = G4ParticleTable::GetParticleTable();
    
    gunDirectory = new G4UIdirectory("/damic/gun/");
    gunDirectory->SetGuidance("Particle Source control commands.");
    
    ListCmd = new G4UIcmdWithoutParameter("/damic/gun/list", this);
    ListCmd->SetGuidance("List available particles.");
    ListCmd->SetGuidance(" Invoke G4ParticleTable.");

    EnergyDirectory = new G4UIdirectory("/damic/gun/energy/");
    EnergyDirectory->SetGuidance("Energy control commands.");
    
    DirectionDirectory = new G4UIdirectory("/damic/gun/direction/");
    DirectionDirectory->SetGuidance("Direction control commands.");
    
    PositionDirectory = new G4UIdirectory("/damic/gun/position/");
    PositionDirectory->SetGuidance("Position control commands.");

   
    // READER
    ShieldReaderCmd = new G4UIcmdWithAString("/damic/gun/file", this);
    ShieldReaderCmd->SetGuidance("Read info from a previous output file. Tree ShieldOut is expected to be in the file.");
    ShieldReaderCmd->SetParameterName("ShieldFileName", true, true);
    ShieldReaderCmd->SetDefaultValue("testShielding.root");
    
    // COSMOGENIC
   
    CosmoDirectory = new G4UIdirectory("/damic/gun/cosmogenic/");
    CosmoDirectory->SetGuidance("Cosmogenic particles control commands. Energy distribution as cosmogenic one");

    DoCosmoPartCmd = new G4UIcmdWithAString("/damic/gun/cosmogenic/particleType", this);
    DoCosmoPartCmd->SetGuidance("Set a cosmogenic particle.");
    DoCosmoPartCmd->SetParameterName("CosmoPart", true, true);
    DoCosmoPartCmd->SetDefaultValue("muon");
    DoCosmoPartCmd->SetCandidates("muon neutron");

    DoCosmoAngCmd = new G4UIcmdWithAString("/damic/gun/cosmogenic/angularDist", this);
    DoCosmoAngCmd->SetGuidance("Set angular distribution for cosmogenic particles. ");
    DoCosmoAngCmd->SetParameterName("CosmoAngDist", true, true);
    DoCosmoAngCmd->SetDefaultValue("cosmo");
    DoCosmoAngCmd->SetCandidates("iso cosmo direction");

    // source position type  
    DoCosmoPosCmd = new G4UIcmdWithAString("/damic/gun/cosmogenic/position", this);
    DoCosmoPosCmd->SetGuidance("Set position for cosmogenic particles.");
    DoCosmoPosCmd->SetParameterName("CosmoPos", true, true);
    DoCosmoPosCmd->SetDefaultValue("cuboid");
    DoCosmoPosCmd->SetCandidates("cuboid"); //add here also wall and disk? 
  
    // source shape - if plane is selected
    DoCosmoPosCentreCmd = new G4UIcmdWith3VectorAndUnit("/damic/gun/cosmogenic/srcCentre",this);
    DoCosmoPosCentreCmd->SetGuidance("Sets source shape centre.");
    DoCosmoPosCentreCmd->SetParameterName("X","Y","Z",true);
    DoCosmoPosCentreCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    DoCosmoPosCentreCmd->SetDefaultUnit("mm");
    DoCosmoPosCentreCmd->SetUnitCandidates("nm um mm cm m km");

    // size of source
    DoCosmoPosSizeCmd = new G4UIcmdWith3VectorAndUnit("/damic/gun/cosmogenic/sizeXYZ",this);
    DoCosmoPosSizeCmd->SetGuidance("Set z (half) thickness of source.");
    DoCosmoPosSizeCmd->SetParameterName("deltaX","deltaY","deltaZ",true,true);
    DoCosmoPosSizeCmd->SetDefaultUnit("mm");
    DoCosmoPosSizeCmd->SetUnitCandidates("nm um mm cm m km"); 

    
    // ENERGY
    DoMonoNRJCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/energy/mono", this);
    DoMonoNRJCmd->SetGuidance("Set a monoenergetic beam.");
    DoMonoNRJCmd->SetParameterName("Energy", true, true);
    DoMonoNRJCmd->SetDefaultUnit("keV");
    DoMonoNRJCmd->SetUnitCategory("Energy");
    DoMonoNRJCmd->SetUnitCandidates("eV keV MeV GeV TeV");
    
    DoDistriNRJCmd = new G4UIcmdWithAString("/damic/gun/energy/distri", this);
    DoDistriNRJCmd->SetGuidance("Set a distribution energy.");
    DoDistriNRJCmd->SetParameterName("DisType", true, true);
    DoDistriNRJCmd->SetDefaultValue("Uniform");
    DoDistriNRJCmd->SetCandidates("Uniform Tritium Histo Neutron_UTh_Rock Neutron_Cosmo_Rock Neutron_Lead Max");
    
    HistoNRJCmd = new G4UIcmdWithAString("/damic/gun/energy/histo",this);
    HistoNRJCmd->SetGuidance("Put an Histogram");
    HistoNRJCmd->SetDefaultValue("");
    
    UniformBotCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/energy/bot", this);
    UniformBotCmd->SetGuidance("Bottom energy value for uniform and Max.");
    UniformBotCmd->SetParameterName("Energy", true, true);
    UniformBotCmd->SetDefaultUnit("keV");
    UniformBotCmd->SetUnitCategory("Energy");
    UniformBotCmd->SetUnitCandidates("eV keV MeV GeV TeV");
    
    UniformTopCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/energy/top", this);
    UniformTopCmd->SetGuidance("Top energy value for uniform and Max.");
    UniformTopCmd->SetParameterName("Energy", true, true);
    UniformTopCmd->SetDefaultUnit("keV");
    UniformTopCmd->SetUnitCategory("Energy");
    UniformTopCmd->SetUnitCandidates("eV keV MeV GeV TeV");
    
    AvgMaxCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/energy/AvgMax", this);
    AvgMaxCmd->SetGuidance("Avg energy value for Max.");
    AvgMaxCmd->SetParameterName("Energy", true, true);
    AvgMaxCmd->SetDefaultUnit("keV");
    AvgMaxCmd->SetUnitCategory("Energy");
    AvgMaxCmd->SetUnitCandidates("eV keV MeV GeV TeV");
    
    // DIRECTION    
    DoOneDCmd = new G4UIcmdWithoutParameter("/damic/gun/direction/oned", this);
    DoOneDCmd->SetGuidance("Set Direction to one direction.");
    
    DoOneDXCmd = new G4UIcmdWithADouble("/damic/gun/direction/onedX", this );
    DoOneDXCmd->SetGuidance("Set X of the direction.");
    DoOneDXCmd->SetParameterName("X",true,true);
    
    DoOneDYCmd = new G4UIcmdWithADouble("/damic/gun/direction/onedY", this);
    DoOneDYCmd->SetGuidance("Set Y of the direction.");
    DoOneDYCmd->SetParameterName("Y",true,true);
    
    DoOneDZCmd = new G4UIcmdWithADouble("/damic/gun/direction/onedZ", this);
    DoOneDZCmd->SetGuidance("Set Z of the direction.");
    DoOneDZCmd->SetParameterName("Z",true,true);
    
    DoDistriDCmd = new G4UIcmdWithAString("/damic/gun/direction/distri", this);
    DoDistriDCmd->SetGuidance("Set Direction to a distribution.");
    DoDistriDCmd->SetParameterName("DisType", true, true);
    DoDistriDCmd->SetDefaultValue("Isotropic");
    DoDistriDCmd->SetCandidates("Isotropic Meietal LSM");
    
    IsoThetaMinCmd = new G4UIcmdWithADouble("/damic/gun/direction/thetamin", this);
    IsoThetaMinCmd->SetGuidance("Theta minimal value.");
    IsoThetaMinCmd->SetParameterName("Angle", true, true);
    IsoThetaMinCmd->SetRange("Angle <= 180 && Angle>= 0");
    
    IsoThetaMaxCmd = new G4UIcmdWithADouble("/damic/gun/direction/thetamax", this);
    IsoThetaMaxCmd->SetGuidance("Theta maximal value.");
    IsoThetaMaxCmd->SetParameterName("Angle", true, true);
    IsoThetaMaxCmd->SetRange("Angle <= 180 && Angle>= 0");
    
    IsoPhiMinCmd = new G4UIcmdWithADouble("/damic/gun/direction/phimin", this);
    IsoPhiMinCmd->SetGuidance("Phi minimal value.");
    IsoPhiMinCmd->SetParameterName("Angle", true, true);
    IsoPhiMinCmd->SetRange("Angle <= 360 && Angle>= -180");
    
    IsoPhiMaxCmd = new G4UIcmdWithADouble("/damic/gun/direction/phimax", this);
    IsoPhiMaxCmd->SetGuidance("Phi maximal value.");
    IsoPhiMaxCmd->SetParameterName("Angle", true, true);
    IsoPhiMaxCmd->SetRange("Angle <= 360 && Angle>= -180");
    
    
    // Position
    DoMaterialCmd = new G4UIcmdWithAString("/damic/gun/position/domaterial", this);
    DoMaterialCmd->SetGuidance(" Set position to material defined.");
    DoMaterialCmd->SetParameterName("Material", true, true);
    DoMaterialCmd->SetDefaultValue("G4_Cu");
    
    MotherVolumeCmd = new G4UIcmdWithAString("/damic/gun/position/mother", this);
    MotherVolumeCmd->SetGuidance("Set Mother Volume of the Material.");
    MotherVolumeCmd->SetParameterName("Mother", true, true);
    MotherVolumeCmd->SetDefaultValue("WorldLV");
    
    DoVolumeCmd = new G4UIcmdWithoutParameter("/damic/gun/position/dovolume", this);
    DoVolumeCmd->SetGuidance("Set position to volume defined.");
    
    DoSurfaceCmd = new G4UIcmdWithoutParameter("/damic/gun/position/dosurface", this);
    DoSurfaceCmd->SetGuidance("Set position to surface defined");
    
    MinEmbedDistCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/setminembeddist", this);
    MinEmbedDistCmd->SetGuidance("Set minimum depth to embed particles generated on the surface.");
    MinEmbedDistCmd->SetParameterName("Distance", true, true);
    MinEmbedDistCmd->SetUnitCategory("Length");
    MinEmbedDistCmd->SetDefaultUnit("nm");
    MinEmbedDistCmd->SetUnitCandidates("nm um mm cm m");
    
    MaxEmbedDistCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/setmaxembeddist", this);
    MaxEmbedDistCmd->SetGuidance("Set maximum depth to embed particles generated on the surface.");
    MaxEmbedDistCmd->SetParameterName("Distance", true, true);
    MaxEmbedDistCmd->SetUnitCategory("Length");
    MaxEmbedDistCmd->SetDefaultUnit("nm");
    MaxEmbedDistCmd->SetUnitCandidates("nm um mm cm m");
    
    SetCCDSurfaceOnlyCmd = new G4UIcmdWithABool("/damic/gun/position/setccdfacing", this);
    SetCCDSurfaceOnlyCmd->SetGuidance("Sets whether to do only CCD facing surfaces or not. Meant mostly for simulations in the lead.");
    SetCCDSurfaceOnlyCmd->SetParameterName("Facing", true);
    SetCCDSurfaceOnlyCmd->SetDefaultValue(true);
    
    SetAdvDiffModelCmd = new G4UIcmdWithABool("/damic/gun/position/setadvdiffmodel", this);
    SetAdvDiffModelCmd->SetGuidance("Sets whether to use the advanced diffusion model or a uniform embeddoing");
    SetAdvDiffModelCmd->SetParameterName("Bool", true);
    SetAdvDiffModelCmd->SetDefaultValue(true);
    
    AddVolumeCmd = new G4UIcommand("/damic/gun/position/addvolume", this);
    AddVolumeCmd->SetGuidance("Add Volume. A weight can be assigned to the volume to mimic the activity or concentration  of the radioactive element. Set multiplicity to run simulations over a family of volumes which PV names are given by Volume_1_PV, Volume_2_PV, ... Volume_N_PV; if multiplicity is set, three extra parameters must be given: min and max copy number to be considered, as well as, the steping value for the iteration algorithm (for instance if only even volumes of the family KCext should be taken into account: KCext 1 true 0 10 2. )");    
    G4UIparameter* param;
    param = new G4UIparameter("Volume", 's', false);
    param->SetDefaultValue("NULL");
    AddVolumeCmd->SetParameter(param);
    param = new G4UIparameter("Concentration", 'd', false);
    param->SetDefaultValue("1");
    param->SetOmittable(true);
    AddVolumeCmd->SetParameter(param);
    // ADDING multiplicity functionality to run several volumes at the same time (the use of 
    // the GDML's loop create volumes called automatically as anyname_1_PV, anyname_2_PV, ... )
    // If we want them considered as the same "volume" for simulations, this parameter must be set
    param = new G4UIparameter("Multiplicity", 'b', false);
    param->SetDefaultValue(false);
    param->SetOmittable(true);
    AddVolumeCmd->SetParameter(param);
    
    param = new G4UIparameter("minCopyNum", 'i', false);
    param->SetGuidance("Starting value of the copyNum to iterate over the defined volume family");
    param->SetDefaultValue(0);
    param->SetOmittable(true);
    AddVolumeCmd->SetParameter(param);
    
    param = new G4UIparameter("maxCopyNum", 'i', false);
    param->SetGuidance("Last copyNum to include in the iteration over the defined volume family");
    param->SetDefaultValue(0);
    param->SetOmittable(true);
    AddVolumeCmd->SetParameter(param);
    
    param = new G4UIparameter("step", 'i', false);
    param->SetGuidance("Step of the iteration algorithm over the volume family");
    param->SetDefaultValue(0);
    param->SetOmittable(true);
    AddVolumeCmd->SetParameter(param);
    
    DoPositionListCmd = new G4UIcmdWithAString("/damic/gun/position/dopositionlist", this);
    DoPositionListCmd->SetGuidance(" Set file to load positions from");
    DoPositionListCmd->SetParameterName("Postion File", false, true);
    
    /*EraseVolumeCmd = new G4UIcmdWithAString("/damic/gun/position/resetvolume", this);
     * EraseVolumeCmd->SetGuidance("Erase volume from the sources volumes.");
     * EraseVolumeCmd->SetParameterName("Volume", true, true);*/
    
    DoSourceCmd = new G4UIcmdWithoutParameter("/damic/gun/position/dosource", this);
    DoSourceCmd->SetGuidance("Set position to source. ");
    
    DoShapeCmd = new G4UIcmdWithAString("/damic/gun/position/doshape", this);
    DoShapeCmd->SetGuidance("Set position to shape");
    DoShapeCmd->SetParameterName("Shape", true, true);
    DoShapeCmd->SetDefaultValue("Para");
    DoShapeCmd->SetCandidates("Para Sphere Cyl");
    
    CenterXCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/centerx", this);
    CenterXCmd->SetGuidance(" Set center coordinate x of the shape or the source");
    CenterXCmd->SetParameterName("X", true, true);
    CenterXCmd->SetUnitCategory("Length");
    CenterXCmd->SetDefaultUnit("mm");
    CenterXCmd->SetUnitCandidates("nm um mm cm m");
    
    CenterYCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/centery", this);
    CenterYCmd->SetGuidance(" Set center coordinate y of the shape or the source");
    CenterYCmd->SetParameterName("Y", true, true);
    CenterYCmd->SetUnitCategory("Length");
    CenterYCmd->SetDefaultUnit("mm");
    CenterYCmd->SetUnitCandidates("nm um mm cm m");
    
    CenterZCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/centerz", this);
    CenterZCmd->SetGuidance(" Set center coordinate z of the shape or the source");  
    CenterZCmd->SetParameterName("Z", true, true);
    CenterZCmd->SetUnitCategory("Length");
    CenterZCmd->SetDefaultUnit("mm");
    CenterZCmd->SetUnitCandidates("nm um mm cm m");
    
    // Para
    ParaXCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/parax", this);
    ParaXCmd->SetGuidance(" Set the Length of the X Side of the Parallelepiped Rectangle");
    ParaXCmd->SetParameterName("X", true, true);
    ParaXCmd->SetUnitCategory("Length");
    ParaXCmd->SetDefaultUnit("mm");
    ParaXCmd->SetUnitCandidates("nm um mm cm m");
    
    ParaYCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/paray", this);
    ParaYCmd->SetGuidance(" Set the Length of the Y Side of the Parallelepiped Rectangle");
    ParaYCmd->SetParameterName("Y", true, true);
    ParaYCmd->SetUnitCategory("Length");
    ParaYCmd->SetDefaultUnit("mm");
    ParaYCmd->SetUnitCandidates("nm um mm cm m");
    
    ParaZCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/paraz", this);
    ParaZCmd->SetGuidance(" Set the Length of the Z Side of the Parallelepiped Rectangle");
    ParaZCmd->SetParameterName("Z", true, true);
    ParaYCmd->SetUnitCategory("Length");
    ParaZCmd->SetDefaultUnit("mm");
    ParaZCmd->SetUnitCandidates("nm um mm cm m");
    
    // Sphere
    SphereRCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/spherer", this);
    SphereRCmd->SetGuidance(" Set the radius of the sphere");
    SphereRCmd->SetParameterName("R", true, true);
    SphereRCmd->SetUnitCategory("Length");
    SphereRCmd->SetDefaultUnit("mm");
    SphereRCmd->SetUnitCandidates("nm um mm cm m");
    SphereRCmd->SetRange(" R >= 0");
    
    //Cyl
    CylRCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/cylr", this);
    CylRCmd->SetGuidance(" Set the radius of the cylinder");
    CylRCmd->SetParameterName("R", true, true);
    CylRCmd->SetUnitCategory("Length");
    CylRCmd->SetDefaultUnit("mm");
    CylRCmd->SetUnitCandidates("nm um mm cm m");
    CylRCmd->SetRange(" R >= 0");
    
    CylHCmd = new G4UIcmdWithADoubleAndUnit("/damic/gun/position/cylh", this);
    CylHCmd->SetGuidance(" Set the height of the cylinder");
    CylHCmd->SetParameterName("H", true, true);
    CylHCmd->SetUnitCategory("Length");
    CylHCmd->SetDefaultUnit("mm");
    CylHCmd->SetUnitCandidates("nm um mm cm m");
    
    // PArticle Specs
    ParticleCmd = new G4UIcmdWithAString("/damic/gun/particle",this);
    ParticleCmd->SetGuidance("Set particle to be generated.");
    ParticleCmd->SetGuidance(" (geantino is default)");
    ParticleCmd->SetGuidance(" (ion can be specified for shooting ions)");
    ParticleCmd->SetParameterName("particleName",true);
    ParticleCmd->SetDefaultValue("geantino");
    G4String candidateList;
    G4int nPtcl = particleTable->entries();
    for(G4int i=0;i<nPtcl;i++)
    {
        candidateList += particleTable->GetParticleName(i);
        candidateList += " ";
    }
    candidateList += "ion ";
    ParticleCmd->SetCandidates(candidateList);
    
    ChargeCmd= new G4UIcmdWithADouble("/damic/gun/charge", this);
    ChargeCmd->SetGuidance("Set charge of the particle.");
    ChargeCmd->SetParameterName("Charge", true);
    
    IonCmd = new G4UIcommand("/damic/gun/ion",this);
    IonCmd->SetGuidance("Set properties of ion to be generated.");
    IonCmd->SetGuidance("[usage] /gun/ion Z A Q E");
    IonCmd->SetGuidance("        Z:(int) AtomicNumber");
    IonCmd->SetGuidance("        A:(int) AtomicMass");
    IonCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
    IonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");
    
    G4UIparameter* parami;
    parami = new G4UIparameter("Z",'i',true);
    parami->SetDefaultValue("1");
    IonCmd->SetParameter(parami);
    parami = new G4UIparameter("A",'i',true);
    parami->SetDefaultValue("1");
    IonCmd->SetParameter(parami);
    parami = new G4UIparameter("Q",'i',false);
    parami->SetDefaultValue("0");
    IonCmd->SetParameter(parami);
    parami = new G4UIparameter("E",'d',false);
    parami->SetDefaultValue("0.0");
    IonCmd->SetParameter(parami);

}

DAMICParticleSourceMessenger::~DAMICParticleSourceMessenger()
{
    delete ChargeCmd;
    delete ParticleCmd;

    delete CylHCmd;
    delete CylRCmd;
    
    delete SphereRCmd;
    delete ParaZCmd;
    delete ParaYCmd;
    delete ParaXCmd;
    //delete IonCmd;
    //zerezerezere
    
    delete CenterZCmd;
    delete CenterYCmd;
    delete CenterXCmd;
    
    delete DoShapeCmd;
    delete DoSourceCmd;
    delete DoSurfaceCmd;
    
    delete MinEmbedDistCmd;
    delete MaxEmbedDistCmd;
    delete SetCCDSurfaceOnlyCmd;
    delete SetAdvDiffModelCmd;
    //delete EraseVolumeCmd;
    delete AddVolumeCmd;
    delete DoVolumeCmd;
    
    delete DoPositionListCmd;
    
    delete DoMaterialCmd;
    delete MotherVolumeCmd;
    
    delete DoOneDCmd;
    delete DoOneDXCmd;
    delete DoOneDYCmd;
    delete DoOneDZCmd;
    
    delete DoDistriDCmd;
    delete IsoThetaMinCmd;
    delete IsoThetaMaxCmd;
    delete IsoPhiMinCmd;
    delete IsoPhiMaxCmd;
    
    delete DoMonoNRJCmd;
    delete DoDistriNRJCmd;
    delete HistoNRJCmd;
    delete UniformBotCmd;
    delete UniformTopCmd;
    delete AvgMaxCmd;
        
    delete ListCmd;
    
    delete EnergyDirectory;
    delete DirectionDirectory;
    delete PositionDirectory;
    
    delete gunDirectory;

    delete CosmoDirectory;
    delete DoCosmoPartCmd;
    delete DoCosmoAngCmd;
    delete DoCosmoPosCmd;
    delete DoCosmoPosCentreCmd;
    delete DoCosmoPosSizeCmd;
}

void DAMICParticleSourceMessenger::SetNewValue(G4UIcommand *command, G4String newValues)
{
    
    //G4cout << "SetNew "  << G4endl;
    if(command == ListCmd)
    {
        particleTable->DumpTable();
    }// SPECS
    else if(command == IonCmd)
    {
        //G4cout << "on teste " << G4endl;
        if(fShootIon)
        {
            G4Tokenizer next( newValues );
            // check argument
            fAtomicNumber = StoI(next());
            fAtomicMass = StoI(next());
            G4String sQ = next();
            if(sQ.empty())
            {
                fIonCharge = fAtomicNumber;
            }
            else
            {
                fIonCharge = StoI(sQ);
                sQ = next();
                if(sQ.empty())
                {
                    fIonExciteEnergy = 0.0;
                }
                else
                {
                    fIonExciteEnergy = StoD(sQ) * keV;
                }
            }
            
            G4ParticleDefinition* ion;
            ion =  G4IonTable::GetIonTable()->GetIon(fAtomicNumber,fAtomicMass,fIonExciteEnergy);
            if(ion==0) 
            {
                G4cout << "Ion with Z=" << fAtomicNumber;
                G4cout << " A=" << fAtomicMass << "is not well defined" << G4endl;
            }
            else 
            {
                fParticleGun->SetParticleDefinition(ion);
                fParticleGun->SetParticleCharge(fIonCharge*eplus);
            }
        }
        else 
        {
            G4cout<<"Set /damic/gun/particle to ion before using /damic/gun/ion command";
            G4cout<<G4endl;
        }
    }
    else if(command == ParticleCmd)
    {
        if(newValues == "ion")
        {
            fShootIon = true;
        }
        else
        {
            fShootIon = false;
            G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
            if(pd != NULL)
            {
                fParticleGun->SetParticleDefinition(pd);
            }
        }
    }
    else if(command == ChargeCmd)
    {
        fParticleGun->SetParticleCharge(ChargeCmd->GetNewDoubleValue(newValues));
    }// ENERGY
    else if(command == DoMonoNRJCmd)
    {
        fParticleGun->DoMonoNRJ();
        fParticleGun->ResetEnergyNum();
        fParticleGun->SetEnergyNumSize(1);
        fParticleGun->SetEnergyNumValue(DoMonoNRJCmd->GetNewDoubleValue(newValues),0);
    }
    else if(command == DoDistriNRJCmd)
    {
        fParticleGun->DoDistriNRJ(newValues);
        fParticleGun->ResetEnergyNum();
        fParticleGun->SetEnergyNumSize(2);
    }
    else if(command == HistoNRJCmd)
    {
        std::ifstream file(newValues, std::ios::in);
        if(file)
        {
            G4cout << " File Exists" << G4endl;
            fParticleGun->SetHistoNRJ(newValues);
            fParticleGun->ResetEnergyNum();
            file.close();
        }
        else
        {
            G4cout << "File doesn't exist" << G4endl;
        }
    }
    else if(command == UniformBotCmd)
    {
        fParticleGun->SetEnergyNumValue(UniformBotCmd->GetNewDoubleValue(newValues),0);
    }
    else if(command == UniformTopCmd)
    {
        fParticleGun->SetEnergyNumValue(UniformTopCmd->GetNewDoubleValue(newValues),1);
    }
    
    else if(command == AvgMaxCmd)
    {
        fParticleGun->SetEnergyNumValue(AvgMaxCmd->GetNewDoubleValue(newValues),2);
    }
    
    else if(command == DoOneDCmd)
    {
        fParticleGun->DoOneD();
        fParticleGun->ResetDirectionNum();
        fParticleGun->SetDirectionNumSize(3);
    }
    else if(command == DoOneDXCmd)
    {
        G4cout << "prob" << G4endl;
        fParticleGun->SetDirectionNumValue(DoOneDXCmd->GetNewDoubleValue(newValues),0);
    }
    else if(command == DoOneDYCmd)
    {
        fParticleGun->SetDirectionNumValue(DoOneDYCmd->GetNewDoubleValue(newValues),1);
    }
    else if(command == DoOneDZCmd)
    {
        fParticleGun->SetDirectionNumValue(DoOneDZCmd->GetNewDoubleValue(newValues),2);
    }
    else if(command == DoDistriDCmd)
    {
        fParticleGun->DoDistriD(newValues);
        fParticleGun->ResetDirectionNum();
        fParticleGun->SetDirectionNumSize(4);
    }
    else if(command == IsoThetaMinCmd)
    {
        fParticleGun->SetDirectionNumValue(IsoThetaMinCmd->GetNewDoubleValue(newValues)*M_PI/180,0);
    }
    else if(command == IsoThetaMaxCmd)
    {
        fParticleGun->SetDirectionNumValue(IsoThetaMaxCmd->GetNewDoubleValue(newValues)*M_PI/180,1);
    }
    else if(command == IsoPhiMinCmd)
    {
        fParticleGun->SetDirectionNumValue(IsoPhiMinCmd->GetNewDoubleValue(newValues)*M_PI/180,2);
    }
    else if(command == IsoPhiMaxCmd)
    {
        fParticleGun->SetDirectionNumValue(IsoPhiMaxCmd->GetNewDoubleValue(newValues)*M_PI/180,3);
    }
    else if(command == DoMaterialCmd)
    {
        fParticleGun->DoMaterial();
        fParticleGun->SetMaterial(newValues);
    }
    else if(command == MotherVolumeCmd)
    {
        fParticleGun->SetMotherVolume(newValues);
    }
    else if(command ==  DoVolumeCmd)
    {
        fParticleGun->DoVolume();
    }
    else if(command == AddVolumeCmd)
    {
        // G4Tokenizer next(newValues);
        // fParticleGun->AddVolume(next(), StoD(next()));
        std::vector<G4String> tokens;
        G4Analysis::Tokenize(newValues, tokens);
        // four new parameters have been added to account for simulations
        // // under several volumes
        if(StoB(tokens[2]))
        {
            G4cout << "\t MESSAGe: Multiple Simulations" << G4endl;
            G4cout << "\t The following volumes will be considered on this simulation " << G4endl;
            for(int copynum= StoI(tokens[3]); copynum <= StoI(tokens[4]); copynum=copynum+StoI(tokens[5]))
            {
                G4cout << "\t" << tokens[0]+"_"+ItoS(copynum)+"_PV" << G4endl;
                G4String Volume = tokens[0]+"_"+ItoS(copynum)+"_PV";
                fParticleGun->AddVolume(Volume, StoD(tokens[1]));
            }
        }
        else
        {
            fParticleGun->AddVolume(tokens[0], StoD(tokens[1]));
        }
    }
    else if(command == DoShapeCmd)
    {
        fParticleGun->DoShape(newValues);
        fParticleGun->ResetPositionNum();
        fParticleGun->SetPositionNumSize(6);
    }
    else if(command == DoSourceCmd)
    {
        fParticleGun->DoSource();
        fParticleGun->ResetPositionNum();
        fParticleGun->SetPositionNumSize(3);
    }
    else if(command == DoSurfaceCmd) 
    {
        fParticleGun->DoSurface();
    }
    else if(command == MinEmbedDistCmd)
    {
        fParticleGun->SetMinEmbeddedDistance(MinEmbedDistCmd->GetNewDoubleValue(newValues));
    }
    else if(command == MaxEmbedDistCmd) 
    {
        fParticleGun->SetMaxEmbeddedDistance(MaxEmbedDistCmd->GetNewDoubleValue(newValues));
    }
    else if(command== SetCCDSurfaceOnlyCmd) 
    {
        fParticleGun->SetCCDFacingOnly(SetCCDSurfaceOnlyCmd->GetNewBoolValue(newValues));
    }
    else if(command== SetAdvDiffModelCmd)
    {
        fParticleGun->SetAdvDiffModel(SetAdvDiffModelCmd->GetNewBoolValue(newValues));
    }
    else if(command == DoPositionListCmd)
    {
        
        fParticleGun->SetPositionList(newValues);
    }
    else if(command == CenterXCmd)
    {
        fParticleGun->SetPositionNumValue(CenterXCmd->GetNewDoubleValue(newValues), 0);
    }
    else if(command == CenterYCmd)
    {
        fParticleGun->SetPositionNumValue(CenterYCmd->GetNewDoubleValue(newValues), 1);
    }
    else if(command == CenterZCmd)
    {
        fParticleGun->SetPositionNumValue(CenterZCmd->GetNewDoubleValue(newValues), 2);
    }
    else if(command == ParaXCmd)
    {
        fParticleGun->SetPositionNumValue(ParaXCmd->GetNewDoubleValue(newValues), 3);
    }
    else if(command == ParaYCmd)
    {
        fParticleGun->SetPositionNumValue(ParaYCmd->GetNewDoubleValue(newValues), 4);
    }
    else if(command == ParaZCmd)
    {
        fParticleGun->SetPositionNumValue(ParaZCmd->GetNewDoubleValue(newValues), 5);
    }
    else if(command == SphereRCmd)
    {
        fParticleGun->SetPositionNumValue(SphereRCmd->GetNewDoubleValue(newValues), 3);
    }
    else if(command == CylRCmd)
    {
        fParticleGun->SetPositionNumValue(CylRCmd->GetNewDoubleValue(newValues), 3);
    }
    else if(command == CylHCmd)
    {
        fParticleGun->SetPositionNumValue(CylHCmd->GetNewDoubleValue(newValues), 4);
    }
    //cosmogenics - shielding
    else if(command == DoCosmoPartCmd) {
      fParticleGun->DoCosmogenic();
      fParticleGun->SetCosmoPart(newValues);
    }
    else if(command == DoCosmoAngCmd)
      fParticleGun->SetCosmoAngDistType(newValues);
        
    else if(command == DoCosmoPosCentreCmd)
      fParticleGun->SetCosmoSrcCenter(DoCosmoPosCentreCmd->GetNew3VectorValue(newValues));

    else if(command == DoCosmoPosSizeCmd)
      fParticleGun->SetCosmoSrcSize(DoCosmoPosSizeCmd->GetNew3VectorValue(newValues));

    else if(command == ShieldReaderCmd) {
      fParticleGun->SetReadingShieldFile(true);
      fParticleGun->ReadParticleFile(newValues);
    }
}
