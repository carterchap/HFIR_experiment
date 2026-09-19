#include "DAMICDetectorConstruction.hh"
#include "DAMICCCDSD.hh"

#include "G4GDMLParser.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include <set>

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

#include "G4SystemOfUnits.hh"

#include "globals.hh"

#include <cmath>

namespace
{
    G4VisAttributes* GetMaterialVisAttributes(const G4Material* mat)
    {
        if (!mat) {
            auto vis = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0)); // fallback white
            vis->SetForceSolid(true);
            return vis;
        }

        const G4String& name = mat->GetName();

        // Bright, distinct colors
        if (name == "G4_AIR" || name == "G4_Galactic") {
            auto vis = new G4VisAttributes(false); // invisible for air/vacuum
            return vis;
        }
        else if (name == "G4_Si") {
            auto vis = new G4VisAttributes(G4Colour(0.1, 1.0, 0.1)); // bright green
            vis->SetForceSolid(true);
            return vis;
        }
        else if (name == "G4_Pb" || name == "Lead") {
            auto vis = new G4VisAttributes(G4Colour(1.0, 0.2, 0.2)); // bright red
            vis->SetForceSolid(true);
            return vis;
        }
        else if (name == "G4_STAINLESS-STEEL" || name == "StainlessSteel") {
            auto vis = new G4VisAttributes(G4Colour(0.2, 0.7, 1.0)); // bright cyan/blue
            vis->SetForceSolid(true);
            return vis;
        }
        else if (name == "G4_Al") {
            auto vis = new G4VisAttributes(G4Colour(1.0, 0.8, 0.1)); // bright yellow
            vis->SetForceSolid(true);
            return vis;
        }
        else if (name == "G4_Cu") {
            auto vis = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0)); // orange
            vis->SetForceSolid(true);
            return vis;
        }
        else if (name == "G4_KAPTON") {
            auto vis = new G4VisAttributes(G4Colour(0.8, 0.0, 1.0)); // magenta
            vis->SetForceSolid(true);
            return vis;
        }

        // Default bright color for anything else
        auto vis = new G4VisAttributes(G4Colour(1.0, 0.4, 1.0)); // bright pink
        vis->SetForceSolid(true);
        return vis;
    }

    void ApplyMaterialColorsRecursive(G4LogicalVolume* lv)
    {
        if (!lv) return;

        lv->SetVisAttributes(GetMaterialVisAttributes(lv->GetMaterial()));

        auto nDaughters = lv->GetNoDaughters();
        for (G4int i = 0; i < nDaughters; ++i) {
            auto phys = lv->GetDaughter(i);
            if (phys && phys->GetLogicalVolume()) {
                ApplyMaterialColorsRecursive(phys->GetLogicalVolume());
            }
        }
    }
}

DAMICDetectorConstruction::DAMICDetectorConstruction(const G4GDMLParser& parser) :
    _gdmlParser(parser)
{
}

DAMICDetectorConstruction::~DAMICDetectorConstruction()
{
}

namespace
{
    void DumpGlobalExtent(G4VPhysicalVolume* pv, const G4AffineTransform& parentToWorld,
                           const std::set<G4String>& targets, int depth=0)
    {
        G4AffineTransform local(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
        G4AffineTransform toWorld = local * parentToWorld;

        G4LogicalVolume* lv = pv->GetLogicalVolume();
        if(targets.count(lv->GetName()))
        {
            G4ThreeVector pMin, pMax;
            lv->GetSolid()->BoundingLimits(pMin, pMax);
            G4ThreeVector corners[8];
            int idx=0;
            for(int ix=0; ix<2; ix++) for(int iy=0; iy<2; iy++) for(int iz=0; iz<2; iz++)
            {
                G4ThreeVector local_c(ix? pMax.x():pMin.x(), iy? pMax.y():pMin.y(), iz? pMax.z():pMin.z());
                corners[idx++] = toWorld.TransformPoint(local_c);
            }
            G4ThreeVector gmin=corners[0], gmax=corners[0];
            for(auto& c: corners)
            {
                gmin.setX(std::min(gmin.x(),c.x())); gmax.setX(std::max(gmax.x(),c.x()));
                gmin.setY(std::min(gmin.y(),c.y())); gmax.setY(std::max(gmax.y(),c.y()));
                gmin.setZ(std::min(gmin.z(),c.z())); gmax.setZ(std::max(gmax.z(),c.z()));
            }
            G4cout << "TEMP DEBUG EXTENT '" << lv->GetName() << "' (pv=" << pv->GetName() << ", depth=" << depth
                   << ") global AABB: "
                   << "x[" << gmin.x() << "," << gmax.x() << "] "
                   << "y[" << gmin.y() << "," << gmax.y() << "] "
                   << "z[" << gmin.z() << "," << gmax.z() << "]" << G4endl;
        }

        for(int i=0;i<lv->GetNoDaughters();++i)
        {
            DumpGlobalExtent(lv->GetDaughter(i), toWorld, targets, depth+1);
        }
    }
}

G4VPhysicalVolume* DAMICDetectorConstruction::Construct()
{
    G4VPhysicalVolume* WorldPV = _gdmlParser.GetWorldVolume();
    // FIXME
    G4cout << "World Physical Volume renamed by the parser to: '" << WorldPV->GetName()
    << "' --> THIS SHOULD BE OBTAINED AND PROPAGATED (used by DAMICVolumeCoordinateStore, for instance)" << G4endl;

    //RegionInitialization();

    // TEMP DEBUG: verify CCD sub-layer orientation, not just the outer box
    {
        std::set<G4String> targets = {"CCD","CCDSensor","CCDdeadLayerTop","CCDdeadLayerBot","SiO2Layer_1","SiO2Layer_0"};
        G4AffineTransform identity;
        DumpGlobalExtent(WorldPV, identity, targets);
    }

    G4cout << " Info -- List of used materials :: " << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    return WorldPV;
}

void DAMICDetectorConstruction::ConstructSDandField()
{
    //------------------------------------------------
    // Sensitive detectors
    //------------------------------------------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    DAMICCCDSD* aTrackerSD = new DAMICCCDSD("CCD");
    SDman->AddNewDetector(aTrackerSD);

    ///////////////////////////////////////////////////////////////////////
    //
    // Example how to retrieve Auxiliary Information for sensitive detector
    //
    const G4GDMLAuxMapType* auxmap = _gdmlParser.GetAuxMap();
    G4cout << "Found " << auxmap->size() << " volume(s) with auxiliary information." << G4endl << G4endl;

    // The same as above, but now we are looking for
    // sensitive detectors setting them for the volumes
    for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++)
    {
        G4cout << "Volume " << ((*iter).first)->GetName() << " has the following list of auxiliary information: " << G4endl << G4endl;

        for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();vit!=(*iter).second.end();vit++)
        {
                if((*vit).type=="SensDet")
                {
                    G4cout << "Attaching sensitive detector " << (*vit).value
                    << " to volume " << ((*iter).first)->GetName()
                    <<  G4endl << G4endl;

                    G4VSensitiveDetector* mydet = SDman->FindSensitiveDetector((*vit).value);
                    if(mydet)
                    {
                        G4LogicalVolume* myvol = (*iter).first;
                        myvol->SetSensitiveDetector(mydet);
                    }
                    else
                    {
                        G4cout << (*vit).value << " detector not found" << G4endl;
                    }
                }
                if((*vit).type=="Color")
                {
                    auto color = _GetColor((*vit).value);
                    auto va= new G4VisAttributes( );
                    va->SetColor( color );
                    va->SetVisibility(true);
                    G4LogicalVolume* myvol = (*iter).first;
                    myvol->SetVisAttributes(va);
                }
                if((*vit).type=="Visibility")
                {
                    auto va= new G4VisAttributes( );
                    if((*vit).value=="0")
                    {
                        va->SetVisibility(false);
                    }
                    else
                    {
                        va->SetVisibility(true);
                    }
                    G4LogicalVolume* myvol = (*iter).first;
                    myvol->SetVisAttributes(va);
                }
                if ((*vit).type=="Region")
                {
                    // From GDML Region attribute will have the region_name, and the energy cut, as
                    // a comma separated values  (region_name,Ecut)
                    auto region_tuple = _GetRegionAndCut( (*vit).value );

                    // create if do not exist the region
                    G4Region* region = G4RegionStore::GetInstance()->FindOrCreateRegion(std::get<0>(region_tuple)); 
                    G4cout << " Info ---- Region " << region->GetName() << " is added to Region Store with length cut " 
                        << std::get<1>(region_tuple) << " mm for e-, " <<std::get<2>(region_tuple) << " mm for e+, " <<std::get<3>(region_tuple) << " mm for gamma, " << std::get<4>(region_tuple) << " mm for prot"<< G4endl;

                    // Add volume to region
                    G4LogicalVolume* theLV = (*iter).first;
                    region->AddRootLogicalVolume(theLV);

                    // Add cut to region
                    G4ProductionCuts* Lcut = new G4ProductionCuts;
                    Lcut->SetProductionCut(std::get<1>(region_tuple)*mm,G4ProductionCuts::GetIndex("e-"));
		    Lcut->SetProductionCut(std::get<2>(region_tuple)*mm,G4ProductionCuts::GetIndex("e+"));
		    Lcut->SetProductionCut(std::get<3>(region_tuple)*mm,G4ProductionCuts::GetIndex("gamma"));
		    Lcut->SetProductionCut(std::get<4>(region_tuple)*mm,G4ProductionCuts::GetIndex("proton"));
                    region->SetProductionCuts(Lcut);
                }
        }
    }

    // Default per-material visualization colors, applied last so they win
    // over the GDML "Color" auxiliary tags just processed above (the
    // CAD-exported files carry hundreds of per-part colors from the
    // original CAD assembly, which would otherwise silently override this).
    ApplyMaterialColorsRecursive(_gdmlParser.GetWorldVolume()->GetLogicalVolume());
}

std::tuple<std::string, double, double, double, double> DAMICDetectorConstruction::_GetRegionAndCut(const std::string& in) 
{
    // Attribute region will have two value: region_name, lenght_cut (in units of mm)
    char sep = ',';
    std::string::size_type b = 0;
    std::vector<std::string> result;

    auto k=0;
    while((b = in.find_first_not_of(sep, b)) != std::string::npos) 
    {
        auto e = in.find_first_of(sep, b);
        result.push_back( in.substr(b, e-b) );
        b = e;
        k = k+1;
    }
    
    // 0.00001 mm as the length cut by default in case no length is given through Region attribute
    if(k == 1)
    {
        result.push_back( "0.00001" );
        result.push_back( "0.00001" );
        result.push_back( "0.00001" );
        result.push_back( "0.00001" );
    }
    if(k ==2){
	result.push_back( result.at(1) );
	result.push_back( result.at(1) );
	result.push_back( result.at(1) );	

    }
     auto rt = std::make_tuple( result.at(0), std::stod(result.at(1)),std::stod(result.at(2)),std::stod(result.at(3)),std::stod(result.at(4)) );
    return rt;
}



G4ThreeVector DAMICDetectorConstruction::_GetColor(const std::string& in)
{
    // Some GDML exporters (e.g. CAD-to-GDML tools) write colors as a
    // "#RRGGBB" or "#AARRGGBB" hex string instead of DAMIC's own
    // "r,g,b" comma-separated float convention.
    if(!in.empty() && in[0]=='#')
    {
        std::string hex = in.substr(1);
        std::string rgb = (hex.size()>=8) ? hex.substr(2,6) : hex;
        unsigned long val = std::stoul(rgb, nullptr, 16);
        double r = ((val>>16)&0xFF)/255.0;
        double g = ((val>>8)&0xFF)/255.0;
        double b = (val&0xFF)/255.0;
        return G4ThreeVector(r,g,b);
    }

    char sep = ',';
    std::string::size_type b = 0;
    std::vector<std::string> result;

    while((b = in.find_first_not_of(sep, b)) != std::string::npos) 
    {
        auto e = in.find_first_of(sep, b);
        result.push_back( in.substr(b, e-b) );
        b = e;
    }

    G4ThreeVector color = G4ThreeVector(std::stod(result.at(0)),std::stod(result.at(1)),std::stod(result.at(2)));
    return color;
}


