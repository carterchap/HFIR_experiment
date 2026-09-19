#include "DAMICVolumeCoordinateStore.hh"
#include "DAMICPositionFileReader.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "Randomize.hh"
#include <cmath>
#include <iostream>
#include <fstream>

#include "G4AutoLock.hh"
namespace {G4Mutex PositionFileReaderMutex = G4MUTEX_INITIALIZER;}

DAMICPositionFileReader* DAMICVolumeCoordinateStore::fileReader=0;

DAMICVolumeCoordinateStore::DAMICVolumeCoordinateStore(G4PhysicalVolumeStore* aStore) : theStore(aStore) {
    NumberOfVolumes=0;
    minEmbeddedDistance=100*nm;
    maxEmbeddedDistance=100*nm;
    CCDFacingOnly=false;
    advDiffModel=false;
    CCDOrigin=G4ThreeVector(0.,0.,0.);
}

DAMICVolumeCoordinateStore::~DAMICVolumeCoordinateStore() {
    for (std::vector<G4VSolid*>::size_type i=0; i < Rotations.size(); i++) {
        delete Rotations[i];
    }
    G4AutoLock lock(&PositionFileReaderMutex);
    if (fileReader) delete fileReader;
}

void DAMICVolumeCoordinateStore::SetEmbeddedDistance(G4double min_distance, G4double max_distance) {
    minEmbeddedDistance=min_distance;
    maxEmbeddedDistance=max_distance;
}


void DAMICVolumeCoordinateStore::SetMinEmbeddedDistance(G4double min_distance) {
    minEmbeddedDistance=min_distance;
}
void DAMICVolumeCoordinateStore::SetMaxEmbeddedDistance(G4double max_distance) {
    maxEmbeddedDistance=max_distance;
}

void DAMICVolumeCoordinateStore::SetCCDFacingOnly(G4bool val) {
    CCDFacingOnly=val;
}

void DAMICVolumeCoordinateStore::SetAdvDiffModel(G4bool val) {
    advDiffModel=val;
}

void DAMICVolumeCoordinateStore::SetPositionList(G4String fname) {
    G4cout << "Loading position file" << G4endl;
    G4AutoLock lock(&PositionFileReaderMutex);
    if (!fileReader) fileReader = new DAMICPositionFileReader(fname);
    else fileReader->LoadFile(fname);
}

G4bool DAMICVolumeCoordinateStore::AddVolume(G4String Volume, G4double Concentration)
{
    histogramBuilt=false;
    G4bool done=false;
    G4bool checked=false;
    G4int volsfound=0;
    G4VPhysicalVolume *WorldPV = theStore->GetVolume("World_PV");
    G4VPhysicalVolume *currentPV = WorldPV;
    G4LogicalVolume *currentLV = nullptr;
    std::vector<G4VPhysicalVolume*> currentMothers;

    //Holds the tree of which daughter we've gone down
    std::vector<G4int> currentDaughter={0};
    G4int depth=0;
    //Descend through all the daughters looking for all instances of the Volume
    //And add it to our list
    while (!done)
    {
        currentLV = currentPV->GetLogicalVolume();
        //First check if our current volume is one of the ones we want
        if((currentPV->GetName()==Volume) && !checked)
        {
            _updateStore(currentPV,currentLV,Concentration,currentMothers);
            volsfound+=1;
            //TODO: add global coordinates of those volumes
            //Crude, but we add the current PV to the end of the Mothers, then remove it later (simpifies the algorithm)
            std::vector<G4RotationMatrix*> rotation_vec;
            rotation_vec.push_back(new G4RotationMatrix(0,0,0));
            std::vector<G4ThreeVector> translation_vec;
            translation_vec.push_back(G4ThreeVector(0,0,0));
            /// --
            currentMothers.push_back(currentPV);
            for (std::vector<G4VPhysicalVolume>::size_type i=0; i < currentMothers.size(); i++)
            {
                if(currentMothers[i]->IsReplicated())
                {
                    // Obtain replica parameters
                    EAxis axis= EAxis::kUndefined;
                    G4int  nreplicas=-1;
                    G4double width = -1;
                    G4double offset=-1;
                    G4bool consuming=false;
                    currentMothers[i]->GetReplicationData(axis,nreplicas,width,offset,consuming);
                    // The global position (i.e. the translation and rotation) of a replicated
                    // volume must be propagated to the translation vector, because any daughter is
                    // expressed only with respect this replicated volume mother (in local
                    // coordinates of the mother).
                    // Therefore, for each element of the current translation_vec, nreplicas element
                    // should born. The total number of elements for the new translation_vec is
                    // the (old)translation_vec.size()*nreplicas.
                    //
                    // Copy a temporary translation and rotation vector
                    auto old_translation_vec = translation_vec;
                    const auto old_size = old_translation_vec.size();
                    translation_vec.clear();

                    auto old_rotation_vec = rotation_vec;
                    if(old_size != old_rotation_vec.size())
                    {
                        // XXX G4 throw exception
                        G4cout << "Unexpected ERROR. Contact developer [E011]" << G4endl;
                    }
                    rotation_vec.clear();
                    // Prepare nreplicas*old_translation_vec.size() elements
                    // and initialize them to 0,0,0
                    translation_vec.resize(old_size*nreplicas);
                    rotation_vec.resize(old_size*nreplicas);
                    for(unsigned int old_i=0; old_i < old_size; ++old_i)
                    {
                        // For each old element, nreplicas elements must be added to it
                        for(unsigned int iR=0; iR < static_cast<unsigned int>(nreplicas); ++iR)
                        {
                            const auto new_i = old_i*nreplicas+iR;
                            // Add to the previous translation, the current translation wrt the
                            // mother
                            // ??? ->> translation_vec[el]+= (*rotation)*currentMothers[i]->GetTranslation();??
                            translation_vec[new_i] = old_translation_vec[old_i];
                            // Is cartesian
                            if(axis == EAxis::kXAxis || axis == EAxis::kYAxis || axis == EAxis::kZAxis)
                            {
                                // Replicate the element along the axis
                                // http://www-geant4.kek.jp/lxr/source/geometry/volumes/include/G4PVReplica.hh#L51
                                /* ??? WHY IS IT NOT NEEDED THIS?
                                 * G4ThreeVector tmp_vect(translation_vec[new_i]);
                                tmp_vect[axis]+= (-width*(nreplicas-1)*0.5+iR*width);
                                if(old_rotation_vec[old_i] != nullptr)
                                {
                                    translation_vec[new_i] = (*(old_rotation_vec[old_i]))*tmp_vect;
                                } --  */
                                translation_vec[new_i][axis] += -width*(nreplicas-1)*0.5+iR*width;
                            }
                            else if(axis == EAxis::kRho)
                            {
                                // http://www-geant4.kek.jp/lxr/source/geometry/volumes/include/G4PVReplica.hh#L57
                                // XXX: TO BE CODIFIED??
                            }
                            else if(axis == EAxis::kPhi)
                            {
                                // http://www-geant4.kek.jp/lxr/source/geometry/volumes/include/G4PVReplica.hh#L62
                                // XXX: TO BE CODIFIED??
                            }
                            else
                            {
                                // XXX -- Use G4 Exception and G4 errors!
                                G4cout << "Undefined EAxis on the replicated volume" << G4endl;
                                exit(-1);
                            }
                            // Rotation vector has to be updated as well
                            // XXX TO BE DELETED!!!???
                            rotation_vec[new_i] = new G4RotationMatrix(*(old_rotation_vec[old_i]));
                            if(currentMothers[i]->GetRotation() != nullptr)
                            {
                                *(rotation_vec[new_i]) *= *(currentMothers[i]->GetObjectRotation());
                            }
                        }
                    }
                    // Add the same amount of volumes in the store (mimicking as they were actually
                    // placed as G4VPhysicalVolume), remove 1 just because was already added before
                    // at the begining of the if(currentLV == Volume blabha) -> L88
                    for(unsigned int inew_store=0; inew_store < (nreplicas*old_size)-1; ++inew_store)
                    {
                        _updateStore(currentPV,currentLV,Concentration,currentMothers);
                    }
                }
                else
                {
                    for(unsigned int el=0; el < translation_vec.size(); ++el)
                    {
                        translation_vec[el]+= (*(rotation_vec[el]))*currentMothers[i]->GetTranslation();
                        if(currentMothers[i]->GetRotation() != nullptr)
                        {
                            *(rotation_vec[el]) *= *(currentMothers[i]->GetObjectRotation());
                        }
                    }
                }
            }
            //G4cout << "Total Translation: " << translation << G4endl;
            if( translation_vec.size() != rotation_vec.size() )
            {
                // XXX throw G4 exception
                G4cout << "Unexplicable ERRROOOOOOR!!!" << G4endl;
                exit(-2);
            }
            for(unsigned int i_final=0; i_final < translation_vec.size(); ++i_final)
            {
                Translations.push_back(translation_vec[i_final]);
                Rotations.push_back(rotation_vec[i_final]);
            }
            currentMothers.pop_back();
        }

        //Next get the daughter, if it has one
        if(currentLV->GetNoDaughters() > currentDaughter[depth] )
        {
            checked=false;
            //Add our mother to the chain
            currentMothers.push_back(currentPV);
            //get the daughter
            currentPV=currentLV->GetDaughter(currentDaughter[depth]);
            //increment our counters
            currentDaughter[depth]+=1;
            currentDaughter.push_back(0);
            depth+=1;
        }
        else if(currentPV->GetName()=="World_PV")
        {
            done=true;
        }
        else
        {
            //Go back up the chain
            checked=true; // We've already checked if want to add this to our Volumes
            currentPV=currentMothers.back();
            currentMothers.pop_back();
            currentDaughter.pop_back();
            depth-=1;
        }
    }
    if (volsfound == 0) {
        G4cout << "Warning: found no volumes" << G4endl;
        return false;
    }
    NumberOfVolumes+=volsfound;
    //G4cout << "Found " << volsfound << " volumes" << G4endl;
    return true;
}

G4bool DAMICVolumeCoordinateStore::RemoveVolume(G4String Volume) {
    histogramBuilt=false;
    G4bool removeall = (Volume=="DELETEALL");
    //Go backwards so if we delete anything we don't have to reindex
    G4int volsdeleted=0;
    for (G4int i=_physicalVolumes.size()-1; i >= 0; i-- ) {
        if ((_physicalVolumes[i]->GetName()==Volume) || removeall) {
            _physicalVolumes.erase(_physicalVolumes.begin()+i);
            _logicalVolumes.erase(_logicalVolumes.begin()+i);
            _solids.erase(_solids.begin()+i);
            Masses.erase(Masses.begin()+i);
            CubicVolumes.erase(CubicVolumes.begin()+i);
            SurfaceAreas.erase(SurfaceAreas.begin()+i);
            Concentrations.erase(Concentrations.begin()+i);
            Mothers.erase(Mothers.begin()+i);
            delete Rotations[i];
            Rotations.erase(Rotations.begin()+i);
            Translations.erase(Translations.begin()+i);
            volsdeleted+=1;
        }
    }
    if (volsdeleted == 0) {
        return false;
    }
    NumberOfVolumes-=volsdeleted;
    return true;
}

void DAMICVolumeCoordinateStore::_updateStore(G4VPhysicalVolume * currentPV, G4LogicalVolume * currentLV,
            G4double concentration, std::vector<G4VPhysicalVolume*> currentMothers)
{
    _physicalVolumes.push_back(currentPV);
    Masses.push_back(currentLV->GetMass());
    _logicalVolumes.push_back(currentLV);
    _solids.push_back(currentLV->GetSolid());
    CubicVolumes.push_back(currentLV->GetSolid()->GetCubicVolume());
    SurfaceAreas.push_back(currentLV->GetSolid()->GetSurfaceArea());
    Concentrations.push_back(concentration);
    Mothers.push_back(currentMothers);
}



void DAMICVolumeCoordinateStore::Clear() {
    RemoveVolume("DELETEALL");
}


G4ThreeVector DAMICVolumeCoordinateStore::GetPointOnSurfaceByName(G4String VolumeName) {
    for (std::vector<G4VPhysicalVolume>::size_type i=0; i < _physicalVolumes.size(); i++) {
        if (_physicalVolumes[i]->GetName()==VolumeName) {
            return GetPointOnSurfaceByID(i);
        }
    }
    G4cout << "Warning: Volume not found in coordinate store" << G4endl;
    return G4ThreeVector(0,0,0);
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointInVolumeByName(G4String VolumeName) {
    for (std::vector<G4VPhysicalVolume>::size_type i=0; i < _physicalVolumes.size(); i++) {
        if (_physicalVolumes[i]->GetName()==VolumeName) {
            return GetPointInVolumeByID(i);
        }
    }
    G4cout << "Warning: Volume not found in coordinate store" << G4endl;
    return G4ThreeVector(0,0,0);
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointOnSurfaceByID(std::vector<G4VSolid>::size_type VolumeID) {
    if (VolumeID >= _solids.size()) {
        G4cout << "Warning: volume not found in coordinate store" << G4endl;
        return G4ThreeVector(0,0,0);
    }
    G4bool accepted=false;
    G4ThreeVector globalPoint=G4ThreeVector(0,0,0);
    while (!accepted)
    {
        G4ThreeVector localPoint=_solids[VolumeID]->GetPointOnSurface();
        G4ThreeVector normal=_solids[VolumeID]->SurfaceNormal(localPoint);
        G4double embed_depth;
        if(advDiffModel) 
        {
            G4bool foundDist=false;
            while (!foundDist)
            {
                //Depth=sqrt(pi)/k*erfc(x/k)
                G4double temp=G4UniformRand()*2;
                if(G4UniformRand()<std::erfc(temp))
                {
                    foundDist=true;
                    embed_depth=temp*maxEmbeddedDistance;
                }
            }
        }
        else
        {
            if(minEmbeddedDistance>=maxEmbeddedDistance)
            {
                embed_depth=maxEmbeddedDistance;
            }
            else
            {
                embed_depth=minEmbeddedDistance+(maxEmbeddedDistance-minEmbeddedDistance)*G4UniformRand();
            }
        }
        G4ThreeVector shift=normal*embed_depth;
        localPoint-=shift;
        globalPoint=((*Rotations[VolumeID])*localPoint)+Translations[VolumeID];
        if(CCDFacingOnly)
        {
            //Should be the global orientation of our normal vector
            G4ThreeVector globalNormal=(*Rotations[VolumeID])*normal;
            G4ThreeVector delta=CCDOrigin-globalPoint;
            if(globalNormal.dot(delta)>0)
            {
                accepted=true;
            }
        }
        else
        {
            accepted=true;
        }
    }
    //  G4ThreeVector globalPoint = localPoint * (*Rotations[VolumeID]) + Translations[VolumeID];
    //  G4ThreeVector globalPoint=((_solids[VolumeID]->GetPointOnSurface())+Translations[VolumeID]);
    return globalPoint;
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointInVolumeByID(std::vector<G4VSolid>::size_type VolumeID)
{
    //Algorithm based on EstimateCubicVolume from G4VSolid source code
    if(VolumeID >= _solids.size())
    {
        G4cout << "Warning: volume not found in coordinate store" << G4endl;
        return G4ThreeVector(0,0,0);
    }
    G4VSolid *solid=_solids[VolumeID];
    EInside in=kOutside;
    G4double minX, maxX,minY,maxY,minZ,maxZ,px,py,pz;
    G4ThreeVector localPoint;
    G4VoxelLimits limit;
    G4AffineTransform origin;
    //Compute the extent of the volume along x,y,z axis
    solid->CalculateExtent(kXAxis, limit, origin,minX,maxX);
    solid->CalculateExtent(kYAxis, limit, origin,minY,maxY);
    solid->CalculateExtent(kZAxis, limit, origin,minZ,maxZ);
/*  G4cout << (*solid).GetEntityType() << " solid dimensions Xmin, Xmax " << minX << " " << maxX << G4endl ;
    G4cout << " solid dimensions Ymin, Ymax " << minY << " " << maxY << G4endl ;
    G4cout << " solid dimensions Zmin, Zmax " << minZ << " " << maxZ << G4endl ;*/
/*  while (in==kOutside) {
        px=-maxX+2*maxX*G4UniformRand();
        py=-maxY+2*maxY*G4UniformRand();
        pz=-maxZ+2*maxZ*G4UniformRand();
        localPoint=G4ThreeVector(px,py,pz);
        in=solid->Inside(localPoint);
    }*/
    G4double Xdim = std::max(std::abs(minX),std::abs(maxX));
    G4double Ydim = std::max(std::abs(minY),std::abs(maxY));
    G4double Zdim = std::max(std::abs(minZ),std::abs(maxZ));
    while (in==kOutside) {
        px=-Xdim+2*Xdim*G4UniformRand();
        py=-Ydim+2*Ydim*G4UniformRand();
        pz=-Zdim+2*Zdim*G4UniformRand();
        localPoint=G4ThreeVector(px,py,pz);
        in=solid->Inside(localPoint);
    }

    G4ThreeVector globalPoint=((*Rotations[VolumeID])*localPoint)+Translations[VolumeID];
    
/*  G4cout << " translation " << Translations[VolumeID] << G4endl;
    G4cout << " rotation " << (*Rotations[VolumeID]).phi() << G4endl;
    G4cout << globalPoint << G4endl;*/
    return globalPoint;
}

void DAMICVolumeCoordinateStore::BuildHistogram() {
    totalSurfaceProbability=0;
    totalVolumeProbability=0;
    surfaceProbabilities.clear();
    volumeProbabilities.clear();
    for (std::vector<G4VPhysicalVolume>::size_type i=0; i <=_physicalVolumes.size(); i++) {
        surfaceProbabilities.push_back(totalSurfaceProbability);
        volumeProbabilities.push_back(totalVolumeProbability);
        if (i < _physicalVolumes.size()) {
            totalSurfaceProbability+=SurfaceAreas[i]*Concentrations[i];
            totalVolumeProbability+=CubicVolumes[i]*Concentrations[i];
        }
    }
    histogramBuilt=true;
}

G4int DAMICVolumeCoordinateStore::GetRandomSurfaceGeneratorID() {
    //Note: this algorithm can be made faster/cleaner
    if (!histogramBuilt) BuildHistogram();
    G4double prob=totalSurfaceProbability*G4UniformRand();
    G4int ID=-1;
    std::vector<G4double>::size_type i=0;
    for (; i < surfaceProbabilities.size()-1; i++ ) {
        if ((surfaceProbabilities[i]<prob) && surfaceProbabilities[i+1]>prob) {
            ID=i;
            break;
        }
    }
    if (ID==-1) ID=i+1;
    return ID;
}

G4int DAMICVolumeCoordinateStore::GetRandomVolumeGeneratorID() {
    if (!histogramBuilt) BuildHistogram();
    G4double prob=totalVolumeProbability*G4UniformRand();
    G4int ID=-1;
    std::vector<G4double>::size_type i=0;
    for (; i < volumeProbabilities.size()-1; i++ ) {
        if ((volumeProbabilities[i]<prob) && volumeProbabilities[i+1]>prob) {
            ID=i;
            break;
        }
    }
    if (ID==-1) ID=i+1;
    return ID;
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointOnSurface() {
    G4int ID = GetRandomSurfaceGeneratorID();
    G4ThreeVector Pos = GetPointOnSurfaceByID(ID);
    //  G4cout << "Spawning at: " << Pos << G4endl;
    //  return GetPointOnSurfaceByID(ID);
    return Pos;
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointInVolume()
{
    G4int ID = GetRandomVolumeGeneratorID();
    G4ThreeVector Pos = GetPointInVolumeByID(ID);
    return Pos;
}

G4ThreeVector DAMICVolumeCoordinateStore::GetPointInList() {
    G4ThreeVector Pos(-1,-1,-1);
    if (!fileReader) {
        G4cout << "Error, fileReader not set, a list of positions was probably not loaded" << G4endl;
    } else {
        G4AutoLock lock(&PositionFileReaderMutex);
        Pos = fileReader->GetNextPosition();
    }
    return Pos;
}

