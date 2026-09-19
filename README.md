DAMICM_G4Sims
=============
This software is intended to be used for DAMIC-M simulations (up to now).  
The code was developed in Geant4 10.4. It was later adapted for Geant4 **v.11.1.1**. See validation plots [here](https://gev.uchicago.edu/cgi-bin/DocDB/ShowDocument?docid=1259).
A Geometry Definition Markup Language (GDML) is used to describe the geometry 
of the detector. 
[GDML](http://lcgapp.cern.ch/project/simu/framework/GDML/doc/GDMLmanual.pdf 
"GDML user's guide") is an application-indepedent geometry description format
based on XML. 

[Wiki](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/wikis/home "Wiki home")


## Writting code Guidelines

- Indentation style [Allman](https://en.wikipedia.org/wiki/Indentation_style) 
- (suggested number of indented spaces 4):
```
while (x == y)
{
    something();
    somethingelse();
}
```
- Naming conventions:
    - **private data member** name starts with an underscore followed by lower case, ex. ***_ccdSensor***
    - **private function** name starts with an underscore followed by upper case, ex. ***_FillTTree()***
    - **public funciton** name starts with upper case, ex. ***ConstructDetector()***

 
## Test before any commit/simulation

To detect some bugs on the implemented algorithms (and/or geometry implementation) 
some test should be done for any change on the GDML geometry file. 


Here a list of some usefull plots to check whether the primary generator has been 
simulated properly. 

#### Check the position of the primary generator

Its validation depends on the way the position of the primary generator was
defined. There is four ways to define its position:
    1.1 from a point source: ***damic/gun/position/dosource***
    1.2 following a shape: ***damic/gun/position/doshape***
    1.3 from a specific material:  ***damic/gun/position/domaterial***
    1.4 from different parts of the damic geometry (different volumes):  ***damic/gun/position/dovolume***

For instance, if the method 1.4 (dovolume) is used to generate thr position of 
the primary generator, the position plot of the generated events should draw 
the shape of the simulated volume (for instance a box in the case of example1, see
line ***damic/gun/position/addvolume CopperBox_PV 1*** in Example 1 below). 

Note that the geometry is replicated on the Z-axes (check the number of CCDs 
simulated on that direction)

Some example using the DAMICG4 root output file,

##### 3D position plot for easy-viewer geometries
    

```
$ root -l output_DAMICG4_root_file.root

root [1] EventOut->Draw("posx:posy:posz")
```

##### or frontal and transversal plot for more complicated ones
    

```
$ root -l output_DAMICG4_root_file.root
    
root [1] EventOut->Draw("posx:posy")
root [2] EventOut->Draw("posx:posz")
```
    
With the distribution of the X, Y and Z axes the radial profile can also be checked.

#### Check the distribution of the CCD number.

```
$ root -l output_DAMICG4_root_file.root

root [1] CCDOut->Draw("CCDid")
```

It will be also useful to plot several output parameters related with the sensitive 
detector, such as

#### Check the energy deposit

#### Check the position of the events generated in the CCD Sensor part

Note the stored positions are the LOCAL POSITIONS, not the global one, so you will
only see the shape of one CCD.


## How to run the simulations
Find instruction on [how to install](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/wikis/HowTo 
"building DAMICM_G4Sims") the DAMICM_G4Sims package.

```
#!sh
export PATH=${PATH}:$DAMICM_RUN_DIR
```

### Interactive mode
Run the command with the 
[main GDML input file](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/blob/DAMICMSims_GDML_geometry_modularized/DAMICMSims/DamicG4/gdml/detector_main.gdml 
"current DAMIC-M geometry") as the only input argument
```
$ DAMICG4 gdml/detector_main.gdml
```
This will open the 
[Geant4 Visulization display](http://geant4.slac.stanford.edu/Presentations/vis/G4OpenGLTutorial/G4OpenGLTutorial.html 
"OpenGL Event Display"). To load the geometry run the following commands in 
the ***session*** command line
```
/vis/open OGL
/control/execute vis.mac
```
[Here](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/blob/DAMICMSims_GDML_geometry/DAMICMSims/DamicG4/vis.mac) 
an example of the ***vis.mac*** file. For extra options check the 
[G4 visualization commands](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/AllResources/Control/UIcommands/_vis_.html). 


### Bash mode
For a batch mode, run ***DAMICG4*** as follows
```
$ DAMICG4 gdml/detector_main.gdml mac.mac
```
where [detector_main.gdml](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/blob/DAMICMSims_GDML_geometry_modularized/DAMICMSims/DamicG4/gdml/detector_main.gdml) 
is the main GDML file that contains the geometry modules to be simulated.  The macro file,
***input.mac***, is a file that contains a list of commands. An example of a macro file
is shown below, ***example 1***, and can be found 
[here](https://gitlab.in2p3.fr/damicm/DAMICM_G4Sims/blob/DAMICMSims_GDML_geometry_modularized/DAMICMSims/DamicG4/macfDAM_example1.mac).
Indeed, you can re-execute **DAMICG4** with different run conditions without recompiling 
anything, just using different macro input files.


### Macrofile Skeleton
For the input macro file there is only two mandatory commands:
```
/run/initialize
/run/beamOn N
```
N is the number of primary particels that you want to simulate. For other 
Geant4 commands take a look at the [Geant4 documentation]. Do not use the Fun Particel 
commands of Geant4 they will not work.


#### Example 1: primary generator defined as a volume of the 
DAMIC-M geometry

```
/control/verbose 0
/run/verbose 0
/run/initialize
/tracking/verbose 0
/event/verbose 0
/random/setSeeds 30991 7750
/damic/gun/particle ion
/damic/gun/ion 27 60 0 0
/damic/gun/energy/mono 0 eV
/damic/gun/direction/oned
/damic/gun/direction/onedX 1
/damic/gun/direction/onedY 0
/damic/gun/direction/onedZ 0
/damic/gun/position/dovolume
/damic/gun/position/addvolume CopperBox_PV 1
/analysis/setFileName CopperBox_PV60a27z_Seeds30991_7750
/run/beamOn 20000
```



