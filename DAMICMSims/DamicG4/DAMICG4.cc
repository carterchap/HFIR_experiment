#ifndef G4VERSION_NUMBER
#include "G4Version.hh"
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"
#include "QBBC.hh"

#include "DAMICDetectorConstruction.hh"
#include "DAMICPhysicsListLivermore.hh"
//#include "DAMICPhysicsList_test.hh"
#include "DAMICActionInitialization.hh"
#include "DAMICAnalysisManager.hh"
#include "G4PhysListFactory.hh"

//To parse command line options
#include "tclap/CmdLine.h"
#include <string>
#include <vector>

#include "G4VisExecutive.hh"

#include "G4UIQt.hh"
#include "G4UIExecutive.hh"


#include "G4GDMLParser.hh"

int main(int argc, char** argv)
{


   	if(argc<2)
   	{
	  G4cout << G4endl  << "   Error! InputFile and/or gdml files are not specified!" << G4endl;
	G4cout << "   Usage: " << argv[0] << " gdml_det input.mac  (<input_gdml_file:mandatory>)" << G4endl;
	G4cout << G4endl;
      	//return -1;
	exit(0);
   	}

   	//G4GDMLParser parser;
   	//parser.Read(argv[1]);

	bool batchMode=false;
	G4double milliCharge=1.0;
	G4double milliMass=0.105;
	G4String macroName="";
	G4double sim_version=1.1;
	G4String physicsListName="Def";

	try
	{
		/*
		For help on adding new commands, see: tclap.sourceforge.net/manual.html
	 	*/
		TCLAP::CmdLine cmd("DAMIC100 Geant simulation", ' ', std::to_string(sim_version).c_str());
		TCLAP::UnlabeledMultiArg<std::string> posArgsArg("args", "[GDMLfile] [macroFile] [physics list {Def, Pen, Liv, Mic, Shielding_LIV}]", false,"string", cmd);

		//You can only specify one optional unrequired positional argument in TCLAP, so this doesn't work, hence above
		//		TCLAP::UnlabeledValueArg<std::string> macroNameArg("macro", "Macro file to run", false, macroName, "string", cmd);
		//		TCLAP::UnlabeledValueArg<int> numThreadsArg("threads", "Number of threads to use", false, nofThreads, "integer",cmd);
		//		TCLAP::UnlabeledValueArg<std::string> physicsListNameArg("physicsList", "Physics list to use, from one of: Def, Pen, Liv, Micro", false, physicsListName, "string", cmd);

		TCLAP::ValueArg<G4double> milliChargeArg("c","charge", "Millicharged particle charge", false, milliCharge, "double", cmd);
		TCLAP::ValueArg<G4double> milliMassArg("m", "mass", "Millicharged particle mass (in MeV)", false, milliMass, "double", cmd);
		cmd.parse(argc, argv);
		//Unlabeled args passed in by position
		std::vector<std::string> posArgs=posArgsArg.getValue();

		if(posArgs.size()>1)
		{
			macroName=posArgs[1];
		}

		if(posArgs.size()>2)
		{
			physicsListName=posArgs[2];
		}

		milliCharge=milliChargeArg.getValue();
		milliMass=milliMassArg.getValue();
	} catch(TCLAP::ArgException &e)
	{
		G4cerr << "Error: " << e.error() << " for arg " << e.argId() << G4endl;
		// exit(-1);
	};

	if(macroName!="")
	{
        batchMode=true;
	}

	// Prepare and initialize Random engine
	G4Random::setTheEngine(new CLHEP::HepJamesRandom);
	
	// Prepare and initialze run manager
	G4RunManager* runManager = new G4RunManager;

   	// Read and parsing geometry from gdml file
   	G4GDMLParser parser;
   	parser.Read(argv[1]);

	// Read some parameter from macro file
	auto dman = DAMICAnalysisManager::GetInstance();
	dman->ParseMacro(macroName);

	// Initilize detector construction
	DAMICDetectorConstruction* det= new DAMICDetectorConstruction(parser);
	runManager->SetUserInitialization(det);

	// get the pointer to the User Interface manager
	G4UImanager* UI = G4UImanager::GetUIpointer();

	/*======================Set up physics list (from command line, if given)===========================*/
	G4PhysListFactory physListFactory; 	
	G4VUserPhysicsList* physicsList = nullptr;

	
	if (physicsListName == G4String("Shielding_LIV")) 
	{
    	G4cout << " you are using Shielding_LIV" << G4endl;
		physicsList = physListFactory.GetReferencePhysList(physicsListName);
	} 
	else if (physicsListName == G4String("Liv")) 
	{
		G4cout << " you are using default DAMIC physics list: Livermore" << G4endl;
		physicsList = new DAMICPhysicsListLivermore();
	}
	else
	{
		G4cout << " you are using default DAMIC physics list: Livermore" << G4endl;
		physicsList = new DAMICPhysicsListLivermore();
	}

	physicsList->SetVerboseLevel(1);
    
	runManager->SetUserInitialization(physicsList);    

	// ---- Set lowest kinetic energy for e-/e+ (10 eV) via EM parameters
	auto em = G4EmParameters::Instance();
	em->SetLowestElectronEnergy(10*CLHEP::eV);

	// ---- Production cut limit in energy
	G4double cutEnergy = 10*CLHEP::eV;
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(cutEnergy, 100*CLHEP::eV);
	
	runManager->SetUserInitialization(new DAMICActionInitialization(det));
	
	//Initialize G4 kernel
	runManager->Initialize();

	if(batchMode)// batch mode
	{
		G4String command = "/control/execute ";
		UI->ApplyCommand(command+macroName);
	}
	else //define visualization and UI terminal for interactive mode
	{

		G4VisManager* visManager = new G4VisExecutive;
		visManager->Initialize();

		G4UIExecutive * ui = new G4UIExecutive(argc,argv);
		//G4UIQt *ui = new G4UIQt(argc,argv);
		//UI->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
		delete ui;

		//Close-out analysis:
		//Save histograms

		delete visManager;

	}

	if(physicsList!=nullptr)
	{
		delete physicsList;
		physicsList=nullptr;
	}
	if(det!=nullptr)
	{
		delete det;
		det=nullptr;
	}
	//if(runManager != nullptr)
	//{
	//	delete runManager;
	//	runManager = nullptr;
	//}
	return 0;
}
