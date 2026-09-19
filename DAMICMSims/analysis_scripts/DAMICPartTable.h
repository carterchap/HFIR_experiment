#ifndef DAMICPartTable_h
#define DAMICPartTable_h

class DAMICPartTable {
 public:
	DAMICPartTable() {};
	~DAMICPartTable() {};

	Double_t GetRadiogenicActivity(TString iso, Double_t Texp, Double_t Trun, Double_t Tcool=51.) {
		double sat[]={230., 1800., 1650., 2100., 215., 455., 53., 125.};
		double conv=0.000001*60*60*24;
		double Thalf[]={77.236, 271, 70.83, 5.271*365, 312.13, 44.495, 83.788, 12.32*365};
		TString names[]={"56a27z", "57a27z", "58a27z", "60a27z", "54a25z", "59a26z", "46a21z", "3a1z"};
		Double_t act=-1.;
		for (int i=0; i < 8; i++ ) {
			if (iso==names[i]) {
				act=sat[i]*conv*(1.-exp(-log(2)/Thalf[i]*Texp))*exp(-log(2)/Thalf[i]*Tcool)*(1./(Trun*log(2)/Thalf[i])*(1-exp(-log(2)/Thalf[i]*Trun)));
				break;
			}
		}
		return act;
	}

	TString GetPartName(TString v) {
		//Copper Modules
		if        (v=="CopperBasePlatePV" || 
							 v=="CopperMountBarPV" ||
							 v=="CopperTopPlate44PV" ||
							 v=="CopperScrewPV") {
			return TString("Copper Modules");			
		}	else if (v=="KaptonCable11PV" ||
							 v=="FlexCable44PV" ||
							 v=="Flex2CablePV" ||
							 v=="KaptonCable2PV") {
			return TString("Flex Cable");
		} else if (v=="BottomPlatePV" || // Copper box parts
							 v=="RearPlatePV" ||
							 v=="UpperPlatePV" ||
							 v=="TopPlatePV" ||
							 v=="SidePlate1PV" ||
							 v=="SidePlate2PV" ||
							 v=="EndCoverBar1PV" ||
							 v=="EndCoverBar2PV" || 
							 v=="EndCoverPlatePV" ||
							 v=="BlockBotSupportPV" ||
							 v=="CopperOnLeadBPV") {
			return TString("Copper Box");
		} else if (v=="SensPartPV" || 
							 v=="CCDSensor44PV" ) {
			return TString("CCD");
		} else if (v=="ModScrewsPV" || 
							 v=="PanelScrewPV" ) {
			return TString("Module Screws");
		} else if (v=="LowerTubePV") {
			return TString("Vessel");
		} else if (v=="BlockBotPV" || 
							 v=="BlockTopPV" ||
							 v=="CornerBlocks1PV" ||
							 v=="CornerBlocks2PV" ||
							 v=="CornerBlocks3PV" ||
							 v=="CornerBlocks4PV" ||
							 v=="InsideLead1PV" ||
							 v=="BoxBottomLeadPV" ) {
			return TString("Chicago Lead");
		}
		return TString("");
	}
			
	TString GetIsotopeSource(TString v) {
		//primeval
		//U238 chain
		if        (v=="234a91z" ||
							 v=="234a90z" ||
							 v=="214a82z" ||
							 v=="214a83z" ||
							 v=="210a82z" ||
							 v=="210a83z" ||	
							 //Th232
							 v=="228a88z" ||
							 v=="228a89z" ||
							 v=="212a82z" ||
							 v=="212a83z" ||
							 v=="208a81z" ||
							 //Potassium / Rubidium
							 v=="40a19z" ||
							 v=="87a37z" ) {
			return TString("Primeval");
		} else if (v=="54a25z" || 
							 v=="56a27z" ||
							 v=="57a27z" ||
							 v=="58a27z" ||
							 v=="59a26z" ||
							 v=="60a27z" || 
							 v=="32a14z" ||
							 v=="32a15z" ||
							 v=="3a1z") {
			return TString("Radiogenic");
		} else {
			std::cout << "Error with: " << v << std::endl;
			return TString("ERROR");
		}
	}
	TString GetIsotopeChain(TString v, Bool_t in_eq=true) {
		//U238 chain
		TString chain;
		if        (v=="234a91z" ||
							 v=="234a90z" ) {
			chain=TString("U238");
		} else if (v=="214a82z" ||
							 v=="214a83z" ) {
			chain=TString("Ra226");
		} else if (v=="210a82z" ||
							 v=="210a83z") {
			chain=TString("Pb210");
			//Th232
		} else if (v=="228a88z" ||
							 v=="228a89z" ||
							 v=="212a82z" ||
							 v=="212a83z" ||
							 v=="208a81z" ) { 
			chain=TString("Th232");
			//Potassium
		} else if (v=="40a19z" ) {
			chain=TString("K40");
		} else if (v=="56a27z" ||
							 v=="57a27z" ||
							 v=="58a27z" ||
							 v=="60a27z" ) {
			chain=TString("Co");
		} else if (v=="54a25z" ||
							 v=="87a37z" ||
							 v=="59a26z" ) {
			chain=TString("other");
		} else if (v=="32a14z" || 
							 v=="32a15z") {
			chain=TString("Si32");
		} else {
			chain=TString("error");
		}
		if (in_eq) {
			if (chain=="Pb210" || chain=="Ra226") {
				chain=TString("U238");
			}
		}
		return chain;
	}
		
	TString GetMaskName(Int_t mask) {
		if (mask==0) {
			return TString("");
		} else if (mask==1) {
			return TString("Screws");
		}
		return TString("");
	}

	Double_t GetMass(TString v) {
		if (v=="FlexCable44PV") return 0.0143;
		if (v=="Flex2CablePV") return 0.00367;
		if (v=="KaptonCable2PV") return 0.0269;
		if (v=="KaptonCable11PV") return 0.00224;
		if (v=="BottomPlatePV") return 0.827;
		if (v=="CopperOnLeadBPV") return 0.087;
		if (v=="EndCoverBar1PV") return 0.152;
		if (v=="EndCoverBar2PV") return 0.152;
		if (v=="EndCoverPlatePV") return 0.552;
		if (v=="RearPlatePV") return 1.052;
		if (v=="SidePlate1PV") return 1.350;
		if (v=="SidePlate2PV") return 1.348;
		if (v=="BlockBotSupportPV") return 0.042;
		if (v=="TopPlatePV") return 1.30;
		if (v=="UpperPlatePV") return 0.528;
		if (v=="CopperMountBarPV") return 0.34;
		if (v=="CopperBasePlatePV") return 1.06;
		if (v=="CopperTopPlate44PV") return 0.66;
		if (v=="SensPartPV") return 0.057;
		if (v=="CCDSensor44PV") return 0.057;
		if (v=="PanelScrewPV") return 0.00035*2*8;
		if (v=="ModScrewsPV") return 0.00017*8*8;
		if (v=="BlockBotPV" || v=="BlockTopPV") return 4.755;
		if (v=="LowerTubePV") return 13.77;

		//		std::cout << "Mass not found for " << v << std::endl;
		return -0.01;
	}
	Double_t GetActivity(TString v, TString iso) {
		TString name=GetPartName(v);
		TString chain=GetIsotopeChain(iso, false);
		Double_t act=-1.0;
		//Copper module exposure time
		Double_t Texp_box= 10.*30;
		Double_t Texp_mod= 8.*30;
		Double_t Texp_pnnl= 5*.30;
		Double_t Texp_CCD=2.45*365.;
		Double_t Texp_vess=1000*365.;
		Double_t Tcool_mod=9*30;
		Double_t Tcool_box=51.;
		Double_t Tcool_vess=6*365;
		Double_t Trun=7*30.;
		if (v=="CopperScrewPV") {
		} else if (name=="Copper Modules") {
			if (chain=="U238") act=11;
			if (chain=="Ra226") act=11;
			if (chain=="Pb210") act=11;
			if (chain=="Th232") act=3.5;
			if (iso=="40a19z") act=2.7;
			if (iso=="87a37z") act=7.4;
			if (GetIsotopeSource(iso)=="Radiogenic") {
				act=GetRadiogenicActivity(iso, Texp_mod, Trun, Tcool_mod);
			}
		}	else if (name=="Copper Box" || name=="Vessel") {
			if (chain=="U238") act=11;
			if (chain=="Ra226") act=11;
			if (chain=="Pb210") act=11;
			if (chain=="Th232") act=3.5;
			if (iso=="40a19z") act=2.7;
			if (iso=="87a37z") act=7.4;
			if (GetIsotopeSource(iso)=="Radiogenic") {
				if (name=="Copper Box") 
					act=GetRadiogenicActivity(iso, Texp_box, Trun, Tcool_box);
				if (name=="Vessel")
					act=GetRadiogenicActivity(iso, Texp_vess,Trun, Tcool_vess);
			}
		} else if (name=="Flex Cable") {
			if (chain=="U238") act=7300.;
			if (chain=="Ra226") act=415.;
			if (chain=="Pb210") act=415.;
			if (chain=="Th232") act=164.;
			if (iso=="40a19z") act=10000;
			if (iso=="87a37z") act=7.4;
			if (GetIsotopeSource(iso)=="Radiogenic") {
				act=GetRadiogenicActivity(iso,Texp_box, Trun);
				//Only half of it is copper, so drop the activity by half (rather than the mass)
				act*=.5;
			}
		} else if (name=="Module Screws") {
			if (chain=="U238") act=1420.;
			if (chain=="Ra226") act=138.;
			if (chain=="Pb210") act=138;
			if (chain=="Th232") act=200.;
			if (iso=="40a19z") act=2430.;
			if (iso=="87a37z") act=7.4;
			if (GetIsotopeSource(iso)=="Radiogenic") {
				act=GetRadiogenicActivity(iso,365*10., Trun);
				//Only half of it is copper, so drop the activity by half (rather than the mass)
				act*=.5;
			}			
		}  else if (name=="CCD") {
			if (chain=="Si32") act=12.;
			else if (chain=="Pb210") act=68.;
			else if (GetIsotopeSource(iso)=="Radiogenic") {
				act=GetRadiogenicActivity(iso, Texp_CCD, Trun);
			}
		} else if (name=="Chicago Lead") {
			if (chain=="U238") act=1.2;
			if (chain=="Ra226") act=2.0;
			if (chain=="Pb210") act=285.;
			if (chain=="Th232") act=0.19;
			if (iso=="40a19z") act=0.54;
			if (iso=="124a51z") act=310*2.44e-4;
			if (iso=="110a47z") act=42*0.13;
			if (iso=="137a55z") act=5.7;
			if (iso=="60a27z") act=0.73;
		}
		if (iso=="208a81z") act*=.36;
		if (iso=="212a83z") act*=.64;
		//		std::cout << "Using activity rate:" << act << " for:" << iso << ":" <<v <<std::endl;
		if (act < 0) {
			std::cout << "Error, found no activity for iso: " << iso << " in " << name << std::endl;
			act=0;
		}
		return act;
	}

	TString GetIsotopeMeasurement(TString v, TString iso) {
		TString name=GetPartName(v);
		TString chain=GetIsotopeChain(iso, false);
		TString type="Unknown";
		if (iso=="87a37z") {
			type=TString("Upper");
		}else if (chain=="Co" || chain=="other") {
			type=TString("Activation");
		}  else if (name=="Copper Modules") {
			type=TString("Upper");
		} else if (name=="Copper Box" || name=="Vessel") {
			type=TString("Upper");
		} else if (name=="Flex Cable" ||
							 name=="Module Screws") {
			if (chain=="U238" || 
					chain=="Th232" || 
					iso=="40a19z") {
				type=TString("Measurement");
			} else if ((chain=="Ra226" || chain=="Pb210") && name=="Flex Cable") {
				type=TString("Measurement");
			}
			else if ((chain=="Ra226" || chain=="Pb210") && name=="Module Screws") {
				type=TString("Upper");
			}

		} else if (name=="CCD") {
			if (iso=="3a1z" || chain=="Si32" || chain=="Pb210") {
				type=TString("Measurement");
			} 

		} else if (name=="Chicago Lead") {
			if (chain=="Th232") {
				type=TString("Measurement");
			} else if (chain=="Pb210") {
				type=TString("Measurement");
			} else if (iso=="124a51z" || iso=="110a47z") {
				type=TString("Measurement");
			} else {
				type=TString("Upper");
			}
		}
		if (type=="Unknown") {
			std::cout << "iso " << iso << " in volume " << v << " unknown measurement" << std::endl;
		}
		return type;
	}

	TString GetProposalCategory(TString v, TString iso) {
		TString name=GetPartName(v);
		if (name=="CCD" && (iso=="210a82" || iso=="210a83" || iso=="32a14z" || iso=="32a15z")) {
			return TString("CCD Surface");
		} else if (name=="CCD") {
			return TString("CCD Other");
		} else if (GetIsotopeChain(iso)=="Co") {
			return TString("Cobalt Activation");
		} else if (GetIsotopeMeasurement(v,iso)=="Upper") {
			return TString("Upper Limit");
		} else if (GetIsotopeMeasurement(v,iso)=="Measurement" || GetIsotopeChain(iso)=="other") {
			return TString("Central Limit");
		}
		else {
			std::cout << "Unknown category, volume: " << v << " iso: " << iso <<std::endl;
			return TString("Unknown");
		}
	}

};
#endif
