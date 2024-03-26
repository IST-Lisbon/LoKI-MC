#ifndef __GeneralDefinitions__
#define __GeneralDefinitions__

#include <string>
#include <map>
#include <vector>
#include "External/eigen-3.4.0/Eigen/Dense"
#include "LoKI-MC/Headers/Constant.h"

/* 
This file defines variables that are used in multiple places along the code, to avoid multiple definitions
*/

namespace GeneralDefinitions{
	// ----- structures used in 'PrescribedEedf', 'BoltzmannMC', 'BoltzmannMCTimeDep', 'Output' ----- //

	struct PowerStruct{
		std::map<std::string,double> Map;
		std::map<std::string, std::map<std::string,double>> gasesMap;
	};

	struct RateCoeffStruct{
		int collID = -1;
		double ineRate = Constant::NON_DEF;
		double supRate = Constant::NON_DEF;
		double ineRateMC = Constant::NON_DEF;
		double supRateMC = Constant::NON_DEF;	
		std::string collDescription;
	};

	// ----- variables used in 'BoltzmannMC' and 'BoltzmannMCTimeDep' ----- //

	// definition of the process types
	static int conservativeType = 0;
	static int ionizationType = 1;
	static int attachmentType = 2;

	// definition of special "process" IDs
	static int nullCollisionID = -1;
	static int partialFreeFlightID = -2;

	// definition of types of energy sharing in ionization events
	static int equalSharingIonizType = 0;
	static int oneTakesAllIonizType = 1;
	static int usingSDCSIonizType = 2;
	static int randomUniformIonizType = 3;

	// definition of types of gas-temperature effect
	static int falseGasTempEffectID = 0;
	static int trueGasTempEffectID = 1;
	static int smartActivationGasTempEffectID = 2;
};	

#endif