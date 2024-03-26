#ifndef __Setup__
#define __Setup__

#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Gas.h"
#include "LoKI-MC/Headers/State.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/Output.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Output.h"
#include "LoKI-MC/Headers/StatePropertyFunctions.h"
#include "LoKI-MC/Headers/GasPropertyFunctions.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/GUI.h"
#include "LoKI-MC/Headers/AngularDistributionFunctions.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include <omp.h>
#include <string>
#include <limits>
#include <numeric>
#include <vector>
#include <map>
#include <typeinfo>
#include <fstream>
#include <algorithm>
#include <time.h>
#include <chrono>

class Collision;

template <class ElectronKineticsType>
class Setup{
public:

	std::string fileName;

	std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

	/************* Configuration properties ****************/

	bool enableGui = false;				   // determines if the gui must be activated
	bool enableOutput = false;			   // determines if the output must be activated
	bool enableElectronKinetics = false;   // determines if the electron kinetics module must be activated

	int numberOfJobs = 1; // total number of jobs to be run
	int currentJobID = 0; // ID of the job currently running

	/************* Objects of the simulation ****************/

	WorkingConditions* workCond;         				   		//  
	GUI<ElectronKineticsType>* gui = NULL;		                //  
	Output<ElectronKineticsType>* output = NULL;	            // ->  General objects of the simulation

	ElectronKineticsType* electronKinetics;				        // (ElectronKineticsType can be BoltzmannMC or PrescribedEedf)
	std::vector<EedfGas*> electronKineticsGasArray;			   	// 
	std::vector<EedfState*> electronKineticsStateArray;		   	// -> objects related with the electron kinetics module
	std::vector<Collision*> electronKineticsCollisionArray;	    // 
	Grid* energyGrid = NULL;							        //

	// ----- class methods ----- //
	
	Setup(std::string fileName1){

		// save the initial time
		start = std::chrono::high_resolution_clock::now();
			
		// save the setup filename
		fileName = fileName1;

		// perform a diagnostic of the correctness of the configuration provided by the user
		selfDiagnostic();

		// check which modules of the simulation must be enabled
		// GUI
		if ( FieldInfo::getField("gui.isOn") ){
			enableGui = Parse::str2bool( FieldInfo::getFieldValue("gui.isOn") );
		}
		// Output
		if ( FieldInfo::getField("output.isOn") ){
			enableOutput = Parse::str2bool( FieldInfo::getFieldValue("output.isOn") );
		}
		// Electron kinetics
		if ( FieldInfo::getField("electronKinetics.isOn") ){
			enableElectronKinetics = Parse::str2bool( FieldInfo::getFieldValue("electronKinetics.isOn") );
		}
	}

	void initializeSimulation(){
		// initializeSimulation creates the objects necessary to run the simulation specified in a particular setup.

		// ----- INITIAL MESSAGE -----
		std::cout<<"Starting simulation..."<<std::endl;

		if ( Parse::str2bool(FieldInfo::getFieldValue("gui.isOn")) ){
			for (auto option: FieldInfo::getFieldChildNames("gui.terminalDisp")){
				if (option == "setup"){		
					// printing the setup info
					FieldInfo::printSetupInfo();
					std::cout<<std::endl;
					break;
				}
			}
		}

		// ----- SETTING UP THE WORKING CONDITIONS OF THE SIMULATION -----
		workCond = new WorkingConditions( FieldInfo::getFieldMap("workingConditions") );


		// ----- SETTING UP THE ELECTRON KINETICS -----
		if (enableElectronKinetics){
			// setting up the gas mixture that is going to be used to solve the electron kinetics
			createElectronKineticsGasMixture();

			/*std::cout<<"----- PRINTING GASES -----\n\n";
			double temp;
			for (auto& gas: electronKineticsGasArray){
				gas->disp();
				std::cout<<"-------------------------------------"<<std::endl;
			}*/

			/*std::cout<<"----- PRINTING STATES -----\n\n";
			for (auto& state: electronKineticsStateArray){
				state->disp();
				std::cout<<"-------------------------------------"<<std::endl;
			}*/

			/*std::cout<<"----- PRINTING COLLISIONS -----\n\n";
			for (auto& collision: electronKineticsCollisionArray){
				std::cout<<collision->description()<<std::endl;
				std::cout<<"  ID: "<<collision->ID<<std::endl;
				std::cout<<"  isExtra: "<<collision->isExtra<<std::endl;
				std::cout<<"  isReverse: "<<collision->isReverse<<std::endl;
				std::cout<<"  threshold: "<<collision->threshold<<std::endl;
				std::cout<<"-------------------------------------"<<std::endl;
				int crossSectionSize = collision->rawCrossSection[0].size();
				for (int i = 0; i < crossSectionSize; ++i){
					std::cout<<collision->rawCrossSection[0][i]<<" "<<collision->rawCrossSection[1][i]<<std::endl;
				}
				std::cout<<"  angularScatteringType: "<<collision->angularScatteringType<<std::endl;
				std::cout<<"  angularScatteringParams: ";
				for (auto param: collision->angularScatteringParams){
					std::cout<<param<<", ";
				}
				std::cout<<"-------------------------------------"<<std::endl;
				int crossSectionSize = collision->rawCrossSection[0].size();
				for (int i = 0; i < crossSectionSize; ++i){
					std::cout<<collision->rawCrossSection[0][i]<<" "<<collision->rawCrossSection[2][i]<<std::endl;
				}				
				std::printf("\n\n");
			}
			std::getchar();*/

			// create energy grid object and save it for later use
			energyGrid = new Grid();

			// adjusting (and linking) the collision's cross sections to the energy grid
			Collision::adjustToEnergyGrid(energyGrid, electronKineticsCollisionArray);

			electronKinetics = new ElectronKineticsType(this);
		}

		// ----- SETTING UP JOBS INFORMATION -----
		if (workCond->variableCondition == "electronTemperature"){
			numberOfJobs = workCond->electronTemperatureArray.size();
		}
		else if (workCond->variableCondition == "reducedElecField"){
			numberOfJobs = workCond->reducedElecFieldArray.size();
		}
		else if (workCond->variableCondition == "reducedMagField"){
			numberOfJobs = workCond->reducedMagFieldArray.size();
		}
		else if (workCond->variableCondition == "elecFieldAngle"){
			numberOfJobs = workCond->elecFieldAngleArray.size();
		}
		else if (workCond->variableCondition == "excitationFrequency"){
			numberOfJobs = workCond->excitationFrequencyArray.size();
		}

		// ----- SETTING UP THE GUI AND OUTPUT OF THE SIMULATION -----
		if (enableGui){
			gui = new GUI<ElectronKineticsType>(this);
		}
		if (enableOutput){
			output = new Output<ElectronKineticsType>(this);
		}
	}

	void nextJob(){
		
		static std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");

		// increase the current job ID
		++currentJobID;

		// set up next job (if there is)
		if (currentJobID < numberOfJobs){
			std::vector<std::string> propertiesToUpdate;
			std::vector<double> newValues;

			if (workCond->variableCondition == "electronTemperature"){
				propertiesToUpdate = {"electronTemperature"};
				newValues = {workCond->electronTemperatureArray[currentJobID]};		
			}
			else if (workCond->variableCondition == "reducedElecField"){
				propertiesToUpdate = {"reducedElecField"};
				newValues = {workCond->reducedElecFieldArray[currentJobID]};
			}
			else if (workCond->variableCondition == "reducedMagField"){
				propertiesToUpdate = {"reducedMagField"};
				newValues = {workCond->reducedMagFieldArray[currentJobID]};	
			}
			else if (workCond->variableCondition == "elecFieldAngle"){
				propertiesToUpdate = {"elecFieldAngle"};
				newValues = {workCond->elecFieldAngleArray[currentJobID]};	
			}
			else if (workCond->variableCondition == "excitationFrequency"){
				propertiesToUpdate = {"excitationFrequency"};
				newValues = {workCond->excitationFrequencyArray[currentJobID]};	
			}				
			// set properties for the next job
			workCond->update(propertiesToUpdate, newValues);
		}
	}

	void createElectronKineticsGasMixture(){
	    // createElectronKineticsGasMixture creates the gas mixture that is going to be considered to solve the electron
	    // kinetics. This includes: a gasArray with all the information related to the gas mixture, a stateArray with all
	    // the states of the different gases and a collisionArray with all the electron collisions of all gases.
	    //
	    // This is done for a particular setup specified in the input file.

	    // create arrays of gases, states and collisions with the information of the LXCat files
	    LXCatData();

	    // fill properties of gases with the information of the input file
	    gasProperties("electronKinetics", electronKineticsGasArray);

	    // fill properties of gases with the information of the input file
	    stateProperties("electronKinetics", electronKineticsStateArray);

	    // evaluate densities
	    for (auto& state: electronKineticsStateArray){
	    	state->evaluateDensity();
	    }

	    // check for non standard populations in effective cross sections

	    // if there are non standard populations, initialize the vectors of all gases
	    std::map<std::string,double> effectiveCrossSectionPopulationsMap = FieldInfo::getFieldNumericMap("electronKinetics.effectiveCrossSectionPopulations");
	    if (!effectiveCrossSectionPopulationsMap.empty()){
	    	for (auto& gas: electronKineticsGasArray){
	    		gas->effectivePopulations.assign(gas->stateArray.size(), 0);
	    	}
	    }
	    for (auto stateString: Parse::getMapKeys(effectiveCrossSectionPopulationsMap)){
	    	double value = effectiveCrossSectionPopulationsMap[stateString];
	    	Parse::rawStateStruct rawState = Parse::getRawState(stateString);
	    	std::vector<int> stateIDs = EedfState::find(rawState.gasName, rawState.ionCharg, rawState.eleLevel, rawState.vibLevel, rawState.rotLevel, electronKineticsStateArray);
	    	if (stateIDs[0] == -1){
	    		continue;
	    	}
	    	for (auto stateID: stateIDs){
	    		electronKineticsStateArray[stateID]->gas->effectivePopulations[stateID] = value;
	    	}
	    }


	    // assign angular scattering types/models
	    assignAngularScatteringTypes();

	    // different checks for each gas in gasArray
	    for (auto& gas: electronKineticsGasArray){
	    	// avoid dummy gases (gases created for extra cross sections or for the sake of a pretty output)
	    	if (gas->collisionArray.empty()){
	    		continue;
	    	}
	    	// check for each (non dummy) gas to have its mass defined
	    	if (gas->mass == Constant::NON_DEF){
	    		Message::error(std::string("Mass of gas ") + gas->name + " not found.\nNeeded for the evaluation of the elastic collision operator (Boltzmann).\nCheck input file.");
	    	}
	    	// check for the distribution of states to be properly normalized
	    	gas->checkPopulationNorms();
	    	// check for an Elastic collision to be defined, for each electronic state with population
	    	gas->checkElasticCollisions(electronKineticsCollisionArray);
	    }
	}

	void LXCatData(){
		// LXCatData parses the LXCat files (regular or extra) included in the setup and creates the corresponding arrays of gases, states and 
		// collisions. Those objects are created only with the information available in the LXCat files, so their properties (like, masses, 
		// populations, etc.) must be filled later on with the information available in the input file.
		
		// parse LXCat files
		std::vector<Parse::LXCatEntryStruct> LXCatEntryArray = Parse::LXCatFiles( FieldInfo::getFieldChildNames("electronKinetics.LXCatFiles") );


		// create gases, states and collisions from the LXCat parsed info
		int gasID, targetID, productID;

		for (auto& LXCatEntry : LXCatEntryArray){
			gasID = EedfGas::add(LXCatEntry.target[0].gasName, electronKineticsGasArray);
			targetID = EedfState::add(electronKineticsGasArray[gasID], LXCatEntry.target[0].ionCharg, LXCatEntry.target[0].eleLevel, LXCatEntry.target[0].vibLevel, 
				       				  LXCatEntry.target[0].rotLevel, electronKineticsStateArray);
			EedfState* target = electronKineticsStateArray[targetID];
			
			std::vector<EedfState*> productArray;
			for (auto& product: LXCatEntry.productArray){
				gasID = EedfGas::add(product.gasName, electronKineticsGasArray);
				productID = EedfState::add(electronKineticsGasArray[gasID], product.ionCharg, product.eleLevel, product.vibLevel, 
				       				  	   product.rotLevel, electronKineticsStateArray);
				productArray.push_back(electronKineticsStateArray[productID]);
			}

			Collision::add(LXCatEntry.type, target, productArray, LXCatEntry.productStoiCoeff, LXCatEntry.isReverse, LXCatEntry.threshold,
			               LXCatEntry.rawIntegralCrossSection, LXCatEntry.rawMomTransfCrossSection, electronKineticsCollisionArray, false);		
		}

		// create "extra" collisions (and corresponding gases/states) in case they are specified in the setup file
		if (FieldInfo::getField("electronKinetics.LXCatFilesExtra")){
			// parse LXCat extra files
			LXCatEntryArray = Parse::LXCatFiles( FieldInfo::getFieldChildNames("electronKinetics.LXCatFilesExtra") );

			// create gases, states and collisions from the LXCat parsed info
			for (auto& LXCatEntry : LXCatEntryArray){
				gasID = EedfGas::add(LXCatEntry.target[0].gasName, electronKineticsGasArray);
				targetID = EedfState::add(electronKineticsGasArray[gasID], LXCatEntry.target[0].ionCharg, LXCatEntry.target[0].eleLevel, LXCatEntry.target[0].vibLevel, 
					       				  LXCatEntry.target[0].rotLevel, electronKineticsStateArray);
				EedfState* target = electronKineticsStateArray[targetID];
				
				std::vector<EedfState*> productArray;
				for (auto& product: LXCatEntry.productArray){
					gasID = EedfGas::add(product.gasName, electronKineticsGasArray);
					productID = EedfState::add(electronKineticsGasArray[gasID], product.ionCharg, product.eleLevel, product.vibLevel, 
					       				  	   product.rotLevel, electronKineticsStateArray);
					productArray.push_back(electronKineticsStateArray[productID]);
				}

				Collision::add(LXCatEntry.type, target, productArray, LXCatEntry.productStoiCoeff, LXCatEntry.isReverse, LXCatEntry.threshold,
				               LXCatEntry.rawIntegralCrossSection, LXCatEntry.rawMomTransfCrossSection, electronKineticsCollisionArray, true);		
			}
		}

		// create states needed to fix orphan states
		EedfState::fixOrphanStates(electronKineticsStateArray);
	}

	void assignAngularScatteringTypes(){

		double angleNumber = 0;

		if (FieldInfo::getField("electronKinetics.anisotropicScattering") && 
			Parse::str2bool(FieldInfo::getFieldValue("electronKinetics.anisotropicScattering.isOn"))){

			angleNumber = FieldInfo::getFieldNumericValue("electronKinetics.anisotropicScattering.angleNumber");

			std::vector<std::string> collisionStrings;
			// join all collision strings
			for (auto str: FieldInfo::getFieldChildNames("electronKinetics.anisotropicScattering.collisions")){
				// if it is a collision string 
				if (str.find(';') != std::string::npos){
					collisionStrings.push_back(str);
				}
				// else, it is a file name for several collision strings
				else{
					collisionStrings = MathFunctions::append(collisionStrings, Parse::readFileStrings(std::string("Input/") + str));
				}			
			}		

			// get all collisions in the setup file
			for (auto collisionString: collisionStrings){

				// separate information
				std::vector<std::string> stringPortions = Parse::tokenizeCharacters(collisionString,(char*)";");

				// get the grouping string: "group" or "single"
				std::string groupingString = stringPortions[0];

				bool collisionFound = false;
				if (groupingString == "group"){
					if (stringPortions.size() < 4){
						Message::error(std::string("Error while reading the following setup line of electronKinetics.anisotropicScattering.collisions:\n") + collisionString);
					}
					// get the different components
					std::string gasName = stringPortions[1];
					std::string collisionType = stringPortions[2];
					std::string angularScatteringType = stringPortions[3];
					std::vector<double> angularScatteringParams;
					if (stringPortions.size() > 4){
						for (auto strParam: Parse::tokenizeCharacters(stringPortions[4],(char*)",")){
							angularScatteringParams.push_back(Parse::str2value(strParam));
						}	
					}	
					// find the gas and the corresponding collisions
					int gasIndex = EedfGas::find(gasName,electronKineticsGasArray);
					if (gasIndex != -1){
						EedfGas* gas = electronKineticsGasArray[gasIndex];
						for (auto& collision: gas->collisionArray){
							if (collision->type == collisionType){
								collision->angularScatteringType = angularScatteringType;
								collision->angularScatteringParams = angularScatteringParams;
								collisionFound = true;
							}
						}
					}	
				}

				else if (groupingString == "single"){
					if (stringPortions.size() < 3){
						Message::error(std::string("Error while reading the following setup line of electronKinetics.anisotropicScattering.collisions:\n") + collisionString);
					}
					// get the different components
					std::string collisionDescription = stringPortions[1];
					std::string angularScatteringType = stringPortions[2];
					std::vector<double> angularScatteringParams;
					if (stringPortions.size() > 3){
						for (auto strParam: Parse::tokenizeCharacters(stringPortions[3],(char*)",")){
							angularScatteringParams.push_back(Parse::str2value(strParam));
						}	
					}
					// find the collision through the description
					for (auto& collision: electronKineticsCollisionArray){
						if (collision->description() == collisionDescription){
							collision->angularScatteringType = angularScatteringType;
							collision->angularScatteringParams = angularScatteringParams;
							collisionFound = true;
							break;							
						}
					}
				}

				else{
					Message::error(std::string("Error while reading the following setup line of electronKinetics.anisotropicScattering.collisions:\n") + collisionString);
				}

				if (!collisionFound){
					Message::error(std::string("Error while reading the following setup line of electronKinetics.anisotropicScattering.collisions:\n") + collisionString + 
						"\nThe collision (group or single) was not found.");
				}				
			}
		}

		// assign the angular scattering functions of all collisions (including isotropic ones)
		// evaluate the raw momentum-transfer cross-section
		// assure that all "Effective" an "Attachment" collisions are isotropic
		for (auto& gas: electronKineticsGasArray){
			for (auto& collision: gas->collisionArray){
				collision->angularDistributionFunction = AngularDistributionFunctions::functionMap(collision->angularScatteringType,collision->angularScatteringParams);
				collision->evaluateRawMomTransfCrossSection(angleNumber);
				if (collision->angularScatteringType != "isotropic" && (collision->type == "Effective" || collision->type == "Attachment")){
					Message::error(std::string("Error in the following collision:\n") + collision->description() + "\n'" + 
					collision->type + "' collisions cannot have a user-prescribed angular scattering model.");
				}
			}
		}	
	}

	template <class GasType>
	void gasProperties(std::string module, std::vector<GasType*> &gasArray){
		// the following properties are set from the setup file: mass, fraction, harmonic frequency, anharmonic frequency,
		// rotational constant, electric quadrupole moment, OPB parameter, Lennard Jones Depth and Distance

		setGasProperty("mass", FieldInfo::getFieldMap(module + ".gasProperties.mass"), gasArray);
		setGasProperty("harmonicFrequency", FieldInfo::getFieldMap(module + ".gasProperties.harmonicFrequency"), gasArray);
		setGasProperty("anharmonicFrequency", FieldInfo::getFieldMap(module + ".gasProperties.anharmonicFrequency"), gasArray);
		setGasProperty("rotationalConstant", FieldInfo::getFieldMap(module + ".gasProperties.rotationalConstant"), gasArray);
		setGasProperty("lennardJonesDistance", FieldInfo::getFieldMap(module + ".gasProperties.lennardJonesDistance"), gasArray);
		setGasProperty("lennardJonesDepth", FieldInfo::getFieldMap(module + ".gasProperties.lennardJonesDepth"), gasArray);
		setGasProperty("electricDipolarMoment", FieldInfo::getFieldMap(module + ".gasProperties.electricDipolarMoment"), gasArray);
		setGasProperty("electricQuadrupoleMoment", FieldInfo::getFieldMap(module + ".gasProperties.electricQuadrupoleMoment"), gasArray);
		setGasProperty("polarizability", FieldInfo::getFieldMap(module + ".gasProperties.polarizability"), gasArray);
		setGasProperty("fraction", FieldInfo::getFieldMap(module + ".gasProperties.fraction"), gasArray);
		setGasProperty("heatCapacity", FieldInfo::getFieldMap(module + ".gasProperties.heatCapacity"), gasArray);
		setGasProperty("thermalConductivity", FieldInfo::getFieldMap(module + ".gasProperties.thermalConductivity"), gasArray);
		setGasProperty("OPBParameter", FieldInfo::getFieldMap(module + ".gasProperties.OPBParameter"), gasArray);

		// check for the proper normalization of the gas fractions
		GasType::checkFractionNorm(gasArray);
	}

	template <class GasType>
	void setGasProperty(std::string property, std::map<std::string,std::string> propertyMap, std::vector<GasType*> &moduleGasArray){
		std::string propertyString;
		std::vector<std::string> gasNames;
		std::vector<double> argumentArray;
		std::vector<std::string> argumentStrArray;

		gasNames = Parse::getMapKeys(propertyMap); 
		for(auto& gasName: gasNames){

			int gasID = GasType::find(gasName, moduleGasArray);
			if (gasID == -1){
				continue;
			}

			propertyString = propertyMap[gasName];
			
			if (propertyString.find("@") != std::string::npos){
				findPropertyFunctionArguments(propertyString, argumentArray, argumentStrArray); 
			}

			GasPropertyFunctions::evaluateGasPropertyFunction(propertyString, moduleGasArray[gasID], property, argumentArray, argumentStrArray, workCond);
			
			argumentArray.clear();
			argumentStrArray.clear();
		}	
	}	

	template <class StateType>
	void stateProperties(std::string module, std::vector<StateType*> &stateArray){

		setStateProperty("energy", FieldInfo::getFieldMap(module + ".stateProperties.energy"), stateArray);
		setStateProperty("statisticalWeight", FieldInfo::getFieldMap(module + ".stateProperties.statisticalWeight"), stateArray);
		setStateProperty("reducedDiffCoeff", FieldInfo::getFieldMap(module + ".stateProperties.reducedDiffCoeff"), stateArray);
		setStateProperty("reducedMobility", FieldInfo::getFieldMap(module + ".stateProperties.reducedMobility"), stateArray);
		setStateProperty("population", FieldInfo::getFieldMap(module + ".stateProperties.population"), stateArray);
	}

	template <class StateType>
	void setStateProperty(std::string property, std::map<std::string,std::string> propertyMap, std::vector<StateType*> &moduleStateArray){
		std::string propertyString;
		std::vector<std::string> stateNames;
		std::vector<StateType*> stateArray;
		Parse::rawStateStruct rawState;
		std::vector<double> argumentArray;
		std::vector<std::string> argumentStrArray;

		stateNames = Parse::getMapKeys(propertyMap); 
		for(auto& stateName: stateNames){
			rawState = Parse::getRawState(stateName);
			stateArray = StateType::findPointer(rawState.gasName, rawState.ionCharg, rawState.eleLevel, rawState.vibLevel, rawState.rotLevel, moduleStateArray);
			if (stateArray.empty()){
				continue;
				//Message::error(std::string("Trying to assign the property ''") + property + "'' to a non-existing state or set of states: " + stateName+ ".");
			}
			
			propertyString = propertyMap[stateName];
			
			if (propertyString.find("@") != std::string::npos){
				findPropertyFunctionArguments(propertyString, argumentArray, argumentStrArray); 
			}

			StatePropertyFunctions::evaluateStatePropertyFunction(propertyString, stateArray, property, argumentArray, argumentStrArray, workCond);
			
			argumentArray.clear();
			argumentStrArray.clear();
		}	
	}

	template <class StateType>
	void convertRawToState(std::vector<Parse::rawStateStruct> rawStateArray, std::vector<StateType*> &stateArray, std::vector<StateType*> &moduleStateArray){
		stateArray.clear();
		for (auto& rawState: rawStateArray){
			stateArray.push_back(StateType::findPointer(rawState.gasName, rawState.ionCharg, rawState.eleLevel, rawState.vibLevel, rawState.rotLevel, moduleStateArray)[0]);
		}
	}

	void findPropertyFunctionArguments(std::string propertyString, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray){
		std::string argumentString = Parse::tokenizeCharacters(propertyString, (char*) "@")[1];
		argumentStrArray = Parse::tokenizeCharacters(argumentString, (char*) ",");
		static std::map<std::string,std::string> workingConditionsMap = FieldInfo::getFieldMap("workingConditions");

		for(auto& argument: argumentStrArray){
			if (workingConditionsMap.find(argument) != workingConditionsMap.end()){
				argumentArray.push_back( workCond->getValue(argument) );
			}
			else{
				argumentArray.push_back( Parse::str2value(argument) );
			}
		}
	}

	void selfDiagnostic(){
		// 'selfDiagnostic' is a function that performs a diagnostic of the values provided by the user in the setup file, checking for the correctness of the simulation configuration.

		std::string fieldValue;
		double numericValue;

		// check configuration of the electron kinetic module (in case it is present in the setup file)
		if (FieldInfo::getField("electronKinetics")){
			// check whether the isOn field is present and logical. Then in case it is true the checking continues
			fieldValue = FieldInfo::getFieldValue("electronKinetics.isOn");
			if (!FieldInfo::getField("electronKinetics.isOn")){
				Message::error("Error found in the configuration of the setup file.\n''isOn'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
			}
			else if (!Parse::isLogical(fieldValue)){
				Message::error("Error found in the configuration of the setup file.\nWrong value for the field electronKinetics>isOn.\nValue should be logical (''true'' or ''false'').\nPlease, fix the problem and run the code again.");
			}
			else if (Parse::str2bool(fieldValue)){
				// check whether the mandatory fields of the electron kinetic module (apart from isOn) are present when the electron kinetic module is activated
				// --- 'eedfType' field
				fieldValue = FieldInfo::getFieldValue("electronKinetics.eedfType");
				if (!FieldInfo::getField("electronKinetics.eedfType")){
					Message::error("Error found in the configuration of the setup file.\n''eedfType'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				else if (fieldValue != "boltzmannMC" && fieldValue != "prescribedEedf"){
					Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>eedfType''.\nValue should be: ''boltzmannMC'' or ''prescribedEedf''.\nPlease, fix the problem and run the code again.");
				}
				else if (fieldValue == "prescribedEedf"){
					// check whether the shapeParameter field is present and the value is correct (between 1 and 2)
					numericValue = FieldInfo::getFieldNumericValue("electronKinetics.shapeParameter");
					if (!FieldInfo::getField("electronKinetics.shapeParameter")){
						Message::error("Error found in the configuration of the setup file.\n''shapeParameter'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (numericValue < 1 || numericValue > 2){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>shapeParameter''.\nValue should a number between 1 and 2 (1 for Maxwellian and 2 for Druyvesteyn).\nPlease, fix the problem and run the code again.");
					}
				}
				else{
					// --- 'ionizationOperator' field
					if (!FieldInfo::getField("electronKinetics.ionizationOperatorType")){
						Message::error("Error found in the configuration of the setup file.\n''ionizationOperatorType'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldValue("electronKinetics.ionizationOperatorType") != "oneTakesAll" &&
						     FieldInfo::getFieldValue("electronKinetics.ionizationOperatorType") != "equalSharing" &&
						     FieldInfo::getFieldValue("electronKinetics.ionizationOperatorType") != "usingSDCS" &&
						     FieldInfo::getFieldValue("electronKinetics.ionizationOperatorType") != "randomUniform"){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>ionizationOperatorType''.\nValue should be either: ''oneTakesAll'', ''equalSharing'', ''usingSDCS'' or ''randomUniform''.\nPlease, fix the problem and run the code again.");
					}
				}			
				// --- 'LXCatFiles' field
				if (!FieldInfo::getField("electronKinetics.LXCatFiles")){
					Message::error("Error found in the configuration of the setup file.\n''LXCatFiles'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				// --- 'gasProperties' field
				if (!FieldInfo::getField("electronKinetics.gasProperties")){
					Message::error("Error found in the configuration of the setup file.\n''gasProperties'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				else if (!FieldInfo::getField("electronKinetics.gasProperties.mass")){
					Message::error("Error found in the configuration of the setup file.\n''mass'' field not found in the ''electronKinetics>gasProperties'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				else if (!FieldInfo::getField("electronKinetics.gasProperties.fraction")){
					Message::error("Error found in the configuration of the setup file.\n''fraction'' field not found in the ''electronKinetics>gasProperties'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				// --- 'stateProperties' field
				if (!FieldInfo::getField("electronKinetics.stateProperties")){
					Message::error("Error found in the configuration of the setup file.\n''stateProperties'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				else if (!FieldInfo::getField("electronKinetics.stateProperties.population")){
					Message::error("Error found in the configuration of the setup file.\n''population'' field not found in the ''electronKinetics>stateProperties'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				// --- 'numerics' field
				if (FieldInfo::getFieldValue("electronKinetics.eedfType") == "prescribedEedf"){
					if (!FieldInfo::getField("electronKinetics.numerics")){
						Message::error("Error found in the configuration of the setup file.\n''numerics'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numerics.energyGrid")){
						Message::error("Error found in the configuration of the setup file.\n''energyGrid'' field not found in the ''electronKinetics>numerics'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numerics.energyGrid.maxEnergy")){
						Message::error("Error found in the configuration of the setup file.\n''maxEnergy'' field not found in the ''electronKinetics>numerics>energyGrid'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.maxEnergy") <= 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numerics>energyGrid>maxEnergy''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numerics.energyGrid.cellNumber")){
						Message::error("Error found in the configuration of the setup file.\n''cellNumber'' field not found in the ''electronKinetics>numerics>energyGrid'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.cellNumber") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.cellNumber") <= 0,1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numerics>energyGrid>cellNumber''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
				}
				// --- 'numericsMC' field
				else if ((FieldInfo::getFieldValue("electronKinetics.eedfType")).find("boltzmannMC") != std::string::npos){
					if (!FieldInfo::getField("electronKinetics.numericsMC")){
						Message::error("Error found in the configuration of the setup file.\n''numericsMC'' field not found in the ''electronKinetics'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numericsMC.nElectrons")){
						Message::error("Error found in the configuration of the setup file.\n''nElectrons'' field not found in the ''electronKinetics>numericsMC'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numericsMC.gasTemperatureEffect")){
						Message::error("Error found in the configuration of the setup file.\n''gasTemperatureEffect'' field not found in the ''electronKinetics>numericsMC'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldValue("electronKinetics.numericsMC.gasTemperatureEffect") != "false" &&
							 FieldInfo::getFieldValue("electronKinetics.numericsMC.gasTemperatureEffect") != "true" &&
							 FieldInfo::getFieldValue("electronKinetics.numericsMC.gasTemperatureEffect") != "smartActivation"){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field electronKinetics>numericsMC.\nValue should be ''false'', ''true'' or ''smartActivation''.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.nEnergyCells") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nEnergyCells") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nEnergyCells"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nEnergyCells''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.nCosAngleCells") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nCosAngleCells") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nCosAngleCells"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nCosAngleCells''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.nAxialVelocityCells") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nAxialVelocityCells") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nAxialVelocityCells"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nAxialVelocityCells''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.nRadialVelocityCells") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nRadialVelocityCells") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nRadialVelocityCells"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nRadialVelocityCells''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.nInterpPoints") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nInterpPoints") <= 0 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nInterpPoints"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nInterpPoints''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.synchronizationTimeXMaxCollisionFrequency") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.synchronizationTimeXMaxCollisionFrequency") <= 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>synchronizationTimeXMaxCollisionFrequency''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getField("electronKinetics.numericsMC.synchronizationOverSampling") &&
							 FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.synchronizationOverSampling") < 1 ||
							 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.synchronizationOverSampling"), 1) != 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>synchronizationOverSampling''.\nValue should be a single integer >= 1.\nPlease, fix the problem and run the code again.");
					}
					if (FieldInfo::getFieldValue("electronKinetics.eedfType") == "boltzmannMC"){
						if (FieldInfo::getField("electronKinetics.numericsMC.initialElecTempOverGasTemp") && FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.initialElecTempOverGasTemp") <= 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>initialElecTempOverGasTemp''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}					
						else if (FieldInfo::getField("electronKinetics.numericsMC.minCollisionsBeforeSteadyState") && FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.minCollisionsBeforeSteadyState") < 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>minCollisionsBeforeSteadyState''.\nValue should be a single non-negative number.\nPlease, fix the problem and run the code again.");
						}					
						else if (FieldInfo::getField("electronKinetics.numericsMC.maxCollisionsBeforeSteadyState") && FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.maxCollisionsBeforeSteadyState") <= 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>maxCollisionsBeforeSteadyState''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.maxCollisionsAfterSteadyState") && FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.maxCollisionsAfterSteadyState") <= 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>maxCollisionsAfterSteadyState''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}																
						else if (FieldInfo::getField("electronKinetics.numericsMC.nIntegrationPoints") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPoints") < 500 ||
								std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPoints"), 1) != 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nIntegrationPoints''.\nValue should be a single integer >= 500.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.nIntegrationPhases") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPhases") < 1 ||
								std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPhases"), 1) != 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nIntegrationPhases''.\nValue should be a single integer > 1.\nPlease, fix the problem and run the code again.");
						}						
						else if (FieldInfo::getField("electronKinetics.numericsMC.nIntegratedSSTimes") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegratedSSTimes") <= 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>nIntegratedSSTimes''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.integratedAbsoluteTime") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.integratedAbsoluteTime") <= 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>integratedAbsoluteTime''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}																			
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.meanEnergy") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.meanEnergy") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>meanEnergy''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.fluxDriftVelocity") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.fluxDriftVelocity") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>fluxDriftVelocity''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>fluxDiffusionCoeffs''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.bulkDriftVelocity") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.bulkDriftVelocity") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>bulkDriftVelocity''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.bulkDiffusionCoeffs") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.bulkDiffusionCoeffs") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>bulkDiffusionCoeffs''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}						
						else if (FieldInfo::getField("electronKinetics.numericsMC.relError.powerBalance") &&
								FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.powerBalance") <= 0 ){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numericsMC>relError>powerBalance''.\nValue should be a single positive number.\nPlease, fix the problem and run the code again.");
						}
					}
				}
				// --- check configuration of the CAR in case it is activated
				if (FieldInfo::getField("electronKinetics.CARgases")){
					if (!FieldInfo::getField("electronKinetics.gasProperties.rotationalConstant")){
						Message::error("Error found in the configuration of the setup file.\n''rotationalConstant'' field not found in the ''electronKinetics>gasProperties'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.gasProperties.electricQuadrupoleMoment")){
						Message::error("Error found in the configuration of the setup file.\n''electricQuadrupoleMoment'' field not found in the ''electronKinetics>gasProperties'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}			
				}
				// check configuration of the smart grid in case it is activated
				if (FieldInfo::getField("electronKinetics.numerics.energyGrid.smartGrid")){
					if (!FieldInfo::getField("electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay")){
						Message::error("Error found in the configuration of the setup file.\n''minEedfDecay'' field not found in the ''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay") <= 0){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numerics>energyGrid>smartGrid>minEedfDecay''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay")){
						Message::error("Error found in the configuration of the setup file.\n''maxEedfDecay'' field not found in the ''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay") <= 0 ||
							 FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay") >= FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay")){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numerics>energyGrid>smartGrid>maxEedfDecay''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
					}
					else if (!FieldInfo::getField("electronKinetics.numerics.energyGrid.smartGrid.updateFactor")){
						Message::error("Error found in the configuration of the setup file.\n''updateFactor'' field not found in the ''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.updateFactor") <= 0 || 
						     FieldInfo::getFieldNumericValue("electronKinetics.numerics.energyGrid.smartGrid.updateFactor") >= 1){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>numerics>energyGrid>smartGrid>updateFactor''.\nValue should be a number larger than 0 and smaller than 1.\nPlease, fix the problem and run the code again.");
					}
				}
				// check configuration of the anisotropicScattering
				if (FieldInfo::getField("electronKinetics.anisotropicScattering")){
					fieldValue = FieldInfo::getFieldValue("electronKinetics.anisotropicScattering.isOn");
					// check whether the isOn field is present and logical. Then in case it is true the checking continues
					if (!FieldInfo::getField("electronKinetics.anisotropicScattering.isOn")){
						Message::error("Error found in the configuration of the setup file.\n''isOn'' field not found in the ''electronKinetics>anisotropicScattering'' section of the setup file.\nPlease, fix the problem and run the code again.");
					}
					else if (!Parse::isLogical(fieldValue)){
						Message::error("Error found in the configuration of the setup file.\nWrong value for the field electronKinetics>anisotropicScattering>isOn.\nValue should be logical (''true'' or ''false'').\nPlease, fix the problem and run the code again.");
					}
					else if (Parse::str2bool(fieldValue)){
						if (!FieldInfo::getField("electronKinetics.anisotropicScattering.collisions")){
							Message::error("Error found in the configuration of the setup file.\n''collisions'' field not found in the ''electronKinetics>anisotropicScattering'' section of the setup file.\nPlease, fix the problem and run the code again.");
						}
						else if (!FieldInfo::getField("electronKinetics.anisotropicScattering.angleNumber")){
							Message::error("Error found in the configuration of the setup file.\n''angleNumber'' field not found in the ''electronKinetics>anisotropicScattering'' section of the setup file.\nPlease, fix the problem and run the code again.");
						}						
						else if (FieldInfo::getFieldNumericValue("electronKinetics.anisotropicScattering.angleNumber") <= 0 ||
								 std::fmod(FieldInfo::getFieldNumericValue("electronKinetics.anisotropicScattering.angleNumber"), 1) != 0){
							Message::error("Error found in the configuration of the setup file.\nWrong value for the field ''electronKinetics>anisotropicScattering>angleNumber''.\nValue should be a single positive integer.\nPlease, fix the problem and run the code again.");
						}												
					}
				}			
			}
		}

		// check for 'empty' simulations (simulations with no module activated)
		if (!Parse::str2bool(FieldInfo::getFieldValue("electronKinetics.isOn"))){
			Message::error("Error found in the configuration of the setup file.\n''electronKinetics'' module is not activated.\nPlease, fix the problem and run the code again.");
		}

		// check configuration of the graphical user interface (in case it is present in the setup file)
		if (FieldInfo::getField("gui")){
			// check whether the isOn field is present and logical. Then in case it is true the checking continues
			fieldValue = FieldInfo::getFieldValue("gui.isOn");
			if (!FieldInfo::getField("gui.isOn")){
				Message::error("Error found in the configuration of the setup file.\n''isOn'' field not found in the ''gui'' section of the setup file.\nPlease, fix the problem and run the code again.");
			}
			else if (!Parse::isLogical(fieldValue)){
				Message::error("Error found in the configuration of the setup file.\nWrong value for the field gui>isOn.\nValue should be logical (''true'' or ''false'').\nPlease, fix the problem and run the code again.");
			}
			else if (Parse::str2bool(fieldValue)){
				// check whether the mandatory fields of the gui module (apart from isOn) are present when the gui module is activated
				// 'plotOptions' field
				if (!FieldInfo::getField("gui.plotOptions")){
					Message::error("Error found in the configuration of the setup file.\n''plotOptions'' field not found in the ''gui'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				std::vector<std::string> plotOptions = FieldInfo::getFieldChildNames("gui.plotOptions");
				std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");
				std::vector<std::string> possibleOptions;
				std::string possibleOptionsString;
				if (eedfType == "boltzmannMC"){
					possibleOptions = {"MCTemporalInfo","MCTemporalInfo_periodic", "distributionFunctions", "swarmParameters", "powerBalance"};
					possibleOptionsString = "MCTemporalInfo, MCTemporalInfo_periodic, distributionFunctions, swarmParameters or powerBalance";
				}
				else if (eedfType == "prescribedEedf"){
					possibleOptions = {"distributionFunctions", "swarmParameters", "powerBalance"};
					possibleOptionsString = "distributionFunctions, swarmParameters or powerBalance";
				}
				for (auto plotOption: plotOptions){
					bool validOption = false;
					for (auto possibleOption: possibleOptions){
						if (plotOption == possibleOption){
							validOption = true;
							break;
						}
					}
					if (!validOption){
						Message::error(std::string("Error found in the configuration of the setup file.\nWrong value for the field ''gui>plotOptions''. Possible options are: ") + possibleOptionsString + ".\nPlease, fix the problem and run the code again.");
					}
				}
				// 'terminalDisp' field
				if (!FieldInfo::getField("gui.terminalDisp")){
					Message::error("Error found in the configuration of the setup file.\n''terminalDisp'' field not found in the ''gui'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				plotOptions = FieldInfo::getFieldChildNames("gui.terminalDisp");
				if (eedfType == "boltzmannMC"){
					possibleOptions = {"setup", "MCStatus"};
					possibleOptionsString = "setup or MCStatus";
				}
				else if (eedfType == "prescribedEedf"){
					possibleOptions = {"setup"};
					possibleOptionsString = "setup";
				}
				for (auto plotOption: plotOptions){
					bool validOption = false;
					for (auto possibleOption: possibleOptions){
						if (plotOption == possibleOption){
							validOption = true;
							break;
						}
					}
					if (!validOption){
						Message::error(std::string("Error found in the configuration of the setup file.\nWrong value for the field ''gui>terminalDisp''. Possible options are: ") + possibleOptionsString + ".\nPlease, fix the problem and run the code again.");
					}
				}
			}			
		}

		// check configuration of the output (in case it is present in the setup file)
		if (FieldInfo::getField("output")){
			// check whether the isOn field is present and logical. Then in case it is true the checking continues
			fieldValue = FieldInfo::getFieldValue("output.isOn");
			if (!FieldInfo::getField("output.isOn")){
				Message::error("Error found in the configuration of the setup file.\n''isOn'' field not found in the ''output'' section of the setup file.\nPlease, fix the problem and run the code again.");
			}
			else if (!Parse::isLogical(fieldValue)){
				Message::error("Error found in the configuration of the setup file.\nWrong value for the field output>isOn.\nValue should be logical (''true'' or ''false'').\nPlease, fix the problem and run the code again.");
			}
			else if (Parse::str2bool(fieldValue)){
				// check whether the mandatory fields of the output module (apart from isOn) are present when the output module is activated
				if (!FieldInfo::getField("output.dataFiles")){
					Message::error("Error found in the configuration of the setup file.\n''dataFiles'' field not found in the ''output'' section of the setup file.\nPlease, fix the problem and run the code again.");
				}
				std::vector<std::string> dataFiles = FieldInfo::getFieldChildNames("output.dataFiles");
				std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");
				std::vector<std::string> possibleOptions;
				std::string possibleOptionsString;
				if (eedfType == "boltzmannMC"){
					possibleOptions = {"eedf", "evdf", "swarmParameters", "rateCoefficients", "powerBalance", "MCTemporalInfo", "MCTemporalInfo_periodic", "MCSimDetails", "lookUpTable"};
					possibleOptionsString = "eedf, evdf, swarmParameters, rateCoefficients, powerBalance, MCTemporalInfo, MCTemporalInfo_periodic, MCSimDetails, lookUpTable";
				}
				else if (eedfType == "prescribedEedf"){
					possibleOptions = {"eedf", "swarmParameters", "rateCoefficients", "powerBalance", "lookUpTable"};
					possibleOptionsString = "eedf, swarmParameters, rateCoefficients, powerBalance, lookUpTable";
				}
				for (auto dataFile: dataFiles){
					bool validOption = false;
					for (auto possibleOption: possibleOptions){
						if (dataFile == possibleOption){
							validOption = true;
							break;
						}
					}
					if (!validOption){
						Message::error(std::string("Error found in the configuration of the setup file.\nWrong value for the field ''gui>dataFiles''. Possible data files are: ") + possibleOptionsString + ".\nPlease, fix the problem and run the code again.");
					}
				}
			}			
		}
	}

	void finishSimulation(){
		// get the end instant and calculate the elapsed time

		end = std::chrono::high_resolution_clock::now();
		auto elap = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		double elapsed = (elap.count())*1E-3;

		std::cout<<"Finished!"<<std::endl;
		std::cout<<"Elapsed time is "<<elapsed<<" seconds.\n";
	}
};
#endif