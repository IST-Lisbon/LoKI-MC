#ifndef __GUI__
#define __GUI__
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/PrescribedEedf.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <boost/signals2.hpp>
#include <map>

#if defined (_MSC_VER)  // Visual studio
    #define SUPRESS_OUTPUT " > nul 2>&1"
#else 
    #define SUPRESS_OUTPUT " >>/dev/null 2>>/dev/null"
#endif

template <class ElectronKineticsType>
class GUI{
public:

	// ----- class attributes ----- //

	char gnuplotTermOptions[500];

	bool MCTemporalInfoIsToBePlotted = false;
	bool distributionFunctionsIsToBePlotted = false;
	bool swarmParamsIsToBePlotted = false;
	bool powerBalanceIsToBePlotted = false;

	std::string eedfType;
	ElectronKineticsType* electronKinetics = NULL;

	double inputParameter;
	int numberOfJobs, currentJobID;

	std::string variableCondition;
	std::string xLabel;	

	// ----- class methods ----- //

	template <class SetupType>
	GUI(SetupType* setup){

		// create a directory for the temporary data files used in the gnuplot scripts
		std::filesystem::create_directories("LoKI-MC/GnuplotFunctions/TempData");

		// save the gnuplot terminal options
		std::string gnuplotTerminal = "qt";
		if (!FieldInfo::getFieldValue("gui.gnuplotTerminal").empty()){
			gnuplotTerminal = FieldInfo::getFieldValue("gui.gnuplotTerminal");
		}
		int fontSize = 9;
		if (!FieldInfo::getFieldValue("gui.fontSize").empty()){
			fontSize = FieldInfo::getFieldNumericValue("gui.fontSize");
		}
		std::sprintf(gnuplotTermOptions, "set term %s size 1150,600 position 100,50 font 'Arial,%d' raise", (char*)gnuplotTerminal.c_str(), fontSize);

		// save which information must be saved
		for (auto option: FieldInfo::getFieldChildNames("gui.plotOptions")){
			if (option == "MCTemporalInfo"){
				MCTemporalInfoIsToBePlotted = true;
			}
			else if (option == "distributionFunctions"){
				distributionFunctionsIsToBePlotted = true;
			}
			else if (option == "swarmParameters"){
				swarmParamsIsToBePlotted = true;
			}
			else if (option == "powerBalance"){
				powerBalanceIsToBePlotted = true;
			}
		}

		// save the eedf type
		eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");

		// save the electron kinetics
		electronKinetics = setup->electronKinetics;

		// save the total number of jobs and initialize the job ID
		numberOfJobs = setup->numberOfJobs;
		currentJobID = -1;

		// assign the variable condition and xLabel
		variableCondition = electronKinetics->workCond->variableCondition;
		if (variableCondition == "electronTemperature"){
			xLabel = "T_e (eV)";
		}
		else if (variableCondition == "reducedElecField"){
			xLabel = "E/N (Td)";
		}
		
		if (setup->enableElectronKinetics){
		    // connect the signal 'obtainedNewEedfSignal' to 'GUI::electronKineticsSolution', in order to update the GUI each time a new eedf is found
		    // idea taken from https://stackoverflow.com/questions/3047381/boost-signals-and-passing-class-method
			electronKinetics->obtainedNewEedfSignal.connect(boost::bind(&GUI::electronKineticsSolution, this));
		}
	}

	void electronKineticsSolution(){

		// update input parameter 
		updateInputParameter();

		// plot selected results of the electron kinetics
		if (MCTemporalInfoIsToBePlotted){
			plotMCTemporalInfo(electronKinetics->samplingTimes, electronKinetics->meanEnergies, electronKinetics->meanPositions, electronKinetics->positionCovariances, electronKinetics->meanVelocities);
		}
		if (distributionFunctionsIsToBePlotted){
			plotDistributionFunctions(electronKinetics->eedf, electronKinetics->energyGrid->cell, electronKinetics->evdf, electronKinetics->radialVelocityCells ,electronKinetics->axialVelocityCells);
		}
		if (powerBalanceIsToBePlotted){
			plotPower(electronKinetics->power);
		}
		if (swarmParamsIsToBePlotted){
			plotSwarm(electronKinetics->swarmParam);
		}
	}

	void updateInputParameter(){
		// updates the input parameter 

		++currentJobID;

		if (variableCondition == "electronTemperature"){
			inputParameter = electronKinetics->workCond->electronTemperature;
		}
		else if (variableCondition == "reducedElecField"){
			inputParameter = electronKinetics->workCond->reducedElecField;
		}
	}	

	void plotMCTemporalInfo(Eigen::ArrayXd &samplingTimes, Eigen::ArrayXd &meanEnergies, Eigen::ArrayXXd &meanPositions, Eigen::ArrayXXd &positionCovariances, Eigen::ArrayXXd &meanVelocities){

		// write the time-dependent data in an auxiliary file
		FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat", "w");
		std::fprintf(fileID, "#Time(s)             MeanEnergy(eV)      xPos(m)             yPos(m)             zPos(m)             xSqWidth(m2)        ySqWidth(m2)        zSqWidth(m2)        xVel(m/s)           yVel(m/s)           zVel(m/s)           \n");
		int nSamplingPoints = electronKinetics->nSamplingPoints;
		for (int i = 0; i < nSamplingPoints; ++i){
			std::fprintf(fileID, "%-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e\n", samplingTimes[i], meanEnergies[i], meanPositions(i, 0), meanPositions(i, 1), meanPositions(i, 2), 
				    positionCovariances(i,0), positionCovariances(i,4), positionCovariances(i,8), meanVelocities(i,0), meanVelocities(i,1), meanVelocities(i,2));
		}
		std::fclose(fileID);

		// call gnuplot
		char commandLine[1000];
		std::sprintf(commandLine, "gnuplot -s -p -e \"steadyStateTime=%.6e; %s title 'MCTemporalInfo: E/N = %.4e Td'; load 'LoKI-MC/GnuplotFunctions/MCTemporalInfo.p'\" %s",
				electronKinetics->steadyStateTime, gnuplotTermOptions, electronKinetics->workCond->reducedElecField, SUPRESS_OUTPUT);
		int systemReturn = std::system(commandLine);
	}

	void plotDistributionFunctions(Eigen::ArrayXd &eedf, Eigen::ArrayXd &eedfEnergyCells, Eigen::ArrayXXd &evdf, Eigen::ArrayXd &radialVelocityCells, Eigen::ArrayXd &axialVelocityCells){

		// write the eedf data in an auxiliary file
		int numberOfEnergyCells = eedfEnergyCells.size();
		FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/eedf.dat", "w");
		std::fprintf(fileID, "#Energy(eV)          EEDF(eV^-(3/2))     \n");
		for (int i = 0; i < numberOfEnergyCells; ++i){
			std::fprintf(fileID, "%-18.10e   %-18.10e\n", eedfEnergyCells[i], eedf[i]);
		}
		std::fclose(fileID);

		char commandLine[1000];
		if (eedfType == "prescribedEedf"){
			std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'EEDF: T_e = %.4e eV'; load 'LoKI-MC/GnuplotFunctions/eedf.p'\" %s", 
				    gnuplotTermOptions, electronKinetics->workCond->electronTemperature, SUPRESS_OUTPUT);
		}
		else if (eedfType == "boltzmannMC" ){
			// write the evdf data in an auxiliary file
			int numberOfAxialCells = radialVelocityCells.size(), numberOfRadialCells = axialVelocityCells.size();
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/evdf.dat", "w");
			std::fprintf(fileID, "#v_r(m/s)          v_z(m/s)           EVDF(m-3s-3)       \n");
			for (int i = 0; i < numberOfRadialCells; ++i){
				for (int j = 0; j < numberOfAxialCells; ++j){
					std::fprintf(fileID, "%-18.10e %-18.10e %-18.10e\n", radialVelocityCells[i], axialVelocityCells[j], evdf(i,j));
				}
			}
			std::fclose(fileID);

			std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'EEDF and EVDF: E/N = %.4e Td'; load 'LoKI-MC/GnuplotFunctions/eedfEvdf.p'\" %s", 
				    gnuplotTermOptions, electronKinetics->workCond->reducedElecField, SUPRESS_OUTPUT);
		}

		// call gnuplot
		int systemReturn = std::system(commandLine);
	}

	void plotPower(PowerStruct power){

		// local copy of the power values
		std::map<std::string,double> powerMap = electronKinetics->power.Map;
		double ref = powerMap["reference"];

		// write the headers of the auxiliary file
		if (currentJobID == 0){
			// write the auxiliary file
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/power.dat", "w");
			std::fprintf(fileID, "#inputParameter    Field              elasticGain        elasticLoss        CARGain            CARLoss            rotationalGain     rotationalLoss     vibrationalGain    vibrationalLoss    electronicGain     electronicLoss     ionization         attachment         growth             \n");
			std::fclose(fileID);
		}

		// write the values of the current job
		FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/power.dat", "a");
		std::fprintf(fileID, "%-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e \n",
				inputParameter, powerMap["field"]/ref, powerMap["elasticGain"]/ref, powerMap["elasticLoss"]/ref, powerMap["carGain"]/ref, powerMap["carLoss"]/ref, powerMap["rotationalSup"]/ref, powerMap["rotationalIne"]/ref,
				powerMap["vibrationalSup"]/ref, powerMap["vibrationalIne"]/ref, powerMap["excitationSup"]/ref, powerMap["excitationIne"]/ref, powerMap["ionizationIne"]/ref, 
				powerMap["attachmentIne"]/ref, powerMap["eDensGrowth"]/ref);
		std::fclose(fileID);

		// if this is the last parameter, call gnuplot
		if (currentJobID == numberOfJobs-1){
			// call gnuplot
			char commandLine[1000]; 
			if (eedfType == "prescribedEedf"){
				std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'Power Balance'; load 'LoKI-MC/GnuplotFunctions/power_PrescribedEedf.p'\" %s", gnuplotTermOptions, SUPRESS_OUTPUT);
			}
			else if (eedfType == "boltzmannMC"){
				std::sprintf(commandLine, "gnuplot -s -p -e \"xLabel='%s';%s title 'Power Balance'; load 'LoKI-MC/GnuplotFunctions/power_BoltzmannMC.p'\" %s", xLabel.c_str(), gnuplotTermOptions, SUPRESS_OUTPUT);
			}
			int systemReturn = std::system(commandLine);
		}
	}

	void plotSwarm(std::map<std::string,double> swarmParam){

		if (eedfType == "prescribedEedf"){

			// write the headers of the auxiliary file
			if (currentJobID == 0){
				FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "w");
				std::fprintf(fileID, "#inputParameter    RedDif(1/(ms))     RedMob(1/(msV))    RedTow(m2)         RedAtt(m2)         MeanE(eV)          CharE(eV)          \n");
				std::fclose(fileID);
			}

			// write the values of the current job
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "a");
			std::fprintf(fileID, "%-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e %-18.10e \n", inputParameter, swarmParam["redDiffCoeff"], swarmParam["redMobCoeff"], swarmParam["redTownsendCoeff"],
					swarmParam["redAttCoeff"], swarmParam["meanEnergy"], swarmParam["characEnergy"]);
			std::fclose(fileID);

			// plot the swarm parameters if this is the last solution
			if (currentJobID == numberOfJobs-1){
				char commandLine[1000];
				std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'Swarm Parameters'; load 'LoKI-MC/GnuplotFunctions/swarmParams_PrescribedEedf.p'\" %s", gnuplotTermOptions, SUPRESS_OUTPUT);
				int systemReturn = std::system(commandLine);
			}
		}


		else if (eedfType == "boltzmannMC"){ 

			// write the headers of the auxiliary file
			if (currentJobID == 0){
				FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "w");
				std::fprintf(fileID, "#inputParameter           FluxRedTransvDif(1/(ms))  FluxRedLongDif(1/(ms))    FluxRedMob(1/(msV))       FluxCharE(eV)             BulkRedTransvDif(1/(ms))  BulkRedLongDif(1/(ms))    BulkRedMob(1/(msV))       BulkCharE(eV)             MeanE(eV)                 ionCoeff(m-3)             attCoeff(m-3)\n");
				std::fclose(fileID);
			}

			// write the values of the current job
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "a");
			std::fprintf(fileID, "%-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e\n", inputParameter, 
					swarmParam["fluxRedTransvDiffCoeff"], swarmParam["fluxRedLongDiffCoeff"], swarmParam["fluxRedMobCoeff"], swarmParam["fluxCharacEnergy"],
					swarmParam["bulkRedTransvDiffCoeff"], swarmParam["bulkRedLongDiffCoeff"], swarmParam["bulkRedMobCoeff"], swarmParam["bulkCharacEnergy"], swarmParam["meanEnergy"], swarmParam["totalIonRateCoeff"], swarmParam["totalAttRateCoeff"]);
			std::fclose(fileID);

		    // plot the swarm parameters if this is the last solution
			if (currentJobID == numberOfJobs-1){
				char commandLine[1000];
				std::sprintf(commandLine, "gnuplot -s -p -e \"xLabel='%s'; %s title 'Swarm Parameters'; load 'LoKI-MC/GnuplotFunctions/swarmParams_BoltzmannMC.p'\" %s", xLabel.c_str(), gnuplotTermOptions, SUPRESS_OUTPUT);
				int systemReturn = std::system(commandLine);
			}
		}
	}
};

#endif