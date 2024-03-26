#ifndef __GUI__
#define __GUI__
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/GeneralDefinitions.h"
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

class BoltzmannMC;

template <class ElectronKineticsType>
class GUI{
public:

	// ----- class attributes ----- //

	char gnuplotTermOptions[500];

	bool MCTemporalInfoIsToBePlotted = false;
	bool MCTemporalInfoPeriodicIsToBePlotted = false;
	bool distributionFunctionsIsToBePlotted = false;
	bool swarmParamsIsToBePlotted = false;
	bool powerBalanceIsToBePlotted = false;

	std::string eedfType;
	ElectronKineticsType* electronKinetics = NULL;

	double inputParameter;
	int numberOfJobs, currentJobID;

	std::string variableCondition;
	bool isCylindricallySymmetric;
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
			else if (option == "MCTemporalInfo_periodic"){
				MCTemporalInfoPeriodicIsToBePlotted = true;
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
		else if (variableCondition == "reducedMagField"){
			xLabel = "B/N (Hx)";
		}
		else if (variableCondition == "elecFieldAngle"){
			xLabel = "E-angle (º)";
		}
		else if (variableCondition == "excitationFrequency"){
			xLabel = "Exc. Frequency (Hz)";
		}

		isCylindricallySymmetric = setup->workCond->isCylindricallySymmetric;

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
		if (MCTemporalInfoIsToBePlotted && eedfType == "boltzmannMC" ){
			plotMCTemporalInfo(electronKinetics->samplingTimes, electronKinetics->meanEnergies, electronKinetics->meanPositions, electronKinetics->positionCovariances, electronKinetics->meanVelocities);
		}
		if (MCTemporalInfoPeriodicIsToBePlotted && eedfType == "boltzmannMC"){
			plotMCTemporalInfoPeriodic(electronKinetics->integrationPhases, electronKinetics->meanEnergies_periodic, electronKinetics->fluxVelocities_periodic, electronKinetics->bulkVelocities_periodic, electronKinetics->fluxDiffusionCoeffs_periodic, electronKinetics->bulkDiffusionCoeffs_periodic);
		}		
		if (distributionFunctionsIsToBePlotted){
			plotDistributionFunctions(electronKinetics->eedf, electronKinetics->energyGrid->cell, electronKinetics->evdf, electronKinetics->radialVelocityCells ,electronKinetics->axialVelocityCells);
		}
		if (powerBalanceIsToBePlotted){
			plotPower(electronKinetics->power.Map);
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
		else if (variableCondition == "reducedMagField"){
			inputParameter = electronKinetics->workCond->reducedMagField;
		}
		else if (variableCondition == "elecFieldAngle"){
			inputParameter = electronKinetics->workCond->elecFieldAngle;
		}
		else if (variableCondition == "excitationFrequency"){
			inputParameter = electronKinetics->workCond->excitationFrequency;
		}
	}	

	void plotMCTemporalInfo(Eigen::ArrayXd &samplingTimes, Eigen::ArrayXd &meanEnergies, Eigen::ArrayXXd &meanPositions, Eigen::ArrayXXd &positionCovariances, Eigen::ArrayXXd &meanVelocities){

		// write the time-dependent data in an auxiliary file
		FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat", "w");
		std::vector<std::string> vars = {"#Time(s)", "MeanEnergy(eV)", "xPos(m)", "yPos(m)", "zPos(m)", "xSqWidth(m2)", "ySqWidth(m2)", "zSqWidth(m2)", "xVel(m/s)", "yVel(m/s)", "zVel(m/s)"};
		std::string header;
		for (auto var: vars){
			header += var + std::string(21-var.size(), ' ');
		}		
		std::fprintf(fileID, "%s\n", header.c_str());
		int nSamplingPoints = electronKinetics->nSamplingPoints;
		for (int i = 0; i < nSamplingPoints; ++i){
			std::fprintf(fileID, "%-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e %-20.10e\n", samplingTimes[i], meanEnergies[i], meanPositions(i, 0), meanPositions(i, 1), meanPositions(i, 2), 
				    positionCovariances(i,0), positionCovariances(i,4), positionCovariances(i,8), meanVelocities(i,0), meanVelocities(i,1), meanVelocities(i,2));
		}
		std::fclose(fileID);

		// call gnuplot
		char commandLine[1000];
		std::sprintf(commandLine, "gnuplot -s -p -e \"steadyStateTime=%.6e; %s title 'MCTemporalInfo: E/N = %.4e Td , E-angle = %.4e º, ExcFreq = %.4e Hz , B/N = %.4e Hx'; load 'LoKI-MC/GnuplotFunctions/MCTemporalInfo.p'\" %s",
				electronKinetics->steadyStateTime, gnuplotTermOptions, electronKinetics->workCond->reducedElecField, electronKinetics->workCond->elecFieldAngle, electronKinetics->workCond->excitationFrequency, electronKinetics->workCond->reducedMagField, SUPRESS_OUTPUT);
		int systemReturn = std::system(commandLine);
	}

	void plotMCTemporalInfoPeriodic(Eigen::ArrayXd &integrationPhases, Eigen::ArrayXd &meanEnergies_periodic, Eigen::ArrayXXd &fluxVelocities_periodic, Eigen::ArrayXXd &bulkVelocities_periodic, Eigen::ArrayXXd &fluxDiffusionCoeffs_periodic, Eigen::ArrayXXd &bulkDiffusionCoeffs_periodic){

		double angularFrequency = electronKinetics->workCond->excitationFrequency*2.0*M_PI;
		if (angularFrequency == 0){
			return;
		}

		FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat", "w");
		std::vector<std::string> vars = {"#Phase(rad)", "Phase(s)","E/N(Td)", "MeanEnergy(eV)", "FluxV_x(m/s)","FluxV_y(m/s)","FluxV_z(m/s)","BulkV_x(m/s)", "BulkV_y(m/s)","BulkV_z(m/s)", "FluxND_xx(1/(ms))", "FluxND_xy(1/(ms))", "FluxND_xz(1/(ms))", "FluxND_yx(1/(ms))", "FluxND_yy(1/(ms))", "FluxND_yz(1/(ms))", "FluxND_zx(1/(ms))", "FluxND_zy(1/(ms))", "FluxND_zz(1/(ms))", "BulkND_xx(1/(ms))", "BulkND_xy(1/(ms))", "BulkND_xz(1/(ms))", "BulkND_yx(1/(ms))", "BulkND_yy(1/(ms))", "BulkND_yz(1/(ms))", "BulkND_zx(1/(ms))", "BulkND_zy(1/(ms))", "BulkND_zz(1/(ms))"};
		std::string header;
		for (auto var: vars){
			header += var + std::string(20-var.size(), ' ');
		}
		std::fprintf(fileID, "%s\n", header.c_str());
		double EN = electronKinetics->workCond->reducedElecField;
		double N = electronKinetics->workCond->gasDensity;		
		int nIntegrationPhases = integrationPhases.size();
		for (int phaseIdx = 0; phaseIdx < nIntegrationPhases; ++phaseIdx){			
			std::fprintf(fileID, "%-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e %-19.10e\n", 
				integrationPhases[phaseIdx], integrationPhases[phaseIdx]/angularFrequency, std::sqrt(2)*EN*std::cos(integrationPhases[phaseIdx]), meanEnergies_periodic[phaseIdx], fluxVelocities_periodic(phaseIdx,0), fluxVelocities_periodic(phaseIdx,1), fluxVelocities_periodic(phaseIdx,2), bulkVelocities_periodic(phaseIdx,0), bulkVelocities_periodic(phaseIdx,1), bulkVelocities_periodic(phaseIdx,2),
				N*fluxDiffusionCoeffs_periodic(phaseIdx,0), N*fluxDiffusionCoeffs_periodic(phaseIdx,1), N*fluxDiffusionCoeffs_periodic(phaseIdx,2), N*fluxDiffusionCoeffs_periodic(phaseIdx,3), N*fluxDiffusionCoeffs_periodic(phaseIdx,4), N*fluxDiffusionCoeffs_periodic(phaseIdx,5), N*fluxDiffusionCoeffs_periodic(phaseIdx,6), N*fluxDiffusionCoeffs_periodic(phaseIdx,7), N*fluxDiffusionCoeffs_periodic(phaseIdx,8), 
				N*bulkDiffusionCoeffs_periodic(phaseIdx,0), N*bulkDiffusionCoeffs_periodic(phaseIdx,1), N*bulkDiffusionCoeffs_periodic(phaseIdx,2), N*bulkDiffusionCoeffs_periodic(phaseIdx,3), N*bulkDiffusionCoeffs_periodic(phaseIdx,4), N*bulkDiffusionCoeffs_periodic(phaseIdx,5), N*bulkDiffusionCoeffs_periodic(phaseIdx,6), N*bulkDiffusionCoeffs_periodic(phaseIdx,7), N*bulkDiffusionCoeffs_periodic(phaseIdx,8));
		}
		std::fclose(fileID);

		// call gnuplot
		char commandLine[1000];
		std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'MCTemporalInfo_periodic: E/N = %.4e Td , E-angle = %.4e º, ExcFreq = %.4e Hz , B/N = %.4e Hx'; load 'LoKI-MC/GnuplotFunctions/MCTemporalInfo_periodic.p'\" %s",
				gnuplotTermOptions, electronKinetics->workCond->reducedElecField, electronKinetics->workCond->elecFieldAngle, electronKinetics->workCond->excitationFrequency, electronKinetics->workCond->reducedMagField, SUPRESS_OUTPUT);
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
		else if (eedfType == "boltzmannMC"  && isCylindricallySymmetric){
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

			std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'EEDF and EVDF: E/N = %.4e Td , E-angle = %.4e º , ExcFreq = %.4e Hz , B/N = %.4e Hx'; load 'LoKI-MC/GnuplotFunctions/eedfEvdf.p'\" %s", 
				    gnuplotTermOptions, electronKinetics->workCond->reducedElecField, electronKinetics->workCond->elecFieldAngle, electronKinetics->workCond->excitationFrequency, electronKinetics->workCond->reducedMagField, SUPRESS_OUTPUT);
		}
		else{
			std::sprintf(commandLine, "gnuplot -s -p -e \"%s title 'EEDF: E/N = %.4e Td , E-angle = %.4e º , ExcFreq = %.4e Hz , B/N = %.4e Hx'; load 'LoKI-MC/GnuplotFunctions/eedf.p'\" %s", 
				    gnuplotTermOptions, electronKinetics->workCond->reducedElecField, electronKinetics->workCond->elecFieldAngle, electronKinetics->workCond->excitationFrequency, electronKinetics->workCond->reducedMagField, SUPRESS_OUTPUT);
		}




		// call gnuplot
		int systemReturn = std::system(commandLine);
	}

	void plotPower(std::map<std::string,double> powerMap){

		double ref = powerMap["reference"];

		// write the headers of the auxiliary file
		if (currentJobID == 0){
			// write the auxiliary file
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/power.dat", "w");
			std::vector<std::string> vars = {"#inputParameter", "Field", "elasticGain", "elasticLoss", "CARGain", "CARLoss", "rotationalGain", "rotationalLoss", "vibrationalGain", "vibrationalLoss", "electronicGain", "electronicLoss", "ionization", "attachment", "growth"};
			std::string header;
			for (auto var: vars){
				header += var + std::string(19-var.size(), ' ');
			}			
			std::fprintf(fileID, "%s\n", header.c_str());
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

		bool BoltzmannMCSymmetricPlot = eedfType == "boltzmannMC" && isCylindricallySymmetric; //&& 
										//electronKinetics->workCond->reducedMagField == 0 && electronKinetics->workCond->reducedMagFieldArray.size() == 1 &&
										//electronKinetics->workCond->excitationFrequency == 0 && electronKinetics->workCond->excitationFrequencyArray.size() == 1;	

		if (eedfType == "prescribedEedf"){

			// write the headers of the auxiliary file
			if (currentJobID == 0){
				FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "w");
				std::vector<std::string> vars = {"#inputParameter", "RedDif(1/(ms))", "RedMob(1/(msV))", "RedTow(m2)", "RedAtt(m2)", "MeanE(eV)", "CharE(eV)"}; 
				std::string header;	
				for (auto var: vars){
					header += var + std::string(19-var.size(), ' ');
				}				
				std::fprintf(fileID, "%s\n", header.c_str());
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


		else if (BoltzmannMCSymmetricPlot){ 

			// write the headers of the auxiliary file
			if (currentJobID == 0){
				FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "w");
				std::vector<std::string> vars = {"#inputParameter", "FluxRedTransvDif(1/(ms))", "FluxRedLongDif(1/(ms))", "FluxRedMob(1/(msV))", "FluxCharE(eV)", "BulkRedTransvDif(1/(ms))", "BulkRedLongDif(1/(ms))", "BulkRedMob(1/(msV))", "BulkCharE(eV)", "MeanE(eV)", "ionCoeff(m-3)", "attCoeff(m-3)","RedDif_eedf(1/(ms))","RedMob_eedf(1/(msV)","charE_eedf(eV)"}; 
				std::string header;	
				for (auto var: vars){
					header += var + std::string(26-var.size(), ' ');
				}				
				std::fprintf(fileID, "%s\n", header.c_str());
				std::fclose(fileID);
			}

			// write the values of the current job
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "a");
			std::fprintf(fileID, "%-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e\n", inputParameter, 
					swarmParam["fluxRedTransvDiffCoeff"], swarmParam["fluxRedLongDiffCoeff"], swarmParam["fluxRedMobCoeff"], swarmParam["fluxCharacEnergy"],
					swarmParam["bulkRedTransvDiffCoeff"], swarmParam["bulkRedLongDiffCoeff"], swarmParam["bulkRedMobCoeff"], swarmParam["bulkCharacEnergy"], swarmParam["meanEnergy"], swarmParam["totalIonRateCoeff"], swarmParam["totalAttRateCoeff"],
					swarmParam["redDiffCoeff_eedf"], swarmParam["redMobCoeff_DC_eedf"], swarmParam["characEnergy_eedf"]);
			std::fclose(fileID);

		    // plot the swarm parameters if this is the last solution
			if (currentJobID == numberOfJobs-1){
				char commandLine[1000];
				std::sprintf(commandLine, "gnuplot -s -p -e \"xLabel='%s'; %s title 'Swarm Parameters'; load 'LoKI-MC/GnuplotFunctions/swarmParams_BoltzmannMC.p'\" %s", xLabel.c_str(), gnuplotTermOptions, SUPRESS_OUTPUT);
				int systemReturn = std::system(commandLine);
			}
		}

		else{
			Eigen::Matrix3d fluxDiff = electronKinetics->totalGasDensity*electronKinetics->averagedFluxDiffusionCoeffs;
		    Eigen::Matrix3d bulkDiff = electronKinetics->totalGasDensity*electronKinetics->averagedBulkDiffusionCoeffs;
			Eigen::Array3d fluxV = Eigen::abs(electronKinetics->averagedFluxDriftVelocity);
			Eigen::Array3d bulkV = Eigen::abs(electronKinetics->averagedBulkDriftVelocity);

		    // write the headers of the auxiliary file
		    if (currentJobID == 0){
				FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "w");
				std::vector<std::string> vars = {"#inputParameter", "fluxNDxx(1/(ms))", "fluxNDyy(1/(ms))", "fluxNDzz(1/(ms))", "bulkNDxx(1/(ms))", "bulkNDyy(1/(ms))", "bulkNDzz(1/(ms))","fluxVx(m/s)","fluxVy(m/s)","fluxVz(m/s)","bulkVx(m/s)","bulkVy(m/s)", "bulkVz(m/s)", "MeanE(eV)", "ionCoeff(m-3)", "attCoeff(m-3)", "RedDif_eedf(1/(ms))","RedMob_eedf(1/(msV))","charE_eedf(eV)"}; 
				std::string header;	
				for (auto var: vars){
					header += var + std::string(26-var.size(), ' ');
				}					
				std::fprintf(fileID, "%s\n", header.c_str());		    	
		    	std::fclose(fileID);
		    }

			// write the values of the current job
			FILE* fileID = std::fopen("LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat", "a");
			std::fprintf(fileID, "%-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e %-25.10e\n", inputParameter, 
				fluxDiff(0,0), fluxDiff(1,1), fluxDiff(2,2), bulkDiff(0,0), bulkDiff(1,1), bulkDiff(2,2),
				fluxV[0], fluxV[1], fluxV[2], bulkV[0], bulkV[1], bulkV[2],
				swarmParam["meanEnergy"], swarmParam["totalIonRateCoeff"], swarmParam["totalAttRateCoeff"],  
				swarmParam["redDiffCoeff_eedf"], swarmParam["redMobCoeff_DC_eedf"], swarmParam["characEnergy_eedf"]);
			std::fclose(fileID);					    

		    // plot the swarm parameters if this is the last solution
			if (currentJobID == numberOfJobs-1){
				char commandLine[1000];
				std::sprintf(commandLine, "gnuplot -s -p -e \"xLabel='%s'; %s title 'Swarm Parameters'; load 'LoKI-MC/GnuplotFunctions/swarmParams_BoltzmannMC_nonSymmetric.p'\" %s", 
					xLabel.c_str(), gnuplotTermOptions, SUPRESS_OUTPUT);
				int systemReturn = std::system(commandLine);
			}
		}
	}
};

#endif