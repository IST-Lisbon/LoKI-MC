#ifndef __Output__
#define __Output__
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/PrescribedEedf.h"
#include "LoKI-MC/Headers/BoltzmannMC.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/Parse.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <map>
#include <boost/signals2.hpp>

template <class ElectronKineticsType>
class Output{
public:

	// ----- class attributes -----

	std::string folder;       // main output folder
	std::string subFolder;    // sub folder for output of different jobs

	bool eedfIsToBeSaved = false;
	bool evdfIsToBeSaved = false;
	bool powerBalanceIsToBeSaved = false;
	bool swarmParamsIsToBeSaved = false;
	bool rateCoeffsIsToBeSaved = false;
	bool lookUpTableIsToBeSaved = false; 
	bool MCTemporalInfoIsToBeSaved = false;
	bool MCSimDetailsIsToBeSaved = false;

	std::string eedfType;
	ElectronKineticsType* electronKinetics = NULL;
	bool enableElectronKinetics = false;
	std::string variableCondition;
	int numberOfJobs, currentJobID;

	// ----- class methods ----- //

	template <class SetupType>
	Output(SetupType* setup){
		// set output folder (if not specified in the setup, a generic folder with a timestamp is created)
		if (FieldInfo::getField("output.folder")){
			folder = std::string("Output/") + FieldInfo::getFieldValue("output.folder");
		}
		else{
			time_t now = time(0);
			std::string timestamp(ctime(&now));
			folder = std::string("Output/") + timestamp;
		}
		
		// create the directory of the output folder
		std::filesystem::create_directories(folder);

		// save the eedf type
		eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");

		// save the variable working condition and the job indices
		variableCondition = setup->workCond->variableCondition;
		numberOfJobs = setup->numberOfJobs;
		currentJobID = -1;

		// save which information must be saved
		for (auto dataFile: FieldInfo::getFieldChildNames("output.dataFiles") ){
			if (dataFile == "eedf"){
				eedfIsToBeSaved = true;
			}
			else if (dataFile == "evdf"){
				evdfIsToBeSaved = true;
			}
			else if (dataFile == "powerBalance"){
				powerBalanceIsToBeSaved = true;
			}
			else if (dataFile == "swarmParameters"){
				swarmParamsIsToBeSaved = true;
			}
			else if (dataFile == "rateCoefficients"){
				rateCoeffsIsToBeSaved = true;
			}
			else if (dataFile == "lookUpTable" && setup->numberOfJobs > 1){
				lookUpTableIsToBeSaved = true;
			}
			else if (dataFile == "MCTemporalInfo" && eedfType == "boltzmannMC"){
				MCTemporalInfoIsToBeSaved = true;
			}
			else if (dataFile == "MCSimDetails"){
				MCSimDetailsIsToBeSaved = true;
			}
		}

		// save the setup information (for reference)
		FieldInfo::saveSetupInfo(folder + "/setup.txt");

		// save the electron kinetics
		electronKinetics = setup->electronKinetics;
		enableElectronKinetics = setup->enableElectronKinetics;

		if (enableElectronKinetics){
		    // connect the signal 'obtainedNewEedfSignal' to 'Output::electronKineticsSolution', in order to update the output each time a new eedf is found
		    // idea taken from https://stackoverflow.com/questions/3047381/boost-signals-and-passing-class-method
			electronKinetics->obtainedNewEedfSignal.connect(boost::bind(&Output::electronKineticsSolution, this));
		}
	}

	void electronKineticsSolution(){

		createSubFolder();

		// save selected results of the electron kinetics
		if (eedfIsToBeSaved){
			saveEedf(electronKinetics->eedf, electronKinetics->efadf, electronKinetics->esadf, electronKinetics->energyGrid->cell);
		}
		if (evdfIsToBeSaved){
			saveEvdf(electronKinetics->evdf, electronKinetics->radialVelocityCells, electronKinetics->axialVelocityCells);
		}
		if (powerBalanceIsToBeSaved){
			savePower(electronKinetics->power);
		}
		if (swarmParamsIsToBeSaved){
			saveSwarm(electronKinetics->swarmParam, electronKinetics->workCond->reducedElecField);
		}
		if (rateCoeffsIsToBeSaved){
			saveRateCoefficients(electronKinetics->rateCoeffAll, electronKinetics->rateCoeffExtra);
		}
		if (lookUpTableIsToBeSaved){
			saveLookUpTable(electronKinetics->workCond, electronKinetics->power, electronKinetics->swarmParam);
		}
		if (MCTemporalInfoIsToBeSaved){
			saveMCTemporalInfo(electronKinetics->samplingTimes, electronKinetics->meanEnergies, electronKinetics->meanPositions, electronKinetics->positionCovariances, electronKinetics->meanVelocities);
		}
		if (MCSimDetailsIsToBeSaved){
			saveMCSimDetails();
		}
	}

	void createSubFolder(){
		++currentJobID;
		if (numberOfJobs > 1){
			// obtain the subFolder name
			char cond[100];
			if (variableCondition == "electronTemperature"){
				std::sprintf(cond, "%g",electronKinetics->workCond->electronTemperatureArray[currentJobID]);
			}
			else if (variableCondition == "reducedElecField"){
				std::sprintf(cond, "%g",electronKinetics->workCond->reducedElecFieldArray[currentJobID]);
			}
			subFolder = std::string("/") + variableCondition + "_" + cond;
			// create the subfolder directory
			std::filesystem::create_directories(folder + subFolder);
		}
		else{
			subFolder = "";
		}
	}

	void saveEedf(Eigen::ArrayXd &eedf, Eigen::ArrayXd &efadf, Eigen::ArrayXd &esadf , Eigen::ArrayXd &energy){
		
		// create file name
		std::string fileName = folder + subFolder + "/eedf.txt";
		
		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");
		int eedfSize = eedf.size();

		if (eedfType == "prescribedEedf"){
			// save information into the file
			std::fprintf(fileID, "Energy(eV)           EEDF(eV^-(3/2))\n");
			for (int i = 0; i < eedfSize; ++i){
				std::fprintf(fileID, "%-20.14e %-20.14e \n", energy[i], eedf[i]);
			}
		}
		else{
			// save information into the file       
			std::fprintf(fileID, "Energy(eV)           EEDF(eV^-(3/2))      First Anisotropy     Second Anisotropy    \n");
			for (int i = 0; i < eedfSize; ++i){
				std::fprintf(fileID, "%-20.14e %-20.14e %-20.14e %-20.14e\n", energy[i], eedf[i], efadf[i], esadf[i]);
			}
		}

		std::fclose(fileID);
	}

	void saveEvdf(Eigen::ArrayXXd &evdf, Eigen::ArrayXd &radialVelocityCells, Eigen::ArrayXd &axialVelocityCells){
		// create file name
		std::string fileName = folder + subFolder + "/evdf.txt";

		int numberOfAxialCells = radialVelocityCells.size(), numberOfRadialCells = axialVelocityCells.size();

		FILE* fileID = std::fopen(fileName.c_str(), "w");
		std::fprintf(fileID, "v_r(m/s)           v_z(m/s)           EVDF(m-3s-3)       \n");
		for (int i = 0; i < numberOfRadialCells; ++i){
			for (int j = 0; j < numberOfAxialCells; ++j){
				std::fprintf(fileID, "%-18.10e %-18.10e %-18.10e\n", radialVelocityCells[i], axialVelocityCells[j], evdf(i,j));
			}
		}
		std::fclose(fileID);
	}

	void saveSwarm(std::map<std::string,double> &swarmParam, double reducedElecField){
		// create file name
		std::string fileName = folder + subFolder + "/swarmParameters.txt";
		
		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");

		if (eedfType == "prescribedEedf"){
			// save information into the file
			std::fprintf(fileID, "        Reduced electric field = %#.14e (Td)\n", reducedElecField);
			std::fprintf(fileID, " Reduced diffusion coefficient = %#.14e (1/(ms))\n", swarmParam["redDiffCoeff"]);
			std::fprintf(fileID, "  Reduced mobility coefficient = %#.14e (1/(msV))\n", swarmParam["redMobCoeff"]);
			std::fprintf(fileID, "  Reduced Townsend coefficient = %#.14e (m2)\n", swarmParam["redTownsendCoeff"]);
			std::fprintf(fileID, "Reduced attachment coefficient = %#.14e (m2)\n", swarmParam["redAttCoeff"]);
			std::fprintf(fileID, "                   Mean energy = %#.14e (eV)\n", swarmParam["meanEnergy"]);
			std::fprintf(fileID, "         Characteristic energy = %#.14e (eV)\n", swarmParam["characEnergy"]);
			std::fprintf(fileID, "          Electron temperature = %#.14e (eV)\n", swarmParam["Te"]);
			std::fprintf(fileID, "                Drift velocity = %#.14e (m/s)\n", swarmParam["driftVelocity"]);
		}
		else if (eedfType == "boltzmannMC"){
			std::fprintf(fileID, "                    Reduced electric field = %#.14e (Td)\n", reducedElecField);
			std::fprintf(fileID, "\n%s\n\n", (std::string(35,'*') + " " + "Flux parameters" + " " + std::string(35,'*')).c_str() );
			std::fprintf(fileID, "  Reduced transverse diffusion coefficient = %#.14e (1/(ms)) ; Rel. std: %-9.3e%%\n", swarmParam["fluxRedTransvDiffCoeff"], swarmParam["fluxRedTransvDiffCoeffError"]/swarmParam["fluxRedTransvDiffCoeff"]*100.0);
			std::fprintf(fileID, "Reduced longitudinal diffusion coefficient = %#.14e (1/(ms)) ; Rel. std: %-9.3e%%\n", swarmParam["fluxRedLongDiffCoeff"], swarmParam["fluxRedLongDiffCoeffError"]/swarmParam["fluxRedLongDiffCoeff"]*100.0);
			std::fprintf(fileID, "              Reduced mobility coefficient = %#.14e (1/(msV)); Rel. std: %-9.3e%%\n", swarmParam["fluxRedMobCoeff"], swarmParam["fluxRedMobCoeffError"]/swarmParam["fluxRedMobCoeff"]*100.0);
			std::fprintf(fileID, "                     Characteristic energy = %#.14e (eV)     ; Rel. std: %-9.3e%%\n", swarmParam["fluxCharacEnergy"], swarmParam["fluxCharacEnergyError"]/swarmParam["fluxCharacEnergy"]*100.0);
			std::fprintf(fileID, "              Reduced Townsend coefficient = %#.14e (m2)\n", swarmParam["fluxRedTownsendCoeff"]);
			std::fprintf(fileID, "            Reduced attachment coefficient = %#.14e (m2)\n", swarmParam["fluxRedAttCoeff"]);
			std::fprintf(fileID, "\n");
			Eigen::ArrayXd velRelErr = electronKinetics->averagedFluxDriftVelocityError/Eigen::abs(electronKinetics->averagedFluxDriftVelocity)*100.0;
			std::fprintf(fileID, "                                   | v_x |   | %-15.8e |                       | %-9.3e%% |\n", electronKinetics->averagedFluxDriftVelocity[0], velRelErr[0]);
			std::fprintf(fileID, "                                   | v_y | = | %-15.8e |  (m/s)    ; Rel. std: | %-9.3e%% |\n", electronKinetics->averagedFluxDriftVelocity[1], velRelErr[1]);
			std::fprintf(fileID, "                                   | v_z |   | %-15.8e |                       | %-9.3e%% |\n\n", electronKinetics->averagedFluxDriftVelocity[2], velRelErr[2]);
			Eigen::ArrayXd redDiff = electronKinetics->averagedFluxDiffusionCoeffs*electronKinetics->totalGasDensity, redDiffRelErr = electronKinetics->averagedFluxDiffusionCoeffsError/Eigen::abs(electronKinetics->averagedFluxDiffusionCoeffs)*100.0;
			std::fprintf(fileID, "                      | ND_xx ND_xy ND_xz|   | %-15.8e %-15.8e %-15.8e |                      | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[0], redDiff[1], redDiff[2], redDiffRelErr[0], redDiffRelErr[1], redDiffRelErr[2]);
			std::fprintf(fileID, "                      | ND_yx ND_yy ND_yz| = | %-15.8e %-15.8e %-15.8e | (1/(ms)) ; Rel. std: | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[3], redDiff[4], redDiff[5], redDiffRelErr[3], redDiffRelErr[4], redDiffRelErr[5]);
			std::fprintf(fileID, "                      | ND_zx ND_zy ND_zz|   | %-15.8e %-15.8e %-15.8e |                      | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[6], redDiff[7], redDiff[8], redDiffRelErr[6], redDiffRelErr[7], redDiffRelErr[8]);
			std::fprintf(fileID, "\n%s\n\n", (std::string(35,'*') + " " + "Bulk parameters" + " " + std::string(35,'*')).c_str() );
			std::fprintf(fileID, "  Reduced transverse diffusion coefficient = %#.14e (1/(ms)) ; Rel. std: %-9.3e%%\n", swarmParam["bulkRedTransvDiffCoeff"], swarmParam["bulkRedTransvDiffCoeffError"]/swarmParam["bulkRedTransvDiffCoeff"]*100.0);
			std::fprintf(fileID, "Reduced longitudinal diffusion coefficient = %#.14e (1/(ms)) ; Rel. std: %-9.3e%%\n", swarmParam["bulkRedLongDiffCoeff"], swarmParam["bulkRedLongDiffCoeffError"]/swarmParam["bulkRedLongDiffCoeff"]*100.0);
			std::fprintf(fileID, "              Reduced mobility coefficient = %#.14e (1/(msV)); Rel. std: %-9.3e%%\n", swarmParam["bulkRedMobCoeff"], swarmParam["bulkRedMobCoeffError"]/swarmParam["bulkRedMobCoeff"]*100.0);
			std::fprintf(fileID, "                     Characteristic energy = %#.14e (eV)     ; Rel. std: %-9.3e%%\n", swarmParam["bulkCharacEnergy"], swarmParam["bulkCharacEnergyError"]/swarmParam["bulkCharacEnergy"]*100.0);
			std::fprintf(fileID, "              Reduced Townsend coefficient = %#.14e (m2)\n", swarmParam["bulkRedTownsendCoeff"]);
			std::fprintf(fileID, "            Reduced attachment coefficient = %#.14e (m2)\n", swarmParam["bulkRedAttCoeff"]);
			std::fprintf(fileID, "\n");
			velRelErr = electronKinetics->averagedBulkDriftVelocityError/Eigen::abs(electronKinetics->averagedBulkDriftVelocity)*100.0;
			std::fprintf(fileID, "                                   | v_x |   | %-15.8e |                       | %-9.3e%% |\n", electronKinetics->averagedBulkDriftVelocity[0], velRelErr[0]);
			std::fprintf(fileID, "                                   | v_y | = | %-15.8e |  (m/s)    ; Rel. std: | %-9.3e%% |\n", electronKinetics->averagedBulkDriftVelocity[1], velRelErr[1]);
			std::fprintf(fileID, "                                   | v_z |   | %-15.8e |                       | %-9.3e%% |\n\n", electronKinetics->averagedBulkDriftVelocity[2], velRelErr[2]);
			redDiff = electronKinetics->averagedBulkDiffusionCoeffs*electronKinetics->totalGasDensity; redDiffRelErr = electronKinetics->averagedBulkDiffusionCoeffsError/Eigen::abs(electronKinetics->averagedBulkDiffusionCoeffs)*100.0;
			std::fprintf(fileID, "                      | ND_xx ND_xy ND_xz|   | %-15.8e %-15.8e %-15.8e |                      | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[0], redDiff[1], redDiff[2], redDiffRelErr[0], redDiffRelErr[1], redDiffRelErr[2]);
			std::fprintf(fileID, "                      | ND_yx ND_yy ND_yz| = | %-15.8e %-15.8e %-15.8e | (1/(ms)) ; Rel. std: | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[3], redDiff[4], redDiff[5], redDiffRelErr[3], redDiffRelErr[4], redDiffRelErr[5]);
			std::fprintf(fileID, "                      | ND_zx ND_zy ND_zz|   | %-15.8e %-15.8e %-15.8e |                      | %-9.3e%% %-9.3e%% %-9.3e%% |\n", redDiff[6], redDiff[7], redDiff[8], redDiffRelErr[6], redDiffRelErr[7], redDiffRelErr[8]);
			std::fprintf(fileID, "\n%s\n\n", (std::string(35,'*') + " " + "Energy parameters" + " " + std::string(35,'*')).c_str() );
			std::fprintf(fileID, "                               Mean energy = %#.14e (eV) ; Rel. std: %-9.3e%%\n", swarmParam["meanEnergy"], swarmParam["meanEnergyError"]/swarmParam["meanEnergy"]*100.0);
			std::fprintf(fileID, "                      Electron temperature = %#.14e (eV) ; Rel. std: %-9.3e%%\n", swarmParam["Te"], swarmParam["TeError"]/swarmParam["Te"]*100.0);
			std::fprintf(fileID, "\n%s\n\n", (std::string(27,'*') + " " + "Parameters obtained from the EEDF" + " " + std::string(27,'*')).c_str() );
			std::fprintf(fileID, "                    Ionization coefficient = %#.14e (m-3)\n", swarmParam["totalIonRateCoeff"]);
			std::fprintf(fileID, "                    Attachment coefficient = %#.14e (m-3)\n", swarmParam["totalAttRateCoeff"]);
			std::fprintf(fileID, "      Reduced energy diffusion coefficient = %#.14e (eV/(ms))\n", swarmParam["redDiffCoeffEnergy_eedf"]);
			std::fprintf(fileID, "                   Reduced energy mobility = %#.14e (eV/(msV))\n", swarmParam["redMobCoeffEnergy_eedf"]);
			std::fprintf(fileID, "             Reduced diffusion coefficient = %#.14e (1/(ms))\n", swarmParam["redDiffCoeff_eedf"]);
			std::fprintf(fileID, "                          Reduced mobility = %#.14e (1/(msV))\n", swarmParam["redMobCoeff_DC_eedf"]);
			std::fprintf(fileID, "                     Characteristic energy = %#.14e (eV)\n\n", swarmParam["characEnergy_eedf"]);
		}

		// close file
		std::fclose(fileID);
	}

	void savePower(PowerStruct &power){
		// create file name
		std::string fileName = folder + subFolder + "/powerBalance.txt";
		
		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");

		// save information into the file
		std::fprintf(fileID, "                               Field = %#+.14e (eVm3/s)\n", power.Map["field"]);
		std::fprintf(fileID, "           Elastic collisions (gain) = %#+.14e (eVm3/s)\n", power.Map["elasticGain"]);
		std::fprintf(fileID, "           Elastic collisions (loss) = %#+.14e (eVm3/s)\n", power.Map["elasticLoss"]);
		std::fprintf(fileID, "                          CAR (gain) = %#+.14e (eVm3/s)\n", power.Map["carGain"]);
		std::fprintf(fileID, "                          CAR (loss) = %#+.14e (eVm3/s)\n", power.Map["carLoss"]);
		std::fprintf(fileID, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["excitationIne"]);
		std::fprintf(fileID, "  Excitation superelastic collisions = %#+.14e (eVm3/s)\n", power.Map["excitationSup"]);
		std::fprintf(fileID, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["vibrationalIne"]);
		std::fprintf(fileID, " Vibrational superelastic collisions = %#+.14e (eVm3/s)\n", power.Map["vibrationalSup"]);
		std::fprintf(fileID, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["rotationalIne"]);
		std::fprintf(fileID, "  Rotational superelastic collisions = %#+.14e (eVm3/s)\n", power.Map["rotationalSup"]);
		std::fprintf(fileID, "               Ionization collisions = %#+.14e (eVm3/s)\n", power.Map["ionizationIne"]);
		std::fprintf(fileID, "               Attachment collisions = %#+.14e (eVm3/s)\n", power.Map["attachmentIne"]);
		std::fprintf(fileID, "             Electron density growth = %#+.14e (eVm3/s) +\n", power.Map["eDensGrowth"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "                       Power Balance = %#+.14e (eVm3/s)\n", power.Map["balance"]);
		std::fprintf(fileID, "              Relative Power Balance = %#.14e%%\n\n", power.Map["relativeBalance"]*100.0);
		std::fprintf(fileID, "           Elastic collisions (gain) = %#+.14e (eVm3/s)\n", power.Map["elasticGain"]);
		std::fprintf(fileID, "           Elastic collisions (loss) = %#+.14e (eVm3/s) +\n", power.Map["elasticLoss"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "            Elastic collisions (net) = %#+.14e (eVm3/s)\n\n", power.Map["elasticNet"]);
		std::fprintf(fileID, "                          CAR (gain) = %#+.14e (eVm3/s)\n", power.Map["carGain"]);
		std::fprintf(fileID, "                          CAR (gain) = %#+.14e (eVm3/s) +\n", power.Map["carLoss"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "                           CAR (net) = %#+.14e (eVm3/s)\n\n", power.Map["carNet"]);
		std::fprintf(fileID, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["excitationIne"]);
		std::fprintf(fileID, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", power.Map["excitationSup"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", power.Map["excitationNet"]);
		std::fprintf(fileID, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["vibrationalIne"]);
		std::fprintf(fileID, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", power.Map["vibrationalSup"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", power.Map["vibrationalNet"]);
		std::fprintf(fileID, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", power.Map["rotationalIne"]);
		std::fprintf(fileID, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", power.Map["rotationalSup"]);
		std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
		std::fprintf(fileID, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n", power.Map["rotationalNet"]);

		// power balance by gases
		for (auto gas: Parse::getMapKeys(power.gasesMap)){
			// assign the power map of the gas
			std::map<std::string,double> gasPowerMap = power.gasesMap[gas];

			std::fprintf(fileID, "\n%s\n\n", (std::string(37,'*') + " " + gas + " " + std::string(39-gas.size(),'*')).c_str() );
			std::fprintf(fileID, "     Excitation inelastic collisions = %#+.14e (eVm3/s)\n", gasPowerMap["excitationIne"]);
			std::fprintf(fileID, "  Excitation superelastic collisions = %#+.14e (eVm3/s) +\n", gasPowerMap["excitationSup"]);
			std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
			std::fprintf(fileID, "         Excitation collisions (net) = %#+.14e (eVm3/s)\n\n", gasPowerMap["excitationNet"]);
			std::fprintf(fileID, "    Vibrational inelastic collisions = %#+.14e (eVm3/s)\n", gasPowerMap["vibrationalIne"]);
			std::fprintf(fileID, " Vibrational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPowerMap["vibrationalSup"]);
			std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
			std::fprintf(fileID, "        Vibrational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPowerMap["vibrationalNet"]);
			std::fprintf(fileID, "     Rotational inelastic collisions = %#+.14e (eVm3/s)\n", gasPowerMap["rotationalIne"]);
			std::fprintf(fileID, "  Rotational superelastic collisions = %#+.14e (eVm3/s) +\n", gasPowerMap["rotationalSup"]);
			std::fprintf(fileID, " %s\n", std::string(73,'-').c_str());
			std::fprintf(fileID, "         Rotational collisions (net) = %#+.14e (eVm3/s)\n\n", gasPowerMap["rotationalNet"]);
			std::fprintf(fileID, "               Ionization collisions = %#+.14e (eVm3/s)\n", gasPowerMap["ionizationIne"]);
			std::fprintf(fileID, "               Attachment collisions = %#+.14e (eVm3/s)\n", gasPowerMap["attachmentIne"]);
		}

		// close file
		std::fclose(fileID);
	}

	void saveRateCoefficients(std::vector<RateCoeffStruct> &rateCoeffAll, std::vector<RateCoeffStruct> &rateCoeffExtra){
		// create file name
		std::string fileName = folder + subFolder + "/rateCoefficients.txt";

		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");

		// save information into the file
		std::fprintf(fileID, " ID  Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)   Description\n");

		for (auto rateCoeff: rateCoeffAll){
			if (rateCoeff.supRate == Constant::NON_DEF){
				// we sum 1 to the ID so as to have the same ID as in Matlab
				std::fprintf(fileID, "%4d %20.14e (N/A)                %s\n", rateCoeff.collID+1, rateCoeff.ineRate, rateCoeff.collDescription.c_str());			
			}
			else{
	          std::fprintf(fileID, "%4d %20.14e %20.14e %s\n", rateCoeff.collID+1, rateCoeff.ineRate, rateCoeff.supRate, rateCoeff.collDescription.c_str());			
			}
		}

		if (!rateCoeffExtra.empty()){
			std::fprintf(fileID, "\n%s\n* Extra Rate Coefficients *\n%s\n\n", std::string(27,'*').c_str(), std::string(27,'*').c_str());
			for (auto rateCoeff: rateCoeffExtra){
				if (rateCoeff.supRate == Constant::NON_DEF){
					// we sum 1 to the ID so as to have the same ID as in Matlab
					std::fprintf(fileID, "%4d %20.14e (N/A)                %s\n", rateCoeff.collID+1, rateCoeff.ineRate, rateCoeff.collDescription.c_str());			
				}
				else{
		          std::fprintf(fileID, "%4d %20.14e %20.14e %s\n", rateCoeff.collID+1, rateCoeff.ineRate, rateCoeff.supRate, rateCoeff.collDescription.c_str());			
				}
			}		
		}

		// close file
		std::fclose(fileID);

		if (eedfType == "boltzmannMC"){
			// create filename
			fileName = folder + subFolder + "/rateCoefficientsMC.txt";

			// open file
			fileID = std::fopen(fileName.c_str(), "w");

			// save information into the file
			std::fprintf(fileID, " ID  Ine.R.Coeff.(m3/s)   Sup.R.Coeff.(m3/s)   Description\n");
			for (auto rateCoeff: rateCoeffAll){
				if (rateCoeff.collDescription.find("Effective") != std::string::npos){
					continue;
				}
				else if (rateCoeff.supRate == Constant::NON_DEF){
					// we sum 1 to the ID so as to have the same ID as in Matlab
					std::fprintf(fileID, "%4d %20.14e (N/A)                %s\n", rateCoeff.collID+1, rateCoeff.ineRateMC, rateCoeff.collDescription.c_str());		
				}
				else{
		          	std::fprintf(fileID, "%4d %20.14e %20.14e %s\n", rateCoeff.collID+1, rateCoeff.ineRateMC, rateCoeff.supRateMC, rateCoeff.collDescription.c_str());		          				
				}
			}

			// close file
			std::fclose(fileID);
		}
	}

	void saveLookUpTable(WorkingConditions* &workCond, PowerStruct &power, std::map<std::string,double> &swarmParam){
		static bool initialized = false;

		if (eedfType == "prescribedEedf"){
			// create file name
			std::string fileName = folder + "/lookUpTable.txt";
			// initialize the file in case it is needed
			if (!initialized){
				// open file
				FILE* fileID = std::fopen(fileName.c_str(), "w");
				// write file headers
		        std::fprintf(fileID, "RedField(Td)         RedDif(1/(ms))       RedMob(1/(msV))      RedTow(m2)           RedAtt(m2)           MeanE(eV)            CharE(eV)            EleTemp(eV)          DriftVelocity(m/s)   RelativePowerBalance\n");
		        // close file
		        std::fclose(fileID);
		        initialized = true;
			}

			// open file
			FILE* fileID = std::fopen(fileName.c_str(), "a");
			// append new line with data
		    std::fprintf(fileID, "%20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %20.14e %19.14e%%\n", workCond->reducedElecField, swarmParam["redDiffCoeff"], swarmParam["redMobCoeff"], 
		        swarmParam["redTownsendCoeff"], swarmParam["redAttCoeff"], swarmParam["meanEnergy"], swarmParam["characEnergy"], swarmParam["Te"], swarmParam["driftVelocity"], power.Map["relativeBalance"]*100.0);
		   	// close file
		    std::fclose(fileID);
		}

		else{
			std::string fileName1 = folder + "/lookUpTableSwarm.txt";
			std::string fileName2 = folder + "/lookUpTablePower.txt";
			std::string fileName3 = folder + "/lookUpTableRateCoeff.txt";
			std::string fileName4 = folder + "/lookUpTableEedf.txt";
			std::string fileName5 = folder + "/lookUpTableEGrid.txt";

			if (!initialized){
				std::string varCond = "RedElecField(Td)      ";				
				// ---- open files ---- //
				FILE* fileID1 = std::fopen(fileName1.c_str(), "w");
				FILE* fileID2 = std::fopen(fileName2.c_str(), "w");
				FILE* fileID3 = std::fopen(fileName3.c_str(), "w");
				FILE* fileID4 = std::fopen(fileName4.c_str(), "w");
				FILE* fileID5 = std::fopen(fileName5.c_str(), "w");

				// ---- write file headers ---- //

				// swarm file (1)
				std::vector<std::string> swarmVars;
				swarmVars = {"FluxND_xx(1/(ms))","FluxND_yy(1/(ms))","FluxND_zz(1/(ms))","FluxV_x(m/s)","FluxV_y(m/s)","FluxV_z(m/s)","BulkND_xx(1/(ms))","BulkND_yy(1/(ms))","BulkND_zz(1/(ms))","BulkV_x(m/s)","BulkV_y(m/s)","BulkV_z(m/s)","MeanE(eV)","EleTemp(eV)","IonCoeff(m3/s)","AttCoeff(m3/s)","RedDiffE_eedf(eV/(ms))","RedMobE_eedf(eV/(msV))","RedDiff_eedf(1/(ms))", "RedMob_DC_eedf(1/(msV))","CharacEnergy_eedf(eV)"};
				std::string swarmHeader = varCond;
				for (auto var: swarmVars){
					swarmHeader += var + std::string(24-var.size(), ' ');
				}
				std::fprintf(fileID1, "%s\n", swarmHeader.c_str());

				// power file (2)
				std::vector<std::string> powerVars = {"PowerField(eVm3/s)","PwrElaGain(eVm3/s)","PwrElaLoss(eVm3/s)","PwrElaNet(eVm3/s)","PwrEleGain(eVm3/s)","PwrEleLoss(eVm3/s)","PwrEleNet(eVm3/s)","PwrVibGain(eVm3/s)","PwrVibLoss(eVm3/s)","PwrVibNet(eVm3/s)","PwrRotGain(eVm3/s)","PwrRotLoss(eVm3/s)","PwrRotNet(eVm3/s)","PwrIon(eVm3/s)","PwrAtt(eVm3/s)","PwrGrowth(eVm3/s)","PwrBalance(eVm3/s)","RelPwrBalance"};
				std::string powerHeader = varCond;
				for (auto var: powerVars){
					powerHeader += var + std::string(22-var.size(), ' ');
				}
				std::fprintf(fileID2, "%s\n", powerHeader.c_str());

				// rate-coeff file (3)
				std::fprintf(fileID3, "%s\n# %-76s #\n", std::string(80, '#').c_str(), "ID   Description");
				std::string rateCoeffHeader = varCond, strAux;
				for (auto rateCoeff: electronKinetics->rateCoeffAll){
					// we sum 1 to the ID so as to have the same ID as in Matlab
					int ID = rateCoeff.collID+1;
					std::fprintf(fileID3, "# %-4d %-71s #\n", ID, rateCoeff.collDescription.c_str());
					strAux = std::string("R") + std::to_string(ID) + "_ine(m3/s)";
					rateCoeffHeader += strAux + std::string(22-strAux.size(), ' ');
					if (rateCoeff.supRate != Constant::NON_DEF){
						strAux = std::string("R") + std::to_string(ID) + "_sup(m3/s)";
						rateCoeffHeader += strAux + std::string(22-strAux.size(), ' ');					
					}
				}
				std::fprintf(fileID3, "#%s#\n# %-76s #\n#%s#\n# %-76s #\n", std::string(78, ' ').c_str(), "*** Extra rate coefficients ***", std::string(78, '#').c_str(), "ID   Description");
				for (auto rateCoeff: electronKinetics->rateCoeffExtra){
					// we sum 1 to the ID so as to have the same ID as in Matlab
					int ID = rateCoeff.collID+1;
					std::fprintf(fileID3, "# %-4d %-71s #\n", ID, rateCoeff.collDescription.c_str());
					strAux = std::string("R") + std::to_string(ID) + "_ine(m3/s)";
					rateCoeffHeader += strAux + std::string(22-strAux.size(), ' ');
					if (rateCoeff.supRate != Constant::NON_DEF){
						strAux = std::string("R") + std::to_string(ID) + "_sup(m3/s)";
						rateCoeffHeader += strAux + std::string(22-strAux.size(), ' ');					
					}
				}
				std::fprintf(fileID3, "%s\n\n%s\n", std::string(80,'#').c_str(), rateCoeffHeader.c_str());				

				std::fclose(fileID1);
				std::fclose(fileID2);
				std::fclose(fileID3);
				std::fclose(fileID4);
				std::fclose(fileID5);
				initialized = true;
			}

			// open files
			FILE* fileID1 = std::fopen(fileName1.c_str(), "a");
			FILE* fileID2 = std::fopen(fileName2.c_str(), "a");
			FILE* fileID3 = std::fopen(fileName3.c_str(), "a");
			FILE* fileID4 = std::fopen(fileName4.c_str(), "a");
			FILE* fileID5 = std::fopen(fileName5.c_str(), "a");

			// check the type of variable condition
			double conditionValue = workCond->reducedElecField;

			// append new line with swarm data
			Eigen::ArrayXd fluxRedDiff = electronKinetics->averagedFluxDiffusionCoeffs*electronKinetics->totalGasDensity;
			Eigen::ArrayXd bulkRedDiff = electronKinetics->averagedBulkDiffusionCoeffs*electronKinetics->totalGasDensity;
		    std::fprintf(fileID1, "%-21.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e %-23.14e\n", conditionValue, 
		    fluxRedDiff[0], fluxRedDiff[4], fluxRedDiff[8], electronKinetics->averagedFluxDriftVelocity[0], electronKinetics->averagedFluxDriftVelocity[1], electronKinetics->averagedFluxDriftVelocity[2] , 
		    bulkRedDiff[0], bulkRedDiff[4], bulkRedDiff[8], electronKinetics->averagedBulkDriftVelocity[0], electronKinetics->averagedBulkDriftVelocity[1], electronKinetics->averagedBulkDriftVelocity[2],
		    swarmParam["meanEnergy"], swarmParam["Te"], swarmParam["totalIonRateCoeff"], swarmParam["totalAttRateCoeff"], swarmParam["redDiffCoeffEnergy_eedf"], swarmParam["redMobCoeffEnergy_eedf"],
		    swarmParam["redDiffCoeff_eedf"], swarmParam["redMobCoeff_DC_eedf"], swarmParam["characEnergy_eedf"]);				

		    // append new line with power data
		    std::fprintf(fileID2, "%-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %-21.14e %19.14e%%\n", conditionValue, 
		    power.Map["field"], power.Map["elasticGain"], power.Map["elasticLoss"], power.Map["elasticNet"], power.Map["excitationSup"], power.Map["excitationIne"], power.Map["excitationNet"],
		    power.Map["vibrationalSup"], power.Map["vibrationalIne"], power.Map["vibrationalNet"], power.Map["rotationalSup"], power.Map["rotationalIne"], power.Map["rotationalNet"],
		    power.Map["ionizationIne"], power.Map["attachmentIne"], power.Map["eDensGrowth"], power.Map["balance"], power.Map["relativeBalance"]*100.0);

		    // append new line with rate-coeff data
		    std::fprintf(fileID3, "%-21.14e ", conditionValue);
			for (auto rateCoeff: electronKinetics->rateCoeffAll){
				std::fprintf(fileID3, "%-21.14e ", rateCoeff.ineRate);
				if (rateCoeff.supRate != Constant::NON_DEF){
					std::fprintf(fileID3, "%-21.14e ", rateCoeff.supRate);				
				}
			}
			for (auto rateCoeff: electronKinetics->rateCoeffExtra){
				std::fprintf(fileID3, "%-21.14e ", rateCoeff.ineRate);
				if (rateCoeff.supRate != Constant::NON_DEF){
					std::fprintf(fileID3, "%-21.14e ", rateCoeff.supRate);				
				}
			}	     
			std::fprintf(fileID3, "\n");

			// append line with the current eedf
			// append line with the correspondent energy grid
			Eigen::ArrayXd eedf = electronKinetics->eedf;
			Eigen::ArrayXd cell = electronKinetics->energyGrid->cell;
			int nCells = cell.size();
			for (int i = 0; i < nCells; ++i){
				std::fprintf(fileID4, "%-20.14e ", eedf[i]);
				std::fprintf(fileID5, "%-20.14e ", cell[i]);	
			}
			std::fprintf(fileID4,"\n");
			std::fprintf(fileID5,"\n");

			// close files
			std::fclose(fileID1);
			std::fclose(fileID2);
			std::fclose(fileID3);
			std::fclose(fileID4);
			std::fclose(fileID5);	
		}
	}


	void saveMCTemporalInfo(Eigen::ArrayXd &samplingTimes, Eigen::ArrayXd &meanEnergies, Eigen::ArrayXXd &meanPositions, Eigen::ArrayXXd &positionCovariances, Eigen::ArrayXXd &meanVelocities){
		std::string fileName = folder + subFolder + "/MCTemporalInfo.txt";

		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");

		// write the time-dependent data in an auxiliary file
		std::fprintf(fileID, "Time(s)              MeanEnergy(eV)      xPos(m)             yPos(m)             zPos(m)             xSqWidth(m2)        ySqWidth(m2)        zSqWidth(m2)        xVel(m/s)           yVel(m/s)           zVel(m/s)           \n");
		int nSamplingPoints = electronKinetics->nSamplingPoints;
		for (int i = 0; i < nSamplingPoints; ++i){
			std::fprintf(fileID, "%#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e %#19.12e\n", samplingTimes[i], meanEnergies[i], meanPositions(i, 0), meanPositions(i, 1), meanPositions(i, 2), 
				    positionCovariances(i,0), positionCovariances(i,4), positionCovariances(i,8), meanVelocities(i,0), meanVelocities(i,1), meanVelocities(i,2));
		}
		std::fclose(fileID);
	}

	void saveMCSimDetails(){

		// create file name
		std::string fileName = folder + subFolder + "/MCSimDetails.txt";

		// open file
		FILE* fileID = std::fopen(fileName.c_str(), "w");

		// save information into the file
		std::fprintf(fileID, "                          number of electrons: %e\n", electronKinetics->nElectrons);
		std::fprintf(fileID, "                        final simulation time: %e s\n", electronKinetics->time);
		std::fprintf(fileID, "                            steady-state time: %e s\n", electronKinetics->steadyStateTime);
		std::fprintf(fileID, "                 number of integration points: %d\n\n", electronKinetics->nIntegrationPoints);
		std::fprintf(fileID, "********************************** Collisions ************************************\n\n");
		std::fprintf(fileID, "number of real collisions before steady-state: %e\n", (double)electronKinetics->collisionCounterAtSS);
		std::fprintf(fileID, " number of real collisions after steady-state: %e\n", (double)electronKinetics->totalCollisionCounter-electronKinetics->collisionCounterAtSS);
		std::fprintf(fileID, "----------------------------------------------------------------------------------\n");
		std::fprintf(fileID, "              total number of real collisions: %e\n\n", (double)electronKinetics->totalCollisionCounter);
		std::fprintf(fileID, "number of null collisions before steady-state: %e\n", (double)electronKinetics->nullCollisionCounterAtSS);
		std::fprintf(fileID, " number of null collisions after steady-state: %e\n", (double)electronKinetics->nullCollisionCounter-electronKinetics->nullCollisionCounterAtSS);
		std::fprintf(fileID, "----------------------------------------------------------------------------------\n");
		std::fprintf(fileID, "              total number of null collisions: %e\n\n", (double)electronKinetics->nullCollisionCounter);
		std::fprintf(fileID, "                  fraction of real collisions: %e%%\n", (double)electronKinetics->totalCollisionCounter/(electronKinetics->totalCollisionCounter+electronKinetics->nullCollisionCounter)*100.0);
		std::fprintf(fileID, "                  fraction of null collisions: %e%%\n", (double)electronKinetics->nullCollisionCounter/(electronKinetics->totalCollisionCounter+electronKinetics->nullCollisionCounter)*100.0);
		std::fprintf(fileID, "----------------------------------------------------------------------------------\n");
		std::fprintf(fileID, "                                        total: %e%%\n\n", 100.0);
		std::fprintf(fileID, "************************** Simulation relative errors ****************************\n\n");
		std::fprintf(fileID, "                                  Mean energy: %.4e\n", electronKinetics->averagedMeanEnergyError/electronKinetics->averagedMeanEnergy);
		Eigen::ArrayXd relError = electronKinetics->averagedFluxDriftVelocityError/Eigen::abs(electronKinetics->averagedFluxDriftVelocity);
		std::fprintf(fileID, "                          Flux drift velocity: %.4e %.4e %.4e\n", relError[0], relError[1], relError[2]);
		relError = electronKinetics->averagedFluxDiffusionCoeffsError/Eigen::abs(electronKinetics->averagedFluxDiffusionCoeffs);
		std::fprintf(fileID, "                  Flux diffusion coefficients: %.4e %.4e %.4e\n", relError[0], relError[4], relError[8]);
		std::fprintf(fileID, "                                Power balance: %.4e\n", electronKinetics->powerBalanceRelError);
		std::fprintf(fileID, "**********************************************************************************\n\n");
		std::fprintf(fileID, "                                 Elapsed time: %e s\n", electronKinetics->elapsedTime);
		std::fclose(fileID);
	}
};

#endif