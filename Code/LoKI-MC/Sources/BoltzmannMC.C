#include "LoKI-MC/Headers/BoltzmannMC.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/GeneralDefinitions.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <cstdio>
#include <vector>
#include <string>
#include <map>
#include <boost/signals2.hpp>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fit.h>
#include <random>
#include <fstream>
#include <ctime>

// the constructor is written in the '.h' file to avoid compilation problems related with the templates

void BoltzmannMC::allocateEvaluateVariablesFirstTime(){
	// 'allocateEvaluateVariablesFirstTime' evaluate and allocate variables when the constructor is called

	// evaluate nProcesses and nGases
	nProcesses = 0;
	nGases = 0;
	for (auto& gas: gasArray){
		if (!gas->collisionArray.empty()){
			++nGases;
			for (auto& collision: gas->collisionArray){
				if (collision->type == "Effective"){ // avoid effective collisions
					continue;
				}
				++nProcesses;
				if (collision->isReverse){
					++nProcesses;
				}
			}
		}
	}

	// reallocate the memory according with the number of processes (remember that here we separate inelastics from superelastics)
	processTypes = new int[nProcesses];
	realCollisionPointers.resize(nProcesses, NULL);
	gasPointers.resize(nProcesses, NULL);
	isElastic = new bool[nProcesses];
	isSuperElastic = new bool[nProcesses];
	isIonization = new bool[nProcesses];
	superElasticStatWeightFactors = new double[nProcesses];
	crossSectionEnergies.resize(nProcesses); crossSectionValues.resize(nProcesses);
	energyMinLimits = new double[nProcesses]; energyMaxLimits = new double[nProcesses];
	crossSectionInterpolations.resize(nProcesses, NULL); crossSectionInterpAccelerators.resize(nProcesses, NULL);
	relDensities = new double[nProcesses]; targetMasses = new double[nProcesses]; reducedMasses = new double[nProcesses]; energyLosses = new double[nProcesses]; thermalStdDeviations = new double[nProcesses];
	crossSectionEnergyGrid = new double[interpolCrossSectionSize];
	interpolCrossSectionsXrelDens = new double*[interpolCrossSectionSize];
	cumulSumInterpolCrossSectionsXrelDens = new double*[interpolCrossSectionSize];
	for (int i = 0; i < interpolCrossSectionSize; ++i){
		interpolCrossSectionsXrelDens[i] = new double[nProcesses];
		cumulSumInterpolCrossSectionsXrelDens[i] = new double[nProcesses];
	}
	totalCollisionFrequencies = new double[interpolCrossSectionSize];
	maxCollisionFrequencies = new double[interpolCrossSectionSize];
	energyGainProcesses = new double[nProcesses]; energyLossProcesses = new double[nProcesses];
	collisionCounters = new unsigned long long int[nProcesses];
	averagedRateCoeffs = new double[nProcesses];
	averagedPowerGainProcesses = new double[nProcesses]; averagedPowerLossProcesses = new double[nProcesses];
	wParameters = new double[nProcesses];
	targetGasIDs = new int[nProcesses];
	isMomentumConservationIonizationScattering = new bool[nProcesses];	
	angularScatteringFunctions.resize(nProcesses, NULL);
	firstProcessIndexPerGas = new double [nGases];
	lastProcessIndexPerGas = new double [nGases];
	gasFractions = new double [nGases];	
	chosenProcessIDs = new int[(int)nElectrons];
	ejectedElectronPositions.resize(nElectrons,3); ejectedElectronVelocities.resize(nElectrons,3); 
	ejectedElectronEnergies.resize(nElectrons);
	electronEnergyChanges.resize(nElectrons);
	electronEnergyChangesOverIncidEnergies.resize(nElectrons);
	energyGainsField.resize(nElectrons);

	// assign the data for each MC process (remember that here we separate inelastics from superelastics)
	int iterProcess = 0;
	int iterGas = 0;

	for (auto& gas: gasArray){

		if (gas->collisionArray.empty()){
			continue;
		}	
		firstProcessIndexPerGas[iterGas] = iterProcess;
		gasFractions[iterGas] = gas->fraction; 

		for (auto& collision: gas->collisionArray){

			// assign process type 
			if (collision->type == "Effective"){ // avoid effective collisions
				continue;
			}
			else if (collision->type == "Ionization"){
				processTypes[iterProcess] = GeneralDefinitions::ionizationType;
				isIonization[iterProcess] = true;
			}
			else if (collision->type == "Attachment"){
				processTypes[iterProcess] = GeneralDefinitions::attachmentType;
				isIonization[iterProcess] = false;
			}
			else{
				processTypes[iterProcess] = GeneralDefinitions::conservativeType;
				isIonization[iterProcess] = false;
				if (collision->type == "Elastic"){
					isElastic[iterProcess] = true;
				}
			}

			// assign 'Collision' pointer
			realCollisionPointers[iterProcess] = collision;
			// assign 'Gas' pointer
			gasPointers[iterProcess] = gas;

			// assign the process direction
			isSuperElastic[iterProcess] = false;

			// assign the superelastic factor (set to zero, since this an inelastic process)
			superElasticStatWeightFactors[iterProcess] = 0;

			// Assign the cross section data
			std::vector<double> tempEnergyVector = collision->rawIntegralCrossSection[0];
			std::vector<double> tempValueVector = collision->rawIntegralCrossSection[1];
			double tempThreshold = collision->threshold;
			// eliminate the points with energy smaller than the threshold
			while (tempEnergyVector[0] < tempThreshold){
				tempEnergyVector.erase(tempEnergyVector.begin());
				tempValueVector.erase(tempValueVector.begin());
			}
			// if there is not a point with the threshold energy, add it
			if (tempEnergyVector[0] != tempThreshold){
				tempEnergyVector.insert(tempEnergyVector.begin(), tempThreshold);
				tempValueVector.insert(tempValueVector.begin(), 0);
			}
			// convert vector to array
			crossSectionEnergies[iterProcess] = MathFunctions::vectorToArray(tempEnergyVector);
			crossSectionValues[iterProcess] = MathFunctions::vectorToArray(tempValueVector);

			// assign the min energy limit
			energyMinLimits[iterProcess] = tempThreshold;
			// assign the max energy limit
			energyMaxLimits[iterProcess] = tempEnergyVector.back();
			// compare with energyMaxElastic
			if (collision->type == "Elastic"){
				energyMaxElastic = std::fmin(energyMaxElastic, energyMaxLimits[iterProcess]);
			}

			// initialize the interpolation objects
			crossSectionInterpolations[iterProcess] = gsl_spline_alloc(gsl_interp_linear, crossSectionEnergies[iterProcess].size());
			gsl_spline_init(crossSectionInterpolations[iterProcess], crossSectionEnergies[iterProcess].data(), crossSectionValues[iterProcess].data(), crossSectionEnergies[iterProcess].size());
			crossSectionInterpAccelerators[iterProcess] = gsl_interp_accel_alloc();

			// save the masses and thresholds
			targetMasses[iterProcess] = collision->target->gas->mass;
			reducedMasses[iterProcess] = Constant::electronMass*targetMasses[iterProcess] / (Constant::electronMass+targetMasses[iterProcess]);	
			energyLosses[iterProcess] = collision->threshold;
			thermalStdDeviations[iterProcess] = std::sqrt(Constant::boltzmann*gasTemperature/targetMasses[iterProcess]);

			// assign the w parameter, if it is an ionization process. Important when using SDCS
			if (collision->type == "Ionization"){
				if (collision->target->gas->OPBParameter == Constant::NON_DEF){
					wParameters[iterProcess] = collision->threshold;
				}
				else{
					wParameters[iterProcess] = collision->target->gas->OPBParameter;
				}
			}
			else{
				wParameters[iterProcess] = 0; // not important in this case
			}

			// assign the target gas ID
			targetGasIDs[iterProcess] = collision->target->gas->ID;

			// assign the angular scattering function
			angularScatteringFunctions[iterProcess] = AngularScatteringFunctions::functionMap(collision->angularScatteringType, collision->angularScatteringParams, this);

			// assign the angular scattering parameters
			angularScatteringParams.push_back(collision->angularScatteringParams);

			// check if the angular scattering function is 'momentumConservationIonization'
			if (collision->angularScatteringType == "momentumConservationIonization"){
				if (!isIonization[iterProcess]){
					Message::error(std::string("Trying to assign the angularScatteringType 'momentumConservationIonization' to the process\n") + collision->description() + "\nwhich is not 'Ionization'");
				}
				isMomentumConservationIonizationScattering[iterProcess] = true;
			}
			else{
				isMomentumConservationIonizationScattering[iterProcess] = false;
			}

			// increment the iterator
			++iterProcess;

			// enter here if there is a superelastic
			if (collision->isReverse){
				// assign process type
				processTypes[iterProcess] = GeneralDefinitions::conservativeType;

				// assign 'Collision' pointer
				realCollisionPointers[iterProcess] = collision;
				// assign 'Gas' pointer
				gasPointers[iterProcess] = gas;	

				// assign the elastic boolean
				isElastic[iterProcess] = false;		
				
				// assign the process direction
				isSuperElastic[iterProcess] = true;

				// assign the ionization boolean
				isIonization[iterProcess] = false;

				// assign the superelastic factor 
				superElasticStatWeightFactors[iterProcess] = collision->target->statisticalWeight/collision->productArray[0]->statisticalWeight;

				// assign the min energy limit
				energyMinLimits[iterProcess] = 0;
				// assign the max energy limit
				energyMaxLimits[iterProcess] = energyMaxLimits[iterProcess-1]-collision->threshold;

				// Assign the cross section data (not relevant for superelastics)
				crossSectionEnergies[iterProcess] = Eigen::ArrayXd::Zero(1);
				// get the superelastic cross section values (not relevant for superelastics)
				crossSectionValues[iterProcess] = Eigen::ArrayXd::Zero(1);

				// initialize the interpolation objects (not relevant for superelastics)
				crossSectionInterpolations[iterProcess] = NULL;
				crossSectionInterpAccelerators[iterProcess] = NULL;

				// save the masses and thresholds
				targetMasses[iterProcess] = collision->productArray[0]->gas->mass;
				reducedMasses[iterProcess] = Constant::electronMass*targetMasses[iterProcess] / (Constant::electronMass+targetMasses[iterProcess]);	
				energyLosses[iterProcess] = -collision->threshold;
				thermalStdDeviations[iterProcess] = std::sqrt(Constant::boltzmann*gasTemperature/targetMasses[iterProcess]);

				// assign the wParameter, although it is only used for ionization
				wParameters[iterProcess] = 0;

				// assign the target gas ID
				targetGasIDs[iterProcess] = collision->productArray[0]->gas->ID;

				// assign the angular scattering function
				angularScatteringFunctions[iterProcess]	= angularScatteringFunctions[iterProcess-1];

				// assign the angular scattering parameters
				angularScatteringParams.push_back(angularScatteringParams[iterProcess-1]);				
				
				// not relevant in superelastics
				isMomentumConservationIonizationScattering[iterProcess] = false;			

				++iterProcess;
			}
		}
		lastProcessIndexPerGas[iterGas] = iterProcess-1;
		++iterGas;
	}
}

void BoltzmannMC::solve(){

	// save the initial time
	const auto start = std::chrono::high_resolution_clock::now();

	// evaluate eedf
	evaluateEEDF();

	// evaluate power balance
	evaluatePower();

	// evaluate rate coefficients
	evaluateRateCoeff();

	// evaluate swarm parameters
	evaluateSwarmParameters();

	// save the final time
	const auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	elapsedTime = (elapsed.count())*1E-3;

	// broadcast obtention of a solution for the EEDF
	obtainedNewEedfSignal();
}

void BoltzmannMC::evaluateEEDF(){

	// ----- evaluate the variables that may be different each time this function is called ----- //

	evaluateNonConstantVariables();

	// calculate the mean data for the initial time
	nSamplingPoints = 1;
	nSynchronizationPoints = 1;
	currentSamplingIndex = 0;
	samplingTimes[currentSamplingIndex] = 0;
	calculateMeanDataForSwarmParams();
	collisionCounterAfterSS = 0;
	nIntegrationPoints = 0;
	if (dispMCStatus){
		dispInfo();
	}

	// ----- perform the Monte Carlo algorithm ----- //

	// the MC calculations are performed until one of these criteria is fulfilled: max collisions reached, required statistical errors, required integration points, required integrated steaty-state times
	while ((!goodStatisticalErrors && errorsToBeChecked) || nIntegrationPoints < requiredIntegrationPoints ||
		totalIntegratedTime/steadyStateTime < requiredIntegratedSSTimes  || totalIntegratedTime < requiredIntegratedAbsoluteTime){

		if (collisionCounterAfterSS >= maxCollisionsAfterSteadyState && nIntegrationPoints > 100){
			break;
		}

		// perform electron dynamics (accellerations and collisions) until synchronization time
		electronDynamicsUntilSynchronization();

		// perform sampling at the synchronization instant
		if (nSynchronizationPoints % synchronizationOverSampling == 0){
	 						
	 		// update the number of sampling points
			++nSamplingPoints;
			currentSamplingIndex = nSamplingPoints-1;

			// update the sampling times
			samplingTimes[currentSamplingIndex] = time;

			// calculate average data for swarm parameters
			calculateMeanDataForSwarmParams();

			// define the number of points between steady-state checks. It increases for larger number of sampling points, to avoid unnecessary checks
			int deucades = std::fmax(std::log(nSamplingPoints)/std::log(2.0)-11, 6);
			nPointsBetweenSteadyStateCheck = std::pow(2,deucades);

			// calculate the swarm parameters if the steady-state was already achieved
			if (steadyStateTime != Constant::NON_DEF){

				// update the number of points used for integration of the swarm parameters
				++nIntegrationPoints;
				collisionCounterAfterSS = totalCollisionCounter - collisionCounterAtSS;

				// update the integration times
				timeIntervalInteg = time - totalIntegratedTime;
				totalIntegratedTime = time - steadyStateTime;

				// get the time-dependent distribution functions
				getTimeDependDistributions();
			}
			// check if the system has reached the steady state, taking into account the temporal evolution of the mean kinetic energy
			else{
				if (nSamplingPoints >= 100 && nSamplingPoints % nPointsBetweenSteadyStateCheck == 0 && totalCollisionCounter > minCollisionsBeforeSteadyState){
					checkSteadyState();
				}
			}		

			// define the number of points between stat-error checks. It increases for larger number of integration points, to avoid unnecessary checks
			deucades = std::fmax(std::log(nIntegrationPoints)/std::log(2.0)-9, 7);
			nPointsBetweenStatErrorsCheck = std::pow(2,deucades);

			// check statistical errors
			if (steadyStateTime != Constant::NON_DEF && nIntegrationPoints > 200 && nIntegrationPoints % nPointsBetweenStatErrorsCheck == 0){
				checkStatisticalErrors();
			}

			nPointsBetweenDispInfo = nPointsBetweenSteadyStateCheck;
			// display info
			if (dispMCStatus && nSamplingPoints % nPointsBetweenDispInfo == 0){
				dispInfo();
			}
		}

	}

	// calculate the time-averaged swarm parameters when the simulation has finished
	totalIntegratedTime = time - steadyStateTime;
	getTimeAverageEnergyParams();
	getTimeAverageDistributions();
	getTimeAverageFluxParams();
	getTimeAverageBulkParams();
	getTimeAverageRateCoeffs();
	getTimeAveragePowerBalance();
	if (excitationFrequencyRadians != 0){
		getAveragedPeriodicParams();
	}

	if (collisionCounterAfterSS > maxCollisionsAfterSteadyState){
		if (dispMCStatus){
			std::cout<<" "<<std::endl;
			for (int i = 0; i < 30; ++i){
				// move up in the terminal
				std::printf("%c[1A", 0x1B);
				// clear terminal line
				std::printf("%c[2K", 0x1B);
			}
			Message::warning("Monte Carlo simulation ended after reaching ''maxCollisionsAfterSteadyState'' indicated in the setup file. [E/N = " + std::to_string(reducedElecField) + " Td, excitFreq = " + std::to_string(excitationFrequency) + " Hz, E-angle = " + std::to_string(elecFieldAngle) + " deg, B/N = " + std::to_string(reducedMagField) + " Hx]\n");
			for (int i = 0; i < 29; ++i){
				std::printf("\n");
			}
		}
		else{
			Message::warning("Monte Carlo simulation ended after reaching ''maxCollisionsAfterSteadyState'' indicated in the setup file. [E/N = " + std::to_string(reducedElecField) + " Td, excitFreq = " + std::to_string(excitationFrequency) + " Hz, E-angle = " + std::to_string(elecFieldAngle) + " deg, B/N = " + std::to_string(reducedMagField) + " Hx]\n");
		}
	}

	// calculate the final statistical errors
	checkStatisticalErrors();	

	// display the final status
	if (dispMCStatus){
		dispInfo();
	}

	energyGrid->updateMaxValue(maxEedfEnergy);
}

void BoltzmannMC::evaluateNonConstantVariables(){
	// 'evaluateNonConstantVariables' evaluates some variables that may be different each time 'evaluateEedf' is called

	// save the gas temperature
	gasTemperature = workCond->gasTemperature;
	gasEnergy = 1.5*Constant::boltzmann*gasTemperature/Constant::electronCharge; 

	// save the total gas density
	totalGasDensity = workCond->gasDensity;
 
	// evaluate the relative densities and the thermal deviations of the molecules
	int iterProcess = 0;
	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){
			// avoid effective processes
			if (collision->type == "Effective"){
				continue;
			}
			relDensities[iterProcess] = collision->target->density;
			thermalStdDeviations[iterProcess] = std::sqrt(Constant::boltzmann*gasTemperature/targetMasses[iterProcess]);
			++iterProcess;
			if (collision->isReverse){
				relDensities[iterProcess] = collision->productArray[0]->density;
				thermalStdDeviations[iterProcess] = std::sqrt(Constant::boltzmann*gasTemperature/targetMasses[iterProcess]);
				++iterProcess;
			}
		}
	}

	// initialize the times
	time = 0;
	steadyStateTime = Constant::NON_DEF;

	// update the electric and magnetic fields
	reducedElecField = workCond->reducedElecField;
	electricFieldValue = workCond->reducedElecFieldSI*totalGasDensity;
	elecFieldAngle = workCond->elecFieldAngle;
	if (elecFieldAngle == 180){
		electricField[0] = 0;
		electricField[1] = 0;
		electricField[2] = -electricFieldValue;
	}
	else{
		electricField[0] = electricFieldValue*std::sin(elecFieldAngle/180.0*M_PI);
		electricField[1] = 0;
		electricField[2] = electricFieldValue*std::cos(elecFieldAngle/180.0*M_PI);
	}

	excitationFrequency = workCond->excitationFrequency;
	excitationFrequencyRadians = excitationFrequency*2.0*M_PI;

	if (excitationFrequency != 0){
		electricField *= std::sqrt(2);
	}

	reducedMagField = workCond->reducedMagField;
	double magneticFieldValue = workCond->reducedMagFieldSI*totalGasDensity;
	cyclotronFrequency = Constant::electronCharge*magneticFieldValue/Constant::electronMass;
	isCylindricallySymmetric = workCond->isCylindricallySymmetric;

	// calculate the electron acceleration, ONLY due to the electric field
	accelerationElecField = -Constant::electronCharge/Constant::electronMass * electricField;

	// initialize to zero the position of the electrons
	electronPositions = Eigen::ArrayXXd::Zero(nElectrons, 3);

	// initialize the velocities, using a maxwell-boltzmann distribution at initialElecTempOverGasTemp*gas temperature
	electronVelocities.resize(nElectrons, 3);
	double electronThermalDeviation = std::sqrt(Constant::boltzmann*initialElecTempOverGasTemp*gasTemperature/Constant::electronMass);
	for (int i = 0; i < nElectrons; ++i){
		electronVelocities.row(i) = MathFunctions::unitNormalRand3() * electronThermalDeviation; 
	}

	// calculate the corresponding energies (in eV)
	electronEnergies = 0.5*Constant::electronMass*electronVelocities.matrix().rowwise().squaredNorm() / Constant::electronCharge;

	// initialize electron times
	electronTimes = Eigen::ArrayXd::Zero(nElectrons);
	
	// initialize collision-free times
	collisionFreeTimes = Eigen::ArrayXd::Constant(nElectrons, Constant::NON_DEF);

	// interpolate the cross sections for the first time, until two times the maximum electron energy
	// the maximum collision frequency for each energy point is also calculated
	interpolateCrossSections(2.0*electronEnergies.maxCoeff());

	// initialize the trial collision frequency
	trialCollisionFrequency = maxCollisionFrequencies[interpolCrossSectionSize-1];
	trialCollisionFrequenciesEachElectron = Eigen::ArrayXd::Constant(nElectrons, trialCollisionFrequency);

	// initialize the matrices with the mean values along time (when the size surpasses '100', we will use 'conservativeResize')
	samplingTimes = Eigen::ArrayXd::Zero(100);
	meanEnergies = Eigen::ArrayXd::Zero(100); meanPositions = Eigen::ArrayXXd::Zero(100,3); meanVelocities = Eigen::ArrayXXd::Zero(100,3); bulkVelocities = Eigen::ArrayXXd::Zero(100,3);
	fluxDiffusionCoeffs = Eigen::ArrayXXd::Zero(100,9); positionCovariances = Eigen::ArrayXXd::Zero(100,9); bulkDiffusionCoeffs = Eigen::ArrayXXd::Zero(100,9);

	// initialize the matrices with the mean values along the period (only when omega != 0)
	if (excitationFrequencyRadians != 0){
		integrationPhaseStep = 2.0*M_PI/nIntegrationPhases;
		integrationPhases = Eigen::ArrayXd::LinSpaced(nIntegrationPhases, 0.5*integrationPhaseStep, 2.0*M_PI-0.5*integrationPhaseStep);
		nIntegrationPointsPerPhase = Eigen::ArrayXd::Zero(nIntegrationPhases);
		meanEnergies_periodic = Eigen::ArrayXd::Zero(nIntegrationPhases);
		fluxVelocities_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, 3);
		bulkVelocities_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, 3);
		fluxDiffusionCoeffs_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, 9);
		bulkDiffusionCoeffs_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, 9);
	}

	// initialize the averagedMeanEnergy to Constant::NON_DEF. Important for the 'dispInfo' function
	averagedMeanEnergy = Constant::NON_DEF;
	maxElecEnergy = electronEnergies.maxCoeff();

	// initialize the time-averaged diffusion coeffs to Constant::NON_DEF
	averagedFluxDiffusionCoeffs = Eigen::Matrix3d::Constant(Constant::NON_DEF); averagedFluxDiffusionCoeffsError = Eigen::Matrix3d::Constant(Constant::NON_DEF);
	averagedBulkDiffusionCoeffs = Eigen::Matrix3d::Constant(Constant::NON_DEF); averagedBulkDiffusionCoeffsError = Eigen::Matrix3d::Constant(Constant::NON_DEF);

	// initialize to zero the data to calculate the power balance
	// initialize to zero the collision counters
	energyGainField = 0;
	energyGrowth = 0;
	totalCollisionCounter = 0; 
	nullCollisionCounter = 0;
	for (int i = 0; i < nProcesses; ++i){
		energyGainProcesses[i] = 0;
		energyLossProcesses[i] = 0;
		collisionCounters[i] = 0;
	}

	// initialize the boolean regarding the statistical errors
	goodStatisticalErrors = false;

	failedCollisionDueToGrid = 0;
}

void BoltzmannMC::interpolateCrossSections(double maxEnergy){
	// 'interpolateCrossSections' interpolates the cross sections until a given maximum energy
	// Besides, it calculates an array with the total collision frequency at each energy point
	// and an array with the maximum total collision frequency in the range [0, energy_i]

	// if the interpolation is already done until this energy, return
	if (maxEnergy == maxInterpolatedEnergy){
		return;
	}

	// calculate the energy step of the grid
	maxInterpolatedEnergy = maxEnergy;
	crossSectionEnergyStep = maxEnergy/(double)(interpolCrossSectionSize-1);
	double maxFrequency = 0;

	for (int i = 0; i < interpolCrossSectionSize; ++i){
		double energy = i*crossSectionEnergyStep;
		crossSectionEnergyGrid[i] = energy;
		double currentFrequency = 0;

		// interpolate the cross sections at this energy
		for (int iterProcess = 0; iterProcess < nProcesses; ++iterProcess){
			double value;
			if (isSuperElastic[iterProcess]){
				if (energy > energyMinLimits[iterProcess] && energy <= energyMaxLimits[iterProcess]){
					// Klein-Rosseland: sigma_{f,i}(u) = (g_i/g_f)*(1+threshold/u)*sigma_{i,f}(u+threshold)
					// note that if u=0, value = 0
					value = superElasticStatWeightFactors[iterProcess] * (1.0 + energyMinLimits[iterProcess-1]/energy)*
					        gsl_spline_eval(crossSectionInterpolations[iterProcess-1], energy+energyMinLimits[iterProcess-1], crossSectionInterpAccelerators[iterProcess-1])*
					        relDensities[iterProcess];
				}
				else{
					value = 0;
				}
			}
			else{
				if (energy >= energyMinLimits[iterProcess] && energy <= energyMaxLimits[iterProcess]){
					value = gsl_spline_eval(crossSectionInterpolations[iterProcess], energy, crossSectionInterpAccelerators[iterProcess]) *
							relDensities[iterProcess];
				}
				else{
					value = 0;
				}				
			}
			interpolCrossSectionsXrelDens[i][iterProcess] = value;
			currentFrequency += value;
			cumulSumInterpolCrossSectionsXrelDens[i][iterProcess] = currentFrequency;
		}

		currentFrequency *= totalGasDensity*std::sqrt(energy*2.0*Constant::electronCharge/Constant::electronMass);
		totalCollisionFrequencies[i] = currentFrequency;
		maxFrequency = std::fmax(currentFrequency, maxFrequency);
		maxCollisionFrequencies[i] = maxFrequency;
	}
}

void BoltzmannMC::electronDynamicsUntilSynchronization(){

	// check if the trial collision frequency used in the null-collision method is appropriate
	checkMaxCollisionFrequency();

	deltaTSynchroniz = synchronizationTimeXMaxCollisionFrequency/trialCollisionFrequency;
	nextSynchronizTime = time + deltaTSynchroniz;
	++nSynchronizationPoints;

	// initialize the array of booleans indicating if the electrons need to be advanced
	Eigen::ArrayXb electronsToBeAdvanced = Eigen::ArrayXb::Ones(nElectrons);

	// perform this cycle until all electrons are advanced to the synchronization time
	while ((electronsToBeAdvanced != false).any()){

		electronsToBeAdvanced.fill(false);

		checkMaxCollisionFrequency();

		#pragma omp parallel for
		for (int elecID = 0; elecID < (int)nElectrons; ++elecID){
			// check if the current electron has reached the synchronization time
			double electronTime = electronTimes[elecID];
			if (electronTime == nextSynchronizTime){
				chosenProcessIDs[elecID] = Constant::NON_DEF;
				continue;
			}
			// if not, advance it by performing an acceleration + collision
			Eigen::Array3d electronPosition = electronPositions.row(elecID);
			Eigen::Array3d electronVelocity = electronVelocities.row(elecID);
			double electronEnergy = electronEnergies[elecID];
			double collisionFreeTime = collisionFreeTimes[elecID];
			// if not defined, calculate a collision-free time
			if (collisionFreeTime == Constant::NON_DEF){
				// random = 0 or 1 are removed to avoid infinite or null time-steps
				collisionFreeTime = -std::log(MathFunctions::unitUniformRand(false,false)) / trialCollisionFrequency;
				// save the trial collision frequency that has been used to calculate the time. Then, this frequency will be used in the collision choice
				trialCollisionFrequenciesEachElectron[elecID] = trialCollisionFrequency;
			}
			// if the collision-free time is enough to surpass the synchronization time, do not perform a collision and advance only until the synch time
			if (electronTime + collisionFreeTime > nextSynchronizTime){
				double deltaTAccell = nextSynchronizTime-electronTime;
				accelerateElectron(elecID, electronTime, deltaTAccell, electronPosition, electronVelocity, electronEnergy);
				electronTime = nextSynchronizTime;
				collisionFreeTime -= deltaTAccell;
				chosenProcessIDs[elecID] = GeneralDefinitions::partialFreeFlightID;
			}
			// accelerate, perform collision and calculate next collision-free time
			else{
				accelerateElectron(elecID, electronTime, collisionFreeTime, electronPosition, electronVelocity, electronEnergy);
				electronTime += collisionFreeTime;
				performCollision(elecID, electronPosition, electronVelocity, electronEnergy);
				// random = 0 or 1 are removed to avoid infinite or null time-steps
				collisionFreeTime = -std::log(MathFunctions::unitUniformRand(false,false)) / trialCollisionFrequency;
				// save the trial collision frequency that has been used to calculate the time. Then, this frequency will be used in the collision choice
				trialCollisionFrequenciesEachElectron[elecID] = trialCollisionFrequency;
				// since the electron has not reached the synch time, activate boolean
				electronsToBeAdvanced[elecID] = true;					
			}
			electronPositions.row(elecID) = electronPosition;
			electronVelocities.row(elecID) = electronVelocity;
			electronEnergies[elecID] = electronEnergy;
			electronTimes[elecID] = electronTime;
			collisionFreeTimes[elecID] = collisionFreeTime;
		}
		// perform all operations that cannot be performed in parallel
		nonParallelCollisionTasks();
	}

	// update the general time
	time = nextSynchronizTime;
}

void BoltzmannMC::getMeanCollisionFrequency(){
	// 'getMeanCollisionFrequency' calculates the mean collision frequency 

	meanCollisionFrequency = 0;
	if (steadyStateTime == Constant::NON_DEF){
		for (int elecID = 0; elecID < nElectrons; ++elecID){
			int energyIndex = std::fmin(electronEnergies[elecID]/crossSectionEnergyStep, interpolCrossSectionSize-1);
			meanCollisionFrequency += totalCollisionFrequencies[energyIndex];
		}
		meanCollisionFrequency /= nElectrons;
	}
	else{
		double normalizer = eehSum.sum();
		for (int i = 0; i < nEnergyCells; ++i){
			int energyIndex = std::fmin(eedfEnergyCells[i]/crossSectionEnergyStep, interpolCrossSectionSize-1);
			if (energyIndex != interpolCrossSectionSize-1){
				meanCollisionFrequency += eehSum[i]*0.5*(totalCollisionFrequencies[energyIndex]+totalCollisionFrequencies[energyIndex+1]);
			}
			else{
				meanCollisionFrequency += eehSum[i]*totalCollisionFrequencies[energyIndex];
			}
		}
		meanCollisionFrequency /= normalizer;
	}
}

void BoltzmannMC::checkMaxCollisionFrequency(){
	// 'checkMaxCollisionFrequency' assures that the trial collision frequency is higher than the maximum collision frequency possible
	// Additionally, it checks if the cross-section energy is appropriate for the current energy

	// maximum energy before acceleration
	double maxEnergyBeforeAccel = electronEnergies.maxCoeff();

	// estimate a maximal contribution for the thermal energy of the molecules in case their motion is being considered
	double thermalMolecContribution = 0;
	if (gasTemperatureEffect == GeneralDefinitions::trueGasTempEffectID || gasTemperatureEffect == GeneralDefinitions::smartActivationGasTempEffectID){
		// note: probability to find a molecule with energy higher than 10.0*1.5*k_B Tg is 1.4E-6
		thermalMolecContribution = 10.0*gasEnergy;
	}

	bool updatedTrial = true;
	while (updatedTrial){
		updatedTrial = false;

		// calculate the maximum energy possible after 10.0/trialCollisionFrequency
		// note: probality of having a collision-free time higher than 10.0/trialCollisionFrequency is ~ 4.5E-5
		// therefore, the time-interval used for the maximization is quite adequate
		double maxEnergy = maximizationAccelerationEnergy(maxEnergyBeforeAccel, 10.0/trialCollisionFrequency) + thermalMolecContribution;

		// update the interpolated cross sections, if the maximum interpolated energy is not high enough or if it is too high, 
		// where the last case would cause physical imprecisions due to poor energy discretization
		if (maxEnergy > maxInterpolatedEnergy || 2.5*maxEnergy < maxInterpolatedEnergy){
			// this condition needs to be here, since in the first iterations the trial collision frequency is very low, 
			// leading to very high energies, even higher than the maximum energy defined in the cross sections
			if (2.0*maxEnergy < energyMaxElastic){
				interpolateCrossSections(2.0*maxEnergy);
			}
			else{
				interpolateCrossSections(energyMaxElastic);
			}
		}

		// calculate the maximum collision frequency in the energy range [0, maxEnergy]
		// fmin is used to prevent segmentation fault when maxEnergy = maxInterpolatedEnergy
		int energyIndex = std::fmin(std::ceil(maxEnergy/crossSectionEnergyStep), interpolCrossSectionSize-1);
		double maxCollisionFrequency = maxCollisionFrequencies[energyIndex];	

		// if the trial collision frequency, used in the null-collision method, is smaller than the max coll freq, increase it
		if (trialCollisionFrequency < maxCollisionFrequency){
			updatedTrial = true;
			trialCollisionFrequency *= 1.1;
		}
	}
}

double BoltzmannMC::maximizationAccelerationEnergy(double initialEnergy, double deltaT){
	double e = Constant::electronCharge, me = Constant::electronMass, e_me = e/me;
	double Ex0 = std::abs(electricField[0]), Ez0 = std::abs(electricField[2]);
	double Ex02 = Ex0*Ex0, Ez02 = Ez0*Ez0, E02 = Ex02+Ez02, E0 = std::sqrt(E02);
	double v0 = std::sqrt(initialEnergy*e*2.0/me);
	double w = excitationFrequencyRadians, W = cyclotronFrequency;

	// start by doing a maximization assuming a constant electric field and null mag field (w = W = 0)
	// this maximization works for any configuration but it may overestimate a lot in cases with w != 0 and/or W != 0
	double energyGain = (E0*v0 + 0.5*e_me*E02*deltaT) *deltaT;

	// for other configurations, find the minimum between the previous maximization and a specific maximization

	// null magnetic field
	if (W == 0){
		if (w != 0){
			energyGain = std::fmin(energyGain, 2.0/w*(e_me*E02/w + v0*(Ex0+Ez0)) );
		}
	}
	// non-null DC magnetic field
	else{
		if (w == 0){
			energyGain = std::fmin(energyGain, 0.5*e_me*Ez02*deltaT*deltaT + (2.0*e_me*Ex02/W + 3.0*v0*Ex0)/W + v0*deltaT*Ez0 );
		}
		else if (std::abs(w-W)/w < 1E-6){
			double W2 = W*W;
			energyGain = std::fmin(energyGain, 2.0*e_me*Ez02/W2 + e_me*Ex02/(8.0*W2)*(4.0+W*deltaT*(2.0+W*deltaT)) + 
									(v0*deltaT+v0/W)*Ex0 + 2.0*v0*Ez0/W);
		}
		else{
			double w2 = w*w, W2 = W*W, w2mW2 = w2-W2;
			energyGain = std::fmin(energyGain, 2.0*e_me*Ez02/w2 + 0.5*e_me*Ex02/(w2mW2*w2mW2)*(5.0*w2+8.0*w*W+5.0*W*W) + 
									3.0*v0*Ex0/std::abs(w-W) + 2.0*v0*Ez0/w);
		}
	}

	return initialEnergy + energyGain;
}

void BoltzmannMC::accelerateElectron(int elecID, double electronTime, double deltaT, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy){
	
	double previousEnergy = electronEnergy;

	// accelerate the electron depending on the field configuration
	double e = Constant::electronCharge, me = Constant::electronMass;

	if (reducedMagField == 0){
		if (excitationFrequency == 0){
			electronPosition +=  electronVelocity*deltaT + (accelerationElecField*(0.5*deltaT*deltaT));
			electronVelocity += (accelerationElecField*deltaT);
		}
		else{
			double w = excitationFrequencyRadians;
			double phi = w*electronTime;
			double sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			double wDeltaT = w*deltaT;
			double phase = wDeltaT + phi;
			double sinPhase = std::sin(phase), cosPhase = std::cos(phase);
			double e_me_w = e/(me*w);
			double aux1 = e_me_w/w*(cosPhase+wDeltaT*sinPhi-cosPhi);
			double aux2 = e_me_w*(sinPhi-sinPhase);
			electronPosition[0] += electronVelocity[0]*deltaT + electricField[0]*aux1;
			electronPosition[1] += electronVelocity[1]*deltaT;
			electronPosition[2] += electronVelocity[2]*deltaT + electricField[2]*aux1;
			electronVelocity[0] += electricField[0]*aux2;
			electronVelocity[2] += electricField[2]*aux2;
		}
	}
	else{
		double vx0 = electronVelocity[0], vy0 = electronVelocity[1], vz0 = electronVelocity[2];
		if (excitationFrequency == 0){
			double W = cyclotronFrequency;
			double WDeltaT = W*deltaT;
			double sinWDeltaT = std::sin(WDeltaT), cosWDeltaT = std::cos(WDeltaT);
			double sinWDeltaT_W = sinWDeltaT/W;
			double aux2 = (cosWDeltaT-1.0)/W;
			double vExCyc = e*electricField[0]/(me*W);
			double az = e*electricField[2]/me;
			electronPosition[0] += vx0*sinWDeltaT_W + (vy0+vExCyc)*aux2;
			electronPosition[1] += -vx0*aux2 + vy0*sinWDeltaT_W + vExCyc*(sinWDeltaT_W-deltaT);
			electronPosition[2] += vz0*deltaT - 0.5*az*deltaT*deltaT;
			electronVelocity[0] = vx0*cosWDeltaT - vy0*sinWDeltaT - vExCyc*sinWDeltaT;
			electronVelocity[1] = vx0*sinWDeltaT + vy0*cosWDeltaT + vExCyc*(cosWDeltaT-1.0);
			electronVelocity[2] = vz0 - az*deltaT;
		}
		else if (std::abs(excitationFrequencyRadians-cyclotronFrequency)/excitationFrequencyRadians < 1E-6){
			// electron cyclotron resonance
			double W = cyclotronFrequency;
			double phi = W*electronTime;
			double sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			double WDeltaT = W*deltaT;
			double sinWDeltaT = std::sin(WDeltaT), cosWDeltaT = std::cos(WDeltaT);
			double phase = WDeltaT + phi;
			double sinPhase = std::sin(phase), cosPhase = std::cos(phase);
			double cosOppositePhase =  std::cos(phi - WDeltaT);
			double e_me_W = e/(me*W);
			double vExCyc = e_me_W*electricField[0];
			double vEzAC = e_me_W*electricField[2];
			double sinWDeltaT_W = sinWDeltaT/W;
			double cosWDeltaTm1_W = (cosWDeltaT-1.0)/W;
			double cosPhase_W = cosPhase/W;
			electronPosition[0] += vx0*sinWDeltaT_W + vy0*cosWDeltaTm1_W - 0.25*vExCyc*(cosPhase_W-cosOppositePhase/W+2.0*deltaT*sinPhase);
			electronPosition[1] += -vx0*cosWDeltaTm1_W + vy0*sinWDeltaT_W + 0.5*vExCyc*(deltaT*cosPhase_W-2.0*cosWDeltaTm1_W*sinPhi-cosPhi*sinWDeltaT_W);
			electronPosition[2] += vz0*deltaT + vEzAC*(cosPhase_W-cosPhi/W+deltaT*sinPhi);
			electronVelocity[0] = vx0*cosWDeltaT - vy0*sinWDeltaT - 0.5*vExCyc*(WDeltaT*cosPhase+cosPhi*sinWDeltaT);
			electronVelocity[1] = vx0*sinWDeltaT + vy0*cosWDeltaT - 0.5*vExCyc*(WDeltaT*cosWDeltaT*sinPhi+(WDeltaT*cosPhi-sinPhi)*sinWDeltaT);
			electronVelocity[2] = vz0 + vEzAC*(sinPhi-sinPhase);
		}
		else{
			double w = excitationFrequencyRadians, W = cyclotronFrequency;
			double phi = w*electronTime;
			double sinPhi = std::sin(phi), cosPhi = std::cos(phi);
			double wDeltaT = w*deltaT;
			double sinwDeltaT = std::sin(wDeltaT), coswDeltaT = std::cos(wDeltaT);
			double phase = wDeltaT + phi;
			double sinPhase = std::sin(phase), cosPhase = std::cos(phase);
			double WDeltaT = W*deltaT;
			double sinWDeltaT = std::sin(WDeltaT), cosWDeltaT = std::cos(WDeltaT);
			double w2 = w*w, W2 = W*W, w2mW2 = w2-W2;
			double vExCyc = e/(me*W)*electricField[0];
			double vEzAC = e/(me*w)*electricField[2];
			double vExCyc_w2mW2 = vExCyc/w2mW2;
			double WVExCyc_w2mW2 = W*vExCyc_w2mW2;
			double sinWDeltaT_W = sinWDeltaT/W;
			double cosWDeltaTm1_W = (cosWDeltaT-1.0)/W;
			double wSinPhi = w*sinPhi;	
			electronPosition[0] += vx0*sinWDeltaT_W + vy0*cosWDeltaTm1_W + WVExCyc_w2mW2*(cosPhi*(coswDeltaT-cosWDeltaT)+(wSinPhi*sinWDeltaT_W-sinPhi*sinwDeltaT));
			electronPosition[1] += -vx0*cosWDeltaTm1_W + vy0*sinWDeltaT_W - vExCyc_w2mW2/w*((W2+w2*cosWDeltaT-w2)*sinPhi+W2*(w*cosPhi*sinWDeltaT_W-sinPhase));
			electronPosition[2] += vz0*deltaT + vEzAC*((cosPhase-cosPhi)/w+deltaT*sinPhi);
			electronVelocity[0] = vx0*cosWDeltaT - vy0*sinWDeltaT + WVExCyc_w2mW2*w*(cosWDeltaT*sinPhi-sinPhase+cosPhi/w*W*sinWDeltaT);
			electronVelocity[1] = vx0*sinWDeltaT + vy0*cosWDeltaT + vExCyc_w2mW2*W2*(cosPhase-cosPhi*cosWDeltaT+wSinPhi*sinWDeltaT_W);
			electronVelocity[2] = vz0 + vEzAC*(sinPhi-sinPhase);
		}
	}

	// update the energy of the electrons (in eV)
	electronEnergy = 0.5*Constant::electronMass*electronVelocity.matrix().squaredNorm() / Constant::electronCharge;
	
	// add the energy gain during the acceleration (in eV)
	energyGainsField[elecID] = (electronEnergy - previousEnergy);	
}

void BoltzmannMC::performCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy){

	Eigen::Array3d targetVelocity = Eigen::Array3d::Zero();

	// ----- select the process ----- // 

	int chosenProcessID = GeneralDefinitions::nullCollisionID;

	// if the gas-temperature effect is to be considered
	if (gasTemperatureEffect == GeneralDefinitions::trueGasTempEffectID || (gasTemperatureEffect == GeneralDefinitions::smartActivationGasTempEffectID && electronEnergy < 20.0*gasEnergy)){
		// calculate the random versor of the target velocity
		Eigen::Array3d targetVelocityVersor = MathFunctions::unitNormalRand3();		
		double randCrossSectionFactor = trialCollisionFrequenciesEachElectron[elecID] * MathFunctions::unitUniformRand(false,true) / totalGasDensity;
		// find the process using the bissection method gas by gas
		double previousSumCrossSectionFactor = 0;

		/*// collision choice for debugging
		double sumCrossSectionFactor = 0;
		int chosenProcessID_debug = GeneralDefinitions::nullCollisionID;
		for (int iterProcess = 0; iterProcess < nProcesses; ++iterProcess){
			// calculate the target velocity for the current gas
			targetVelocity = targetVelocityVersor*thermalStdDeviations[iterProcess];
			// calculate the relative speed and energy index
			double relativeSpeed = (electronVelocity-targetVelocity).matrix().norm();
			double relEnergyOverStep = 0.5*reducedMasses[iterProcess]*relativeSpeed*relativeSpeed/Constant::electronCharge/crossSectionEnergyStep;
			int relEnergyIndex = std::fmin(relEnergyOverStep, interpolCrossSectionSize-1);
			sumCrossSectionFactor += interpolCrossSectionsXrelDens[relEnergyIndex][iterProcess] * relativeSpeed;
			// check if the cumulative factor surpassed the random factor, i.e., if this process is to be chosen
			if (sumCrossSectionFactor >= randCrossSectionFactor){
				chosenProcessID_debug = iterProcess;
				break;
			}			
		}*/

		// choose the collision type using a bissection method. However, it needs to be separated by gas, since each gas will have a different relative energy
		for (int iGas = 0; iGas < nGases && chosenProcessID == GeneralDefinitions::nullCollisionID; ++iGas){
			if (gasFractions[iGas] == 0){
				continue;
			}
			// define the limits for the bissection search
			int leftLimit = firstProcessIndexPerGas[iGas];
			int rightLimit = lastProcessIndexPerGas[iGas];
			// calculate the target velocity for the current gas
			targetVelocity = targetVelocityVersor*thermalStdDeviations[leftLimit];
			// calculate the relative speed and energy index
			double relativeSpeed = (electronVelocity-targetVelocity).matrix().norm();
			double relEnergyOverStep = 0.5*reducedMasses[leftLimit]*relativeSpeed*relativeSpeed/Constant::electronCharge/crossSectionEnergyStep;
			// index of the first energy row
			int relEnergyIndex1 = std::fmin(relEnergyOverStep, interpolCrossSectionSize-1); // conversion to int rounds the number downwards!
			// index of the second energy row
			int relEnergyIndex2 = std::fmin(relEnergyIndex1+1, interpolCrossSectionSize-1);
			// assign weight fractions to each row, based on linear interpolation
			double weight1 = (double)relEnergyIndex2 - relEnergyOverStep;
			// in case both indexes correspond to the limit
			if (weight1 < 0){
				weight1 = 0.0;
			}
			else{
				weight1 = 1.0;
			}
			double weight2 = 1.0-weight1;			
			// save the pointers to the two energy rows
			double* cumulSumInterpolCrossSectionsXrelDens_rowPtr1 = cumulSumInterpolCrossSectionsXrelDens[relEnergyIndex1];
			double* cumulSumInterpolCrossSectionsXrelDens_rowPtr2 = cumulSumInterpolCrossSectionsXrelDens[relEnergyIndex2];
			double* interpolCrossSectionsXrelDens_rowPtr1 = interpolCrossSectionsXrelDens[relEnergyIndex1];
			double* interpolCrossSectionsXrelDens_rowPtr2 = interpolCrossSectionsXrelDens[relEnergyIndex2];
			// consider as reference value the cumulative sum immediately before this gas
			double referenceValue = 0;
			if (leftLimit > 0){
				referenceValue = weight1*cumulSumInterpolCrossSectionsXrelDens_rowPtr1[leftLimit-1]+weight2*cumulSumInterpolCrossSectionsXrelDens_rowPtr2[leftLimit-1];
			}
			// if the randCrossSectionFactor is higher than the sumed cross-section factor at the right limit, continue to the next gas
			double limitCrossSectionFactor = previousSumCrossSectionFactor + (weight1*cumulSumInterpolCrossSectionsXrelDens_rowPtr1[rightLimit] +
				weight2*cumulSumInterpolCrossSectionsXrelDens_rowPtr2[rightLimit] - referenceValue)*relativeSpeed;
			if (randCrossSectionFactor > limitCrossSectionFactor){
				// save the previousSumCrossSectionFactor, to be used in the next iteration
				previousSumCrossSectionFactor = limitCrossSectionFactor;
				continue;
			}

			// perform bissection method (taking into account what was already summed from the other gases)
			while(leftLimit != rightLimit){

				// get the tentative index. Note that in this operation, the result is rounded downwards (see arithmetic rules of integers in c++)
				int tentativeIndex = (leftLimit + rightLimit)/2;
				double tentativeValue = previousSumCrossSectionFactor + (weight1*cumulSumInterpolCrossSectionsXrelDens_rowPtr1[tentativeIndex] + 
					weight2*cumulSumInterpolCrossSectionsXrelDens_rowPtr2[tentativeIndex]-referenceValue)*relativeSpeed;
				// if the cumulative value of the tentative index is greater than the random value, the solution is in the left side
				if (randCrossSectionFactor < tentativeValue){
					rightLimit = tentativeIndex;
				}
				// if the cumulative value of the tentative index is smaller than the random value, the solution is in the right side
				else if (randCrossSectionFactor > tentativeValue){
					leftLimit = tentativeIndex + 1;
				}
				else{
					chosenProcessID = tentativeIndex;
					break;
				}
			}

			if (leftLimit == rightLimit){
				chosenProcessID = leftLimit;
			}

			// avoid collisions with null rates
			while (weight1*interpolCrossSectionsXrelDens_rowPtr1[chosenProcessID] == 0 && weight2*interpolCrossSectionsXrelDens_rowPtr2[chosenProcessID] == 0){
				--chosenProcessID;
			}
		}

		// save the chosen process ID
		chosenProcessIDs[elecID] = chosenProcessID;

		/*if (chosenProcessID != chosenProcessID_debug){
			std::cout<<"chosenProcessID_debug: "<<chosenProcessID_debug<<"\n"<<realCollisionPointers[chosenProcessID_debug]->description()<<"\n";
			std::cout<<"chosenProcessID: "<<chosenProcessID<<"\n"<<realCollisionPointers[chosenProcessID]->description()<<"\n";

		}*/						

		// check if this is a null collision
		if (chosenProcessID == GeneralDefinitions::nullCollisionID){
			return;
		}
	}

	// if the gas-temperature effect is to be neglected
	else{
		double randCollisionFrequency = trialCollisionFrequenciesEachElectron[elecID] * MathFunctions::unitUniformRand(false,true);
		double energyOverStep = electronEnergy/crossSectionEnergyStep;
		// index of the first energy row
		int incEnergyIndex1 = std::fmin(energyOverStep, interpolCrossSectionSize-1); // conversion to int rounds the number downwards!
		// index of the second energy row
		int incEnergyIndex2 = std::fmin(incEnergyIndex1+1, interpolCrossSectionSize-1);
		// assign weight fractions to each row, based on linear interpolation
		double weight1 = (double)incEnergyIndex2 - energyOverStep;
		// in case both indexes correspond to the limit
		if (weight1 < 0){
			weight1 = 0.0;
		}				
		double weight2 = 1.0-weight1;

		// check if this is a null collision 
		if (randCollisionFrequency > weight1*totalCollisionFrequencies[incEnergyIndex1]+weight2*totalCollisionFrequencies[incEnergyIndex2]){
			chosenProcessIDs[elecID] = GeneralDefinitions::nullCollisionID;
			return;				
		}
		// save the pointers to the two energy rows
		double* cumulSumInterpolCrossSectionsXrelDens_rowPtr1 = cumulSumInterpolCrossSectionsXrelDens[incEnergyIndex1];
		double* cumulSumInterpolCrossSectionsXrelDens_rowPtr2 = cumulSumInterpolCrossSectionsXrelDens[incEnergyIndex2];
		double* interpolCrossSectionsXrelDens_rowPtr1 = interpolCrossSectionsXrelDens[incEnergyIndex1];
		double* interpolCrossSectionsXrelDens_rowPtr2 = interpolCrossSectionsXrelDens[incEnergyIndex2];

		// use a bissection method to find the process ID. Only possible when the gas-temperature effect is not considered
		double randCrossSectionFactor = randCollisionFrequency / totalGasDensity / electronVelocity.matrix().norm();
		int leftLimit = 0;
		int rightLimit = nProcesses-1;

		int iterations = 0;
		while(leftLimit != rightLimit){

			// get the tentative index. Note that in this operation, the result is rounded downwards (see arithmetic rules of integers in c++)
			int tentativeIndex = (leftLimit + rightLimit)/2;

			// if the cumulative value of the tentative index is greater than the random value, the solution is in the left side
			double tentativeValue = weight1*cumulSumInterpolCrossSectionsXrelDens_rowPtr1[tentativeIndex]+weight2*cumulSumInterpolCrossSectionsXrelDens_rowPtr2[tentativeIndex];
			if (randCrossSectionFactor < tentativeValue){
				rightLimit = tentativeIndex;
			}
			// if the cumulative value of the tentative index is smaller than the random value, the solution is in the right side
			else if (randCrossSectionFactor > tentativeValue){
				leftLimit = tentativeIndex + 1;
			}
			else{
				chosenProcessID = tentativeIndex;
				break;
			}
		}

		if (leftLimit == rightLimit){
			chosenProcessID = leftLimit;
		}	

		// avoid collisions with null rates
		while (weight1*interpolCrossSectionsXrelDens_rowPtr1[chosenProcessID] == 0 && weight2*interpolCrossSectionsXrelDens_rowPtr2[chosenProcessID] == 0){
			--chosenProcessID;
		}

		// save the chosen process ID
		chosenProcessIDs[elecID] = chosenProcessID;		
	}

	// ----- calculate the velocity after the collision ----- //

	if (processTypes[chosenProcessID] == GeneralDefinitions::conservativeType){
		conservativeCollision(elecID, electronPosition, electronVelocity, targetVelocity, electronEnergy, chosenProcessID);
	}

	else if (processTypes[chosenProcessID] == GeneralDefinitions::ionizationType){
		ionizationCollision(elecID, electronPosition, electronVelocity, electronEnergy, chosenProcessID);
	}

	else if (processTypes[chosenProcessID] == GeneralDefinitions::attachmentType){
		attachmentCollision(elecID, electronEnergy);
	}

}

void BoltzmannMC::conservativeCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, Eigen::Array3d &targetVelocity, double &electronEnergy, int chosenProcessID){
	// 'convervativeCollision' solves the dynamics of a conservative collision, i.e., a collision without changes in the electron number
	// Important note: for gasTemperatureEffect = trueGasTempEffectID, it was verified that the energy of the system [electron + target] is conserved!
	// Performance can be improved by not calculating angles, but only sinAngle and cosAngle. This avoids multiple uses of sinus and cossinus

	double targetMass = targetMasses[chosenProcessID];
	double reducedMass = reducedMasses[chosenProcessID];
	double energyLoss = energyLosses[chosenProcessID];
	double incidentEnergy = electronEnergy;

	if (gasTemperatureEffect == GeneralDefinitions::trueGasTempEffectID || (gasTemperatureEffect == GeneralDefinitions::smartActivationGasTempEffectID && incidentEnergy < 20.0*gasEnergy)){
		// calculate the relative velocity before the collision and convert it to spherical coordinates
		Eigen::Array3d relVelocityBefore = electronVelocity - targetVelocity;
		double relSpeedBefore, sinThetaRelVel, cosThetaRelVel, sinPhiRelVel, cosPhiRelVel;
		MathFunctions::cart2sph(relVelocityBefore, relSpeedBefore, sinThetaRelVel, cosThetaRelVel, sinPhiRelVel, cosPhiRelVel);
		double relativeEnergy = 0.5*reducedMass*relSpeedBefore*relSpeedBefore/Constant::electronCharge;
		double relEnergyAfter = relativeEnergy-energyLoss;

		// this occurs very rarely (typically less than once per 1E9 [real + null] collisions)
		// and is caused by the linear interpolation between two energy nodes made in the "performCollision"
		// When the threshold of the collision is in-between the energy nodes, very very rarely an electron with relativeEnergy<energyLoss may be chosen
		// This has ZERO influence in the results
		if (relEnergyAfter <= 0){
			chosenProcessIDs[elecID] = GeneralDefinitions::nullCollisionID;
			return;			
		}		

		// generate the cos(scattering angle,) using the functions defined in 'AngularScatteringFunctions.h'
		double cosChiScattered = angularScatteringFunctions[chosenProcessID](relativeEnergy, relEnergyAfter, chosenProcessID, this);
		double sinChiScattered = std::sqrt(1.0-cosChiScattered*cosChiScattered);

		// generate azimutal angle in ]0,2pi]
		double etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		double sinEtaScattered = std::sin(etaScattered), cosEtaScattered = std::cos(etaScattered);

		// calculate the relative speed after the collision, using energy conservation (note that the energy must be converted from eV to SI)
		// use Euler relations to calculate the relative velocity after the collision in the laboratory, knowing chi and eta in the CM frame
		Eigen::Array3d relVelocityAfter = std::sqrt(relSpeedBefore*relSpeedBefore-2.0/reducedMass*energyLoss*Constant::electronCharge)* 
								   MathFunctions::eulerTransformation(sinChiScattered, cosChiScattered, sinEtaScattered, cosEtaScattered, sinThetaRelVel, cosThetaRelVel, sinPhiRelVel, cosPhiRelVel);

		// get the electron velocity after the collision using the conservation of momentum transfer
		electronVelocity = targetMass/(Constant::electronMass+targetMass)*relVelocityAfter +
						   (Constant::electronMass*electronVelocity + targetMass*targetVelocity)/(Constant::electronMass+targetMass);
	}
	else{
		// convert the velocity of the incident electron to spherical coordinates
		double electronSpeed, sinTheta, cosTheta, sinPhi, cosPhi;
		MathFunctions::cart2sph(electronVelocity, electronSpeed, sinTheta, cosTheta, sinPhi, cosPhi);
		double energyAfter = incidentEnergy-energyLoss;

		// this occurs very rarely (typically less than once per 1E9 [real + null] collisions)
		// and is caused by the linear interpolation between two energy nodes made in the "performCollision"
		// When the threshold of the collision is in-between the energy nodes, very very rarely an electron with electronEnergy<energyLoss may be chosen
		// This has ZERO influence in the results
		if (energyAfter <= 0){
			chosenProcessIDs[elecID] = GeneralDefinitions::nullCollisionID;		
			return;			
		}

		// generate the scattering angle, using the functions defined in 'AngularScatteringFunctions.h'
		double cosChiScattered = angularScatteringFunctions[chosenProcessID](incidentEnergy, energyAfter, chosenProcessID, this);
		double sinChiScattered = std::sqrt(1.0-cosChiScattered*cosChiScattered);

		// generate azimutal angle in ]0,2pi]
		double etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		double sinEtaScattered = std::sin(etaScattered), cosEtaScattered = std::cos(etaScattered);		

		// calculate the electron speed after the collision, using energy conservation and assuming a target molecule at rest (due to the high mass) [Reid 1979]
		// use Euler relations to determine the velocity after the collision
		electronVelocity = std::sqrt((electronSpeed*electronSpeed-2.0/Constant::electronMass*energyLoss*Constant::electronCharge) * 
							    (1.0-2.0*reducedMass/(Constant::electronMass+targetMass)*(1.0-cosChiScattered)))*
						   MathFunctions::eulerTransformation(sinChiScattered, cosChiScattered, sinEtaScattered, cosEtaScattered, sinTheta, cosTheta, sinPhi, cosPhi);
	}

	// calculate the energy change (in eV) due to this collision
	electronEnergy = 0.5*Constant::electronMass*electronVelocity.matrix().squaredNorm() / Constant::electronCharge;
	electronEnergyChanges[elecID] = electronEnergy - incidentEnergy;
	electronEnergyChangesOverIncidEnergies[elecID] = electronEnergyChanges[elecID]/incidentEnergy;
}

void BoltzmannMC::ionizationCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy, int chosenProcessID){
	// 'ionizationCollision' solves the dynamics of an ionization collision

	double incidentEnergy = electronEnergy;
	// check if the incident energy is smaller than the threshold. This happens only when the gas-temperature effect is considered
	// This occurs much less than once per simulation and does not affect the results
	double ionizationThreshold = energyLosses[chosenProcessID];
	if (incidentEnergy < ionizationThreshold){
		chosenProcessIDs[elecID] = GeneralDefinitions::nullCollisionID;
		return;
	}

	// convert the velocity of the incident electron to spherical coordinates
	double electronSpeed, sinTheta, cosTheta, sinPhi, cosPhi;
	MathFunctions::cart2sph(electronVelocity, electronSpeed, sinTheta, cosTheta, sinPhi, cosPhi);

	// calculate the speeds of the scattered and ejected electrons, assuming a background of molecules at rest
	double netEnergy = incidentEnergy - ionizationThreshold;
	double electronEnergyEjected;
	if (energySharingIonizType == GeneralDefinitions::usingSDCSIonizType){
		electronEnergyEjected = wParameters[chosenProcessID]*std::tan(MathFunctions::unitUniformRand(true,true) * std::atan(netEnergy/(2.0*wParameters[chosenProcessID])) );
	}
	else if (energySharingIonizType == GeneralDefinitions::randomUniformIonizType){
		electronEnergyEjected = MathFunctions::unitUniformRand(true,true)*netEnergy;
	}
	else{
		electronEnergyEjected = energySharingFactor*netEnergy;
	}
	electronEnergy = netEnergy - electronEnergyEjected;

	// generate scattering angles depending on the model considered
	double sinChiScattered, cosChiScattered, sinChiEjected, cosChiEjected, sinEtaScattered, cosEtaScattered, sinEtaEjected, cosEtaEjected;

	if (isMomentumConservationIonizationScattering[chosenProcessID]){
		// momentum-conservation scattering	
		// generate the scattering angles, assuming that the incident, scattered and ejected electron velocities are coplanar
		// and that the scattered and ejected electron velocities are perpendicular in the CM frame (according with Boeuf 1982)				
		cosChiScattered = std::sqrt(electronEnergy/netEnergy);
		sinChiScattered = std::sqrt(1.0-cosChiScattered*cosChiScattered);
		double etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		sinEtaScattered = std::sin(etaScattered);
		cosEtaScattered = std::cos(etaScattered);
		cosChiEjected = std::sqrt(electronEnergyEjected/netEnergy);
		sinChiEjected = std::sqrt(1.0-cosChiEjected*cosChiEjected);
		sinEtaEjected = -sinEtaScattered;
		cosEtaEjected = -cosEtaScattered;	
	}
	else{
		cosChiScattered = angularScatteringFunctions[chosenProcessID](incidentEnergy, electronEnergy, chosenProcessID, this);
		sinChiScattered = std::sqrt(1.0-cosChiScattered*cosChiScattered);
		double etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		sinEtaScattered = std::sin(etaScattered);
		cosEtaScattered = std::cos(etaScattered);
		cosChiEjected = angularScatteringFunctions[chosenProcessID](incidentEnergy, electronEnergyEjected, chosenProcessID, this);
		sinChiEjected = std::sqrt(1.0-cosChiEjected*cosChiEjected);
		double etaEjected = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		sinEtaEjected = std::sin(etaEjected);
		cosEtaEjected = std::cos(etaEjected);		
	}	

	// use Euler relations to calculate the electron velocity after the collision in the laboratory, knowing chi_scattered and eta in the CM frame (according with Yousfi 1994)
	electronVelocity = std::sqrt(2.0*electronEnergy*Constant::electronCharge/Constant::electronMass) * 
					   MathFunctions::eulerTransformation(sinChiScattered, cosChiScattered, sinEtaScattered, cosEtaScattered, sinTheta, cosTheta, sinPhi, cosPhi);

	// use Euler relations to calculate the ejected electron velocity in the laboratory, knowing chi_ejected and eta in the CM frame
	ejectedElectronVelocities.row(elecID) = std::sqrt(2.0*electronEnergyEjected*Constant::electronCharge/Constant::electronMass) *
											MathFunctions::eulerTransformation(sinChiEjected, cosChiEjected, sinEtaEjected, cosEtaEjected, sinTheta, cosTheta, sinPhi, cosPhi);

	// save the ejected energy to avoid recalculation
	ejectedElectronEnergies[elecID] = electronEnergyEjected;

	// create the electron in the same position as the scattered one
	ejectedElectronPositions.row(elecID) = electronPosition;

	// update the energy change (in eV) due to this collision
	electronEnergyChanges[elecID] = -ionizationThreshold;
	electronEnergyChangesOverIncidEnergies[elecID] = -ionizationThreshold/incidentEnergy;
}

void BoltzmannMC::attachmentCollision(int elecID, double incidentEnergy){
	// 'attachmentCollision' solves the dynamics of an attachment collision

	// update the energy change (in eV) due to this collision
	electronEnergyChanges[elecID] = -incidentEnergy;
	electronEnergyChangesOverIncidEnergies[elecID] = -1;
}

void BoltzmannMC::nonParallelCollisionTasks(){
	// 'nonParallelCollisionTasks' performs all the tasks involving collisions that are not parallelized:
	// increment collision counters, update the energy gain/loss in processes, create ejected electrons, destroy attached electrons
	// IMPORTANT: now we conserve correctly the distribution of electrons. Attached electrons no longer participate in the duplication of electrons
	// and the ejected electrons are treated in the same manner as the electrons that are in the "normal" ensemble. Previously, we were keeping always
	// the ejected electrons which is NOT CORRECT. 

	std::vector<bool> isAttached(nElectrons,false);
	std::vector<int> ejectElectronIDs;
	std::vector<int> attachElectronIDs;

	for (int elecID = 0; elecID < nElectrons; ++elecID){

		int chosenProcessID = chosenProcessIDs[elecID];

		// if the electron was not advanced since it was already at the synchronization time
		if (chosenProcessID == Constant::NON_DEF){
			continue;
		}

		// update the energy gain due to the field
		energyGainField += energyGainsField[elecID];

		// update the collision counters

		// null collision
		if (chosenProcessID == GeneralDefinitions::nullCollisionID){
			++nullCollisionCounter;
			continue;
		}
		// did not reach the end of the free flight
		else if (chosenProcessID == GeneralDefinitions::partialFreeFlightID){
			continue;
		}
		// collision happened
		else{
			++collisionCounters[chosenProcessID];
			++totalCollisionCounter;
		}

		// update the energy gain due to collisions
		if (electronEnergyChanges[elecID] >= 0){
			energyGainProcesses[chosenProcessID] += electronEnergyChanges[elecID];
		}
		else{
			energyLossProcesses[chosenProcessID] += electronEnergyChanges[elecID];
		}

		// collect IDs that have ejected electrons
		if (processTypes[chosenProcessID] == GeneralDefinitions::ionizationType){
			ejectElectronIDs.push_back(elecID);
		}
		// collect IDs that have attached electrons
		else if (processTypes[chosenProcessID] == GeneralDefinitions::attachmentType){
			isAttached[elecID] = true;
			attachElectronIDs.push_back(elecID);
		}
	}

	int nEjected = ejectElectronIDs.size();

	// replace the attached electrons
	for (auto attachID: attachElectronIDs){
		// if there are still ejected electrons, put the last one in the place of the attached one
		if (nEjected > 0){
			int ejElecID = ejectElectronIDs.back();
			electronPositions.row(attachID) = ejectedElectronPositions.row(ejElecID);
			electronVelocities.row(attachID) = ejectedElectronVelocities.row(ejElecID);
			electronEnergies[attachID] = ejectedElectronEnergies[ejElecID];
			electronTimes[attachID] = electronTimes[ejElecID];			
			collisionFreeTimes[attachID] = Constant::NON_DEF;
			trialCollisionFrequenciesEachElectron[attachID] = trialCollisionFrequency;
			// the electron in this ID is no longer attached and can be used for the ensemble used for copying (done in the next "else")
			isAttached[attachID] = false;
			// eliminate the ID of the ejected electron
			ejectElectronIDs.pop_back();
			--nEjected;
		}
		// if there are no more ejected electrons, copy one electron that was NOT attached (this is important when attachment is strong)
		else{
			int randSelecID = std::fmin(MathFunctions::unitUniformRand(true,false)*nElectrons, nElectrons-1);
			while (isAttached[randSelecID]){
				randSelecID = std::fmin(MathFunctions::unitUniformRand(true,false)*nElectrons, nElectrons-1);
			}
			// add the energy change due to the growth profile
			energyGrowth += electronEnergies[randSelecID];
			// replace electron			
			electronPositions.row(attachID) = electronPositions.row(randSelecID);
			electronVelocities.row(attachID) = electronVelocities.row(randSelecID);
			electronEnergies[attachID] = electronEnergies[randSelecID];
			electronTimes[attachID] = electronTimes[randSelecID];			
			collisionFreeTimes[attachID] = collisionFreeTimes[randSelecID];
			trialCollisionFrequenciesEachElectron[attachID] = trialCollisionFrequenciesEachElectron[randSelecID];
			// note: here we do not update isAttached, since we want to copy from the original distribution
		}
	}

	// if there are more ejected electrons than attached ones
	while (nEjected > 0){
		// select an uniformly over the (nElectrons + nEjected) electrons
		int randSelecID = std::fmin(MathFunctions::unitUniformRand(true,false)*(nElectrons+nEjected), nElectrons+nEjected-1);
		// remove an electron from the "normal" ensemble and put there one of the ejected electrons
		if (randSelecID < nElectrons){
			int ejElecID = ejectElectronIDs.back();
			// add the energy change due to the growth profile
			energyGrowth -= electronEnergies[randSelecID];
			// replace electron
			electronPositions.row(randSelecID) = ejectedElectronPositions.row(ejElecID);
			electronVelocities.row(randSelecID) = ejectedElectronVelocities.row(ejElecID);
			electronEnergies[randSelecID] = ejectedElectronEnergies[ejElecID];
			electronTimes[randSelecID] = electronTimes[ejElecID];		
			collisionFreeTimes[randSelecID] = Constant::NON_DEF;
			trialCollisionFrequenciesEachElectron[randSelecID] = trialCollisionFrequency;
		}
		// remove from the nEjected electrons
		else{
			randSelecID -= nElectrons;
			// add the energy change due to the growth profile
			energyGrowth -= ejectedElectronEnergies[ejectElectronIDs[randSelecID]];
			// put the ID of the ejected electron in the last position, to be later removed
			std::swap(ejectElectronIDs[randSelecID],ejectElectronIDs[nEjected-1]);
		}
		// eliminate the ID of the ejected electron
		ejectElectronIDs.pop_back();
		--nEjected;		
	}
}

void BoltzmannMC::calculateMeanDataForSwarmParams(){
	// 'calculateMeanData' calculates the mean values over the electron ensemble at each sampling time to use later in the calculation of the time-averaged swarm parameters

	// allocate memory when the number of lines is equal to the sampling points. This is done to avoid many memory allocations, which are time-consuming
	if (meanEnergies.size() <= nSamplingPoints){
		samplingTimes.conservativeResize(2*nSamplingPoints);
		meanEnergies.conservativeResize(2*nSamplingPoints);
		meanPositions.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
		meanVelocities.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
		fluxDiffusionCoeffs.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
		positionCovariances.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
		bulkVelocities.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
		bulkDiffusionCoeffs.conservativeResize(2*nSamplingPoints, Eigen::NoChange);
	}

	// check if the maximum energy is higher than the maximum energy limit
	maxElecEnergy = std::fmax(maxElecEnergy, electronEnergies.maxCoeff());
	if (maxElecEnergy > energyMaxElastic){
		Message::error("The energy of the electrons reached " + std::to_string(maxElecEnergy) + " eV, while at least one of the elastic cross sections are defined only until "
			           + std::to_string(energyMaxElastic) + " eV. Please reduce the electric field or change the cross sections.");
	}

	meanEnergies[currentSamplingIndex] = electronEnergies.mean();
	meanPositions.row(currentSamplingIndex) = electronPositions.colwise().mean();
	meanVelocities.row(currentSamplingIndex) = electronVelocities.colwise().mean();

	// calculate the position covariances and the flux diffusion coeffs of the ensemble (at this time-instant)
	Eigen::Matrix3d posCovariance =  Eigen::Matrix3d::Zero(), fluxDiffCoeff =  Eigen::Matrix3d::Zero();
	Eigen::Vector3d electronPosition, electronVelocity;
	for (int i = 0; i < nElectrons; ++i){
		electronPosition = electronPositions.row(i);
		electronVelocity = electronVelocities.row(i);
		posCovariance += electronPosition*electronPosition.transpose();
		fluxDiffCoeff += electronPosition*electronVelocity.transpose();
	}
	posCovariance = posCovariance/nElectrons - meanPositions.row(currentSamplingIndex).matrix().transpose()*meanPositions.row(currentSamplingIndex).matrix();
	fluxDiffCoeff = fluxDiffCoeff/nElectrons - meanPositions.row(currentSamplingIndex).matrix().transpose()*meanVelocities.row(currentSamplingIndex).matrix();

	// save the matrices in a row: xx, xy, xz, yx, yy, yz, zx, zy, zz
	positionCovariances.block(currentSamplingIndex, 0, 1, 3) = posCovariance.row(0); 
	positionCovariances.block(currentSamplingIndex, 3, 1, 3) = posCovariance.row(1);
	positionCovariances.block(currentSamplingIndex, 6, 1, 3) = posCovariance.row(2);
	fluxDiffusionCoeffs.block(currentSamplingIndex, 0, 1, 3) = fluxDiffCoeff.row(0); 
	fluxDiffusionCoeffs.block(currentSamplingIndex, 3, 1, 3) = fluxDiffCoeff.row(1);
	fluxDiffusionCoeffs.block(currentSamplingIndex, 6, 1, 3) = fluxDiffCoeff.row(2);

	if (currentSamplingIndex != 0){
		bulkVelocities.row(currentSamplingIndex) = (meanPositions.row(currentSamplingIndex)-meanPositions.row(currentSamplingIndex-1))/(samplingTimes[currentSamplingIndex]-samplingTimes[currentSamplingIndex-1]);		
		bulkDiffusionCoeffs.row(currentSamplingIndex) = 0.5*(positionCovariances.row(currentSamplingIndex)-positionCovariances.row(currentSamplingIndex-1))/(samplingTimes[currentSamplingIndex]-samplingTimes[currentSamplingIndex-1]);
	}	
	else{
		// just for consistency, not relevant since it is not used for the temporal average
		bulkVelocities.row(0) = meanVelocities.row(0);
		bulkDiffusionCoeffs.row(0) = fluxDiffusionCoeffs.row(0);
	}

	// save results as a function of the period, when there is an AC frequency
	// Do this only when steady-state has been found
	if (excitationFrequencyRadians != 0 && steadyStateTime != Constant::NON_DEF){
		// reduce the phase to 2*pi
		double phase = std::fmod(excitationFrequencyRadians*time, 2.0*M_PI);
		// find the corresponding idx of the phase, avoiding problems if phase is exactly 2*pi
		int phaseIdx = std::fmin(std::floor(phase/integrationPhaseStep), nIntegrationPhases-1);
		// increase the number of points for the corresponding idx 
		nIntegrationPointsPerPhase[phaseIdx] += 1.0;
		// sum the other quantities for the corresponding idx
		meanEnergies_periodic[phaseIdx] += meanEnergies[currentSamplingIndex];
		fluxVelocities_periodic.row(phaseIdx) += meanVelocities.row(currentSamplingIndex);
		bulkVelocities_periodic.row(phaseIdx) += bulkVelocities.row(currentSamplingIndex);
		fluxDiffusionCoeffs_periodic.row(phaseIdx) += fluxDiffusionCoeffs.row(currentSamplingIndex);
		bulkDiffusionCoeffs_periodic.row(phaseIdx) += bulkDiffusionCoeffs.row(currentSamplingIndex);
	}
}

void BoltzmannMC::getTimeAverageEnergyParams(){
	// 'getTimeAverageEnergyParams' calculates the time-averaged mean energy

	// calculate the time-averaged mean energy
	averagedMeanEnergy = meanEnergies.segment(firstIntegrationIndex,nIntegrationPoints).mean();
	averagedMeanEnergyError = MathFunctions::statisticalError(meanEnergies.segment(firstIntegrationIndex,nIntegrationPoints), 50);
}

void BoltzmannMC::getTimeDependDistributions(){
	// 'getTimeDependDistributions' calculates the time-dependent quantities necessary to get the time-averaged distribution functions

	// check if the grid for the calculation of eedf and anisotropies is to be updated
	// we avoid to update the grid for evdf, since it is less relevant
	if (maxElecEnergy > maxEedfEnergy){
		// save temporarily the old grid variables
		double oldEedfEnergyStep = eedfEnergyStep;
		Eigen::ArrayXd oldEedfEnergyNodes = eedfEnergyNodes;
		Eigen::ArrayXd oldEehSum = eehSum;
		Eigen::ArrayXXd oldEahSum = eahSum;
		Eigen::ArrayXXd oldEehSum_periodic = eehSum_periodic;
		// update with an increasing factor of 1.2 to avoid constant updates
		maxEedfEnergy = 1.2*maxElecEnergy;
		// create the new grid and reset histograms to zero
		eedfEnergyNodes = Eigen::ArrayXd::LinSpaced(nEnergyCells+1, 0.0, maxEedfEnergy);
		eedfEnergyStep = eedfEnergyNodes[1];
		eedfEnergyCells = Eigen::ArrayXd::LinSpaced(nEnergyCells, eedfEnergyStep/2.0, maxEedfEnergy-eedfEnergyStep/2.0);
		eehSum = Eigen::ArrayXd::Zero(nEnergyCells);		
		eahSum = Eigen::ArrayXXd::Zero(nEnergyCells, nCosAngleCells);
		// map the values of the old grid into the new grid
		for (int posOldGrid = 0; posOldGrid < nEnergyCells; ++posOldGrid){
			double leftOldNodeEnergy = oldEedfEnergyNodes[posOldGrid];
			double rightOldNodeEnergy = oldEedfEnergyNodes[posOldGrid+1];
			// find the position of the nodes in the new grid (integer is rounded downwards)
			int newLeftPos = leftOldNodeEnergy/eedfEnergyStep;
			int newRightPos = rightOldNodeEnergy/eedfEnergyStep;
			// if the two limits are contained in the same cell
			if (newLeftPos == newRightPos){
				eehSum[newLeftPos] += oldEehSum[posOldGrid];
				eahSum.row(newLeftPos) += oldEahSum.row(posOldGrid);
				if (excitationFrequencyRadians != 0){
					for (int i = 0; i < nIntegrationPhases; ++i){
						eehSum_periodic.col(newLeftPos) += oldEehSum_periodic.col(posOldGrid);
					}
				}
			}
			// else, the right limit is already in the next cell	
			else{
				double leftFraction = (eedfEnergyNodes[newRightPos]-leftOldNodeEnergy)/oldEedfEnergyStep;
				eehSum[newLeftPos] += leftFraction*oldEehSum[posOldGrid];
				eehSum[newRightPos] += (1.0-leftFraction)*oldEehSum[posOldGrid];
				eahSum.row(newLeftPos) += leftFraction*oldEahSum.row(posOldGrid);
				eahSum.row(newRightPos) += (1.0-leftFraction)*oldEahSum.row(posOldGrid);
			}
		}

		if ((oldEehSum.sum()-eehSum.sum())/eehSum.sum() > 1E-9 || (oldEehSum_periodic.sum()-eehSum_periodic.sum())/eehSum_periodic.sum() > 1E-9){
			Message::error(std::string("old normalizer: ") + std::to_string(oldEehSum.sum()) + "\nnew normalizer: " + std::to_string(eehSum.sum()));
		} 
	}

	// create the current electron energy histogram (eeh) and sum it to the previous ones
	eeh = MathFunctions::histogramCount(electronEnergies, eedfEnergyNodes); 
	eehSum += eeh;
	// for conditions with AC E-field, do sampling along period
	if (excitationFrequencyRadians != 0){
		// reduce the phase to 2*pi
		double phase = std::fmod(excitationFrequencyRadians*time, 2.0*M_PI);
		// find the corresponding idx of the phase, avoiding problems if phase is exactly 2*pi
		int phaseIdx = std::fmin(std::floor(phase/integrationPhaseStep), nIntegrationPhases-1); 		
		eehSum_periodic.row(phaseIdx) += eeh;
	}	

	if (isCylindricallySymmetric){
		// calculate the cos(angles) of the velocities relatively to the z axis
		electronCosAngles = electronVelocities.col(2) / electronVelocities.matrix().rowwise().norm().array();

		// create the current electron angular histogram (eah) and sum it to the previous ones
		eahSum += MathFunctions::histogram2DCount(electronEnergies, electronCosAngles, eedfEnergyNodes, cosAngleNodes);

		// create the current electron velocity histogram (evh) and sum it to the previous ones
		evhSum += MathFunctions::histogram2DCount(Eigen::sqrt(electronVelocities.col(0).square() + electronVelocities.col(1).square()), electronVelocities.col(2), radialVelocityNodes, axialVelocityNodes);
	}
}

void BoltzmannMC::getTimeAverageDistributions(){
	// 'getTimeAverageDistributions' calculates the time-averaged distribution functions (performed only at end of the simulation, no impact in performance)

	double normalizer = eehSum.sum();

	// get the electron energy distribution function (eedf)
	eedf = eehSum/(Eigen::sqrt(eedfEnergyCells)*normalizer*eedfEnergyStep);

	if (isCylindricallySymmetric){
		// get the electron angular distribution function (eadf)
		for (int i = 0; i < nEnergyCells; ++i){
			eadf.row(i) = 2.0*eahSum.row(i) / (std::sqrt(eedfEnergyCells[i]) * normalizer*eedfEnergyStep*cosAngleStep);
		}

		// get the electron first-anisotropy distribution function (efadf)
		Eigen::ArrayXd legendreWeights = cosAngleCells;
		double anisotropyDegree = 1;
		for (int i = 0; i < nEnergyCells; ++i){
			efadf[i] = (anisotropyDegree+0.5) * (eadf.row(i).transpose()*legendreWeights).sum()*cosAngleStep;
		}

		// get the electron second-anisotropy distribution function (efadf)
		legendreWeights = 0.5*(3.0*cosAngleCells.square() -1);
		anisotropyDegree = 2;
		for (int i = 0; i < nEnergyCells; ++i){
			esadf[i] = (anisotropyDegree+0.5) * (eadf.row(i).transpose()*legendreWeights).sum()*cosAngleStep;
		}

		// get the electron velocity distribution function (evdf)
		evdf = evhSum / (evhSum.rowwise().sum()*axialVelocityStep*2.0*M_PI*radialVelocityCells*radialVelocityStep).sum();
	}
}

void BoltzmannMC::getTimeAverageFluxParams(){
	// 'getTimeAverageFluxParams' calculates the time-average of the flux drift velocity and the flux diffusion coefficients

	// calculate the time-averaged flux drift velocity
	averagedFluxDriftVelocity = meanVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3).colwise().mean();
	averagedFluxDriftVelocityError = MathFunctions::statisticalErrorColwise(meanVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3), 50);

	// calculate the time-averaged flux diffusion coeffs (9 components)
	Eigen::ArrayXd averagedFluxDiffusionCoeffsLine = fluxDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9).colwise().mean();
	Eigen::ArrayXd averagedFluxDiffusionCoeffsErrorLine = MathFunctions::statisticalErrorColwise(fluxDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9), 50);

	// convert them into matrix
	averagedFluxDiffusionCoeffs.row(0) = averagedFluxDiffusionCoeffsLine.segment(0,3);
	averagedFluxDiffusionCoeffs.row(1) = averagedFluxDiffusionCoeffsLine.segment(3,3);
	averagedFluxDiffusionCoeffs.row(2) = averagedFluxDiffusionCoeffsLine.segment(6,3);
	averagedFluxDiffusionCoeffsError.row(0) = averagedFluxDiffusionCoeffsErrorLine.segment(0,3);
	averagedFluxDiffusionCoeffsError.row(1) = averagedFluxDiffusionCoeffsErrorLine.segment(3,3);
	averagedFluxDiffusionCoeffsError.row(2) = averagedFluxDiffusionCoeffsErrorLine.segment(6,3);	
}

void BoltzmannMC::getTimeAverageRateCoeffs(){
	// 'getTimeAverageRateCoeffs' calculates the time-averaged rate-coefficients (performed only at end of the simulation, no impact in performance)

	for (int i = 0; i < nProcesses; ++i){
		if (relDensities[i] != 0){
			averagedRateCoeffs[i] = (double)collisionCounters[i]/(totalIntegratedTime*nElectrons*relDensities[i]*totalGasDensity);
		}
		else{
			averagedRateCoeffs[i] = 0;
		}
	}
}

void BoltzmannMC::getTimeAveragePowerBalance(){ 
	// 'getTimeAveragePowerBalance' calculates the time-averaged power balance in [eV/electron/s] (performed only at end of the simulation, no impact in performance)

	averagedPowerGainField = energyGainField/(nElectrons*totalIntegratedTime*totalGasDensity);
	averagedPowerGrowth = energyGrowth/(nElectrons*totalIntegratedTime*totalGasDensity);
	for (int i = 0; i < nProcesses; ++i){
		averagedPowerGainProcesses[i] = energyGainProcesses[i]/(nElectrons*totalIntegratedTime*totalGasDensity);
		averagedPowerLossProcesses[i] = energyLossProcesses[i]/(nElectrons*totalIntegratedTime*totalGasDensity);
	}
}

void BoltzmannMC::getTimeAverageBulkParams(){
	// 'getTimeAverageBulkParams' calculates the bulk swarm parameters from the swarm's center of mass and of the position covariances

	// calculate the time-averaged bulk drift velocity
	averagedBulkDriftVelocity = bulkVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3).colwise().mean();
	averagedBulkDriftVelocityError = MathFunctions::statisticalErrorColwise(bulkVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3), 50);

	// calculate the time-averaged bulk diffusion coeffs (9 components)
	Eigen::ArrayXd averagedBulkDiffusionCoeffsLine = bulkDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9).colwise().mean();
	Eigen::ArrayXd averagedBulkDiffusionCoeffsErrorLine = MathFunctions::statisticalErrorColwise(bulkDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9), 50);

	// convert them into matrix
	averagedBulkDiffusionCoeffs.row(0) = averagedBulkDiffusionCoeffsLine.segment(0,3);
	averagedBulkDiffusionCoeffs.row(1) = averagedBulkDiffusionCoeffsLine.segment(3,3);
	averagedBulkDiffusionCoeffs.row(2) = averagedBulkDiffusionCoeffsLine.segment(6,3);
	averagedBulkDiffusionCoeffsError.row(0) = averagedBulkDiffusionCoeffsErrorLine.segment(0,3);
	averagedBulkDiffusionCoeffsError.row(1) = averagedBulkDiffusionCoeffsErrorLine.segment(3,3);
	averagedBulkDiffusionCoeffsError.row(2) = averagedBulkDiffusionCoeffsErrorLine.segment(6,3);	

	/*double c0, c1, cov00, cov01, cov11, sumsq;

	int nSlices = 20;
	int nSlicePoints = nIntegrationPoints/nSlices;

	// ----- perform fits 'f(t) = c0 + c1 * t', where f = center of mass and t = integration time ----- //

	Eigen::ArrayXd integrationTimes = samplingTimes.segment(firstIntegrationIndex,nIntegrationPoints);
	Eigen::ArrayXXd centerOfMass = meanPositions.block(firstIntegrationIndex,0,nIntegrationPoints,3);

	for (int i = 0; i < 3; ++i){
		gsl_fit_linear(integrationTimes.data(), 1, centerOfMass.col(i).data(), 1, nIntegrationPoints, &c0, &c1, &cov00, &cov01 , &cov11, &sumsq);
		averagedBulkDriftVelocity[i] = c1;
	}

	// divide the integration time in 20 slices, in order to estimate the error
	// this must be done, since the errors from the fit are not realistic
	Eigen::ArrayXd sum = Eigen::ArrayXd::Zero(3), squaredSum = Eigen::ArrayXd::Zero(3);
	Eigen::ArrayXd sliceTimes;
	for (int i = 0; i < nSlices; ++i){
		centerOfMass = meanPositions.block(i*nSlicePoints,0,nSlicePoints,3);
		sliceTimes = integrationTimes.segment(i*nSlicePoints,nSlicePoints);
		for (int j = 0; j < 3; ++j){
			gsl_fit_linear(sliceTimes.data(), 1, centerOfMass.col(j).data(), 1, nSlicePoints, &c0, &c1, &cov00, &cov01 , &cov11, &sumsq);
			sum[j] += c1;
			squaredSum[j] += c1*c1;
		}
	}
	Eigen::ArrayXd mean = sum / (double) nSlices;
	averagedBulkDriftVelocityError = Eigen::sqrt((squaredSum/(double)nSlices - mean.square())/nSlices);	

	// ----- perform fits 'f(t) = c0 + c1 * x', where f = half of the covariance and x = integration time ----- //
	// the diffusion coeffs are organized in a row: xx xy xz yx yy yz zx zy zz

	Eigen::ArrayXXd halfCovariances = 0.5*positionCovariances;
	Eigen::ArrayXd averagedBulkDiffusionCoeffsLine(9);

	for (int i = 0; i < 9; ++i){
		gsl_fit_linear(integrationTimes.data(), 1, halfCovariances.col(i).data(), 1, nIntegrationPoints, &c0, &c1, &cov00, &cov01 , &cov11, &sumsq);
		averagedBulkDiffusionCoeffsLine[i] = c1;
	}

	// divide the integration time in 20 slices, in order to estimate the error
	// this must be done, since the errors from the fit are not realistic
	sum = Eigen::ArrayXd::Zero(9), squaredSum = Eigen::ArrayXd::Zero(9);
	for (int i = 0; i < nSlices; ++i){
		halfCovariances = 0.5*positionCovariances.block(i*nSlicePoints,0,nSlicePoints,9);
		sliceTimes = integrationTimes.segment(i*nSlicePoints,nSlicePoints);
		for (int j = 0; j < 9; ++j){
			gsl_fit_linear(sliceTimes.data(), 1, halfCovariances.col(j).data(), 1, nSlicePoints, &c0, &c1, &cov00, &cov01 , &cov11, &sumsq);
			sum[j] += c1;
			squaredSum[j] += c1*c1;
		}
	}
	mean = sum / (double) nSlices;
	Eigen::ArrayXd averagedBulkDiffusionCoeffsErrorLine = Eigen::sqrt((squaredSum/(double)nSlices - mean.square())/nSlices);*/	
}

void BoltzmannMC::getAveragedPeriodicParams(){
	// obtain the average values per phase: divide the summed quantities by the number of integration points per phase
	meanEnergies_periodic /= nIntegrationPointsPerPhase;
	fluxVelocities_periodic.colwise() /= nIntegrationPointsPerPhase;
	bulkVelocities_periodic.colwise() /= nIntegrationPointsPerPhase;
	fluxDiffusionCoeffs_periodic.colwise() /= nIntegrationPointsPerPhase;
	bulkDiffusionCoeffs_periodic.colwise() /= nIntegrationPointsPerPhase;
	// get the electron energy distribution functions (eedf) along the period
	Eigen::ArrayXd normalizer = eehSum_periodic.rowwise().sum();
	Eigen::ArrayXd aux = Eigen::sqrt(eedfEnergyCells);
	for (int i = 0; i < nIntegrationPhases; ++i){
		eedf_periodic.row(i) = eehSum_periodic.row(i)/(aux.matrix().transpose().array()*normalizer[i]*eedfEnergyStep);
	}
}

void BoltzmannMC::checkStatisticalErrors(){
	// 'checkStatisticalErrors' checks if the minimum errors specified by the user on the setup are already verified

	getTimeAverageEnergyParams();
	getTimeAverageFluxParams();
	getTimeAverageBulkParams();
	checkPowerBalance();

	double relErrorAbsFluxVelocity = (std::abs(averagedFluxDriftVelocity[0])*averagedFluxDriftVelocityError[0] + std::abs(averagedFluxDriftVelocity[1])*averagedFluxDriftVelocityError[1] + 
								  std::abs(averagedFluxDriftVelocity[2])*averagedFluxDriftVelocityError[2]) / averagedFluxDriftVelocity.matrix().squaredNorm();
	double relErrorAbsBulkVelocity = (std::abs(averagedBulkDriftVelocity[0])*averagedBulkDriftVelocityError[0] + std::abs(averagedBulkDriftVelocity[1])*averagedBulkDriftVelocityError[1] + 
								  std::abs(averagedBulkDriftVelocity[2])*averagedBulkDriftVelocityError[2]) / averagedBulkDriftVelocity.matrix().squaredNorm();								  
	if (averagedMeanEnergyError/averagedMeanEnergy <= requiredMeanEnergyRelError &&
		relErrorAbsFluxVelocity <= requiredFluxDriftVelocityRelError &&
		averagedFluxDiffusionCoeffsError(0,0)/averagedFluxDiffusionCoeffs(0,0) <= requiredFluxDiffusionCoeffsRelError &&
		averagedFluxDiffusionCoeffsError(1,1)/averagedFluxDiffusionCoeffs(1,1) <= requiredFluxDiffusionCoeffsRelError &&
		averagedFluxDiffusionCoeffsError(2,2)/averagedFluxDiffusionCoeffs(2,2) <= requiredFluxDiffusionCoeffsRelError &&
		relErrorAbsBulkVelocity <= requiredBulkDriftVelocityRelError &&
		averagedBulkDiffusionCoeffsError(0,0)/averagedBulkDiffusionCoeffs(0,0) <= requiredBulkDiffusionCoeffsRelError &&
		averagedBulkDiffusionCoeffsError(1,1)/averagedBulkDiffusionCoeffs(1,1) <= requiredBulkDiffusionCoeffsRelError &&
		averagedBulkDiffusionCoeffsError(2,2)/averagedBulkDiffusionCoeffs(2,2) <= requiredBulkDiffusionCoeffsRelError &&		
		powerBalanceRelError <= requiredPowerBalanceRelError){
		goodStatisticalErrors = true;
	}
}

void BoltzmannMC::checkPowerBalance(){
	// 'checkPowerBalance' calculates the relative error of the power balance

	double totalGain = 0, totalLoss = 0;
	totalGain += energyGainField;
	if (energyGrowth > 0){
		totalGain += energyGrowth;
	}
	else{
		totalLoss += energyGrowth;
	}
	for (int i = 0; i < nProcesses; ++i){
		totalGain += energyGainProcesses[i];
		totalLoss += energyLossProcesses[i];
	}
	powerBalanceRelError = std::abs(totalGain + totalLoss)/totalGain;
}

void BoltzmannMC::checkSteadyState(){
	// 'checkSteadyState' verifies if the steady-state has been reached, by checking if 
	// the mean energy of the interval 50-75% is larger than the one of the interval 75-100%

	// determine the mean value of [0.5*time,0.75*time] and [0.75*time,time], and the std of [0.75*time,time]
	double t1 = 0.5*time, t2 = 0.75*time;
	double averageFirstPart = 0;
	double averageSecondPart = 0;
	double relStdSecondPart = 0;
	int nFirst = 0, nSecond = 0;
	for (int i = 0; i < nSamplingPoints; ++i){
		double sampTime = samplingTimes[i];
		if (sampTime >= t1 && sampTime <= t2){
			averageFirstPart += meanEnergies[i];
			++nFirst;
		}
		else if(sampTime > t2){
			double val = meanEnergies[i];
			averageSecondPart += val;
			relStdSecondPart += val*val;
			++nSecond;
		}	
	}
	averageFirstPart /= nFirst;
	averageSecondPart /= nSecond;
	relStdSecondPart = std::sqrt((relStdSecondPart/nSecond-averageSecondPart*averageSecondPart)/nSecond);

	// check both averages and the std of the second part
	if ((averageFirstPart >= averageSecondPart && relStdSecondPart < 0.01)  || totalCollisionCounter >= maxCollisionsBeforeSteadyState){ 

		if (totalCollisionCounter >= maxCollisionsBeforeSteadyState){
			if (dispMCStatus){
				std::cout<<" "<<std::endl;
				for (int i = 0; i < 30; ++i){
					// move up in the terminal
					std::printf("%c[1A", 0x1B);
					// clear terminal line
					std::printf("%c[2K", 0x1B);
				}
				Message::warning("Steady-state imposed by the parameter 'maxCollisionsBeforeSteadyState'. The electron energy may not be stabilized yet. [E/N = " + std::to_string(reducedElecField) + " Td, excitFreq = " + std::to_string(excitationFrequency) + " Hz, E-angle = " + std::to_string(elecFieldAngle) + " deg, B/N = " + std::to_string(reducedMagField) + " Hx]\n");
				for (int i = 0; i < 29; ++i){
					std::printf("\n");
				}
				dispInfo();
			}
			else{
				Message::warning("Steady-state imposed by the parameter 'maxCollisionsBeforeSteadyState'. The electron energy may not be stabilized yet. [E/N = " + std::to_string(reducedElecField) + " Td, excitFreq = " + std::to_string(excitationFrequency) + " Hz, E-angle = " + std::to_string(elecFieldAngle) + " deg, B/N = " + std::to_string(reducedMagField) + " Hx]\n");
			}
		}

		// save the first integration index
		firstIntegrationIndex = currentSamplingIndex;
		nIntegrationPoints = 1;

		// save the steady state time
		steadyStateTime = time;

		// save the number of collisions at steady state
		collisionCounterAtSS = totalCollisionCounter;
		nullCollisionCounterAtSS = nullCollisionCounter;

		// reset to zero the collision counters
		// reset to zero the net energies per collision
		energyGainField = 0;
		for (int i = 0; i < nProcesses; ++i){
			collisionCounters[i] = 0;
			energyLossProcesses[i] = 0;
			energyGainProcesses[i] = 0;
		}
		energyGrowth = 0;

		// initialize the integration time
		totalIntegratedTime = 0;

		// initialize all variables related with the energy swarm parameters
		maxEedfEnergy = 1.2*maxElecEnergy;
		eedfEnergyNodes = Eigen::ArrayXd::LinSpaced(nEnergyCells+1, 0.0, maxEedfEnergy);
		eedfEnergyStep = eedfEnergyNodes[1];
		eedfEnergyCells = Eigen::ArrayXd::LinSpaced(nEnergyCells, eedfEnergyStep/2.0, maxEedfEnergy-eedfEnergyStep/2.0);
		eehSum = Eigen::ArrayXd::Zero(nEnergyCells); eedf = Eigen::ArrayXd::Zero(nEnergyCells);

		cosAngleNodes = Eigen::ArrayXd::LinSpaced(nCosAngleCells+1, -1.0, 1.0);
		cosAngleStep = cosAngleNodes[1] + 1.0;
		cosAngleCells = Eigen::ArrayXd::LinSpaced(nCosAngleCells, -1.0+cosAngleStep/2.0, 1.0-cosAngleStep/2.0);
		eahSum = Eigen::ArrayXXd::Zero(nEnergyCells, nCosAngleCells); eadf = Eigen::ArrayXXd::Zero(nEnergyCells, nCosAngleCells);
		efadf = Eigen::ArrayXd::Zero(nEnergyCells);
		esadf = Eigen::ArrayXd::Zero(nEnergyCells);

		// create the velocity nodes between [-maxSpeed, maxSpeed]
		double maxSpeed = std::sqrt(2.0*maxEedfEnergy*Constant::electronCharge/Constant::electronMass);
		radialVelocityNodes = Eigen::ArrayXd::LinSpaced(nRadialVelocityCells+1, 0, maxSpeed); 
		radialVelocityStep = radialVelocityNodes[1];
		radialVelocityCells = Eigen::ArrayXd::LinSpaced(nRadialVelocityCells, radialVelocityStep/2.0, maxSpeed-radialVelocityStep/2.0);
		axialVelocityNodes = Eigen::ArrayXd::LinSpaced(nAxialVelocityCells+1, -maxSpeed, maxSpeed); 
		axialVelocityStep = axialVelocityNodes[1]+maxSpeed;
		axialVelocityCells = Eigen::ArrayXd::LinSpaced(nAxialVelocityCells, -maxSpeed+axialVelocityStep/2.0, maxSpeed-axialVelocityStep/2.0);		
		evhSum = Eigen::ArrayXXd::Zero(nRadialVelocityCells, nAxialVelocityCells); evdf = Eigen::ArrayXXd::Zero(nRadialVelocityCells, nAxialVelocityCells);

		// for AC E-field conditions (integration along the period)
		if (excitationFrequencyRadians != 0){
			eehSum_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, nEnergyCells);
			eedf_periodic = Eigen::ArrayXXd::Zero(nIntegrationPhases, nEnergyCells);
		}

		getTimeDependDistributions();
	}
}

void BoltzmannMC::dispInfo(){
	static bool firstEvaluation = true;
	static char terminal_moveup[10], terminal_clearline[10];

	if (firstEvaluation){
		firstEvaluation = false;
		std::sprintf(terminal_moveup, "%c[1A", 0x1B);
		std::sprintf(terminal_clearline, "%c[2K", 0x1B);
		std::cout<<"\n*********Simulation status*********\n\n";
	}
	else{
		std::cout<<" "<<std::endl;
		for (int i = 0; i < 30 ; ++i){
			std::printf("%s",terminal_moveup);
			std::printf("%s",terminal_clearline);
		}
	}

	std::cout<<"E/N: "<<reducedElecField<<" Td\n";
	std::cout<<"excitFreq: "<<excitationFrequency<<" Hz\n";
	std::cout<<"E-field angle: "<<elecFieldAngle<<" degrees\n";
	std::cout<<"B/N: "<<reducedMagField<<" Hx\n\n";
	std::cout<<"Number of real collisions: "<<(double)totalCollisionCounter<<"\n";
	std::cout<<"Number of null collisions: "<<(double)nullCollisionCounter<<"\n\n";
	std::cout<<"Current time: "<<time<<" s\n";
	if (steadyStateTime == Constant::NON_DEF){
		std::cout<<"Mean energy: "<<meanEnergies[currentSamplingIndex]<<"\n";
		for (int i = 0; i < 19; ++i){
			std::printf("\n");
		}	
	}
	else{
		std::cout<<"Steady-state time: "<<steadyStateTime<<" s\n\n";
		if (nIntegrationPoints >= 3*nPointsBetweenStatErrorsCheck && averagedMeanEnergy != Constant::NON_DEF){
			std::cout<<"Number of integration points: "<<(double)nIntegrationPoints<<"\n\n";
			std::cout<<"Mean energy [eV]: "<<averagedMeanEnergy<<std::endl;
			std::cout<<"Relative error: "<<averagedMeanEnergyError/averagedMeanEnergy<<"\n\n";
			std::cout<<"Flux drift velocity [m/s]: "<<averagedFluxDriftVelocity.transpose()<<std::endl;
			std::cout<<"Relative error: "<<(averagedFluxDriftVelocityError/Eigen::abs(averagedFluxDriftVelocity)).transpose()<<"\n\n";
			std::cout<<"Bulk drift velocity [m/s]: "<<averagedBulkDriftVelocity.transpose()<<std::endl;
			std::cout<<"Relative error: "<<(averagedBulkDriftVelocityError/Eigen::abs(averagedBulkDriftVelocity)).transpose()<<"\n\n";			
			std::cout<<"Flux diffusion coefficients [m^2 s^-1]: "<<averagedFluxDiffusionCoeffs(0,0)<<" "<<averagedFluxDiffusionCoeffs(1,1)<<" "<<averagedFluxDiffusionCoeffs(2,2)<<std::endl;
			std::cout<<"Relative error: "<<averagedFluxDiffusionCoeffsError(0,0)/averagedFluxDiffusionCoeffs(0,0)<<" "<<averagedFluxDiffusionCoeffsError(1,1)/averagedFluxDiffusionCoeffs(1,1)<<" "<<averagedFluxDiffusionCoeffsError(2,2)/averagedFluxDiffusionCoeffs(2,2)<<"\n\n";
			std::cout<<"Bulk diffusion coefficients [m^2 s^-1]: "<<averagedBulkDiffusionCoeffs(0,0)<<" "<<averagedBulkDiffusionCoeffs(1,1)<<" "<<averagedBulkDiffusionCoeffs(2,2)<<std::endl;
			std::cout<<"Relative error: "<<averagedBulkDiffusionCoeffsError(0,0)/averagedBulkDiffusionCoeffs(0,0)<<" "<<averagedBulkDiffusionCoeffsError(1,1)/averagedBulkDiffusionCoeffs(1,1)<<" "<<averagedBulkDiffusionCoeffsError(2,2)/averagedBulkDiffusionCoeffs(2,2)<<"\n\n";			
			std::cout<<"Power balance relative error: "<<powerBalanceRelError<<"\n";
		}
		else{
			for (int i = 0; i < 18; ++i){
				std::printf("\n");
			}	
		}
	}
}

void BoltzmannMC::evaluatePower(){
	// initialize power structure;
	std::map<std::string,double> powerMap;
	powerMap["field"] = 0; powerMap["elasticNet"] = 0; powerMap["elasticGain"] = 0; powerMap["elasticLoss"] = 0; powerMap["carNet"] = 0; powerMap["carGain"] = 0; powerMap["carLoss"] = 0;
	powerMap["excitationIne"] = 0; powerMap["excitationSup"] = 0; powerMap["excitationNet"] = 0; powerMap["vibrationalIne"] = 0; powerMap["vibrationalSup"] = 0; powerMap["vibrationalNet"] = 0;
	powerMap["rotationalIne"] = 0; powerMap["rotationalSup"] = 0; powerMap["rotationalNet"] = 0; powerMap["ionizationIne"] = 0; powerMap["attachmentIne"] = 0; powerMap["inelastic"] = 0;
	powerMap["superelastic"] = 0; powerMap["eDensGrowth"] = 0; powerMap["electronElectron"] = 0;

	// evaluate power absorved per electron at unit gas density due to the electric field
	powerMap["field"] = averagedPowerGainField;

	// evaluate power absorved per electron at unit gas density due to elastic processes
	for (int i = 0; i < nProcesses; ++i){
		if (realCollisionPointers[i]->type == "Elastic"){
			powerMap["elasticGain"] += averagedPowerGainProcesses[i];
			powerMap["elasticLoss"] += averagedPowerLossProcesses[i];
		}
	}
	powerMap["elasticNet"] = powerMap["elasticGain"] + powerMap["elasticLoss"];

	// evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic processes
	// loop over each gas in the mixture
	for (auto gas: gasArray){
		std::string gasName = gas->name;
		std::map<std::string,double> gasPower;
		// initialize power balance information of this gas
		gasPower["excitationIne"] = 0; gasPower["excitationSup"] = 0; gasPower["excitationNet"] = 0;
		gasPower["vibrationalIne"] = 0; gasPower["vibrationalSup"] = 0; gasPower["vibrationalNet"] = 0;
		gasPower["rotationalIne"] = 0; gasPower["rotationalSup"] = 0; gasPower["rotationalNet"] = 0;
		gasPower["ionizationIne"] = 0; gasPower["attachmentIne"] = 0;

		// loop over all processes
		for (int i = 0; i < nProcesses; ++i){
			// check if this gas is the target in the present process
			if (gas == gasPointers[i]){
				if (realCollisionPointers[i]->type == "Excitation"){
					if (!isSuperElastic[i]){
						gasPower["excitationIne"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
					else{
						gasPower["excitationSup"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
				}
				else if (realCollisionPointers[i]->type == "Vibrational"){
					if (!isSuperElastic[i]){
						gasPower["vibrationalIne"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
					else{
						gasPower["vibrationalSup"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
				}
				else if (realCollisionPointers[i]->type == "Rotational"){
					if (!isSuperElastic[i]){
						gasPower["rotationalIne"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
					else{
						gasPower["rotationalSup"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
					}
				}
				else if (realCollisionPointers[i]->type == "Ionization"){
					gasPower["ionizationIne"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
				}
				else if (realCollisionPointers[i]->type == "Attachment"){
					gasPower["attachmentIne"] += averagedPowerGainProcesses[i] + averagedPowerLossProcesses[i];
				}
			}
		}
		// evaluate net values (for each gas)
		gasPower["excitationNet"] = gasPower["excitationIne"] + gasPower["excitationSup"];
		gasPower["vibrationalNet"] = gasPower["vibrationalIne"] + gasPower["vibrationalSup"];
		gasPower["rotationalNet"] = gasPower["rotationalIne"] + gasPower["rotationalSup"];
		gasPower["inelastic"] = gasPower["excitationIne"] + gasPower["vibrationalIne"] + gasPower["rotationalIne"] + gasPower["ionizationIne"] + gasPower["attachmentIne"];
		gasPower["superelastic"] = gasPower["excitationSup"] + gasPower["vibrationalSup"] + gasPower["rotationalSup"];

		// evaluate net values (for the gas mixture)
		powerMap["excitationIne"] += gasPower["excitationIne"];
		powerMap["excitationSup"] += gasPower["excitationSup"];
		powerMap["vibrationalIne"] += gasPower["vibrationalIne"];
		powerMap["vibrationalSup"] += gasPower["vibrationalSup"];
		powerMap["rotationalIne"] += gasPower["rotationalIne"];
		powerMap["rotationalSup"] += gasPower["rotationalSup"];
		powerMap["ionizationIne"] += gasPower["ionizationIne"];
		powerMap["attachmentIne"] += gasPower["attachmentIne"];

		// store power balance information of this gas in the main structure
		power.gasesMap[gasName]	= gasPower;
	}

	powerMap["excitationNet"] = powerMap["excitationIne"] + powerMap["excitationSup"];
	powerMap["vibrationalNet"] = powerMap["vibrationalIne"] + powerMap["vibrationalSup"];
	powerMap["rotationalNet"] = powerMap["rotationalIne"] + powerMap["rotationalSup"];
	powerMap["inelastic"] = powerMap["excitationIne"] + powerMap["vibrationalIne"] + powerMap["rotationalIne"] + powerMap["ionizationIne"] + powerMap["attachmentIne"];
	powerMap["superelastic"] = powerMap["excitationSup"] + powerMap["vibrationalSup"] + powerMap["rotationalSup"];
	powerMap["eDensGrowth"] = averagedPowerGrowth;

	// calculate the total power balance
	powerMap["balance"] = powerMap["field"] + powerMap["elasticNet"] + powerMap["inelastic"] + powerMap["superelastic"] + powerMap["eDensGrowth"];

	// calculate the relative power balance
	std::vector<double> powerValues({powerMap["field"], powerMap["elasticGain"], powerMap["elasticLoss"],powerMap["excitationSup"], powerMap["excitationIne"], 
		                        powerMap["vibrationalSup"], powerMap["vibrationalIne"], powerMap["rotationalSup"], powerMap["rotationalIne"], powerMap["eDensGrowth"]});
	double totalGain = 0;
	for (auto powerValue: powerValues){
		if (powerValue > 0){
			totalGain += powerValue;
		}
	}	
	powerMap["relativeBalance"] = std::abs(powerMap["balance"]) / totalGain;
	powerMap["reference"] = totalGain;

	// save power balance information
	power.Map = powerMap;
}

void BoltzmannMC::evaluateSwarmParameters(){

	double reducedElecFieldSI = workCond->reducedElecFieldSI;

	// create the rotation matrix to put the parameters in a frame where the E-field is along z (rotation in turn of y, https://en.wikipedia.org/wiki/Rotation_matrix)
	double rotAngle = -elecFieldAngle/180.0*M_PI;
	Eigen::Matrix3d rotationMatrix;
	rotationMatrix<<std::cos(rotAngle), 0, std::sin(rotAngle),
					0, 1, 0,
					-std::sin(rotAngle), 0, std::cos(rotAngle);	
	Eigen::Matrix3d absRotationMatrix = rotationMatrix.array().abs(); // for the maximization of error of rotated quantities

	// ----- evaluate flux parameters ----- // 	

	// rotate the drift velocity
	rotatedAveragedFluxDriftVelocity = rotationMatrix*averagedFluxDriftVelocity.matrix();
	rotatedAveragedFluxDriftVelocityError = absRotationMatrix*averagedFluxDriftVelocityError.matrix(); // maximization of the error

	// rotate the diffusion matrix		
	rotatedAveragedFluxDiffusionCoeffs = rotationMatrix*averagedFluxDiffusionCoeffs*rotationMatrix.transpose();	
	rotatedAveragedFluxDiffusionCoeffsError = absRotationMatrix*averagedFluxDiffusionCoeffsError*absRotationMatrix.transpose(); // maximization of the error

	// reduced transverse and longitudinal diffusion coeffients (relative to E-field)
	swarmParam["fluxRedTransvDiffCoeff"] = totalGasDensity * (rotatedAveragedFluxDiffusionCoeffs(0,0) + rotatedAveragedFluxDiffusionCoeffs(1,1))/2.0;
	swarmParam["fluxRedTransvDiffCoeffError"] = totalGasDensity * (rotatedAveragedFluxDiffusionCoeffsError(0,0) + rotatedAveragedFluxDiffusionCoeffsError(1,1))/2.0;
	swarmParam["fluxRedLongDiffCoeff"] = totalGasDensity * rotatedAveragedFluxDiffusionCoeffs(2,2);
	swarmParam["fluxRedLongDiffCoeffError"] = totalGasDensity * rotatedAveragedFluxDiffusionCoeffsError(2,2);

	// total ionization and attachment rate-coefficients
	swarmParam["totalIonRateCoeff"] = 0;
	swarmParam["totalAttRateCoeff"] = 0;
	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){
			if (collision->type == "Ionization"){
				swarmParam["totalIonRateCoeff"] += collision->target->density * collision->ineRateCoeff;
			}
			else if (collision->type == "Attachment"){
				swarmParam["totalAttRateCoeff"] += collision->target->density * collision->ineRateCoeff;
			}
		}
	}

	if (excitationFrequency == 0){
		// reduced mobility
		double fluxVelocAlongEfield = std::abs(rotatedAveragedFluxDriftVelocity[2]);
		swarmParam["fluxRedMobCoeff"] = fluxVelocAlongEfield/reducedElecFieldSI;
		swarmParam["fluxRedMobCoeffError"] = std::abs(rotatedAveragedFluxDriftVelocityError[2])/reducedElecFieldSI;

		// reduced Townsend coefficient and reduced attachment coefficient
		swarmParam["fluxRedTownsendCoeff"] = swarmParam["totalIonRateCoeff"]/fluxVelocAlongEfield;
		swarmParam["fluxRedAttCoeff"] = swarmParam["totalAttRateCoeff"]/fluxVelocAlongEfield;

		// characteristic energy
		swarmParam["fluxCharacEnergy"] = swarmParam["fluxRedTransvDiffCoeff"] / swarmParam["fluxRedMobCoeff"];
		swarmParam["fluxCharacEnergyError"] = swarmParam["fluxRedTransvDiffCoeffError"]/swarmParam["fluxRedMobCoeff"] +
											swarmParam["fluxRedMobCoeffError"]*swarmParam["fluxRedTransvDiffCoeff"]/std::pow(swarmParam["fluxRedMobCoeff"], 2);
	}

	// ----- evaluate bulk parameters ----- // 

	// rotate the drift velocity
	rotatedAveragedBulkDriftVelocity = rotationMatrix*averagedBulkDriftVelocity.matrix();
	rotatedAveragedBulkDriftVelocityError = absRotationMatrix*averagedBulkDriftVelocityError.matrix(); // maximization of the error

	// rotate the diffusion matrix		
	rotatedAveragedBulkDiffusionCoeffs = rotationMatrix*averagedBulkDiffusionCoeffs*rotationMatrix.transpose();	
	rotatedAveragedBulkDiffusionCoeffsError = absRotationMatrix*averagedBulkDiffusionCoeffsError*absRotationMatrix.transpose(); // maximization of the error

	// reduced transverse and longitudinal diffusion coeffients (relative to E-field)
	swarmParam["bulkRedTransvDiffCoeff"] = totalGasDensity * (rotatedAveragedBulkDiffusionCoeffs(0,0) + rotatedAveragedBulkDiffusionCoeffs(1,1))/2.0;
	swarmParam["bulkRedTransvDiffCoeffError"] = totalGasDensity * (rotatedAveragedBulkDiffusionCoeffsError(0,0) + rotatedAveragedBulkDiffusionCoeffsError(1,1))/2.0;
	swarmParam["bulkRedLongDiffCoeff"] = totalGasDensity * rotatedAveragedBulkDiffusionCoeffs(2,2);
	swarmParam["bulkRedLongDiffCoeffError"] = totalGasDensity * rotatedAveragedBulkDiffusionCoeffsError(2,2);

	// total ionization and attachment rate-coefficients
	swarmParam["totalIonRateCoeff"] = 0;
	swarmParam["totalAttRateCoeff"] = 0;
	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){
			if (collision->type == "Ionization"){
				swarmParam["totalIonRateCoeff"] += collision->target->density * collision->ineRateCoeff;
			}
			else if (collision->type == "Attachment"){
				swarmParam["totalAttRateCoeff"] += collision->target->density * collision->ineRateCoeff;
			}
		}
	}

	if (excitationFrequency == 0){
		// reduced mobility
		double bulkVelocAlongEField = std::abs(rotatedAveragedBulkDriftVelocity[2]);
		swarmParam["bulkRedMobCoeff"] = bulkVelocAlongEField/reducedElecFieldSI;
		swarmParam["bulkRedMobCoeffError"] = std::abs(rotatedAveragedBulkDriftVelocityError[2])/reducedElecFieldSI;

		// reduced Townsend coefficient and reduced attachment coefficient
		swarmParam["bulkRedTownsendCoeff"] = swarmParam["totalIonRateCoeff"]/bulkVelocAlongEField;
		swarmParam["bulkRedAttCoeff"] = swarmParam["totalAttRateCoeff"]/bulkVelocAlongEField;

		// characteristic energy
		swarmParam["bulkCharacEnergy"] = swarmParam["bulkRedTransvDiffCoeff"] / swarmParam["bulkRedMobCoeff"];
		swarmParam["bulkCharacEnergyError"] = swarmParam["bulkRedTransvDiffCoeffError"]/swarmParam["bulkRedMobCoeff"] +
											swarmParam["bulkRedMobCoeffError"]*swarmParam["bulkRedTransvDiffCoeff"]/std::pow(swarmParam["bulkRedMobCoeff"], 2);

		// ----- evaluate effective SST parameters ---- //

		// effective SST townsend ionization coefficient
		swarmParam["effSSTAverageVelocity"] = 0.5*bulkVelocAlongEField + std::sqrt(0.25*bulkVelocAlongEField*bulkVelocAlongEField - totalGasDensity*rotatedAveragedBulkDiffusionCoeffs(2,2)*(swarmParam["totalIonRateCoeff"]-swarmParam["totalAttRateCoeff"]));
		swarmParam["effSSTRedTownsendCoeff"] = swarmParam["totalIonRateCoeff"]/swarmParam["effSSTAverageVelocity"];
		swarmParam["effSSTRedAttCoeff"] = swarmParam["totalAttRateCoeff"]/swarmParam["effSSTAverageVelocity"];
	}

	// ----- evaluate energy parameters ----- //

	// mean energy
	swarmParam["meanEnergy"] = averagedMeanEnergy;
	swarmParam["meanEnergyError"] = averagedMeanEnergyError;

	// electron temperature
	swarmParam["Te"] = 2.0/3.0 * swarmParam["meanEnergy"];
	swarmParam["TeError"] = 2.0/3.0 * swarmParam["meanEnergyError"];
	workCond->update({"electronTemperature"}, {swarmParam["Te"]});

	// ----- evaluate parameters from the EEDF, using two-term expressions ----- //

	double factor = std::sqrt(2.0*Constant::electronCharge/Constant::electronMass)/3.0;
	double energyStep = energyGrid->step;
	Eigen::ArrayXd energyNode = energyGrid->node;
	Eigen::ArrayXd energyCell = energyGrid->cell;
	int Nnodes = energyNode.size();
	int Ncells = energyCell.size();

	// evaluate total cross sections, which are necessary for the calculations
	totalMomTransfCrossSection = Eigen::ArrayXd::Zero(Nnodes);
	Eigen::ArrayXd totalEnergyRelaxCrossSection = Eigen::ArrayXd::Zero(Nnodes);
	Eigen::ArrayXd totalEnergyRelaxCrossSectionClassic = Eigen::ArrayXd::Zero(Nnodes);
	for (auto& gas: gasArray){
		if (gas->collisionArray.empty()){
			continue;
		}
		// evaluation of the mass ratio
		double massRatio = Constant::electronMass/gas->mass;		
		// loop over each collision with the gas
		for (auto& collision: gas->collisionArray){
			// avoid effective collisions
			if (collision->type == "Effective"){
				continue;
			}
			// add collision cross section to the total cross section (also superelastic)
			totalMomTransfCrossSection += collision->momTransfCrossSection * collision->target->density;
			if (collision->type == "Elastic"){
				totalEnergyRelaxCrossSection += 2.0*massRatio*collision->integralCrossSection * collision->target->density;
			}
			else{
				Eigen::ArrayXd weightedICS = collision->integralCrossSection * collision->target->density;
				Eigen::ArrayXd aux = weightedICS * collision->threshold / energyNode;
				aux[0] = 0;
				totalEnergyRelaxCrossSection += aux;
				totalEnergyRelaxCrossSectionClassic += weightedICS;
			}
			if (collision->isReverse){
				EedfState* prod = collision->productArray[0];
				Eigen::ArrayXd weightedICS = collision->superElasticCrossSection("integral", Eigen::ArrayXd(0)) * prod->density;
				Eigen::ArrayXd aux =  weightedICS * collision->threshold / energyNode;
				aux[0] = 0;				
				totalMomTransfCrossSection += collision->superElasticCrossSection("momTransf", Eigen::ArrayXd(0)) * prod->density;				
				totalEnergyRelaxCrossSection += aux;
				totalEnergyRelaxCrossSectionClassic += weightedICS;
			}
		}
	}

	// momentum-transfer frequency (excluding the growth contribution)
	Eigen::ArrayXd cellMTCS = (totalMomTransfCrossSection.segment(1,Ncells) + totalMomTransfCrossSection.segment(0,Ncells))*0.5;
	swarmParam["momTransfFreq"] = totalGasDensity*3.0*factor*energyStep * (eedf * energyCell * cellMTCS).sum();

	// energy-relaxation frequency (excluding the growth contribution)
	Eigen::ArrayXd cellERCS = (totalEnergyRelaxCrossSection.segment(1,Ncells) + totalEnergyRelaxCrossSection.segment(0,Ncells))*0.5;
	swarmParam["energyRelaxFreq"] = totalGasDensity*3.0*factor*energyStep * (eedf * energyCell * cellERCS).sum();

	// classical formula for energy-relaxation frequency (excluding the growth contribution)
	Eigen::ArrayXd cellERCSClassic = (totalEnergyRelaxCrossSectionClassic.segment(1,Ncells) + totalEnergyRelaxCrossSectionClassic.segment(0,Ncells))*0.5;
	swarmParam["energyRelaxFreqClassic"] = totalGasDensity*3.0*factor*energyStep * (eedf * energyCell * cellERCSClassic).sum();		
	
	// add the growth contribution
	Eigen::ArrayXd totalMTCrossSectionAux = totalMomTransfCrossSection;
	totalMTCrossSectionAux.segment(1,Nnodes-1) += 1.0/(3.0*factor*Eigen::sqrt(energyNode.segment(1,Nnodes-1)))*(swarmParam["totalIonRateCoeff"]-swarmParam["totalAttRateCoeff"]);
	// convert to frequency
	Eigen::ArrayXd totalMTFrequencyAux = totalGasDensity*totalMTCrossSectionAux*3.0*factor*Eigen::sqrt(energyNode);

	// reduced energy diffusion 
	swarmParam["redDiffCoeffEnergy_eedf"] = 2.0*factor*energyStep*(energyCell.square()*eedf/(totalMTCrossSectionAux.segment(0,Nnodes-1)+totalMTCrossSectionAux.segment(1,Nnodes-1))).sum();

	// reduced energy mobility 
	swarmParam["redMobCoeffEnergy_eedf"] = -factor*(energyNode.segment(1,Nnodes-2).square()*(eedf.segment(1,Nnodes-2)-eedf.segment(0,Nnodes-2))/totalMTCrossSectionAux.segment(1,Nnodes-2)).sum();

	// reduced diffusion
	swarmParam["redDiffCoeff_eedf"] = 2.0*factor*energyStep*(energyCell*eedf/(totalMTCrossSectionAux.segment(0,Nnodes-1)+totalMTCrossSectionAux.segment(1,Nnodes-1))).sum();

	// DC non-magnetized reduced mobility
	swarmParam["redMobCoeff_DC_eedf"] = -factor*(energyNode.segment(1,Nnodes-2)*(eedf.segment(1,Nnodes-2)-eedf.segment(0,Nnodes-2))/totalMTCrossSectionAux.segment(1,Nnodes-2)).sum();

	// characteristic energy
	swarmParam["characEnergy_eedf"] = swarmParam["redDiffCoeff_eedf"]/swarmParam["redMobCoeff_DC_eedf"];

	// reduced mobility matrix
	Eigen::ArrayXd aux1 = -2.0/3.0*Constant::electronCharge/Constant::electronMass*Eigen::pow(energyNode.segment(1,Nnodes-2),1.5)*((eedf.segment(1,Nnodes-2)-eedf.segment(0,Nnodes-2)));
	Eigen::ArrayXd totalMTFrequencyAuxSegm = totalMTFrequencyAux.segment(1,Nnodes-2);
	Eigen::ArrayXd aux2 = (totalMTFrequencyAuxSegm.square() + std::pow(excitationFrequencyRadians-cyclotronFrequency, 2))*(totalMTFrequencyAuxSegm.square() + std::pow(excitationFrequencyRadians+cyclotronFrequency, 2));
	Eigen::ArrayXd aux3 = totalMTFrequencyAuxSegm.square()+excitationFrequencyRadians*excitationFrequencyRadians;

	swarmParam["redMobCoeff_xx_real_eedf"] = totalGasDensity*(totalMTFrequencyAuxSegm*(totalMTFrequencyAuxSegm.square()+excitationFrequencyRadians*excitationFrequencyRadians+cyclotronFrequency*cyclotronFrequency)/aux2*aux1).sum();
	swarmParam["redMobCoeff_xx_imag_eedf"] = -totalGasDensity*(excitationFrequencyRadians*(totalMTFrequencyAuxSegm.square()+excitationFrequencyRadians*excitationFrequencyRadians-cyclotronFrequency*cyclotronFrequency)/aux2*aux1).sum();
	swarmParam["redMobCoeff_xy_real_eedf"] = -totalGasDensity*cyclotronFrequency*((totalMTFrequencyAuxSegm.square()-excitationFrequencyRadians*excitationFrequencyRadians+cyclotronFrequency*cyclotronFrequency)/aux2*aux1).sum();
	swarmParam["redMobCoeff_xy_imag_eedf"] = totalGasDensity*2.0*excitationFrequencyRadians*cyclotronFrequency*(totalMTFrequencyAuxSegm/aux2*aux1).sum();
	swarmParam["redMobCoeff_zz_real_eedf"] = totalGasDensity*(totalMTFrequencyAuxSegm/aux3*aux1).sum();
	swarmParam["redMobCoeff_zz_imag_eedf"] = -totalGasDensity*(excitationFrequencyRadians/aux3*aux1).sum();
}

void BoltzmannMC::evaluateRateCoeff(){
	// clear rate coeff vectors
	rateCoeffAll.clear();
	rateCoeffExtra.clear();
	rateCoeffAll_periodic.clear();
	rateCoeffExtra_periodic.clear();

	// evaluate rate coefficients for all MC processes
	for (int i = 0; i < nProcesses; ++i){
		GeneralDefinitions::RateCoeffStruct  auxRateCoeff;
		// if it is inelastic
		if (!isSuperElastic[i]){
			// get the correspondent 'Collision' object
			Collision* realCollision = realCollisionPointers[i];
			// save the ID of the 'Collision' object
			auxRateCoeff.collID = realCollision->ID;
			// save the typical rate coefficients, from the convolution of the eedf with the cross section
			realCollision->evaluateRateCoeff(eedf);
			auxRateCoeff.ineRate = realCollision->ineRateCoeff;
			auxRateCoeff.supRate = realCollision->supRateCoeff;
			// save the MC inelastic rate coefficient 
			auxRateCoeff.ineRateMC = averagedRateCoeffs[i];
			// check if the 'Collision' object has a superelastic
			if (auxRateCoeff.supRate != Constant::NON_DEF){
				// increment i, since the corresponding MC process is immediately after
				++i;
				// save the MC superelastic rate coefficient 
				auxRateCoeff.supRateMC = averagedRateCoeffs[i];				
			}
			// save the collision description
			auxRateCoeff.collDescription = realCollision->description();
		}
		rateCoeffAll.push_back(auxRateCoeff);
	}

	// evaluate rate coefficients of effective collisions, if they exist.
	// evaluate rate coefficients for all 'extra' 'Collision' objects
	for (auto& gas: gasArray){
		GeneralDefinitions::RateCoeffStruct  auxRateCoeff;
		for (auto& collision: gas->collisionArray){
			if (collision->type == "Effective"){
				auxRateCoeff.collID = collision->ID;
				collision->evaluateRateCoeff(eedf);
				auxRateCoeff.ineRate = collision->ineRateCoeff;
				auxRateCoeff.supRate = collision->supRateCoeff;
				auxRateCoeff.collDescription = collision->description();
				rateCoeffAll.push_back(auxRateCoeff);				
			}
		}
		for (auto& collision: gas->collisionArrayExtra){
			auxRateCoeff.collID = collision->ID;
			collision->evaluateRateCoeff(eedf);
			auxRateCoeff.ineRate = collision->ineRateCoeff;
			auxRateCoeff.supRate = collision->supRateCoeff;
			auxRateCoeff.collDescription = collision->description();
			rateCoeffExtra.push_back(auxRateCoeff);
		}
	}

	// evaluate rate-coefficients for along the different phases of the period
	if (excitationFrequencyRadians != 0){
		for (int phaseIdx = 0; phaseIdx < nIntegrationPhases; ++phaseIdx){
			std::vector<GeneralDefinitions::RateCoeffStruct> rateCoeffAll_onePhase;
			std::vector<GeneralDefinitions::RateCoeffStruct> rateCoeffExtra_onePhase;
			Eigen::ArrayXd eedf_onePhase = eedf_periodic.row(phaseIdx);			
			// evaluate rate coefficients for all MC processes
			for (int i = 0; i < nProcesses; ++i){
				GeneralDefinitions::RateCoeffStruct  auxRateCoeff;
				// if it is inelastic
				if (!isSuperElastic[i]){
					// get the correspondent 'Collision' object
					Collision* realCollision = realCollisionPointers[i];
					// save the ID of the 'Collision' object
					auxRateCoeff.collID = realCollision->ID;
					// save the typical rate coefficients, from the convolution of the eedf with the cross section
					realCollision->evaluateRateCoeff(eedf_onePhase);
					auxRateCoeff.ineRate = realCollision->ineRateCoeff;
					auxRateCoeff.supRate = realCollision->supRateCoeff;
					// check if the 'Collision' object has a superelastic
					if (auxRateCoeff.supRate != Constant::NON_DEF){
						// increment i, since the corresponding MC process is immediately after
						++i;			
					}
					// save the collision description
					auxRateCoeff.collDescription = realCollision->description();
				}
				rateCoeffAll_onePhase.push_back(auxRateCoeff);
			}

			// evaluate rate coefficients of effective collisions, if they exist.
			// evaluate rate coefficients for all 'extra' 'Collision' objects
			for (auto& gas: gasArray){
				GeneralDefinitions::RateCoeffStruct  auxRateCoeff;
				for (auto& collision: gas->collisionArray){
					if (collision->type == "Effective"){
						auxRateCoeff.collID = collision->ID;
						collision->evaluateRateCoeff(eedf_onePhase);
						auxRateCoeff.ineRate = collision->ineRateCoeff;
						auxRateCoeff.supRate = collision->supRateCoeff;
						auxRateCoeff.collDescription = collision->description();
						rateCoeffAll_onePhase.push_back(auxRateCoeff);				
					}
				}
				for (auto& collision: gas->collisionArrayExtra){
					auxRateCoeff.collID = collision->ID;
					collision->evaluateRateCoeff(eedf_onePhase);
					auxRateCoeff.ineRate = collision->ineRateCoeff;
					auxRateCoeff.supRate = collision->supRateCoeff;
					auxRateCoeff.collDescription = collision->description();
					rateCoeffExtra_onePhase.push_back(auxRateCoeff);
				}
			}
			rateCoeffAll_periodic.push_back(rateCoeffAll_onePhase);
			rateCoeffExtra_periodic.push_back(rateCoeffExtra_onePhase);		
		}
	}		
}

