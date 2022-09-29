#include "LoKI-MC/Headers/BoltzmannMC.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/PrescribedEedf.h"
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

#define bold "\e[1m"
#define nonBold "\e[0m"
#define blue "\033[0;34m"
#define resetColor "\033[0m"

// definition of the process types
const int conservativeType = 0;
const int ionizationType = 1;
const int attachmentType = 2;

// the constructor is written in the '.h' file to avoid compilation problems related with the templates

void BoltzmannMC::allocateEvaluateVariablesFirstTime(){
	// 'allocateEvaluateVariablesFirstTime' evaluate and allocate variables when the constructor is called

	// evaluate nProcesses
	nProcesses = 0;
	for (auto& gas: gasArray){
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
	interpolCrossSectionsXrelDens = new double*[nProcesses];
	cumulSumInterpolCrossSectionsXrelDens = new double*[nProcesses];
	for (int i = 0; i < nProcesses; ++i){
		interpolCrossSectionsXrelDens[i] = new double[interpolCrossSectionSize];
		cumulSumInterpolCrossSectionsXrelDens[i] = new double[interpolCrossSectionSize];
	}
	totalCollisionFrequencies = new double[interpolCrossSectionSize];
	maxCollisionFrequencies = new double[interpolCrossSectionSize];
	energyGainProcesses = new double[nProcesses]; energyLossProcesses = new double[nProcesses];
	collisionCounters = new unsigned long long int[nProcesses];
	averagedRateCoeffs = new double[nProcesses];
	averagedPowerGainProcesses = new double[nProcesses]; averagedPowerLossProcesses = new double[nProcesses];
	wParameters = new double[nProcesses];
	targetGasIDs = new int[nProcesses];
	chosenProcessIDs = new int[(int)nElectrons];
	ejectedElectronPositions.resize(nElectrons,3); ejectedElectronVelocities.resize(nElectrons,3);
	electronEnergyChanges = new double [(int)nElectrons];

	// assign the data for each MC process (remember that here we separate inelastics from superelastics)
	int iterProcess = 0;

	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){

			// assign process type
			if (collision->type == "Effective"){ // avoid effective collisions
				continue;
			}
			else if (collision->type == "Ionization"){
				processTypes[iterProcess] = ionizationType;
				isIonization[iterProcess] = true;
			}
			else if (collision->type == "Attachment"){
				processTypes[iterProcess] = attachmentType;
				isIonization[iterProcess] = false;
			}
			else{
				processTypes[iterProcess] = conservativeType;
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
			std::vector<double> tempEnergyVector = collision->rawCrossSection[0];
			std::vector<double> tempValueVector = collision->rawCrossSection[1];
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
				energyMaxElastic = fmin(energyMaxElastic, energyMaxLimits[iterProcess]);
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

			// increment the iterator
			++iterProcess;

			// enter here if there is a superelastic
			if (collision->isReverse){
				// assign process type
				processTypes[iterProcess] = conservativeType;

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

				++iterProcess;
			}
		}
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
	while ((collisionCounterAfterSS < maxCollisionsAfterSteadyState && (!goodStatisticalErrors || !errorsToBeChecked) 
		   && nIntegrationPoints < requiredIntegrationPoints && totalIntegratedTime/steadyStateTime < requiredIntegratedSSTimes) || nIntegrationPoints < 500){

		// calculate the time step and perform the accelerated motion of the electrons
		freeFlight();
 						
 		// update the number of sampling points
		++nSamplingPoints;
		currentSamplingIndex = nSamplingPoints-1;

		// update the sampling times
		samplingTimes[currentSamplingIndex] = time;

		// calculate average data for swarm parameters
		calculateMeanDataForSwarmParams();

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

		// select and perform collisions for each electron in the ensemble
		performCollisions();

		if (steadyStateTime != Constant::NON_DEF && nIntegrationPoints > 200 && nIntegrationPoints % nPointsBetweenStatErrorsCheck == 0){
			checkStatisticalErrors();
		}

		if (dispMCStatus && nSamplingPoints % nPointsBetweenDispInfo == 0){
			dispInfo();
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

	if (collisionCounterAfterSS > maxCollisionsAfterSteadyState){
		if (dispMCStatus){
			std::cout<<" "<<std::endl;
			for (int i = 0; i < 21; ++i){
				// move up in the terminal
				std::printf("%c[1A", 0x1B);
				// clear terminal line
				std::printf("%c[2K", 0x1B);
			}
			Message::warning("Monte Carlo simulation ended after reaching ''maxCollisionsAfterSteadyState'' indicated in the setup file. [E/N = " + std::to_string(reducedElecField) + " Td\n");
			for (int i = 0; i < 20; ++i){
				std::printf("\n");
			}
		}
		else{
			Message::warning("Monte Carlo simulation ended after reaching ''maxCollisionsAfterSteadyState'' indicated in the setup file. [E/N = " + std::to_string(reducedElecField) + " Td\n");
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

	// update the electric field
	reducedElecField = workCond->reducedElecField;
	electricField[0] = 0;
	electricField[1] = 0;
	electricField[2] = -workCond->reducedElecFieldSI*totalGasDensity;;

	// calculate the electron acceleration, ONLY due to the electric field
	accelerationElecField = -Constant::electronCharge/Constant::electronMass * electricField;

	// initialize to zero the position of the electrons
	electronPositions = Eigen::ArrayXXd::Zero(nElectrons, 3);

	// initialize the velocities, using a maxwell-boltzmann distribution at gas temperature
	electronVelocities.resize(nElectrons, 3);
	double electronThermalDeviation = std::sqrt(Constant::boltzmann*gasTemperature/Constant::electronMass);
	for (int i = 0; i < nElectrons; ++i){
		electronVelocities(i,0) = MathFunctions::unitNormalRand() * electronThermalDeviation; 
		electronVelocities(i,1) = MathFunctions::unitNormalRand() * electronThermalDeviation;
		electronVelocities(i,2) = MathFunctions::unitNormalRand() * electronThermalDeviation;
	}

	// calculate the corresponding energies (in eV)
	electronEnergies = 0.5*Constant::electronMass*electronVelocities.matrix().rowwise().squaredNorm() / Constant::electronCharge;

	// interpolate the cross sections for the first time, until two times the maximum electron energy
	// the maximum collision frequency for each energy point is also calculated
	interpolateCrossSections(2.0*electronEnergies.maxCoeff());

	// initialize the trial collision frequency
	trialCollisionFrequency = maxCollisionFrequencies[interpolCrossSectionSize-1];

	// initialize the matrices with the mean values along time (when the size surpasses '100', we will use 'conservativeResize')
	samplingTimes = Eigen::ArrayXd::Zero(100);
	meanEnergies = Eigen::ArrayXd::Zero(100); meanPositions = Eigen::ArrayXXd::Zero(100,3); meanVelocities = Eigen::ArrayXXd::Zero(100,3);
	fluxDiffusionCoeffs = Eigen::ArrayXXd::Zero(100,9); positionCovariances = Eigen::ArrayXXd::Zero(100,9);

	// initialize the averagedMeanEnergy to Constant::NON_DEF. Important for the 'dispInfo' function
	averagedMeanEnergy = Constant::NON_DEF;
	maxElecEnergy = electronEnergies.maxCoeff();

	// initialize the time-averaged diffusion coeffs to Constant::NON_DEF
	averagedFluxDiffusionCoeffs = Eigen::ArrayXd::Constant(9, Constant::NON_DEF); averagedFluxDiffusionCoeffsError = Eigen::ArrayXd::Constant(9, Constant::NON_DEF);
	averagedBulkDiffusionCoeffs = Eigen::ArrayXd::Constant(9, Constant::NON_DEF); averagedBulkDiffusionCoeffsError = Eigen::ArrayXd::Constant(9, Constant::NON_DEF);

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
			interpolCrossSectionsXrelDens[iterProcess][i] = value;
			currentFrequency += value;
			cumulSumInterpolCrossSectionsXrelDens[iterProcess][i] = currentFrequency;
		}

		currentFrequency *= totalGasDensity*std::sqrt(energy*2.0*Constant::electronCharge/Constant::electronMass);
		totalCollisionFrequencies[i] = currentFrequency;
		maxFrequency = std::fmax(currentFrequency, maxFrequency);
		maxCollisionFrequencies[i] = maxFrequency;
	}
}

void BoltzmannMC::freeFlight(){
	// 'freeFlight' calculates the random time step and performs the accelerated motion of the electrons

	// ----- save the previous energies (in eV) ----- //
	Eigen::ArrayXd previousEnergies = 0.5*Constant::electronMass*electronVelocities.matrix().rowwise().squaredNorm() / Constant::electronCharge;

	
	// ----- calculate the first time-step ----- //
	// random = 0 or 1 are removed to avoid infinite or null time-steps
	deltaT = -std::log(MathFunctions::unitUniformRand(false,false)) / trialCollisionFrequency;

	// ----- assure that the trial collision frequency is higher than the maximum collision frequency possible -----//

	// get the maximum speed before acceleration 
	double maxSpeedBeforeAccel = electronVelocities.matrix().rowwise().norm().maxCoeff();
	double accelerationNorm = accelerationElecField.matrix().norm();

	// worst-case scenario: the electron with maximum speed has a velocity alligned with the acceleration
	double maxSpeed = maxSpeedBeforeAccel + accelerationNorm*deltaT;
	// convert the maximum speed possible (after acceleration) to energy in eV
	// avoid very small energies which would lead to instabilities in the code
	double maxEnergy = std::fmax(0.5*Constant::electronMass*maxSpeed*maxSpeed/Constant::electronCharge, 0.01);

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
	int energyIndex = std::fmin(std::floor(maxEnergy/crossSectionEnergyStep), interpolCrossSectionSize-1);
	double maxCollisionFrequency = maxCollisionFrequencies[energyIndex];

	while (maxCollisionFrequency > trialCollisionFrequency){

		// increment the trial collision frequency and recalculate the time-step (random = 0 or 1 are removed to avoid infinite or null time-steps)
		trialCollisionFrequency *= 1.1;
		deltaT = -std::log(MathFunctions::unitUniformRand(false,false)) / trialCollisionFrequency;

		// worst-case scenario: the electron with maximum speed has a velocity alligned with the acceleration
		maxSpeed = maxSpeedBeforeAccel + accelerationNorm*deltaT;
		// convert the maximum speed possible (after acceleration) to energy in eV
		// avoid very small energies which would lead to instabilities in the code
		maxEnergy = std::fmax(0.5*Constant::electronMass*maxSpeed*maxSpeed/Constant::electronCharge, 0.01);

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
		energyIndex = std::fmin(std::floor(maxEnergy/crossSectionEnergyStep), interpolCrossSectionSize-1);
		maxCollisionFrequency = maxCollisionFrequencies[energyIndex];
	}

	// ----- solve the accelerated motion of the electrons -----//

	// accelerate the electrons
	electronPositions += electronVelocities*deltaT;
	electronPositions.rowwise() += (accelerationElecField*0.5*deltaT*deltaT).transpose();
	electronVelocities.rowwise() += (accelerationElecField*deltaT).transpose();

	// save the new time
	time += deltaT;

	// update the energy of the electrons (in eV)
	electronEnergies = 0.5*Constant::electronMass*electronVelocities.matrix().rowwise().squaredNorm() / Constant::electronCharge;
	
	// add the energy gain during the acceleration (in eV)
	energyGainField += (electronEnergies - previousEnergies).sum();
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
}

void BoltzmannMC::getTimeAverageEnergyParams(){
	// 'getTimeAverageEnergyParams' calculates the time-averaged mean energy

	// calculate the time-averaged mean energy
	averagedMeanEnergy = meanEnergies.segment(firstIntegrationIndex,nIntegrationPoints).mean();
	averagedMeanEnergyError = MathFunctions::statisticalError(meanEnergies.segment(firstIntegrationIndex,nIntegrationPoints), 50);
}

void BoltzmannMC::getTimeDependDistributions(){
	// 'getTimeDependDistributions' calculates the time-dependent quantities necessary to get the time-averaged distribution functions

	// create the current electron energy histogram (eeh) and sum it to the previous ones 
	eehSum += MathFunctions::histogramCount(electronEnergies, eedfEnergyNodes);	

	// calculate the cos(angles) of the velocities relatively to the z axis
	electronCosAngles = electronVelocities.col(2) / electronVelocities.matrix().rowwise().norm().array();

	// create the current electron angular histogram (eah) and sum it to the previous ones
	eahSum += MathFunctions::histogram2DCount(electronEnergies, electronCosAngles, eedfEnergyNodes, cosAngleNodes);

	// create the current electron velocity histogram (evh) and sum it to the previous ones
	evhSum += MathFunctions::histogram2DCount(Eigen::sqrt(electronVelocities.col(0).square() + electronVelocities.col(1).square()), electronVelocities.col(2), radialVelocityNodes, axialVelocityNodes);
}

void BoltzmannMC::getTimeAverageDistributions(){
	// 'getTimeAverageDistributions' calculates the time-averaged distribution functions (performed only at end of the simulation, no impact in performance)

	double normalizer = eehSum.sum();

	// get the electron energy distribution function (eedf)
	eedf = eehSum/(Eigen::sqrt(eedfEnergyCells)*normalizer*eedfEnergyStep);

	// get the electron angular distribution function (eadf)
	// Note: I had to multiply the eadf by two to obtain the correct values in the anisotropies (it is probably related with the contribution of the azymuthal angle)
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

void BoltzmannMC::getTimeAverageFluxParams(){
	// 'getTimeAverageFluxParams' calculates the time-average of the flux drift velocity and the flux diffusion coefficients

	// calculate the time-averaged flux drift velocity
	averagedFluxDriftVelocity = meanVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3).colwise().mean();
	averagedFluxDriftVelocityError = MathFunctions::statisticalErrorColwise(meanVelocities.block(firstIntegrationIndex,0,nIntegrationPoints,3), 50);

	// calculate the time-averaged flux diffusion coeffs (9 components)
	averagedFluxDiffusionCoeffs = fluxDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9).colwise().mean();
	averagedFluxDiffusionCoeffsError = MathFunctions::statisticalErrorColwise(fluxDiffusionCoeffs.block(firstIntegrationIndex,0,nIntegrationPoints,9), 50);
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
	// 'getTimeAverageBulkParams' calculates the bulk swarm parameters using linear regressions of the swarm's center of mass and of the position covariances (performed only at end of the simulation, no impact in performance)

	double c0, c1, cov00, cov01, cov11, sumsq;

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

	for (int i = 0; i < 9; ++i){
		gsl_fit_linear(integrationTimes.data(), 1, halfCovariances.col(i).data(), 1, nIntegrationPoints, &c0, &c1, &cov00, &cov01 , &cov11, &sumsq);
		averagedBulkDiffusionCoeffs[i] = c1;
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
	averagedBulkDiffusionCoeffsError = Eigen::sqrt((squaredSum/(double)nSlices - mean.square())/nSlices);
}

void BoltzmannMC::performCollisions(){
	// 'performCollisions' selects the process that each electron suffers and solves the associated dynamics

	// for each electron
	#pragma omp parallel for
	for (int elecID = 0; elecID < (int)nElectrons; ++elecID){

		// get the temporary values of the position and the velocity
		Eigen::Array3d electronPosition = electronPositions.row(elecID);
		Eigen::Array3d electronVelocity = electronVelocities.row(elecID);
		Eigen::Array3d targetVelocity = Eigen::Array3d::Zero();

		// save the energy before the collision (in eV)
		double incidentEnergy = electronEnergies[elecID];

		// ----- select the process ----- // 

		int chosenProcessID = -1;
		if (gasTemperatureEffect == 1 || (gasTemperatureEffect == 2 && incidentEnergy < 20.0*gasEnergy)){
			double sumCrossSectionFactor = 0;
			double randCrossSectionFactor = trialCollisionFrequency * MathFunctions::unitUniformRand(false,true) / totalGasDensity;
			double incidentSpeed = electronVelocity.matrix().norm(), relativeSpeed;
			int incEnergyIndex = std::fmin(std::floor(incidentEnergy/crossSectionEnergyStep), interpolCrossSectionSize-1), relEnergyIndex;
			int previousTargetGasID = -1;
			for (int iterProcess = 0; iterProcess < nProcesses; ++iterProcess){
				// calculate the target random velocity and the relative energy index only if the gas differs from the previous iteration
				if (targetGasIDs[iterProcess] != previousTargetGasID){
					previousTargetGasID = targetGasIDs[iterProcess];
					targetVelocity[0] = MathFunctions::unitNormalRand(); targetVelocity[1] = MathFunctions::unitNormalRand(); targetVelocity[2] = MathFunctions::unitNormalRand();
					targetVelocity *= thermalStdDeviations[iterProcess];
					// calculate the relative speed and the respective energy index of the grid
					relativeSpeed = (electronVelocity-targetVelocity).matrix().norm();
					relEnergyIndex = std::fmin(std::floor(0.5*reducedMasses[iterProcess]*relativeSpeed*relativeSpeed/Constant::electronCharge/crossSectionEnergyStep), 
											  interpolCrossSectionSize-1);
				}

				if (isIonization[iterProcess]){
					sumCrossSectionFactor += interpolCrossSectionsXrelDens[iterProcess][incEnergyIndex] * incidentSpeed;
				}
				else{
					sumCrossSectionFactor += interpolCrossSectionsXrelDens[iterProcess][relEnergyIndex] * relativeSpeed;
				}

				// check if the cumulative factor surpassed the random factor, i.e., if this process is to be chosen
				if (sumCrossSectionFactor >= randCrossSectionFactor){
					chosenProcessID = iterProcess;
					chosenProcessIDs[elecID] = chosenProcessID;
					break;
				}
			}
			// check if this is a null collision
			if (chosenProcessID == -1){
				chosenProcessIDs[elecID] = -1;
				continue;
			}
		}
		else{
			double randCollisionFrequency = trialCollisionFrequency * MathFunctions::unitUniformRand(false,true);
			int incEnergyIndex = std::fmin(std::floor(incidentEnergy/crossSectionEnergyStep), interpolCrossSectionSize-1);
			// check if this is a null collision
			if (randCollisionFrequency > totalCollisionFrequencies[incEnergyIndex]){
				chosenProcessIDs[elecID] = -1;
				continue;				
			}

			// use a bissection method to find the process ID. Only possible when the gas-temperature effect is not considered
			double randCrossSectionFactor = randCollisionFrequency / totalGasDensity / electronVelocity.matrix().norm();
			int leftLimit = 0;
			int rightLimit = nProcesses-1;

			while(leftLimit != rightLimit){

				// get the tentative index. Note that in this operation, the result is rounded downwards (see arithmetic rules of integers in c++)
				int tentativeIndex = (leftLimit + rightLimit)/2;

				// if the cumulative value of the tentative index is greater than the random value, the solution is in the left side
				double tentativeValue = cumulSumInterpolCrossSectionsXrelDens[tentativeIndex][incEnergyIndex];
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
			while (interpolCrossSectionsXrelDens[chosenProcessID][incEnergyIndex] == 0){
				--chosenProcessID;
			}

			// save the chosen process ID
			chosenProcessIDs[elecID] = chosenProcessID;		
		}

		// ----- calculate the velocity after the collision ----- //

		if (processTypes[chosenProcessID] == conservativeType){
			conservativeCollision(elecID, electronPosition, electronVelocity, targetVelocity, incidentEnergy, chosenProcessID);
		}

		else if (processTypes[chosenProcessID] == ionizationType){
			ionizationCollision(elecID, electronPosition, electronVelocity, incidentEnergy, chosenProcessID);
		}

		else if (processTypes[chosenProcessID] == attachmentType){
			attachmentCollision(elecID, incidentEnergy);
		}

		// update the electron position and velocity
		electronPositions.row(elecID) = electronPosition;
		electronVelocities.row(elecID) = electronVelocity;
	}

	// perform all the tasks involving collisions that cannot be parallelized
	nonParallelCollisionTasks();
}

void BoltzmannMC::conservativeCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, Eigen::Array3d &targetVelocity, double incidentEnergy, int chosenProcessID){
	// 'convervativeCollision' solves the dynamics of a conservative collision, i.e., a collision without changes in the electron number
	// Important note: for gasTemperatureEffect = 1, it was verified that the energy of the system [electron + target] is conserved!

	double targetMass = targetMasses[chosenProcessID];

	if (gasTemperatureEffect == 1 || (gasTemperatureEffect == 2 && incidentEnergy < 20.0*gasEnergy)){
		// calculate the relative velocity before the collision and convert it to spherical coordinates
		Eigen::Array3d relVelocityBefore = electronVelocity - targetVelocity;
		double relSpeedBefore, thetaRelVel, phiRelVel;
		MathFunctions::cart2sph(relVelocityBefore, relSpeedBefore, thetaRelVel, phiRelVel);

		// calculate the relative speed after the collision, using energy conservation (note that the energy must be converted from eV to SI)
		// use Euler relations to calculate the relative velocity after the collision in the laboratory, knowing chi and eta in the CM frame
		// the scattering angles are generated assuming isotropic scattering in the CM frame
		Eigen::Array3d relVelocityAfter = std::sqrt(relSpeedBefore*relSpeedBefore-2.0/reducedMasses[chosenProcessID]*energyLosses[chosenProcessID]*Constant::electronCharge)* 
								   MathFunctions::eulerTransformation(std::acos(1.0-2.0*MathFunctions::unitUniformRand(true,true)), 2.0*M_PI*MathFunctions::unitUniformRand(false,true), thetaRelVel, phiRelVel);

		// get the electron velocity after the collision using the conservation of momentum transfer
		electronVelocity = targetMass/(Constant::electronMass+targetMass)*relVelocityAfter +
						   (Constant::electronMass*electronVelocity + targetMass*targetVelocity)/(Constant::electronMass+targetMass);
	}
	else{
		// convert the velocity of the incident electron to spherical coordinates
		double electronSpeed, theta, phi;
		MathFunctions::cart2sph(electronVelocity, electronSpeed, theta, phi);

		double chiScattered = std::acos(1.0-2.0*MathFunctions::unitUniformRand(true,true));

		// calculate the electron speed after the collision, using energy conservation and assuming a target molecule at rest (due to the high mass) [Reid 1979]
		// use Euler relations to determine the velocity after the collision
		electronVelocity = std::sqrt((electronSpeed*electronSpeed-2.0/Constant::electronMass*energyLosses[chosenProcessID]*Constant::electronCharge) * 
							    (1.0-2.0*Constant::electronMass*targetMass/(Constant::electronMass+targetMass)/(Constant::electronMass+targetMass)*(1.0-std::cos(chiScattered))))*
						   MathFunctions::eulerTransformation(chiScattered, 2.0*M_PI*MathFunctions::unitUniformRand(false,true), theta, phi);
	}

	// calculate the energy change (in eV) due to this collision
	electronEnergyChanges[elecID] = 0.5*Constant::electronMass*electronVelocity.matrix().squaredNorm() / Constant::electronCharge - incidentEnergy;
}

void BoltzmannMC::ionizationCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double incidentEnergy, int chosenProcessID){
	// 'ionizationCollision' solves the dynamics of an ionization collision

	// convert the velocity of the incident electron to spherical coordinates
	double electronSpeed, theta, phi;
	MathFunctions::cart2sph(electronVelocity, electronSpeed, theta, phi);

	// calculate the speeds of the scattered and ejected electrons, assuming a background of molecules at rest
	double netEnergy = incidentEnergy - energyLosses[chosenProcessID];
	double electronEnergyEjected;
	if (usingSDCS){
		electronEnergyEjected = wParameters[chosenProcessID]*std::tan(MathFunctions::unitUniformRand(true,true) * std::atan(netEnergy/(2.0*wParameters[chosenProcessID])) );
	}
	else{
		electronEnergyEjected = energySharingFactor*netEnergy;
	}
	double electronEnergyScattered = netEnergy - electronEnergyEjected;

	double chiScattered, chiEjected, etaScattered, etaEjected;
	if (isoScatteringIonization){
		chiScattered = std::acos(1.0-2.0*MathFunctions::unitUniformRand(true,true));
		etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		chiEjected = std::acos(1.0-2.0*MathFunctions::unitUniformRand(true,true));
		etaEjected = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
	}
	else{
		// generate the scattering angles, assuming that the incident, scattered and ejected electron velocities are coplanar
		// and that the scattered and ejected electron velocities are perpendicular in the CM frame (according with Boeuf 1982)
		chiScattered = std::acos(std::sqrt(electronEnergyScattered/netEnergy));
		etaScattered = 2.0*M_PI*MathFunctions::unitUniformRand(false,true);
		chiEjected = M_PI/2.0-chiScattered;
		etaEjected = etaScattered+M_PI;
	}

	// use Euler relations to calculate the electron velocity after the collision in the laboratory, knowing chi_scattered and eta in the CM frame (according with Yousfi 1994)
	electronVelocity = std::sqrt(2.0*electronEnergyScattered*Constant::electronCharge/Constant::electronMass) * 
					   MathFunctions::eulerTransformation(chiScattered, etaScattered, theta, phi);

	// use Euler relations to calculate the ejected electron velocity in the laboratory, knowing chi_ejected and eta in the CM frame
	ejectedElectronVelocities.row(elecID) = std::sqrt(2.0*electronEnergyEjected*Constant::electronCharge/Constant::electronMass) *
											MathFunctions::eulerTransformation(chiEjected, etaEjected, theta, phi);

	// create the electron in the same position as the scattered one
	ejectedElectronPositions.row(elecID) = electronPosition;

	// update the energy change (in eV) due to this collision
	electronEnergyChanges[elecID] = -energyLosses[chosenProcessID];
}

void BoltzmannMC::attachmentCollision(int elecID, double incidentEnergy){
	// 'attachmentCollision' solves the dynamics of an attachment collision

	// update the energy change (in eV) due to this collision
	electronEnergyChanges[elecID] = -incidentEnergy;
}

void BoltzmannMC::nonParallelCollisionTasks(){
	// 'nonParallelCollisionTasks' performs all the tasks involving collisions that are not parallelized:
	// increment collision counters, update the energy gain/loss in processes, create ejected electrons, destroy attached electrons

	for (int elecID = 0; elecID < nElectrons; ++elecID){

		int chosenProcessID = chosenProcessIDs[elecID];

		// update the collision counters
		if (chosenProcessID == -1){
			++nullCollisionCounter;
			continue;
		}
		else{
			++collisionCounters[chosenProcessID];
			++totalCollisionCounter;
		}

		// update the energy gain
		if (electronEnergyChanges[elecID] >= 0){
			energyGainProcesses[chosenProcessID] += electronEnergyChanges[elecID];
		}
		else{
			energyLossProcesses[chosenProcessID] += electronEnergyChanges[elecID];
		}

		// if it is an ionization process, create the ejected electron and eliminate a randomly chosen electron
		if (processTypes[chosenProcessID] == ionizationType){
			// select the electron that will be replaced
			int selectedID = std::fmin(std::floor(MathFunctions::unitUniformRand(true,false)*nElectrons), nElectrons-1);
			// calculate the energy change due to the growth profile
			energyGrowth -= 0.5*Constant::electronMass/Constant::electronCharge*electronVelocities.row(selectedID).matrix().squaredNorm();
			// update the position and velocity of the selected electron
			electronPositions.row(selectedID) = ejectedElectronPositions.row(elecID);
			electronVelocities.row(selectedID) = ejectedElectronVelocities.row(elecID);
		}
		// if it is an attachment process, eliminate the attached electron and copy a randomly chosen electron
		else if (processTypes[chosenProcessID] == attachmentType){
			// select the electron that will be copied
			int selectedID = std::fmin(std::floor(MathFunctions::unitUniformRand(true,false)*nElectrons), nElectrons-1);
			// calculate the energy change due to the growth profile
			energyGrowth += 0.5*Constant::electronMass/Constant::electronCharge*electronVelocities.row(selectedID).matrix().squaredNorm();
			// copy the position and velocity to the ID of the 'destroyed' electron
			electronPositions.row(elecID) = electronPositions.row(selectedID);
			electronVelocities.row(elecID) = electronVelocities.row(selectedID);
		}
	}
}

void BoltzmannMC::checkStatisticalErrors(){
	// 'checkStatisticalErrors' checks if the minimum errors specified by the user on the setup are already verified

	getTimeAverageEnergyParams();
	getTimeAverageFluxParams();
	checkPowerBalance();

	double relErrorAbsVelocity = (std::abs(averagedFluxDriftVelocity[0])*averagedFluxDriftVelocityError[0] + std::abs(averagedFluxDriftVelocity[1])*averagedFluxDriftVelocityError[1] + 
								  std::abs(averagedFluxDriftVelocity[2])*averagedFluxDriftVelocityError[2]) / averagedFluxDriftVelocity.matrix().squaredNorm();
	if (averagedMeanEnergyError/averagedMeanEnergy <= requiredMeanEnergyRelError &&
		relErrorAbsVelocity <= requiredFluxDriftVelocityRelError &&
		averagedFluxDiffusionCoeffsError[0]/averagedFluxDiffusionCoeffs[0] <= requiredFluxDiffusionCoeffsRelError &&
		averagedFluxDiffusionCoeffsError[4]/averagedFluxDiffusionCoeffs[4] <= requiredFluxDiffusionCoeffsRelError &&
		averagedFluxDiffusionCoeffsError[8]/averagedFluxDiffusionCoeffs[8] <= requiredFluxDiffusionCoeffsRelError &&
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

	int checkPoints = nSamplingPoints*0.25;
	double averageFirstPart = meanEnergies.segment(nSamplingPoints-1-2*checkPoints, checkPoints).mean();
	double averageSecondPart = meanEnergies.segment(nSamplingPoints-1-checkPoints, checkPoints).mean();
	double relStdSecondPart = MathFunctions::standardDeviation(meanEnergies.segment(nSamplingPoints-1-checkPoints, checkPoints))/std::sqrt(checkPoints-1) / averageSecondPart;

	if ((averageFirstPart >= averageSecondPart && relStdSecondPart < 0.01)  || totalCollisionCounter >= maxCollisionsBeforeSteadyState){

		if (totalCollisionCounter >= maxCollisionsBeforeSteadyState){
			if (dispMCStatus){
				std::cout<<" "<<std::endl;
				for (int i = 0; i < 23; ++i){
					// move up in the terminal
					std::printf("%c[1A", 0x1B);
					// clear terminal line
					std::printf("%c[2K", 0x1B);
				}
				Message::warning("Steady-state imposed by the parameter 'maxCollisionsBeforeSteadyState'. The electron energy may not be stabilized yet. [E/N = " + std::to_string(reducedElecField) + " Td]\n");
				for (int i = 0; i < 22; ++i){
					std::printf("\n");
				}
				dispInfo();
			}
			else{
				Message::warning("Steady-state imposed by the parameter 'maxCollisionsBeforeSteadyState'. The electron energy may not be stabilized yet. [E/N = " + std::to_string(reducedElecField) + " Td]\n");
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
		maxEedfEnergy = 1.75*maxElecEnergy;
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
		for (int i = 0; i < 21 ; ++i){
			std::printf("%s",terminal_moveup);
			std::printf("%s",terminal_clearline);
		}
	}

	std::cout<<"E/N: "<<reducedElecField<<" Td\n\n";
	std::cout<<"Number of real collisions: "<<(double)totalCollisionCounter<<std::endl;
	std::cout<<"Number of null collisions: "<<(double)nullCollisionCounter<<std::endl<<std::endl;
	std::cout<<"Current time: "<<time<<" s\n";
	if (steadyStateTime == Constant::NON_DEF){
		std::cout<<"Steady-state time: "<<"non-defined\n\n";
		std::cout<<"\n\n\n\n\n\n\n\n\n\n\n\n";
	}
	else{
		std::cout<<"Steady-state time: "<<steadyStateTime<<" s\n\n";
		if (nIntegrationPoints >= 3*nPointsBetweenStatErrorsCheck && averagedMeanEnergy != Constant::NON_DEF){
			std::cout<<"Number of integration points: "<<(double)nIntegrationPoints<<"\n\n";
			std::cout<<"Mean energy [eV]: "<<averagedMeanEnergy<<std::endl;
			std::cout<<"Relative error: "<<averagedMeanEnergyError/averagedMeanEnergy<<"\n\n";
			std::cout<<"Flux drift velocity [m/s]: "<<averagedFluxDriftVelocity.transpose()<<std::endl;
			std::cout<<"Relative error: "<<(averagedFluxDriftVelocityError/Eigen::abs(averagedFluxDriftVelocity)).transpose()<<"\n\n";
			std::cout<<"Flux diffusion coefficients [m^2 s^-1]: "<<averagedFluxDiffusionCoeffs[0]<<" "<<averagedFluxDiffusionCoeffs[4]<<" "<<averagedFluxDiffusionCoeffs[8]<<std::endl;
			Eigen::ArrayXd temp = averagedFluxDiffusionCoeffsError/averagedFluxDiffusionCoeffs;
			std::cout<<"Relative error: "<<temp[0]<<" "<<temp[4]<<" "<<temp[8]<<"\n\n";
			std::cout<<"Power balance relative error: "<<powerBalanceRelError<<"\n";
		}
		else{
			std::cout<<"\n\n\n\n\n\n\n\n\n\n\n\n";
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

	// ----- initialize transport parameters structure -----/
	swarmParam["fluxRedTransvDiffCoeff"] = 0; swarmParam["fluxRedTransvDiffCoeffError"] = 0; swarmParam["fluxRedLongDiffCoeff"] = 0; swarmParam["fluxRedLongDiffCoeffError"] = 0;
	swarmParam["fluxRedMobCoeff"] = 0; swarmParam["fluxRedMobCoeffError"] = 0;
	swarmParam["fluxRedTownsendCoeff"] = 0; swarmParam["fluxRedTownsendCoeffError"] = 0;
	swarmParam["fluxRedAttCoeff"] = 0; swarmParam["fluxRedAttCoeffError"] = 0;
	swarmParam["fluxCharacEnergy"] = 0; swarmParam["fluxCharacEnergyError"] = 0;
	swarmParam["bulkRedTransvDiffCoeff"] = 0; swarmParam["bulkRedTransvDiffCoeffError"] = 0; swarmParam["bulkRedLongDiffCoeff"] = 0; swarmParam["bulkRedLongDiffCoeffError"] = 0;
	swarmParam["bulkRedMobCoeff"] = 0; swarmParam["bulkRedMobCoeffError"] = 0;
	swarmParam["bulkRedTownsendCoeff"] = 0; swarmParam["bulkRedTownsendCoeffError"] = 0;
	swarmParam["bulkRedAttCoeff"] = 0; swarmParam["bulkRedAttCoeffError"] = 0;
	swarmParam["bulkCharacEnergy"] = 0; swarmParam["bulkCharacEnergyError"] = 0;
	swarmParam["meanEnergy"] = 0; swarmParam["meanEnergyError"] = 0;
	swarmParam["Te"] = 0; swarmParam["TeError"] = 0;
	swarmParam["totalIonRateCoeff"] = 0; swarmParam["totalAttRateCoeff"] = 0;

	// ----- evaluate flux parameters ----- // 

	// reduced diffusion coeffients
	swarmParam["fluxRedTransvDiffCoeff"] = totalGasDensity * (averagedFluxDiffusionCoeffs[0] + averagedFluxDiffusionCoeffs[4])/2.0;
	swarmParam["fluxRedTransvDiffCoeffError"] = totalGasDensity * (averagedFluxDiffusionCoeffsError[0] + averagedFluxDiffusionCoeffsError[4])/2.0;
	swarmParam["fluxRedLongDiffCoeff"] = totalGasDensity * averagedFluxDiffusionCoeffs[8];
	swarmParam["fluxRedLongDiffCoeffError"] = totalGasDensity * averagedFluxDiffusionCoeffsError[8];


	// total ionization and attachment rate-coefficients
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

	// reduced mobility
	double absFluxVz = std::abs(averagedFluxDriftVelocity[2]);
	swarmParam["fluxRedMobCoeff"] = absFluxVz/reducedElecFieldSI;
	swarmParam["fluxRedMobCoeffError"] = absFluxVz/reducedElecFieldSI;
	

	// reduced Townsend coefficient and reduced attachment coefficient
	swarmParam["fluxRedTownsendCoeff"] = swarmParam["totalIonRateCoeff"]/absFluxVz;
	swarmParam["fluxRedAttCoeff"] = swarmParam["totalAttRateCoeff"]/absFluxVz;

	// characteristic energy
	swarmParam["fluxCharacEnergy"] = swarmParam["fluxRedTransvDiffCoeff"] / swarmParam["fluxRedMobCoeff"];
	swarmParam["fluxCharacEnergyError"] = swarmParam["fluxRedTransvDiffCoeffError"]/swarmParam["fluxRedMobCoeff"] +
										  swarmParam["fluxRedMobCoeffError"]*swarmParam["fluxRedTransvDiffCoeff"]/std::pow(swarmParam["fluxRedMobCoeff"], 2);

	// ----- evaluate bulk parameters ----- // 

	// reduced diffusion coeffients
	swarmParam["bulkRedTransvDiffCoeff"] = totalGasDensity * (averagedBulkDiffusionCoeffs[0] + averagedBulkDiffusionCoeffs[4])/2.0;
	swarmParam["bulkRedTransvDiffCoeffError"] = totalGasDensity * (averagedBulkDiffusionCoeffsError[0] + averagedBulkDiffusionCoeffsError[4])/2.0;
	swarmParam["bulkRedLongDiffCoeff"] = totalGasDensity * averagedBulkDiffusionCoeffs[8];
	swarmParam["bulkRedLongDiffCoeffError"] = totalGasDensity * averagedBulkDiffusionCoeffsError[8];

	// reduced mobility
	double absBulkVz = std::abs(averagedBulkDriftVelocity[2]);
	swarmParam["bulkRedMobCoeff"] = absBulkVz/reducedElecFieldSI;
	swarmParam["bulkRedMobCoeffError"] = absBulkVz/reducedElecFieldSI;

	// reduced Townsend coefficient and reduced attachment coefficient
	swarmParam["bulkRedTownsendCoeff"] = swarmParam["totalIonRateCoeff"]/absBulkVz;
	swarmParam["bulkRedAttCoeff"] = swarmParam["totalAttRateCoeff"]/absBulkVz;

	// characteristic energy
	swarmParam["bulkCharacEnergy"] = swarmParam["bulkRedTransvDiffCoeff"] / swarmParam["bulkRedMobCoeff"];
	swarmParam["bulkCharacEnergyError"] = swarmParam["bulkRedTransvDiffCoeffError"]/swarmParam["bulkRedMobCoeff"] +
										  swarmParam["bulkRedMobCoeffError"]*swarmParam["bulkRedTransvDiffCoeff"]/std::pow(swarmParam["bulkRedMobCoeff"], 2);

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

	// evaluate total cross section, which is necessary to the calculations
	totalCrossSection = Eigen::ArrayXd::Zero(Nnodes);
	for (auto& gas: gasArray){
		if (gas->collisionArray.empty()){
			continue;
		}
		// loop over each collision with the gas
		for (auto& collision: gas->collisionArray){
			// avoid effective collisions
			if (collision->type == "Effective"){
				continue;
			}
			// add collision cross section to the total momentum transfer cross section (also superelastic)
			totalCrossSection += collision->crossSection * collision->target->density;
			if (collision->isReverse){
				totalCrossSection += collision->superElasticCrossSection(Eigen::ArrayXd(0)) * collision->productArray[0]->density;
			}
		}
	}
	// add the growth contribution
	Eigen::ArrayXd totalCrossSectionAux = totalCrossSection;
	totalCrossSectionAux.segment(1,Nnodes-1) += 1.0/(3.0*factor*Eigen::sqrt(energyNode.segment(1,Nnodes-1)))*(swarmParam["totalIonRateCoeff"]-swarmParam["totalAttRateCoeff"]);
	// convert to frequency
	Eigen::ArrayXd totalFrequencyAux = totalGasDensity*totalCrossSectionAux*3.0*factor*Eigen::sqrt(energyNode);

	// reduced energy diffusion 
	swarmParam["redDiffCoeffEnergy_eedf"] = 2.0*factor*energyStep*(energyCell.square()*eedf/(totalCrossSectionAux.segment(0,Nnodes-1)+totalCrossSectionAux.segment(1,Nnodes-1))).sum();

	// reduced energy mobility 
	swarmParam["redMobCoeffEnergy_eedf"] = -factor*(energyNode.segment(1,Nnodes-2).square()*(eedf.segment(1,Nnodes-2)-eedf.segment(0,Nnodes-2))/totalCrossSectionAux.segment(1,Nnodes-2)).sum();

	// reduced diffusion
	swarmParam["redDiffCoeff_eedf"] = 2.0*factor*energyStep*(energyCell*eedf/(totalCrossSectionAux.segment(0,Nnodes-1)+totalCrossSectionAux.segment(1,Nnodes-1))).sum();

	// reduced mobility
	swarmParam["redMobCoeff_DC_eedf"] = -factor*(energyNode.segment(1,Nnodes-2)*(eedf.segment(1,Nnodes-2)-eedf.segment(0,Nnodes-2))/totalCrossSectionAux.segment(1,Nnodes-2)).sum();

	// characteristic energy
	swarmParam["characEnergy_eedf"] = swarmParam["redDiffCoeff_eedf"]/swarmParam["redMobCoeff_DC_eedf"];
}

void BoltzmannMC::evaluateRateCoeff(){
	// clear rate coeff vectors
	rateCoeffAll.clear();
	rateCoeffExtra.clear();

	// evaluate rate coefficients for all MC processes
	for (int i = 0; i < nProcesses; ++i){
		RateCoeffStruct  auxRateCoeff;
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
		RateCoeffStruct  auxRateCoeff;
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
}

