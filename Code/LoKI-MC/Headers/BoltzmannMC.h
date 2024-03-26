#ifndef __BoltzmannMC__
#define __BoltzmannMC__

#include "LoKI-MC/Headers/GeneralDefinitions.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Setup.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/AngularScatteringFunctions.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <string>
#include <map>
#include <boost/signals2.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <ctime>

class BoltzmannMC{

public:
	// ----- general class attributes -----//

	std::vector<EedfState*> stateArray;      // handle to the electron kinetics state mixture
	std::vector<EedfGas*> gasArray;		     // handle to the electron kinetics gas mixture
	Grid* energyGrid;				         // handle to the energy grid where the eedf is solved
	WorkingConditions* workCond;	         // handle to the working conditions of the simulation

	std::vector<std::string> CARgases;		 // this is not important for this class. Just for PrescribedEedf

	Eigen::ArrayXd totalMomTransfCrossSection;   // total momentum transfer cross section
	Eigen::ArrayXd elasticCrossSection;	         // total elastic cross section

	GeneralDefinitions::PowerStruct power; 		 // power balance
	std::map<std::string,double> swarmParam;     // swarm parameters obtained with eedf
	std::vector<GeneralDefinitions::RateCoeffStruct> rateCoeffAll;   // rate coefficients obtained with the eedf and collisions in gasArray
	std::vector<GeneralDefinitions::RateCoeffStruct> rateCoeffExtra; // extra rate coefficients for collisions not taken into account to obtain the eedf
	std::vector<std::vector<GeneralDefinitions::RateCoeffStruct>> rateCoeffAll_periodic;   // when there is an AC E-field, rateCoeffAll for each sampling phase (between 0 and 2pi)
	std::vector<std::vector<GeneralDefinitions::RateCoeffStruct>> rateCoeffExtra_periodic; // when there is an AC E-field, rateCoeffExtra for each sampling phase (between 0 and 2pi)

	boost::signals2::signal<void ()> obtainedNewEedfSignal; // signal to activate methods of other classes when a new EEDF is obtained 

	// ----- class attributes directly related with the MC method -----//

	// simulation parameters provided by the user
	int energySharingIonizType;								// integer identifying the ionization operator type, concerning energy sharing between primary and secondary electrons
	double nElectrons;                   	    		    // number of electrons used in the simulation
	int gasTemperatureEffect;					    		// integer concerning the gas temperature effect (0 = "false", 1 = "true", 2 = "smartActivation")
	double energySharingFactor = 0.5;				  	    // energy sharing factor between the scattered and ejected electrons
	double minCollisionsBeforeSteadyState = 0;				// minimum number of electron collisions before we start the steady-state check
	double maxCollisionsBeforeSteadyState = 1E100;		    // number of electron collisions at which we consider that the steady-state has been achieved, independently of the temporal evolution of the energy
	double maxCollisionsAfterSteadyState = 1E100;  	        // maximum number of electron collisions that will be simulated after the steady-state
	int nPointsBetweenSteadyStateCheck = 100; 				// number of sampling points between steady-state checks (fixed at this moment)
	int nPointsBetweenStatErrorsCheck = 200;			    // number of sampling points between statistical error checks (fixed at this moment)
	double nEnergyCells = 1000; 				  			// number of energy cells used in the discretization of the electron distribution functions
	double nCosAngleCells = 100;    					    // number of cos(angle) cells used in the discretization of the electron angular distribution function
	double nAxialVelocityCells = 200;			            // number of axial velocity cells used in the discretization of the electron velocity distribution function
	double nRadialVelocityCells = 200;				        // number of radial velocity cells used in the discretization of the electron velocity distribution function
	double requiredMeanEnergyRelError = 1E100;    	        // maximum error for the mean energy
	double requiredFluxDriftVelocityRelError = 1E100;       // maximum error for the flux drift velocity 
	double requiredFluxDiffusionCoeffsRelError = 1E100;     // maximum error for the flux diffusion coefficients
	double requiredBulkDriftVelocityRelError = 1E100;       // maximum error for the bulk drift velocity 
	double requiredBulkDiffusionCoeffsRelError = 1E100;     // maximum error for the bulk diffusion coefficients	
	double requiredPowerBalanceRelError = 1E100;		    // maximum error for the power balance
	double requiredIntegrationPoints = 200;			        // number of integration points after which the simulation is stopped
	double requiredIntegratedSSTimes = 0;				    // number of steady-state times to be integrated
	double requiredIntegratedAbsoluteTime = 0;				// required absolute time to be integrated
	double synchronizationTimeXMaxCollisionFrequency = 1;   // syncTime/minCollTime = syncTime*maxCollFrequency
	int synchronizationOverSampling = 1;                    // number of synchronization points between sampling (>= 1)
	double initialElecTempOverGasTemp = 0.01;               // initial temperature of the electron ensemble, given as multiples of the gas temperature. Should be much smaller than 1, so as to determine a realistic steady-state/relaxation time under low E/N conditions
	double nIntegrationPhases = 100;                        // number of phases to be discretized along the period (for AC E and/or B)

	bool dispMCStatus;										// boolean to activate the terminal display of the MC status
	int nPointsBetweenDispInfo = 100;						// number of sampling points between updates of the MC status

	// quantities from the working conditions
	double gasTemperature;							    	// gas temperature from the working conditions [K]
	double totalGasDensity;							    	// gas density from the working conditions [m-3]
	double gasEnergy;										// gas energy [eV]

	// information regarding the processes considered in the BoltzmannMC calculations. Inelastic and superelastic collisions are considered as separate processes!
	int nProcesses;							  		             		 // total number of processes, considering inelastic and superelastic collisions as separate processes
	int* processTypes;					                             	 // reduced process types. They can be 'conservativeType', 'ionizationType' or 'attachmentType'
	std::vector<Collision*> realCollisionPointers;		                 // pointers to the 'Collision' objects
	std::vector<EedfGas*> gasPointers;							         // pointers to the target 'Gas' object
	bool* isElastic;							                     	 // booleans indicating the elastic processes
	bool* isSuperElastic;							                     // booleans regarding the process direction
	bool* isIonization;							                         // booleans indicating the ionization processes
	double* superElasticStatWeightFactors;								 // ratios of statistical weights used to calculate the superelastic cross sections (set to zero for the inelastic processes)
	std::vector<Eigen::ArrayXd> crossSectionEnergies, crossSectionValues;// arrays with the cross sections of each process
	double* energyMinLimits;                                             // first value of the cross section energy of each process [eV]
	double* energyMaxLimits;		                                     // last value of the cross section energy of each process [eV]
	double energyMaxElastic = 1E100;   						             // maximum energy (in eV) that will be used in the simulation: the minimum value of the maximum energies of all elastic cross sections
	std::vector<gsl_spline*> crossSectionInterpolations;                 // objects with the interpolations of the cross sections
	std::vector<gsl_interp_accel*> crossSectionInterpAccelerators;       // objects with the accelerators of the interpolations
	double *relDensities, *targetMasses, *reducedMasses, *energyLosses, *thermalStdDeviations;   // constants of each process (energyLosses are in eV!)
	double *crossSectionEnergyGrid;									     // energy grid used to interpolate the cross sections
	double crossSectionEnergyStep;    								     // energy step used in the energy grid
	double maxInterpolatedEnergy;										 // maximum energy of the energy grid
	int interpolCrossSectionSize = 1E4;									 // size of the energy grid
	double **interpolCrossSectionsXrelDens;                              // arrays with the interpolated cross sections (times relDensities) for a given energy grid. Size = {nProcesses, interpolCrossSectionSize}
	double **cumulSumInterpolCrossSectionsXrelDens;                      // arrays with the cumulative sum of the interpolated cross sections (times relDensities) for a given energy grid. Size = {nProcesses, interpolCrossSectionSize}
	double *totalCollisionFrequencies;									 // array with the total collision frequency for each value of the energy grid. Size = nProcesses
	double *maxCollisionFrequencies;									 // array with the maximum collision frequency in the region [0,energy i]. Size = nProcesses
	double *wParameters;                                                 // array with the w parameters (OPB or ionization potential). This is only relevant for the ionization processes
	int* targetGasIDs;													 // array with the IDs of the target gas in each process
	std::vector<AngularScatteringFunctions::functionPointer<BoltzmannMC>> angularScatteringFunctions; // array with the angular scattering function of each process
	std::vector<std::vector<double>> angularScatteringParams;			 // arrays with the parameters the angular scattering of each process
	bool* isMomentumConservationIonizationScattering;					 // booleans indicating the processes where momentumConservationIonization is used
	
	int nGases;							// number of target gases
	double* firstProcessIndexPerGas;    // first process index of each gas (used for the bissection method in performCollisions)
	double* lastProcessIndexPerGas;		// last process index of each gas (used for the bissection method in performCollisions)
	double* gasFractions;

	double time;  							// current GLOBAl simulation time [s]
	double deltaTSynchroniz;				// time-interval after which the electrons will be synchronized [s]
	double nextSynchronizTime;				// absolute time at which the electrons will be synchronized [s]

	Eigen::ArrayXd samplingTimes;			// all sampling times [s]
	double steadyStateTime;					// time instant at which the steady state was achieved [s]
	Eigen::ArrayXd integrationTimes;		// sampling times that are used in the integration [s]
	double timeIntervalInteg;
	double totalIntegratedTime;				// cumulative integrated time [s]
	Eigen::ArrayXd integrationPhases;		// all integration phases [rad], discretized in nIntegrationPhases
	double integrationPhaseStep;            // step between integration phases [rad]
	Eigen::ArrayXd nIntegrationPointsPerPhase;  // number of points that were summed per integration phase. Used for calculate the mean values in the end of the simulation

	double trialCollisionFrequency;							// current trial collision frequency [s^-1]
	Eigen::ArrayXd trialCollisionFrequenciesEachElectron;   // trial collision frequency [s^-1] that is being considered for each electron. Can be different than the general 'trialCollisionFrequency' immediately after an update
	double meanCollisionFrequency;							// mean collision frequency of the ensemble [s^-1]

	double reducedElecField;				// reduced electric field [Td]
	double electricFieldValue;
	Eigen::Array3d electricField;           // electric field vector [V/m]
	double elecFieldAngle;					// angle of the e-field relatively to z [degrees]
	double excitationFrequency;				// excitation frequency of the electric field [Hz]
	double excitationFrequencyRadians;		// excitation frequency in radians/s
	double reducedMagField;				    // reduced magnetic field [Hx]
	double cyclotronFrequency;				// cyclotron frequency [s-1]
	Eigen::Array3d accelerationElecField;   // array with the electron acceleration due to the electric field [m.s-2]
	bool isCylindricallySymmetric;

	Eigen::ArrayXXd electronPositions;				// matrix with electron positions with size {nElectrons, 3} [m]
	Eigen::ArrayXXd electronVelocities;				// matrix with electron velocities with size {nElectrons, 3} [m.s-1]
	Eigen::ArrayXd electronEnergies;				// array with the electron energies [eV]
	Eigen::ArrayXd electronCosAngles;				// array with the cos(angles that the electron velocities make with the z axis)
	Eigen::ArrayXd electronTimes;					// array with the non-synchronized times of each electron 
	Eigen::ArrayXd collisionFreeTimes;				// array with the collision-free times of each electron

	// number of sampling/integration points
	int nSamplingPoints, nIntegrationPoints, nSynchronizationPoints;
	// sampling/integration indices
	int currentSamplingIndex, firstIntegrationIndex;

	// arrays with the mean values for each sampling time
	Eigen::ArrayXd meanEnergies;         // [eV]
	Eigen::ArrayXXd meanPositions;	     // [m]
	Eigen::ArrayXXd meanVelocities;		 // [m.s-1]
	Eigen::ArrayXXd fluxDiffusionCoeffs; // [m2.s-1] each line has a flat matrix with the 9 coefficients: xx, xy, xz, yx, yy, yz, zx, zy, zz
	Eigen::ArrayXXd positionCovariances; // [m2] ensemble averages for each time instant // each line has a flat matrix with the 9 coefficients 
	Eigen::ArrayXXd bulkVelocities;      // [m.s-1]
	Eigen::ArrayXXd bulkDiffusionCoeffs; // [m2.s-1]

	// arrays with the mean values for each integration phase (for this part only the points after SS are considered)
	// only used for AC E-field
	Eigen::ArrayXd meanEnergies_periodic;
	Eigen::ArrayXXd fluxVelocities_periodic;
	Eigen::ArrayXXd bulkVelocities_periodic;
	Eigen::ArrayXXd fluxDiffusionCoeffs_periodic; // each line has a flat matrix with the 9 coefficients: xx, xy, xz, yx, yy, yz, zx, zy, zz
	Eigen::ArrayXXd bulkDiffusionCoeffs_periodic; // each line has a flat matrix with the 9 coefficients: xx, xy, xz, yx, yy, yz, zx, zy, zz	

	// data to calculate the power balance
	double energyGainField;									// total energy gain from the field [eV]
	double *energyGainProcesses, *energyLossProcesses;		// energy gain/loss from the processes [eV]
	double energyGrowth;									// energy difference due to the electron density profile [eV]

	// collision counters
	unsigned long long int totalCollisionCounter, nullCollisionCounter, collisionCounterAtSS, nullCollisionCounterAtSS, collisionCounterAfterSS;
	unsigned long long int *collisionCounters; // size = {nProcesses}

	// variables related with the energy parameters
	double averagedMeanEnergy, averagedMeanEnergyError;
	double maxElecEnergy;

	// variables related with the distribution functions
	double eedfEnergyStep, maxEedfEnergy;		                  // energy step and max energy used in the discretization of the distribution functions
	Eigen::ArrayXd eedfEnergyNodes, eedfEnergyCells;              // arrays with the energies used in the discretization of the distribution functions
	Eigen::ArrayXd eeh, eehSum, eedf;                     		  // arrays used to calculate the time-averaged electron energy distribution function
	Eigen::ArrayXXd eehSum_periodic, eedf_periodic;               // arrays used to calculate the eedfs for each integration phase (for AC E-field conditions)
	Eigen::ArrayXd cosAngleNodes, cosAngleCells; 			      // arrays with the angles used in the discretization of the electron angular distribution function
	double cosAngleStep;							              // cos(angle) step used in the discretization of the electron angular distribution function
	Eigen::ArrayXXd eahSum, eadf;							      // electron angular distribution function
	Eigen::ArrayXd efadf;				            		      // electron first-anisotropy distribution function		
	Eigen::ArrayXd esadf;				            		      // electron second-anisotropy distribution function
	Eigen::ArrayXd	radialVelocityNodes, radialVelocityCells,
			axialVelocityNodes, axialVelocityCells;               // arrays used in the discretization of the electron velocity distribution function
	double axialVelocityStep, radialVelocityStep;                 // velocity steps used in the discretization of the electron velocity distribution function
	Eigen::ArrayXXd evhSum, evdf;							      // electron velocity distribution function

	// variables related with the flux swarm parameters
	Eigen::Array3d averagedFluxDriftVelocity, averagedFluxDriftVelocityError;      // [m.s-1]
	Eigen::Matrix3d averagedFluxDiffusionCoeffs, averagedFluxDiffusionCoeffsError;  // [m2.s-1] // matrix 3x3
	Eigen::Array3d rotatedAveragedFluxDriftVelocity, rotatedAveragedFluxDriftVelocityError;      // [m.s-1], in the ref-frame where E is along z
	Eigen::Matrix3d rotatedAveragedFluxDiffusionCoeffs, rotatedAveragedFluxDiffusionCoeffsError;  // [m2.s-1]. matrix 3x3, in the ref-frame where E is along z

	// variables related with the rate coefficients
	double* averagedRateCoeffs;

	// variables related with the power balance
	double averagedPowerGainField, averagedPowerGrowth;
	double *averagedPowerGainProcesses, *averagedPowerLossProcesses;
	double powerBalanceRelError;

	// variables related with the bulk swarm parameters
	Eigen::Array3d averagedBulkDriftVelocity, averagedBulkDriftVelocityError;      // [m.s-1]
	Eigen::Matrix3d averagedBulkDiffusionCoeffs, averagedBulkDiffusionCoeffsError;  // [m2.s-1] // matrix 3x3
	Eigen::Array3d rotatedAveragedBulkDriftVelocity, rotatedAveragedBulkDriftVelocityError;      // [m.s-1], in the ref-frame where E is along z
	Eigen::Matrix3d rotatedAveragedBulkDiffusionCoeffs, rotatedAveragedBulkDiffusionCoeffsError;  // [m2.s-1]. matrix 3x3, in the ref-frame where E is along z

	// arrays with the information of each electron. Added ultimately to enable parallelization
	int* chosenProcessIDs;                                                  // size = nElectrons
	Eigen::ArrayXXd ejectedElectronPositions, ejectedElectronVelocities;    // only used when an ionization process is chosen. Size = {nElectrons,3}
	Eigen::ArrayXd  ejectedElectronEnergies;								// only used when an ionization process is chosen. Size = nElectrons
	Eigen::ArrayXd electronEnergyChanges;                                   // electron energy changes due to the chosen process. Size = nElectrons 
	Eigen::ArrayXd electronEnergyChangesOverIncidEnergies;					// electron energy changes due to the chosen process divided by the incident energy. Used for the collision counters that are energy weighted. Size = nElectrons
	Eigen::ArrayXd energyGainsField; 										// energy gains from the field of each electron [eV]

	// boolean indicating if the desired statistical errors are satisfied
	bool goodStatisticalErrors = false;
	bool errorsToBeChecked = true;

	// calculation time (reset at each 'solve' evaluation)
	double elapsedTime;

	int failedCollisionDueToGrid;

	// ----- general class methods -----//

	template <class SetupType>
	BoltzmannMC(SetupType* setup){
		// this method is written in the '.h' file to avoid compilation errors associated with the templates

		// store the state array
		stateArray = setup->electronKineticsStateArray;

		// store the gas array
		gasArray = setup->electronKineticsGasArray;

		// store the energy grid
		energyGrid = setup->energyGrid;

	    // store working conditions
	    workCond = setup->workCond;

	    // ----- store the simulation parameters parsed from the setup file ----- //

	    std::string ionizationOperatorType = FieldInfo::getFieldValue("electronKinetics.ionizationOperatorType");
	    if(ionizationOperatorType == "equalSharing"){
	    	energySharingIonizType = GeneralDefinitions::equalSharingIonizType;
	    	energySharingFactor = 0.5;
	    }
	    else if(ionizationOperatorType == "oneTakesAll"){
	    	energySharingIonizType = GeneralDefinitions::oneTakesAllIonizType;
	    	energySharingFactor = 0;
	    }
	    else if (ionizationOperatorType == "usingSDCS"){
	    	energySharingIonizType = GeneralDefinitions::usingSDCSIonizType;
	    }
	    else if (ionizationOperatorType == "randomUniform"){
	    	energySharingIonizType = GeneralDefinitions::randomUniformIonizType;
	    }	    

	    // get the mandatory 'numericsMC' fields (nElectrons and gasTemperatureEffect)
		nElectrons = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nElectrons");
		std::string gasTempString = FieldInfo::getFieldValue("electronKinetics.numericsMC.gasTemperatureEffect");
		if (gasTempString == "false"){
			gasTemperatureEffect = GeneralDefinitions::falseGasTempEffectID;
		}
		else if (gasTempString == "true"){
			gasTemperatureEffect = GeneralDefinitions::trueGasTempEffectID;
		}
		else if (gasTempString == "smartActivation"){
			gasTemperatureEffect = GeneralDefinitions::smartActivationGasTempEffectID;
		}

		// get the optional fields. If not defined, the declaration values are used
		if (FieldInfo::getField("electronKinetics.numericsMC.initialElecTempOverGasTemp")){
			initialElecTempOverGasTemp = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.initialElecTempOverGasTemp");
		}		
		if (FieldInfo::getField("electronKinetics.numericsMC.minCollisionsBeforeSteadyState")){
			minCollisionsBeforeSteadyState = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.minCollisionsBeforeSteadyState") * nElectrons;
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.maxCollisionsBeforeSteadyState")){
			maxCollisionsBeforeSteadyState = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.maxCollisionsBeforeSteadyState") * nElectrons;
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.maxCollisionsAfterSteadyState")){
			maxCollisionsAfterSteadyState = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.maxCollisionsAfterSteadyState") * nElectrons;		
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nEnergyCells")){
			nEnergyCells = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nEnergyCells");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nCosAngleCells")){
			nCosAngleCells = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nCosAngleCells");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nRadialVelocityCells")){
			nRadialVelocityCells = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nRadialVelocityCells");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nAxialVelocityCells")){
			nAxialVelocityCells = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nAxialVelocityCells");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.relError")){
			errorsToBeChecked = true;
		}
		else{
			errorsToBeChecked = false;
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.meanEnergy")){
			requiredMeanEnergyRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.meanEnergy");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.fluxDriftVelocity")){
			requiredFluxDriftVelocityRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.fluxDriftVelocity");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.bulkDriftVelocity")){
			requiredBulkDriftVelocityRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.bulkDriftVelocity");
		}		
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs")){
			requiredFluxDiffusionCoeffsRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.bulkDiffusionCoeffs")){
			requiredBulkDiffusionCoeffsRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.bulkDiffusionCoeffs");
		}		
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.powerBalance")){
			requiredPowerBalanceRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.powerBalance");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nIntegrationPoints")){
			requiredIntegrationPoints = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPoints");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nIntegratedSSTimes")){
			requiredIntegratedSSTimes = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegratedSSTimes");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.integratedAbsoluteTime")){
			requiredIntegratedAbsoluteTime = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.integratedAbsoluteTime");
		}
		if (FieldInfo::getField("electronKinetics.numericsMC.nIntegrationPhases")){
			nIntegrationPhases = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nIntegrationPhases");
		}				

		// set the number of points used for the interpolation of the cross sections
		if (FieldInfo::getField("electronKinetics.numericsMC.nInterpPoints")){
			interpolCrossSectionSize = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nInterpPoints");
		}

		// synchronization time for the electron advance
		if (FieldInfo::getField("electronKinetics.numericsMC.synchronizationTimeXMaxCollisionFrequency")){
			synchronizationTimeXMaxCollisionFrequency = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.synchronizationTimeXMaxCollisionFrequency");
		}

		if (FieldInfo::getField("electronKinetics.numericsMC.synchronizationOverSampling")){
			synchronizationOverSampling = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.synchronizationOverSampling");
		}			

		// check if the display of MC status is on
		dispMCStatus = false;
		if ( Parse::str2bool(FieldInfo::getFieldValue("gui.isOn")) ){
			for (auto option: FieldInfo::getFieldChildNames("gui.terminalDisp")){
				if (option == "MCStatus"){
					dispMCStatus = true;
					break;
				}
			}
		}

	    // Allocate and evaluate some variables that will be used in the MC method
	    allocateEvaluateVariablesFirstTime();
	}
	void allocateEvaluateVariablesFirstTime();

	void updateDensityDependencies(){;}

	void solve();

	// ----- class methods directly related with the MC method -----//

	void evaluateEEDF(); 
	void evaluateNonConstantVariables();
	void interpolateCrossSections(double maxEnergy);

	void electronDynamicsUntilSynchronization();
	void getMeanCollisionFrequency();
	void checkMaxCollisionFrequency();
	double maximizationAccelerationEnergy(double initialEnergy, double deltaT);
	void accelerateElectron(int elecID, double electronTime, double deltaT, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy);

	void performCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy);
	void conservativeCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, Eigen::Array3d &targetVelocity, double &electronEnergy, int chosenProcessID);
	void ionizationCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double &electronEnergy, int chosenProcessID);
	void attachmentCollision(int elecID, double incidentEnergy);
	void nonParallelCollisionTasks();

	void calculateMeanDataForSwarmParams();

	void getTimeAverageEnergyParams();
	void getTimeDependDistributions();
	void getTimeAverageDistributions();
	void getTimeAverageFluxParams();
	void getTimeAverageRateCoeffs();
	void getTimeAveragePowerBalance();
	void getTimeAverageBulkParams();
	void getAveragedPeriodicParams();

	void checkStatisticalErrors();
	void checkPowerBalance();
	void checkSteadyState();

	void dispInfo();

	// ----- continuation of general class methods ---- //

	void evaluatePower();
	void evaluateSwarmParameters();
	void evaluateRateCoeff();

	// variables added to avoid compilation problems
	int nFreeFlights, updateFrequencyCounter, failedTrialCounter;
	Eigen::ArrayXd nElectronsArray;
};

#endif