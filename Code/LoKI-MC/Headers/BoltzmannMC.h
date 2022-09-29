#ifndef __BoltzmannMC__
#define __BoltzmannMC__

#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Setup.h"
#include "LoKI-MC/Headers/PrescribedEedf.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/Parse.h"
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

	std::vector<EedfState*> stateArray;
	std::vector<EedfGas*> gasArray;		     // handle to the electron kinetics gas mixture
	Grid* energyGrid;				         // handle to the energy grid where the eedf is solved
	WorkingConditions* workCond;	         // handle to the working conditions of the simulation

	std::vector<std::string> CARgases;		 // this is not important for this class. Just for PrescribedEedf

	Eigen::ArrayXd totalCrossSection;		     // total momentum transfer cross section
	Eigen::ArrayXd elasticCrossSection;	         // total elastic cross section

	PowerStruct power; 				             // power balance
	std::map<std::string,double> swarmParam;     // swarm parameters obtained with eedf
	std::vector<RateCoeffStruct> rateCoeffAll;   // rate coefficients obtained with the eedf and collisions in gasArray
	std::vector<RateCoeffStruct> rateCoeffExtra; // extra rate coefficients for collisions not taken into account to obtain the eedf

	boost::signals2::signal<void ()> obtainedNewEedfSignal;

	// ----- class attributes directly related with the MC method -----//

	// simulation parameters provided by the user
	bool usingSDCS = false;								    // boolean to know if SDCS is to be used in the ionization
	double nElectrons;                   	    		    // number of electrons used in the simulation
	int gasTemperatureEffect;					    		// integer concerning the gas temperature effect (0 = "false", 1 = "true", 2 = "smartActivation")
	double energySharingFactor;				  	            // energy sharing factor between the scattered and ejected electrons
	bool isoScatteringIonization;
	double minCollisionsBeforeSteadyState = 0;				// minimum number of electron collisions before we start the steady-state check
	double maxCollisionsBeforeSteadyState = 1E100;		    // number of electron collisions at which we consider that the steady-state has been achieved, independently of the temporal evolution of the energy
	double maxCollisionsAfterSteadyState = 1E100;  	        // maximum number of electron collisions that will be simulated after the steady-state
	int nPointsBetweenSteadyStateCheck = 10; 				// number of sampling points between steady-state checks (fixed at this moment)
	int nPointsBetweenStatErrorsCheck = 100;			    // number of sampling points between statistical error checks (fixed at this moment)
	double nEnergyCells = 500; 				  			    // number of energy cells used in the discretization of the electron distribution functions
	double nCosAngleCells = 50;    				            // number of cos(angle) cells used in the discretization of the electron angular distribution function
	double nAxialVelocityCells = 200;			            // number of axial velocity cells used in the discretization of the electron velocity distribution function
	double nRadialVelocityCells = 200;				        // number of radial velocity cells used in the discretization of the electron velocity distribution function
	double requiredMeanEnergyRelError = 1E100;    	        // maximum error for the mean energy
	double requiredFluxDriftVelocityRelError = 1E100;       // maximum error for the flux drift velocity (along z)
	double requiredFluxDiffusionCoeffsRelError = 1E100;     // maximum error for the flux diffusion coefficients
	double requiredPowerBalanceRelError = 1E100;		    // maximum error for the power balance
	double requiredIntegrationPoints = 1E100;				// number of integration points after which the simulation is stopped
	double requiredIntegratedSSTimes = 1E100;				// number of steady-state times to be integrated

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
	int interpolCrossSectionSize;										 // size of the energy grid
	double **interpolCrossSectionsXrelDens;                              // arrays with the interpolated cross sections (times relDensities) for a given energy grid. Size = {nProcesses, interpolCrossSectionSize}
	double **cumulSumInterpolCrossSectionsXrelDens;                      // arrays with the cumulative sum of the interpolated cross sections (times relDensities) for a given energy grid. Size = {nProcesses, interpolCrossSectionSize}
	double *totalCollisionFrequencies;									 // array with the total collision frequency for each value of the energy grid. Size = nProcesses
	double *maxCollisionFrequencies;									 // array with the maximum collision frequency in the region [0,energy i]. Size = nProcesses
	double *wParameters;                                                 // array with the w parameters (OPB or ionization potential). This is only relevant for the ionization processes
	int* targetGasIDs;													 // array with the IDs of the target gas in each process. This is used for the gas temperature effect

	double time;  							// current simulation time [s]
	Eigen::ArrayXd samplingTimes;			// all sampling times [s]
	double deltaT;							// current time step [s]
	double steadyStateTime;					// time instant at which the steady state was achieved [s]
	Eigen::ArrayXd integrationTimes;		// sampling times that are used in the integration [s]
	double timeIntervalInteg;
	double totalIntegratedTime;				// cumulative integrated time [s]

	double trialCollisionFrequency;         // trial collision frequency [s^-1]

	double reducedElecField;				// reduced electric field [Td]
	Eigen::Array3d electricField;           // electric field vector [V/m]
	Eigen::Array3d accelerationElecField;   // array with the electron acceleration due to the electric field [m.s-2]

	Eigen::ArrayXXd electronPositions;				// matrix with electron positions with size {nElectrons, 3} [m]
	Eigen::ArrayXXd electronVelocities;				// matrix with electron velocities with size {nElectrons, 3} [m.s-1]
	Eigen::ArrayXd electronEnergies;				// array with the electron energies [eV]
	Eigen::ArrayXd electronCosAngles;				// array with the cos(angles that the electron velocities make with the z axis)

	// number of sampling/integration points
	int nSamplingPoints, nIntegrationPoints;
	// sampling/integration indices
	int currentSamplingIndex, firstIntegrationIndex;

	// arrays with the mean values for each sampling time
	Eigen::ArrayXd meanEnergies;         // [eV]
	Eigen::ArrayXXd meanPositions;	     // [m]
	Eigen::ArrayXXd meanVelocities;		 // [m.s-1]
	Eigen::ArrayXXd fluxDiffusionCoeffs; // each line has a flat matrix with the 9 coefficients: xx, xy, xz, yx, yy, yz, zx, zy, zz
	Eigen::ArrayXXd positionCovariances; // ensemble averages for each time instant // each line has a flat matrix with the 9 coefficients 

	// data to calculate the power balance
	double energyGainField;									// energy gain from the field [eV]
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
	Eigen::ArrayXd eehSum, eedf;                     			  // arrays used to calculate the time-averaged electron energy distribution function
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
	Eigen::ArrayXd averagedFluxDiffusionCoeffs, averagedFluxDiffusionCoeffsError;  // [m2.s-1] // flat matrix with 9 coeffs: xx, xy, xz, yx, yy, yz, zx, zy, zz

	// variables related with the rate coefficients
	double* averagedRateCoeffs;

	// variables related with the power balance
	double averagedPowerGainField, averagedPowerGrowth;
	double *averagedPowerGainProcesses, *averagedPowerLossProcesses;
	double powerBalanceRelError;

	// variables related with the bulk swarm parameters
	Eigen::Array3d averagedBulkDriftVelocity, averagedBulkDriftVelocityError;     // [m.s-1]
	Eigen::ArrayXd averagedBulkDiffusionCoeffs, averagedBulkDiffusionCoeffsError; // [m2.s-1] // flat matrix with 9 coeffs: xx, xy, xz, yx, yy, yz, zx, zy, zz

	// arrays with the information of each electron. Added ultimately to enable parallelization
	int* chosenProcessIDs;                                                 // size = nElectrons
	Eigen::ArrayXXd ejectedElectronPositions, ejectedElectronVelocities;   // only used when an ionization process is chosen. Size = {nElectrons,3}
	double* electronEnergyChanges;                                         // electron energy changes due to the chosen process. Size = nElectrons

	// boolean indicating if the desired statistical errors are satisfied
	bool goodStatisticalErrors = false;
	bool errorsToBeChecked = true;

	// calculation time (reset at each 'solve' evaluation)
	double elapsedTime;

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
	    if (ionizationOperatorType == "usingSDCS"){
	    	usingSDCS = true;
	    }
	    else if(ionizationOperatorType == "equalSharing"){
	    	usingSDCS = false;
	    	energySharingFactor = 0.5;
	    }
	    else if(ionizationOperatorType == "oneTakesAll"){
	    	usingSDCS = false;
	    	energySharingFactor = 0;
	    }

	    if (FieldInfo::getFieldValue("electronKinetics.ionizationScattering") == "isotropic"){
	    	isoScatteringIonization = true;
	    }
	    else{
	    	isoScatteringIonization = false;
	    }

	    // get the mandatory 'numericsMC' fields (nElectrons and gasTemperatureEffect)
		nElectrons = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nElectrons");
		std::string gasTempString = FieldInfo::getFieldValue("electronKinetics.numericsMC.gasTemperatureEffect");
		if (gasTempString == "false"){
			gasTemperatureEffect = 0;
		}
		else if (gasTempString == "true"){
			gasTemperatureEffect = 1;
		}
		else if (gasTempString == "smartActivation"){
			gasTemperatureEffect = 2;
		}

		// get the optional fields. If not defined, the declaration values are used
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
		if (FieldInfo::getField("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs")){
			requiredFluxDiffusionCoeffsRelError = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.relError.fluxDiffusionCoeffs");
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

		// set the number of points used for the interpolation of the cross sections
		if (FieldInfo::getField("electronKinetics.numericsMC.nInterpPoints")){
			interpolCrossSectionSize = FieldInfo::getFieldNumericValue("electronKinetics.numericsMC.nInterpPoints");
		}
		else{
			interpolCrossSectionSize = 1E4;
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
	double getMaxCollisionFrequency(double maxEnergy);
	void interpolateCrossSections(double maxEnergy);

	void freeFlight();

	void calculateMeanData();
	void calculateMeanDataForSwarmParams();

	void getTimeAverageEnergyParams();
	void getTimeDependDistributions();
	void getTimeAverageDistributions();
	void getTimeAverageFluxParams();
	void getTimeAverageRateCoeffs();
	void getTimeAveragePowerBalance();
	void getTimeAverageBulkParams();

	void performCollisions();
	void conservativeCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, Eigen::Array3d &targetVelocity, double incidentEnergy, int chosenCollisionID);
	void ionizationCollision(int elecID, Eigen::Array3d &electronPosition, Eigen::Array3d &electronVelocity, double incidentEnergy, int chosenCollisionID);
	void attachmentCollision(int elecID, double incidentEnergy);
	void nonParallelCollisionTasks();

	void checkStatisticalErrors();
	void checkPowerBalance();
	void checkSteadyState();

	void dispInfo();

	// ----- continuation of general class methods ---- //

	void evaluatePower();
	void evaluateSwarmParameters();
	void evaluateRateCoeff();
};

#endif