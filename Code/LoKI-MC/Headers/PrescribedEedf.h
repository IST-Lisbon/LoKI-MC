#ifndef __PrescribedEedf__
#define __PrescribedEedf__

#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <string>
#include <map>
#include <boost/signals2.hpp>

// ----- structures used in 'PrescribedEedf' and 'BoltzmannMC' -----

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
// -------

class PrescribedEedf{
  //  PrescribedEedf evaluates a generalized EEDF, corresponding to a particular electron temperature, as well as the
  //  corresponding swarm parameters and electron impact rate coefficients. The generalized expression, that has as 
  //  limiting cases the Maxwellian and Druyvesteyn distribution functions, is as follows:
  // 
  //  f(u) = (gamma(5/(2g))^3/2)/(gamma(3/(2g))^5/2)*(2/(3KbTe))^3/2*exp(-(2u*gamma(5/(2g)/(3KbTe*gamma(3/(2g))))^g)
  // 
  //  where gamma is the gamma function (https://en.wikipedia.org/wiki/Gamma_function), Kb the boltzmann constant, 
  //  Te the electron temperature and g is the parameter thar controls the shape of the distribution function (g=1
  //   for Maxwellian and g=2 for Druyvesteyn).

public:
	// ----- class attributes -----

	std::vector<EedfState*> stateArray;
	std::vector<EedfGas*> gasArray;			 // handle to the electron kinetics gas mixture
	Grid* energyGrid;				             // handle to the energy grid where the eedf is solved
	WorkingConditions* workCond;	       // handle to the working conditions of the simulation

	std::vector<std::string> CARgases;

	Eigen::ArrayXd totalCrossSection;		     // total momentum transfer cross section
	Eigen::ArrayXd elasticCrossSection;	     // total elastic cross section

	double shapeParameter = Constant::NON_DEF;  // value of the parameter controling the shape of the eedf

	Eigen::ArrayXd eedf;					         // eedf

	PowerStruct power; 				         // power balance
	std::map<std::string,double> swarmParam;	    // swarm parameters obtained with eedf
	std::vector<RateCoeffStruct> rateCoeffAll;    // rate coefficients obtained with the eedf and collisions in gasArray
	std::vector<RateCoeffStruct> rateCoeffExtra;  // extra rate coefficients for collisions not taken into account to obtain the eedf

	boost::signals2::signal<void ()> obtainedNewEedfSignal;

	// ----- class methods -----

	template <class SetupType>
	PrescribedEedf(SetupType* setup){
		// this metod is declared in the '.h' file to avoid compilation errors associated with the templates

		// store the state array
		stateArray = setup->electronKineticsStateArray;

		// store the gas array
		gasArray = setup->electronKineticsGasArray;

		// store gases for which the CAR is activated (in case there is any)
		if (FieldInfo::getField("electronKinetics.CARgases")){
			CARgases = FieldInfo::getFieldChildNames("electronKinetics.CARgases");
		}

		// store the energy grid
		energyGrid = setup->energyGrid;

    // connect the signal 'updatedMaxEnergy2Signal' to evaluateTotalAndElasticCrossSections, in order to update these cross sections each time the max energy is changed
    // idea taken from https://stackoverflow.com/questions/3047381/boost-signals-and-passing-class-method
    energyGrid->updatedMaxEnergy2Signal.connect(boost::bind(&PrescribedEedf::evaluateTotalAndElasticCrossSections, this));

    // store working conditions
    workCond = setup->workCond;

    // store the value of the parameter that controls the shape of the eedf (1 Maxwellian, 2 Druyvesteyn)
    shapeParameter = FieldInfo::getFieldNumericValue("electronKinetics.shapeParameter");

    // evaluate total momentum transfer and total elastic cross sections
    evaluateTotalAndElasticCrossSections();
	}
	void solve();
	void updateDensityDependencies();
	void evaluateTotalAndElasticCrossSections();
	void evaluateEEDF();
	void evaluatePower();
	void evaluateSwarmParameters();
	void evaluateRateCoeff();

	// variables added here to avoid compilation problems
	Eigen::ArrayXd samplingTimes, meanEnergies, efadf, esadf, eedfEnergyCells, cosAngleCells, radialVelocityCells, axialVelocityCells, nElectronsArray;
	Eigen::ArrayXd averagedFluxDriftVelocityError, averagedFluxDriftVelocity, averagedFluxDiffusionCoeffsError, averagedFluxDiffusionCoeffs, averagedBulkDriftVelocityError, averagedBulkDriftVelocity, averagedBulkDiffusionCoeffsError, averagedBulkDiffusionCoeffs;					 
	Eigen::ArrayXXd meanPositions, positionCovariances, meanVelocities, eadf, evdf;
	double nElectrons, elapsedTime, steadyStateTime, time, totalGasDensity;
	double averagedMeanEnergyError, averagedMeanEnergy, requiredMeanEnergyRelError, requiredFluxDriftVelocityRelError, requiredFluxDiffusionCoeffsRelError, powerBalanceRelError, requiredPowerBalanceRelError, eedfMeanRelError, requiredEedfMeanRelError;
	unsigned long long int totalCollisionCounter, nullCollisionCounter, collisionCounterAtSS, nullCollisionCounterAtSS;
	int nSamplingPoints, nIntegrationPoints, nFreeFlights, updateFrequencyCounter, failedTrialCounter;
	bool isCylindricallySymmetric;
};

#endif