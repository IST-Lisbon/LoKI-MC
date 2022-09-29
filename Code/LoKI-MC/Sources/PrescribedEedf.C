#include "LoKI-MC/Headers/PrescribedEedf.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/WorkingConditions.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Parse.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/signals2.hpp>
#include <cmath>
#include <boost/algorithm/string.hpp>

// the constructor is declared in the '.h' file

void PrescribedEedf::solve(){

	// when the smart grid is activated the minimum number of decades of decay in the eedf is ensured
	if (energyGrid->isSmart){
		double g = shapeParameter;
		double Te = workCond->electronTemperature;
		double decades = 0.5*(energyGrid->minEedfDecay + energyGrid->maxEedfDecay);
		double maxEnergy = std::pow((decades/std::log10(std::exp(1))), 1.0/g) * 1.5 * Te * std::tgamma(3.0/(2.0*g)) / std::tgamma(5.0/(2.0*g));
		energyGrid->updateMaxValue(maxEnergy);
	}

	// evaluate EEDF
	evaluateEEDF();

	// evaluate power balance
	evaluatePower();

	// evaluate rate coefficients
	evaluateRateCoeff();

	// evaluate transport parameters
	evaluateSwarmParameters();

	// broadcast obtention of a solution for the EEDF
	obtainedNewEedfSignal();
}

void PrescribedEedf::updateDensityDependencies(){
	// updateDensityDependencies is a function that evaluates all the species density dependencies of the prescribedEedf object
	evaluateTotalAndElasticCrossSections();
}

void PrescribedEedf::evaluateTotalAndElasticCrossSections(){

	// reset values to zero
	totalCrossSection = Eigen::ArrayXd::Zero(energyGrid->node.size());
	elasticCrossSection = Eigen::ArrayXd::Zero(energyGrid->node.size());

	// loop over each gas in the mixture
	for (auto& gas: gasArray){

		// avoid dummy gases
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
			double temp;
			// add collision cross section to the total momentum transfer cross section (also superelastic)
			totalCrossSection += collision->crossSection * collision->target->density;
			if (collision->isReverse){
				totalCrossSection += collision->superElasticCrossSection(Eigen::ArrayXd(0)) * collision->productArray[0]->density;
			}
			// add elastic collision cross section to the total elastic cross section (weighted by the mass ratio)
			if (collision->type == "Elastic"){
				elasticCrossSection += collision->crossSection * massRatio * collision->target->density;
				continue;
			}
		}
	}
}

void PrescribedEedf::evaluateEEDF(){

	// evaluate the distribution without normalization
	double g = shapeParameter;
	double Te = workCond->electronTemperature;
	Eigen::ArrayXd energy = energyGrid->cell;
	eedf = Eigen::exp( -Eigen::pow(energy * std::tgamma(5.0/(2.0*g)) * 2.0 / ( std::tgamma(3.0/(2.0*g)) * 3.0 * Te ), g) );

	// renormalizing the distribution
	eedf /= (Eigen::sqrt(energy)*eedf).sum()*energyGrid->step;
}

void PrescribedEedf::evaluatePower(){
	// initialize power structure;
	std::map<std::string,double> powerMap;
	powerMap["field"] = 0; powerMap["elasticNet"] = 0; powerMap["elasticGain"] = 0; powerMap["elasticLoss"] = 0; powerMap["carNet"] = 0; powerMap["carGain"] = 0; powerMap["carLoss"] = 0;
	powerMap["excitationIne"] = 0; powerMap["excitationSup"] = 0; powerMap["excitationNet"] = 0; powerMap["vibrationalIne"] = 0; powerMap["vibrationalSup"] = 0; powerMap["vibrationalNet"] = 0;
	powerMap["rotationalIne"] = 0; powerMap["rotationalSup"] = 0; powerMap["rotationalNet"] = 0; powerMap["ionizationIne"] = 0; powerMap["attachmentIne"] = 0; powerMap["inelastic"] = 0;
	powerMap["superelastic"] = 0; powerMap["eDensGrowth"] = 0; powerMap["electronElectron"] = 0;

	// save a local copy of the energy grid information
	int N = energyGrid->cellNumber;
	double energyStep = energyGrid->step;
	Eigen::ArrayXd energyCell = energyGrid->cell;
	Eigen::ArrayXd energyNode = energyGrid->node;

	// multiplicative constant to obtain the right units
	double factor = std::sqrt(2.0*Constant::electronCharge/Constant::electronMass);

	// evaluate power absorved per electron at unit gas density due to the electric field (auxiliary value)
	// evaluation of the elements of the electric field operator (Boltzmann)
	Eigen::ArrayXd g_E = energyNode/(totalCrossSection*3.0);
	double auxPowerField = factor * (eedf * ( g_E.segment(1, g_E.size()-1) - g_E.segment(0, g_E.size()-1)) ).sum();
	powerMap["field"] = Constant::NON_DEF;

	// auxiliary quantities needed to evaluate the elastic and CAR powers
	double kTg = Constant::boltzmannInEV * workCond->gasTemperature;
	double aux1 = kTg + energyStep*0.5;
	double aux2 = kTg - energyStep*0.5;

	// evaluate power absorved per electron at unit gas density due to elastic collisions
	// evaluation of the elements of the elastic operator (Boltzmann)
	Eigen::ArrayXd g_c = Eigen::pow(energyNode, 2.0) * 2.0 * elasticCrossSection;
	g_c(0) = 0;
	g_c(g_c.size()-1) = 0;
	powerMap["elasticNet"] = factor * (eedf * ( g_c.segment(1, g_c.size()-1)*aux2 - g_c.segment(0, g_c.size()-1)*aux1) ).sum();
	powerMap["elasticGain"] = factor * kTg * (eedf * ( g_c.segment(1, g_c.size()-1) - g_c.segment(0, g_c.size()-1)) ).sum();
	powerMap["elasticLoss"] = powerMap["elasticNet"] - powerMap["elasticGain"];

	// evaluate power absorbed per electron at unit gas density due to rotations CAR
	if (!CARgases.empty()){
		// evaluation of the elements of the CAR operator (Boltzmann)
		double sigma0B = 0;
		for (auto gasName: CARgases){
			int gasID = EedfGas::find(gasName, gasArray);
			if (gasID != -1){
				EedfGas* gas = gasArray[gasID];
				sigma0B += gas->fraction*gas->electricQuadrupoleMoment*gas->rotationalConstant;
			}
		}
		Eigen::ArrayXd g_car = energyNode*4.0*(8.0*M_PI*sigma0B/(15.0*Constant::electronCharge));
		powerMap["carNet"] = factor * (eedf * ( g_car.segment(1, g_car.size()-1)*aux2 - g_car.segment(0, g_car.size()-1)*aux1) ).sum();
		powerMap["carGain"] = factor * kTg * (eedf * ( g_car.segment(1, g_car.size()-1) - g_car.segment(0, g_car.size()-1)) ).sum();
		powerMap["carLoss"] = powerMap["carNet"] - powerMap["carGain"];
	}

	// evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic collisions
	// loop over each gas in the mixture
	for (auto gas: gasArray){
		std::string gasName = gas->name;
		std::map<std::string,double> gasPower;
		// initialize power balance information of this gas
		gasPower["excitationIne"] = 0; gasPower["excitationSup"] = 0; gasPower["excitationNet"] = 0;
		gasPower["vibrationalIne"] = 0; gasPower["vibrationalSup"] = 0; gasPower["vibrationalNet"] = 0;
		gasPower["rotationalIne"] = 0; gasPower["rotationalSup"] = 0; gasPower["rotationalNet"] = 0;
		gasPower["ionizationIne"] = 0; gasPower["attachmentIne"] = 0;

		// loop over each collision with the gas
		for (auto& collision: gas->collisionArray){
			// collision type
			std::string collType = collision->type;
			// avoid Effective or Elastic collisions and collisions whose threshold is larger than the maximum energy
			if (collType == "Effective" || collType == "Elastic" || collision->threshold > energyNode(N)){
				continue;
			}
			else if (collType == "Attachment"){
				// evaluate cross section at cell positions
				Eigen::ArrayXd cellCrossSection = (collision->crossSection.segment(1,N) + collision->crossSection.segment(0,N))*0.5;
				// evaluate departure cell
				int lmin = std::floor(collision->threshold/energyStep)-1;
				gasPower["attachmentIne"] -= factor * collision->target->density * energyStep * 
											 (eedf.segment(lmin+1, N-lmin-1) * Eigen::pow(energyCell.segment(lmin+1, N-lmin-1), 2) * cellCrossSection.segment(lmin+1, N-lmin-1)).sum();
				continue;
			}
			// switch to lower case because of esthetical reasons
			boost::algorithm::to_lower(collType);
			// evaluate cross section at cell positions
			Eigen::ArrayXd cellCrossSection = (collision->crossSection.segment(1,N) + collision->crossSection.segment(0,N))*0.5;
			// evaluate departure cell
			int lmin = std::floor(collision->threshold/energyStep)-1;
			// add contribution to the power due to the inelastic collisions
			gasPower[collType + "Ine"] += -factor * collision->target->density * energyStep * energyNode[lmin+1] * 
											 (eedf.segment(lmin+1, N-lmin-1) * energyCell.segment(lmin+1, N-lmin-1) * cellCrossSection.segment(lmin+1, N-lmin-1)).sum();							 
			// add contribution to the power due to the superelastic collisions
			if (collision->isReverse){
				double statWeightRatio = collision->target->statisticalWeight/collision->productArray[0]->statisticalWeight;
				gasPower[collType + "Sup"] += factor * statWeightRatio * collision->productArray[0]->density * energyStep * energyNode[lmin+1] *
											  (eedf.segment(0,N-lmin-1) * energyCell.segment(lmin+1,N-lmin-1) * cellCrossSection.segment(lmin+1,N-lmin-1) ).sum();							  									  
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

	// evaluate the electric field by ensuring perfect power balance
	powerMap["field"] = -(powerMap["elasticNet"] + powerMap["carNet"] + powerMap["inelastic"] + powerMap["superelastic"]);
	powerMap["balance"] = 0;
	powerMap["relativeBalance"] = 0;
	workCond->update({"reducedElecField"}, {std::sqrt(powerMap["field"]/auxPowerField)/1e-21});

	// evaluate reference power (totalGain)
	std::vector<double> powerValues({powerMap["field"], powerMap["elasticGain"], powerMap["elasticLoss"], powerMap["carGain"], powerMap["carLoss"],
								powerMap["excitationSup"], powerMap["excitationIne"], powerMap["vibrationalSup"], powerMap["vibrationalIne"], 
								powerMap["rotationalSup"], powerMap["rotationalIne"]});
	double totalGain = 0;
	for (auto powerValue: powerValues){
		if (powerValue > 0){
			totalGain += powerValue;
		}
	}	
	powerMap["reference"] = totalGain;

	// save power balance information
	power.Map = powerMap;
}

void PrescribedEedf::evaluateSwarmParameters(){

	int N = energyGrid->cellNumber;

	// initialize transport parameters structure
	swarmParam["redDiffCoeff"] = 0; swarmParam["redMobCoeff"] = 0; swarmParam["redTownsendCoeff"] = 0; swarmParam["redAttCoeff"] = 0;
	swarmParam["meanEnergy"] = 0; swarmParam["characEnergy"] = 0; swarmParam["Te"] = 0; swarmParam["driftVelocity"] = 0;

	// evaluate reduced diffusion coefficient
	swarmParam["redDiffCoeff"] = (2.0/3.0)*std::sqrt(2.0*Constant::electronCharge/Constant::electronMass) * 
								 (energyGrid->cell * eedf / (totalCrossSection.segment(0,N) + totalCrossSection.segment(1,N) )).sum() * energyGrid->step;

	// evaluate reduced mobility coefficient
	swarmParam["redMobCoeff"] = -std::sqrt(2.0*Constant::electronCharge/Constant::electronMass)/3.0 * 
								(energyGrid->node.segment(1,N-1) * (eedf.segment(1,N-1)-eedf.segment(0,N-1)) / totalCrossSection.segment(1,N-1)).sum();

	// evaluate drift velocity
	swarmParam["driftVelocity"] = swarmParam["redMobCoeff"] * workCond->reducedElecFieldSI;

	// evaluate reduced Townsed coefficient 
	double totalIonRateCoeff = 0;
	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){
			if (collision->type == "Ionization"){
				totalIonRateCoeff += collision->target->density * collision->ineRateCoeff;
				//break;
			}
		}
	}
	swarmParam["redTownsendCoeff"] = totalIonRateCoeff / swarmParam["driftVelocity"];

	// evaluate reduced attachment coefficient
	double totalAttRateCoeff = 0;
	for (auto& gas: gasArray){
		for (auto& collision: gas->collisionArray){
			if (collision->type == "Attachment"){
				totalAttRateCoeff += collision->target->density * collision->ineRateCoeff;
				//break;
			}
		}
	}
	swarmParam["redAttCoeff"] = totalAttRateCoeff / swarmParam["driftVelocity"];

	// evaluate mean energy
	swarmParam["meanEnergy"] = (Eigen::pow(energyGrid->cell, 1.5) * eedf * energyGrid->step).sum();

	// evaluate characteristic energy
	swarmParam["characEnergy"] = swarmParam["redDiffCoeff"]/swarmParam["redMobCoeff"];

	// evaluate electron temperature
	swarmParam["Te"] = workCond->electronTemperature; // swarmParam.meanEnergy*2/3;
}

void PrescribedEedf::evaluateRateCoeff(){

	// clear rate coeff vectors
	rateCoeffAll.clear();
	rateCoeffExtra.clear();

	RateCoeffStruct auxRateCoeff;

	// evaluate rate coefficient for all collisions
	for (auto& gas: gasArray){

		// collisions taken into account for solving the eedf
		for (auto& collision: gas->collisionArray){
			auxRateCoeff.collID = collision->ID;
			collision->evaluateRateCoeff(eedf);
			auxRateCoeff.ineRate = collision->ineRateCoeff;
			auxRateCoeff.supRate = collision->supRateCoeff;
			auxRateCoeff.collDescription = collision->description();
			rateCoeffAll.push_back(auxRateCoeff);
		}

		// extra collisions
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

