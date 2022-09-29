#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Grid.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <cstdio>
#include <gsl/gsl_spline.h>
#include <boost/signals2.hpp>

Collision::Collision(std::string type1, EedfState* &target1, std::vector<EedfState*> &productArray1, std::vector<double> &productStoiCoeff1, bool isReverse1, double threshold1, std::vector<std::vector<double>> rawCrossSection1, bool isExtra1){
	static int lastID = -1;
	++lastID;
	ID = lastID;
	type = type1;
	target = target1;
	productArray = productArray1;
	productStoiCoeff = productStoiCoeff1;
    isExtra = isExtra1;
	isReverse = isReverse1;
	threshold = threshold1;
	rawCrossSection = rawCrossSection1;
    // check if the energy ('x') column is strictly increasing
    int rawCrossSectionSize = rawCrossSection[0].size();
    for (int i = 0; i < rawCrossSectionSize-1; ++i){
        if (rawCrossSection[0][i] >= rawCrossSection[0][i+1]){
            Message::error(std::string("The energy column of the cross section corresponding to the collision shown below is not strictly increasing.\n") + description() + "\nPlease fix this in the LXCat files.");
        }
    }
	// update the collision arrays
    if (isExtra){
    	target->gas->collisionArrayExtra.push_back(this);
    	target->collisionArrayExtra.push_back(this);
    }
    else{
        target->gas->collisionArray.push_back(this);
        target->collisionArray.push_back(this);
        target->isTarget = true;       
    }
	if (isReverse){
		if (productArray.size() == 1 && productStoiCoeff[0] == 1){
            if (isExtra){
    			productArray[0]->collisionArrayExtra.push_back(this);
            }
            else{
                productArray[0]->collisionArray.push_back(this);
                productArray[0]->isTarget = true;
            }
		}
		else{
			Message::error(std::string("Error while creating collision '") + description() + "'. Klein-Rosseland microreversibility relation valid only for binary collisions.\n");
		}
	}
    if (type == "Effective" || type == "Elastic"){
        if (target->type != "ele"){
            Message::error(std::string("Found ''") + type + "'' collision with ''" + target->name + "'' as target [" + description() + "].\n" + type + " collisions are only allowed for electronic states. Please check LXCat files.\n");
        }
    }
}

std::string Collision::descriptionShort(){
    std::string collisionStr = "e+" + target->name;
    if (isReverse){
        collisionStr += "<->";
    }
    else{
        collisionStr += "->";
    }
    if (type == "Ionization"){
        collisionStr += "e+e+";
    }
    else if(type != "Attachment"){
        collisionStr += "e+";
    }
    double numberOfProducts = productArray.size();
    EedfState* product;
    std::string stoiCoeffStr;
    for (int i = 0; i<numberOfProducts; ++i){
        product = productArray[i];
        if (productStoiCoeff[i] > 1){
            stoiCoeffStr = std::to_string((int) productStoiCoeff[i]);
        }
        else{
            stoiCoeffStr = "";
        }
        collisionStr += stoiCoeffStr + product->name;
        if (i<numberOfProducts-1){
            collisionStr += "+";
        }
    }
    return collisionStr;
}

std::string Collision::description(){
	std::string collisionStr = descriptionShort();
	collisionStr += "," + type;
	return collisionStr;
}

void Collision::adjustCrossSection(Grid* energyGrid1){
	// adjustCrossSection interpolates the values of the cross section of a certain collision and stores the adjusted values in the 'crossSection'
	// property of the object 'collision'. The function also stores the energy grid used for the interpolation in the 'energyGrid' property.

    // check that 'Elastic' cross section is available for the maximum value of the energy grid
    if (type == "Elastic"){
        double endEnergyGrid = energyGrid1->node[energyGrid1->node.size()-1];
        double endEnergyCross = rawCrossSection[0].back();
        if (endEnergyGrid > endEnergyCross){
            Message::error(std::string("''") + description() + "'' cross section data is not available for the maximum energy of the simulation (" + std::to_string(endEnergyGrid) + " eV).\nSimulation is not reliable under these conditions.");
        }
    }

	// save energy grid
	if (!energyGrid){
		energyGrid = energyGrid1;
        // connect the signal 'updatedMaxEnergy1Signal' to reAdjustCrossSection, in order to update the cross section each time the max energy is changed
        // idea taken from https://stackoverflow.com/questions/3047381/boost-signals-and-passing-class-method
        energyGrid->updatedMaxEnergy1Signal.connect(boost::bind(&Collision::reAdjustCrossSection, this));
	}

	// save interpolated values of the cross section in crossSection property
	crossSection = interpolatedCrossSection(energyGrid->node);
}

Eigen::ArrayXd Collision::interpolatedCrossSection(Eigen::ArrayXd energyValues){
    // interpolatedCrossSection returns the values of the cross section of a
    // "collision" interpolated at certain "energyValues". The interpolation
    // is performed with the gsl_spline_eval function, with linear
    // interpolation and an extrapolation value of 0.0.
    int energyValuesSize = energyValues.size();
    Eigen::ArrayXd crossSection1 = Eigen::ArrayXd::Zero(energyValuesSize);


    int minIndex = -1;
    if (type == "Effective" || type == "Elastic"){
    	minIndex = 0;
    }
    else{
    	for (int i = 0; i < energyValuesSize; ++i){
    		if (energyValues[i] > threshold){
    			minIndex = i;
    			break;
    		}
    	}
    	if (minIndex == -1){
    		return crossSection1;
    	}
    }

    int rawCrossSectionSize = rawCrossSection[0].size();
    double firstEnergy = rawCrossSection[0][0];
    double lastEnergy = rawCrossSection[0][rawCrossSectionSize-1];

    // initialize and allocate the gsl objects
    gsl_spline* interpolation = gsl_spline_alloc(gsl_interp_linear, rawCrossSectionSize);
    gsl_spline_init(interpolation, rawCrossSection[0].data(), rawCrossSection[1].data(), rawCrossSectionSize);
    gsl_interp_accel* accelerator = gsl_interp_accel_alloc();

    for (int i = minIndex; i < energyValuesSize; ++i){
        double energyValue = energyValues[i];
    	if (energyValue >= firstEnergy && energyValue <= lastEnergy){
    		crossSection1[i] = gsl_spline_eval(interpolation, energyValue, accelerator);
    	}
    	else{
    		// set the cross section to zero, when the energy values are higher than the largest energy of the cross section
    		crossSection1[i] = 0;
    	}
    }

    return crossSection1;
}

Eigen::ArrayXd Collision::superElasticCrossSection(Eigen::ArrayXd energyValue){
    // superElasticCrossSection returns the super elastic cross section of a certain collision at the specified 
    // energyValues by using the Klein-Rosseland microreversibility relation:
    // 
    //    sigma_{f,i}(u) = (g_i/g_f)*(1+threshold/u)*sigma_{i,f}(u+threshold)
    //
    // Note: in case energyValue is not provided the function uses the energy values of the collision raw cross section.

    // error checking
    if (!isReverse){
    	Message::error("Collision '" + description() + " is not defined as bidirectional.\n");
    }
    else if (target->statisticalWeight == Constant::NON_DEF){
    	Message::error(std::string("The statistical weight of the state '") + target->name + "' is not defined.\nThe super elastic cross section of the collision '" + description() + " cannot be evaluated.\n");
    }
    else if (productArray[0]->statisticalWeight == Constant::NON_DEF){
    	Message::error(std::string("The statistical weight of the state '") + productArray[0]->name + "' is not defined.\nThe super elastic cross section of the collision '" + description() + " cannot be evaluated.\n");
    }

    // when no energyValue is given, the function returns the super elastic cross section at:
    // energyGrid.node (if available) or rawCrossSection energy values

    if (energyValue.size() == 0){
    	if (!energyGrid){
            double rawCrossSectionSize = rawCrossSection[0].size();
            energyValue.resize(rawCrossSectionSize);
            for (int i = 0; i < rawCrossSectionSize; ++i){
                energyValue[i] = rawCrossSection[0][i];
            }
    	}
    	else{
    		energyValue = energyGrid->node;
    	}
    }

    // initialization of the super elastic cross section
    int energyValueSize = energyValue.size();
    Eigen::ArrayXd crossSection1 = Eigen::ArrayXd::Zero(energyValueSize);
    int minIndex;

    // Klein-Rosseland microreversibility relation (super elastic cross section is zero at zero energy)
    if (energyValue[0] == 0){
    	if (energyValueSize == 1){
    		return crossSection1;
    	}
    	else{
    		minIndex = 1;
    	}
    }
    else{
    	minIndex = 0;
    }

    // Array with the values of the interpolation: constant * sigma_{i,f}(u+threshold)
   	Eigen::ArrayXd interpolation = interpolatedCrossSection(energyValue + threshold) * (target->statisticalWeight/productArray[0]->statisticalWeight);

    for (int i = minIndex; i < energyValueSize; ++i){
    	crossSection1[i] =  (1.0+threshold/energyValue[i]) * interpolation[i];
    }

    return crossSection1;
}

void Collision::evaluateRateCoeff(Eigen::ArrayXd eedf){
    // evaluateRateCoeff evaluates the rate coefficient(s) of a certain collision provided an eedf. The function
    // returns the value(s) of the rate coefficient for the inelastic (and superelastic) collision. The value(s) 
    // is (are) also stored in the collision properties.

    // evaluate auxiliary variables
    double factor = std::sqrt( 2.0*Constant::electronCharge/Constant::electronMass );
    int lmin = std::floor(threshold / energyGrid->step);

    int crossSectionSize = crossSection.size();
    int partialSize = crossSectionSize-1-lmin;

    // preventing subarrays with negative size
    if (partialSize <= 0){
        ineRateCoeff = 0;
        if (isReverse){
            supRateCoeff = 0;
        }
        return;
    }
    
    Eigen::ArrayXd cellCrossSection = (crossSection.segment(lmin,partialSize) + crossSection.segment(lmin+1,partialSize))/2.0;

    Eigen::ArrayXd aux = cellCrossSection * energyGrid->cell.segment(lmin,partialSize);

    // evaluate inelastic rate coefficient
    ineRateCoeff = factor * (aux*eedf.segment(lmin,partialSize)).sum() * energyGrid->step;

    // evaluate superelastic rate coefficient (if collision is reverse)
    if (isReverse){
    	double targetStatWeight = target->statisticalWeight;
    	double productStatWeight = productArray[0]->statisticalWeight;
	    if (targetStatWeight == Constant::NON_DEF){
	    	Message::error(std::string("Unable to find '") + target->name + "' statistical weight for the evaluation of superelastic rate coefficient of " + description() + "\n");
	    }
	    else if (productStatWeight == Constant::NON_DEF){
	    	Message::error(std::string("Unable to find '") + productArray[0]->name + "' statistical weight for the evaluation of superelastic rate coefficient of " + description() + "\n");
	    }
    	double statWeightRatio = targetStatWeight/productStatWeight;

    	supRateCoeff = factor*statWeightRatio * (aux*eedf.segment(0,partialSize)).sum() * energyGrid->step;
    }

}

void Collision::reAdjustCrossSection(){
    // reAdjustCrossSection readjusts the interpolated cross section of the collision whenever the 'updatedMaxEnergy1' event is trigered by the energyGrid object

    // check that 'Elastic' cross section is available for the maximum value of the energy grid
    if (type == "Elastic"){
        double endEnergyGrid = energyGrid->node[energyGrid->node.size()-1];
        double endEnergyCross = rawCrossSection[0].back();
        if (endEnergyGrid > endEnergyCross){
            Message::error(std::string("''") + description() + "'' cross section data is not available for the maximum energy of the simulation (" + std::to_string(endEnergyGrid) + " eV).\nSimulation is not reliable under these conditions.");
        }
    }

    // save new interpolated values of the cross section in crossSection property
    crossSection = interpolatedCrossSection(energyGrid->node);
}

int Collision::add(std::string type, EedfState* &target, std::vector<EedfState*> &productArray, std::vector<double> productStoiCoeff, bool isReverse, double threshold, std::vector<std::vector<double>> rawCrossSection, std::vector<Collision*> &collisionArray, bool isExtra){
	int collisionID = find(type, target, productArray, productStoiCoeff, isReverse, threshold);
	if (collisionID == -1){
		collisionArray.push_back(new Collision(type, target, productArray, productStoiCoeff, isReverse, threshold, rawCrossSection, isExtra));
		collisionID = collisionArray.back()->ID;
	}
	else{
		Message::warning(std::string("Warning! Avoiding duplicated electron impact collision:\n\t ") + collisionArray[collisionID]->description());
	}
	return collisionID;
}

int Collision::find(std::string type, EedfState* target, std::vector<EedfState*> productArray, std::vector<double> productStoiCoeff, bool isReverse, double threshold){
	bool equalProducts;
	int numProducts, collisionID;
	for (auto collision: target->collisionArray){
		if (collision->threshold == threshold && collision->type == type && collision->isReverse == isReverse){
			numProducts = collision->productArray.size();
			if (numProducts == productArray.size()){
				equalProducts = true;
				for (int j = 0; j < numProducts; ++j){
					for (int k = 0; k < numProducts; ++k){
						if (collision->productArray[j] == productArray[k] && collision->productStoiCoeff[j] == productStoiCoeff[k]){
							break;
						}
						else if (k == numProducts-1){
							equalProducts = false;
						}
					}
					if (!equalProducts){
						break;
					}
				}
				if (equalProducts){
					collisionID = collision->ID;
					return collisionID;
				}
			}
		}
	}
	collisionID = -1;
	return collisionID;
}

Collision* Collision::findEquivalent(EedfState* &target, std::vector<EedfState*> &productArray, std::vector<double> productStoiCoeff, bool isReverse){
    // findEquivalent finds an electron collision from the electron kinetics collisionArray with a prescribed directionality, target, productArray and productStoiCoeff
    Collision* eqColl = NULL;
    int numProducts;
    bool equalProducts;

    for (auto& collision: MathFunctions::append(target->collisionArray, target->collisionArrayExtra)){
        if (collision->isReverse == isReverse){
            numProducts = collision->productArray.size();
            if (numProducts == productArray.size()){
                equalProducts = true;
                for (int j = 0; j < numProducts; ++j){
                    for (int k = 0; k < numProducts; ++k){
                        if (collision->productArray[j] == productArray[k] && collision->productStoiCoeff[j] == productStoiCoeff[k]){
                            break;
                        }
                        else if (k == numProducts-1){
                            equalProducts = false;
                        }
                    }
                    if (!equalProducts){
                        break;
                    }
                }
                if (equalProducts){
                    eqColl = collision;
                    break;
                }
            }
        }
    }
    return eqColl;
}

void Collision::adjustToEnergyGrid(Grid* energyGrid, std::vector<Collision*> &collisionArray){
	// adjustic cross sections of each collision to the energy grid
	for (auto& collision: collisionArray){
		collision->adjustCrossSection(energyGrid);
	}
}