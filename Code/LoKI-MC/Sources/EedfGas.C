#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Collision.h"
#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/MathFunctions.h"
#include "LoKI-MC/Headers/AngularDistributionFunctions.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

EedfGas::EedfGas(std::string gasName){
	static int lastID = -1;
	++lastID;
	ID = lastID;
	name = gasName;
}

void EedfGas::dispCollisions(){
	std::printf("%s collisions:\n", name.c_str());
	for (int i = 0; i < collisionArray.size(); ++i){
		std::printf("\t %d.- %s\n", i, collisionArray[i]->description().c_str());
	}
}

void EedfGas::checkPopulationNorms(){
    // checkPopulationNorms checks for the population of the different states of the gas to be properly normalised,
    // i. e. the populations of all sibling states should add to one.

    // avoid dummy gases
    if (fraction == 0){
        return;
    }

    // check norm of electronic/ionic states
    double gasNorm=0,vibNorm=0,rotNorm=0;
    bool electronicStatesToBeChecked = true;
    bool ionStatesToBeChecked = true;
    std::vector<EedfState*> auxArray;


    for (auto& state: stateArray){
        // norm of electronic states
        if (state->type=="ele" && electronicStatesToBeChecked){
            auxArray = state->siblingArray;
            auxArray.insert(auxArray.begin(),state);
            for (auto& eleState: auxArray){
                if (eleState->isTarget && eleState->population!=0){
                    gasNorm = gasNorm + eleState->population;
                    // check norm of vibrational states (if they exist)
                    if (!eleState->childArray.empty()){
                        vibNorm = 0;
                        bool nonDummyVibStates = false;
                        for (auto& vibState: eleState->childArray){
                            if (vibState->isTarget && vibState->population!=0){
                                nonDummyVibStates = true;
                                vibNorm = vibNorm + vibState->population;
                                // check norm of rotational states (if they exist)
                                if (!vibState->childArray.empty()){
                                    rotNorm = 0;
                                    bool nonDummyRotStates = false;
                                    for (auto& rotState: vibState->childArray){
                                        nonDummyRotStates = true;
                                        if (rotState->isTarget && rotState->population!=0){
                                            rotNorm = rotNorm + rotState->population;
                                        }
                                    }
                                    if (nonDummyRotStates && std::abs(rotNorm-1)>10*std::numeric_limits<double>::epsilon()){
                                        Message::error(std::string("Rotational distribution ") + vibState->name.substr(0,vibState->name.size()-1) + ",J=*) is not properly normalized. (Error = " + std::to_string(rotNorm-1) + ")\n");
                                    }
                                }
                            }
                        }
                        if (nonDummyVibStates && std::abs(vibNorm-1)>10*std::numeric_limits<double>::epsilon()){
                            Message::error(std::string("Vibrational distribution ") + eleState->name.substr(0,eleState->name.size()-1) + ",v=*) is not properly normalized. (Error = " + std::to_string(vibNorm-1) + ")\n");
                        }
                    }
                }
            }
            electronicStatesToBeChecked = false;
        }
        if (state->type=="ion" && ionStatesToBeChecked){
            auxArray = state->siblingArray;
            auxArray.insert(auxArray.begin(),state);
            for (auto& ionState: auxArray){
                if (ionState->population!=0){
                    gasNorm = gasNorm + ionState->population;
                }
            }
            ionStatesToBeChecked = false;
        }
    }
    if (std::abs(gasNorm-1)>10*std::numeric_limits<double>::epsilon()){
        Message::error(std::string("Electronic/ionic distribution ") + name + "(*) is not properly normalized. (Error = " + std::to_string(gasNorm-1) + ")\n");
    }       
}

void EedfGas::checkCARconditions(){
    // checkCARconditions checks that the particular gas fulfils the conditions to use the continuous approximation for 
    // rotation: not to have electron collisions defined that causes a rotational transition and to have defined a 
    // value for the electric quadrupole moment and the rotational constant. This function is called whenever the 
    // continuous approximation for rotations is considered for this particular gas.

    // check for the definition of the electric quadrupole moment
    if (electricQuadrupoleMoment == Constant::NON_DEF){
    	Message::error(std::string("A value for the electric quadrupole moment is not found for gas ") + name + ".\nCAR can not be taken into account without it.\nPlease check your setup file.");
    }
    // check for the definition of the rotational constant
    if (rotationalConstant == Constant::NON_DEF){
    	Message::error(std::string("A value for the rotational constant is not found for gas ") + name + ".\nCAR can not be taken into account without it.\nPlease check your setup file.");
    }
    // check for the absence of 'Rotational' collisions
    for (auto& collision: collisionArray){
    	if (collision->type == "Rotational"){
    		Message::error(std::string("Collision '')") + collision->description() + "'' was found while CAR is activated for the gas " + name + ".\nPlease check your setup file.");
    	}
    }
}

void EedfGas::checkElasticCollisions(std::vector<Collision*> &collisionArray){
    // checkElasticCollisions checks for each electronic state (with population) of the gas to have an elastic 
    // collision defined. In case the Elastic collision is not found, the function tries to obtain it from an 
    // Effective one.
    //
    // Note: the function avoid dummy gases (gases without collisions, created for the sake of a pretty output)

    // avoid dummy gases
    if (collisionArray.empty()){
    	return;
    }

    std::vector<EedfState*> auxArray;
    bool elasticCollisionNeedsToBeCreated;
    std::vector<std::vector<double>> rawElasticCrossSection;

    // loop over all the states of the gas
    for (auto& state: stateArray){
    	if (state->type == "ele"){
    		// loop over every electronic states of the gas
    		auxArray = state->siblingArray;
    		auxArray.insert(auxArray.begin(), state);
    		for (auto& eleState: auxArray){
    			// find electronic states that are target (direct or reverse)
    			if (eleState->isTarget){ 
    				// look for an Elastic collision associated with the state
    				elasticCollisionNeedsToBeCreated = true;
    				for (auto& collision: eleState->collisionArray){
    					if (collision->type == "Elastic"){
    						elasticCollisionNeedsToBeCreated = false;
    						break;
    					}
    				}
    				// create Elastic collision in case it is needed
    				if (elasticCollisionNeedsToBeCreated){
    					rawElasticCrossSection = elasticFromEffectiveCrossSection();
    					std::vector<EedfState*> productArray(1, eleState);
    					std::vector<double> productStoiCoeff(1,1);
    					Collision::add(std::string("Elastic"), eleState, productArray, productStoiCoeff, false, 0.0, rawElasticCrossSection, rawElasticCrossSection, collisionArray, false);
                        collisionArray.back()->angularDistributionFunction = AngularDistributionFunctions::functionMap("isotropic",{});
    				}
    			}
    		}
    		break;
    	}
    }
}

std::vector<std::vector<double>> EedfGas::elasticFromEffectiveCrossSection(){
    // elasticFromEffectiveCrossSection evaluates an elastic cross section from an Effective one. In order to do that, 
    // it subtracts from the Effective cross section all the Excitation, Vibrational, Rotational, Ionization and 
    // Attachment cross sections of the gas weighted by the populations of the corresponding targets. 
    //
    // Note: in case that the user do not specify some particular populations for effective cross section in the setup
    // file (exceptional case), they are assumed to be those of equilibrium at 300k (regular case), which is our best 
    // guess for the conditions of measurement of an 'Effective' cross section. That is:
    //     - electronic distribution collapsed to the ground state
    //     - vibrational distribution assumed to be boltzmann at 300k
    //     - rotational distribution assumed to be boltzmann at 300k

    // look for Effective collision
    Collision* effectiveCollision = NULL;
    for (auto& collision: collisionArray){
    	if (collision->type == "Effective"){
    		effectiveCollision = collision;
    		break;
    	}
    }
    if (!effectiveCollision){
    	Message::error(std::string("Gas ''") + name + "'' does not have an ''Effective'' collision defined, so ''Elastic'' collisions cannot be evaluated from it.\nPlease, check the corresponding LXCat file.");
    }

    // initialize Elastic cross section as the Effective one
    std::vector<std::vector<double>> rawElasticCrossSection = effectiveCollision->rawMomTransfCrossSection;

    // find maximum ID of the states of the gas
    int maxID = stateArray[0]->ID;
    for (auto& state: stateArray){
    	if (maxID < state->ID){
    		maxID = state->ID;
    	}
    }

    double norm;

    // calculate populations for the evaluation of the elastic cross section (in case it is needed)
    if (effectivePopulations.empty()){
    	
    	// initialize vector of populations
    	effectivePopulations.assign(maxID+1,0.0);

    	// look for electronic ground state (target of the 'Effective' collision) and assign it a population of 1
    	EedfState* electronicGroundState = effectiveCollision->target;
    	effectivePopulations[electronicGroundState->ID] = 1;
    	EedfState* vibrationalGroundState = NULL;

    	// look for vibrational states of the electronic ground state
    	if (!electronicGroundState->childArray.empty()){
    		// evaluate vibrational distribution (Boltzmann at 300 K)
    		norm = 0;
    		vibrationalGroundState = electronicGroundState->childArray[0];
    		for (auto& state: electronicGroundState->childArray){
	    		if (state->energy == Constant::NON_DEF){
	    			Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of ''Elastic'' cross section of " + state->gas->name + ".\nCheck input file");
	    		}
	    		else if (state->statisticalWeight == Constant::NON_DEF){
	    			Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of ''Elastic'' cross section of " + state->gas->name + ".\nCheck input file");
	    		}
	    		else if (state->energy < vibrationalGroundState->energy){
	    			vibrationalGroundState = state;
	    		}
	    		effectivePopulations[state->ID] = state->statisticalWeight * std::exp(-state->energy/(Constant::boltzmannInEV*300.0));
	    		norm += effectivePopulations[state->ID];
	    	}   		
    	}
    	for (auto& state: electronicGroundState->childArray){
    		effectivePopulations[state->ID] = effectivePopulations[state->ID]/norm;
    	}

    	// look for rotational states of the vibrational ground state of the electronic ground state
        if (vibrationalGroundState){
        	if (!vibrationalGroundState->childArray.empty()){
        		// evaluate rotational distribution (Boltzmann at 300 K)
        		norm = 0;
        		for (auto& state: vibrationalGroundState->childArray){
    	    		if (state->energy == Constant::NON_DEF){
    	    			Message::error(std::string("Unable to find ") + state->name + " energy for the evaluation of ''Elastic'' cross section of " + state->gas->name + ".\nCheck input file");
    	    		}
    	    		else if (state->statisticalWeight == Constant::NON_DEF){
    	    			Message::error(std::string("Unable to find ") + state->name + " statistical weight for the evaluation of ''Elastic'' cross section of " + state->gas->name + ".\nCheck input file");
    	    		}
    	    		effectivePopulations[state->ID] = state->statisticalWeight * std::exp(-state->energy/(Constant::boltzmannInEV*300.0));
    	    		norm += effectivePopulations[state->ID];  	    		
        		}
    	    	for (auto& state: vibrationalGroundState->childArray){
    	    		effectivePopulations[state->ID] = effectivePopulations[vibrationalGroundState->ID]*effectivePopulations[state->ID]/norm;
    	    	}    		
        	}
        }
    }

    while (effectivePopulations.size() < maxID+1){
    	effectivePopulations.push_back(0.0);
    }


    // remove contributions to the effective due to the different collisional mechanisms
    for (auto& collision: collisionArray){
    	if (collision->type == "Effective" || collision->type == "Elastic"){
    		continue;
    	}
    	int rawElasticCrossSectionSize = rawElasticCrossSection[0].size();
    	Eigen::ArrayXd collisionInterpolatedCrossSection = collision->interpolatedCrossSection("momTransf", MathFunctions::vectorToArray(rawElasticCrossSection[0]));
    	for (int i = 0; i < rawElasticCrossSectionSize; ++i){
    		rawElasticCrossSection[1][i] -= effectivePopulations[collision->target->ID] * collisionInterpolatedCrossSection[i];
    	}

    	//remove contributions to the effective due to the super-elastic mechanism (in case it is defined)
    	if (collision->isReverse){
    		collisionInterpolatedCrossSection = collision->superElasticCrossSection("momTransf", MathFunctions::vectorToArray(rawElasticCrossSection[0]));
	    	for (int i = 0; i < rawElasticCrossSectionSize; ++i){
	    		rawElasticCrossSection[1][i] -= effectivePopulations[collision->productArray[0]->ID] * collisionInterpolatedCrossSection[i];
	    	}
    	}
    }

    // check for negative values in the Elastic cross section
    bool negativeValues = false;
    for (auto& value: rawElasticCrossSection[1]){
        if (value < 0){
            negativeValues = true;
            value = 0;
        }
    }
    if (negativeValues){
        Message::warning(std::string("Negative values obtained when evaluating an Elastic cross section from an Effective one (") + name + ").\nNegative values have been clipped to 0 and unreliable results may be obtained.\nPlease, carefully check inputs and outputs of your simulation");
    }
    
    return rawElasticCrossSection;
}