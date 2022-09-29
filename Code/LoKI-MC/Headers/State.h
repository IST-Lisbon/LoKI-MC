#ifndef __State__
#define __State__

#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/StatePropertyFunctions.h"
#include <iostream>
#include <vector>
#include <string>

class Collision;
class WorkingConditions;

template <class StateType>
using StatePropertyFunctionPointer = void (*) (std::vector<StateType*> &stateArray, std::vector<double> &argumentArray, std::vector<std::string> &argumentStrArray, WorkingConditions* workCond);

template <class GasType, class StateType>
class State{
	// Class that stores the information of a certain state of a gas. The information stored here 
	// (in particular the population of the state) is to be used in a Boltzmann solver to obtain the EEDF.

public:
	// ----- class attributes -----
	// missing heatCapacityFunc and thermalConductivityFunc
	int ID=-1;  						     // ID that identifies the state in the state array of a gas
	std::string type; 						 // type of state. Options: "ele", "vib", "rot" and "ion"
	std::string ionCharg; 					 // ionic level of the state
	std::string eleLevel; 					 // electronic level of the state
	std::string vibLevel;					 // vibrational level of the state
	std::string rotLevel; 					 // rotational level of the state
	std::string name; 						 // name of the state gas(ionCharg,eleLevel,vibLevel,rotLevel)

	GasType* gas = NULL; 				     // handle to the gas that the state belongs to (initialized in the subclass)
	StateType* parent = NULL; 			     // handle to the parent state (e.g. the parent of N2(X,v=0) is -> N2(X)). Electronic states have no parent 
	std::vector<StateType*> siblingArray; 	 // handle array to the states with the same parent (e.g. siblings of N2(X) are N2(A), N2(B),...) 
	std::vector<StateType*> childArray;      // handle array to the child states (e.g. children of N2(X) are N2(X,v=0), N2(v=1),...)

	double energy = Constant::NON_DEF; 			 // energy of the state
	StatePropertyFunctionPointer<StateType> energyFunc = NULL;
	std::vector<double> energyParams;
	std::vector<std::string> energyStrParams;

	double statisticalWeight = Constant::NON_DEF; // statistical weight of the state
	StatePropertyFunctionPointer<StateType> statisticalWeightFunc = NULL;
	std::vector<double> statisticalWeightParams;
	std::vector<std::string> statisticalWeightStrParams;

	double population = 0; 				          // population of the state relative to its siblings
	StatePropertyFunctionPointer<StateType> populationFunc = NULL;
	std::vector<double> populationParams;
	std::vector<std::string> populationStrParams;
	bool populationIsUpdated = false;

	double density = 0; 				          // absolute density (relative to the gas density)	

	double reducedDiffCoeff = Constant::NON_DEF;  // reduced free diffusion coefficient
	StatePropertyFunctionPointer<StateType> reducedDiffCoeffFunc = NULL;
	std::vector<double> reducedDiffCoeffParams;
	std::vector<std::string> reducedDiffCoeffStrParams;

	double reducedMobility = Constant::NON_DEF;   // reduced free mobility coefficient
	StatePropertyFunctionPointer<StateType> reducedMobilityFunc = NULL;
	std::vector<double> reducedMobilityParams;
	std::vector<std::string> reducedMobilityStrParams;

	// concerning electron kinetics
	bool isTarget = false;                          // true if the state is the target of a collision, false otherwise
	std::vector<Collision*> collisionArray;         // handle array to the collisions of which the target is state
	std::vector<Collision*> collisionArrayExtra;	

	State(){;}
	// note that 'addFamily' is at 'chemState' and 'eedfState'

	void evaluateDensity(){
		if (type=="rot"){
			density = population * parent->population * parent->parent->population * gas->fraction;
		}
		else if (type=="vib"){
			density = population * parent->population * gas->fraction;
		}
		else if (type=="ele"){
			density = population * gas->fraction;
		}
		else if (type=="ion"){
			density = population * gas->fraction;
		}
	}
	void disp(){
		std::cout<<"ID: "<<ID<<std::endl;
		std::cout<<"name: "<<name<<std::endl;
		if (energy != Constant::NON_DEF){
			std::cout<<"energy: "<<energy<<std::endl;
		}
		if (statisticalWeight != Constant::NON_DEF){
			std::cout<<"statisticalWeight: "<<statisticalWeight<<std::endl;
		}
		if (population != Constant::NON_DEF){
			std::cout<<"population: "<<population<<std::endl;
		}
		if (density != Constant::NON_DEF){
			std::cout<<"density: "<<density<<std::endl;
		}
		if (reducedDiffCoeff != Constant::NON_DEF){
			std::cout<<"reducedDiffCoeff: "<<reducedDiffCoeff<<std::endl;
		}
		if (reducedMobility != Constant::NON_DEF){
			std::cout<<"reducedMobility: "<<reducedMobility<<std::endl;
		}
		if ( parent ){
			std::cout<<"parent: "<<parent->name<<std::endl;
		}
		if (!siblingArray.empty()){
			std::cout<<"siblings:"<<std::endl;
			for(auto& sibling: siblingArray){
				std::cout<<"\t "<<sibling->name<<std::endl;
			}
		}
		if (!childArray.empty()){
			std::cout<<"children:"<<std::endl;
			for (auto& child: childArray){
				std::cout<<"\t "<<child->name<<std::endl;
			}
		}		
	}
	void evaluateName(){
		std::string stateName = gas->name+"(";
		if (!ionCharg.empty()){
			stateName = stateName + ionCharg + ",";
		}
		stateName = stateName + eleLevel;
		if (!vibLevel.empty()){
			stateName = stateName + ",v=" + vibLevel;
			if(!rotLevel.empty()){
				stateName = stateName + ",J=" + rotLevel;
			}
		}
		stateName = stateName + ")";
		name = stateName;		
	}
	double mass(){
		// mass is an alias function to obtain the mass of a certain state which is the mass of the parent gas
		return gas->mass;
	}
	static int add(GasType* &gas1, std::string ionCharg1, std::string eleLevel1, std::string vibLevel1, std::string rotLevel1, std::vector<StateType*> &stateArray1){
		// adds a new state and its position if it is not already defined
		int stateID = StateType::find(gas1->name,ionCharg1,eleLevel1,vibLevel1,rotLevel1,stateArray1)[0];
		if (stateID != -1){
			return stateID;
		}
		else{
			stateArray1.push_back( new StateType(gas1, ionCharg1, eleLevel1, vibLevel1, rotLevel1) );
			return stateArray1.back()->ID;
		}		
	}
	static std::vector<int> find(std::string gasName1,std::string ionCharg1,std::string eleLevel1,std::string vibLevel1,std::string rotLevel1,std::vector<StateType*> stateArray1){
		std::vector<int> stateID;

		if (eleLevel1 == "*"){ // * is the wildCardChar
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && "ele" == state->type){
					stateID.push_back(state->ID);
					for (auto& state1 : state->siblingArray){
						stateID.push_back(state1->ID);
					}
					break;
				}
			}
		}
		else if (vibLevel1 == "*"){
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && eleLevel1 == state->eleLevel && "ele" == state->type){
					for (auto& state1 : state->childArray){
						stateID.push_back(state1->ID);
					}
					break;				
				}
			}
		}
		else if (rotLevel1 == "*"){
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && eleLevel1 == state->eleLevel && vibLevel1 == state->vibLevel && "vib" == state->type){
					for (auto& state1 : state->childArray){
						stateID.push_back(state1->ID);
					}
					break;
				}
			}
		}
		else{
			for (auto& state: stateArray1){
				if (gasName1==state->gas->name && ionCharg1==state->ionCharg && eleLevel1==state->eleLevel && vibLevel1==state->vibLevel && rotLevel1==state->rotLevel){
					stateID.push_back(state->ID);
					break;
				}
			}
		}
		if (stateID.empty()){
			stateID.push_back(-1);
		}
		return stateID;		
	}
	static std::vector<StateType*> findPointer(std::string gasName1,std::string ionCharg1,std::string eleLevel1,std::string vibLevel1,std::string rotLevel1,std::vector<StateType*> stateArray1){
		std::vector<StateType*> statePointer;

		if (eleLevel1 == "*"){ // * is the wildCardChar
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && "ele" == state->type){
					statePointer.push_back(state);
					for (auto& state1 : state->siblingArray){
						statePointer.push_back(state1);
					}
					break;
				}
			}
		}
		else if (vibLevel1 == "*"){
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && eleLevel1 == state->eleLevel && "ele" == state->type){
					for (auto& state1 : state->childArray){
						statePointer.push_back(state1);
					}
					break;				
				}
			}
		}
		else if (rotLevel1 == "*"){
			for (auto& state : stateArray1){
				if (gasName1 == state->gas->name && eleLevel1 == state->eleLevel && vibLevel1 == state->vibLevel && "vib" == state->type){
					for (auto& state1 : state->childArray){
						statePointer.push_back(state1);
					}
					break;
				}
			}
		}
		else{
			for (auto& state: stateArray1){
				if (gasName1==state->gas->name && ionCharg1==state->ionCharg && eleLevel1==state->eleLevel && vibLevel1==state->vibLevel && rotLevel1==state->rotLevel){
					statePointer.push_back(state);
					break;
				}
			}
		}

		return statePointer;
	}
	static void fixOrphanStates(std::vector<StateType*> &stateArray1){
		// fixOrphanStates looks for orphan states and create the needed parent states
		
		for (int i=0; i<stateArray1.size(); ++i){
			if (stateArray1[i]->type=="rot" &&  !stateArray1[i]->parent){
				add(stateArray1[i]->gas,stateArray1[i]->ionCharg,stateArray1[i]->eleLevel,stateArray1[i]->vibLevel,std::string(""),stateArray1) ;
			}
		}
		for (int i=0; i<stateArray1.size(); ++i){
			if (stateArray1[i]->type=="vib" &&  !stateArray1[i]->parent){
				add(stateArray1[i]->gas,stateArray1[i]->ionCharg,stateArray1[i]->eleLevel,std::string(""),std::string(""),stateArray1) ;
			}
		}
	}
	static bool compareArrays(std::vector<StateType*> stateArray1, std::vector<StateType*> stateArray2){
	
		int numberOfStates = stateArray1.size();
		if (numberOfStates == stateArray2.size()){
			if (numberOfStates==0){
				return true;
			}
			else if (numberOfStates==1){
				if (stateArray1[0]->ID == stateArray2[0]->ID){
					return true;
				}
			}
			else{
				for (int i=0; i<numberOfStates; ++i){
					for (int j=0; j<numberOfStates; ++j){
						if (stateArray1[i]->ID==stateArray2[j]->ID){
							break;
						}
						else if (j==numberOfStates-1){
							return false;
						}
					}
				}
				return true;
			}
		}

		return false;		
	}
};

#endif