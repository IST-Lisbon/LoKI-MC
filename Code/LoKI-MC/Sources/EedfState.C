#include "LoKI-MC/Headers/EedfState.h"
#include "LoKI-MC/Headers/EedfGas.h"
#include "LoKI-MC/Headers/Message.h"
#include <iostream>
#include <string>
#include <vector>

EedfState::EedfState(EedfGas* &gas1,std::string ionCharg1,std::string eleLevel1,std::string vibLevel1,std::string rotLevel1){
	static int lastID=-1;
	++lastID;
	ID = lastID;
	gas = gas1;
	ionCharg = ionCharg1;
	eleLevel = eleLevel1;
	vibLevel = vibLevel1;
	rotLevel = rotLevel1;
	if (ionCharg.empty()){
		if (rotLevel.empty()){
			if (vibLevel.empty()){
				type.assign("ele");
			}
			else{
				type.assign("vib");
			}
		}
		else{
			type.assign("rot");
		}
	}
	else{
		type.assign("ion");
	}
	addFamily();
	evaluateName();
	gas->stateArray.push_back(this);
}

void EedfState::addFamily(){
	// addFamily find all the relatives of state and stores the
    // corresponding information in their properties parent, sibling and child.

	int stateArraySize = gas->stateArray.size();

	if (type=="rot"){
		for (int i=0;i<stateArraySize;++i){
			if (gas->stateArray[i]->type=="vib" && gas->stateArray[i]->eleLevel==eleLevel && gas->stateArray[i]->vibLevel==vibLevel){
				parent = gas->stateArray[i];
				gas->stateArray[i]->childArray.push_back(this);
			}
			else if (gas->stateArray[i]->type=="rot" && gas->stateArray[i]->eleLevel==eleLevel && gas->stateArray[i]->vibLevel==vibLevel){
				siblingArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->siblingArray.push_back(this);
			}			
		}
	}

	else if (type=="vib"){
		for (int i=0;i<stateArraySize;++i){
			if (gas->stateArray[i]->type=="ele" && gas->stateArray[i]->eleLevel==eleLevel){
				parent = gas->stateArray[i];
				gas->stateArray[i]->childArray.push_back(this);
			}
			else if (gas->stateArray[i]->type=="vib" && gas->stateArray[i]->eleLevel==eleLevel){
				siblingArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->siblingArray.push_back(this);
			}
			else if (gas->stateArray[i]->type=="rot" && gas->stateArray[i]->eleLevel==eleLevel && gas->stateArray[i]->vibLevel==vibLevel){
				childArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->parent = this;
			}		
		}
	}

	else if (type=="ele"){
		for (int i=0;i<stateArraySize;++i){
			if (gas->stateArray[i]->type=="ele"){
				siblingArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->siblingArray.push_back(this);
			}
			else if (gas->stateArray[i]->type=="vib" && gas->stateArray[i]->eleLevel==eleLevel){
				childArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->parent = this;
			}		
		}
	}

	else if (type=="ion"){
		for (int i=0;i<stateArraySize;++i){
			if (gas->stateArray[i]->type=="ion"){
				siblingArray.push_back(gas->stateArray[i]);
				gas->stateArray[i]->siblingArray.push_back(this);
			}
		}
	}
}
