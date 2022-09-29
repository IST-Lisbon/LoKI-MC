#ifndef __EedfState__
#define __EedfState__

#include "LoKI-MC/Headers/State.h"
#include "LoKI-MC/Headers/Constant.h"
#include <vector>
#include <string>

class EedfGas;
class Collision;

class EedfState: public State<EedfGas,EedfState>{
public:
	// ----- class attributes -----

	bool isTarget = false; 			   // true if the state is the target of a collision, false otherwise
	std::vector<Collision*> collisionArray; // handle array to the collisions of which the target is state
	std::vector<Collision*> collisionArrayExtra;

	// ----- class methods -----
	
	EedfState(EedfGas* &gas1,std::string ionCharg1,std::string eleLevel1,std::string vibLevel1,std::string rotLevel1);
	void addFamily();
};

#endif