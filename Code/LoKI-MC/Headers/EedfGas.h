#ifndef __EedfGas__
#define __EedfGas__

#include "LoKI-MC/Headers/Gas.h"
#include "LoKI-MC/Headers/Constant.h"
#include <string>
#include <vector>

class Collision;
class EedfState;

class EedfGas : public Gas<EedfGas, EedfState>{
public:
	// ----- class attributes -----

	std::vector<Collision*> collisionArray;       // electron collisions array with all the collisions of the gas
	std::vector<Collision*> collisionArrayExtra;
	std::vector<double> effectivePopulations;     // populations to evaluate an elastic cross section from an effective one

	// ----- class methods -----

	EedfGas(std::string gasName);
	void dispCollisions();
	void checkPopulationNorms();
	void checkCARconditions();
	void checkElasticCollisions(std::vector<Collision*> &collisionArray);
	std::vector<std::vector<double>> elasticFromEffectiveCrossSection();
	void renormalizeWithDensities();
};

#endif