#ifndef __Collision__
#define __Collision__

#include "LoKI-MC/Headers/Constant.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <string>

class EedfState;
class Grid;

class Collision{
	//  Class that stores the information of an electron collision read from an LXCat file. The information stored here (in particular the cross 
	//  section) is to be used in a Boltzmann solver to obtain the EEDF.
	
public:
	// ----- class attributes -----

	int ID = -1;        // ID that identifies the collision in the collision array of the gas
	std::string type;   // type of collision as defined in LXCat file posible values are: 
				        // 'Elastic', 'Effective', 'Excitation', 'Vibrational', 'Rotational', 'Ionization' and 'Attachment'

	EedfState* target;             		     // handle to the target of the collision
    std::vector<EedfState*> productArray;    // handle to the products of the collision
    std::vector<double> productStoiCoeff;  	 // array of stoichiometric coefficients for the products
    
    bool isExtra = false;            		 // true when the collision is not meant to be used in the electron kinetics calculations
    bool isReverse = false;          		 // true when super elastic collision is defined
    double threshold = 0.0;          		 // energy threshold of the collision (eV)
    std::vector<std::vector<double>> rawCrossSection;  // as read from LXCat file, 2 rows (eV m^2)
    
    Grid* energyGrid = NULL;    		     // handle to the energy grid of the simulation
    Eigen::ArrayXd crossSection;          	 // interpolated into the grid, 1 row (m^2)
    
    double ineRateCoeff = Constant::NON_DEF; // inelastic rate coefficient of the collision obtained once the eedf is known (->)
    double supRateCoeff = Constant::NON_DEF; // superelastic rate coefficient of the collision obtained once the eedf is known (if <->)

    // ----- class methods -----

    Collision(std::string, EedfState*&, std::vector<EedfState*>&, std::vector<double>&, bool, double, std::vector<std::vector<double>>, bool);
    void disp();
    std::string descriptionShort();
    std::string description();
    void adjustCrossSection(Grid*);
    Eigen::ArrayXd interpolatedCrossSection(Eigen::ArrayXd);
    Eigen::ArrayXd superElasticCrossSection(Eigen::ArrayXd);
    void evaluateRateCoeff(Eigen::ArrayXd);
    void reAdjustCrossSection();
    
    static int add(std::string, EedfState*&, std::vector<EedfState*>&, std::vector<double>, bool , double, std::vector<std::vector<double>>, std::vector<Collision*>&, bool isExtra);
    static int find(std::string, EedfState*, std::vector<EedfState*>, std::vector<double>, bool, double);
	static Collision* findEquivalent(EedfState*&, std::vector<EedfState*>&, std::vector<double> productStoiCoeff, bool isReverse);
    static void adjustToEnergyGrid(Grid*, std::vector<Collision*>&);
};

#endif