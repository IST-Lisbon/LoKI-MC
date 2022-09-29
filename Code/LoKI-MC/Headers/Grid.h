#ifndef __Grid__
#define __Grid__
 
#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <vector>
#include <string>
#include <boost/signals2.hpp>

class Grid{
  	//   Class that defines a grid and its interactions with other objects. A
    //   grid object stores three properties: the values of the grid at node 
    //   positions, the values of the grid at cell positions (middle point) and 
    //   the step of the grid. 
    //
    //   Node values -> |     |     |     |     |     |     |     |
    //   Grid        -> o-----o-----o-----o-----o-----o-----o-----o
    //   Cell values ->    |     |     |     |     |     |     |

public:
	// ----- class attributes -----

    Eigen::ArrayXd node;              	// values of the grid at node positions
    Eigen::ArrayXd cell;              	// values of the grid at cell positions (i.e. between two consecutive nodes)
    double step;              			// difference between the values at two consecutive nodes
    double cellNumber;        			// number of cells in the energy grid
    bool isSmart = false;        		// smart properties of the energy grid (deactivated by default)
    double minEedfDecay = Constant::NON_DEF;      // minimum number of decades of decay for the EEDF
    double maxEedfDecay = Constant::NON_DEF;      // maximum number of decades of decay for the EEDF
    double updateFactor = Constant::NON_DEF;      // percentage factor to update the maximum value of the energy grid

    boost::signals2::signal<void ()> updatedMaxEnergy1Signal;				  
    boost::signals2::signal<void ()> updatedMaxEnergy2Signal;

    // ----- class methods -----

    Grid();
    void updateMaxValue(double);
    Eigen::ArrayXd adjustToGrid(std::vector<std::vector<double>>, double, std::string);
    double integrate(Eigen::ArrayXd);
    double maxwellianRateCoeff(Eigen::ArrayXd, double);
};
#endif