#ifndef __WorkingConditions__
#define __WorkingConditions__

#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <boost/signals2.hpp>

class WorkingConditions{
    //  WorkingConditions class contains the information about the working conditions of the simulation
    //  This class contains the information about the working conditins of the simulation. It is also the class
    //  responsible for handling changes in the working conditions by broadcasting any change in them to the
    //  corresponding objects (e.g. when the reduced electric field is updated the boltzmann object is notified, so it can
    //  reevaluate the field operator).

public:
	double gasPressure = Constant::NON_DEF;           // in Pa
	double gasTemperature = Constant::NON_DEF;        // in K
    double wallTemperature = Constant::NON_DEF;       // in K
	double gasDensity = Constant::NON_DEF;            // in m-3
	double electronDensity = Constant::NON_DEF;       // in m-3
	double electronTemperature = Constant::NON_DEF;   // in eV
	double chamberLength = Constant::NON_DEF;         // in m
	double chamberRadius = Constant::NON_DEF;         // in m
    double reducedElecField = Constant::NON_DEF;      // in Td
    double reducedElecFieldSI = Constant::NON_DEF;    // in V.cm2

    std::vector<double> electronTemperatureArray;
    std::vector<double> reducedElecFieldArray;

    // signals that are equivalent to the matlab events
    boost::signals2::signal<void ()> updatedGasPressure;
    boost::signals2::signal<void ()> updatedGasDensity;  
    boost::signals2::signal<void ()> updatedElectronDensity;
    boost::signals2::signal<void ()> updatedElectronTemperature;
    boost::signals2::signal<void ()> updatedReducedElecField;

    std::string variableCondition; // can be reducedElecField", "electronTemperature"

    WorkingConditions(std::map<std::string,std::string> workingConditionsMap){

        gasPressure = Parse::str2value(workingConditionsMap["gasPressure"]);
        gasTemperature = Parse::str2value(workingConditionsMap["gasTemperature"]);
        gasDensity = gasPressure/(Constant::boltzmann*gasTemperature);
        if (workingConditionsMap.find("wallTemperature") != workingConditionsMap.end()){
            wallTemperature = Parse::str2value(workingConditionsMap["wallTemperature"]);
        }
        if (workingConditionsMap.find("electronDensity") != workingConditionsMap.end()){
            electronDensity = Parse::str2value(workingConditionsMap["electronDensity"]);
        }
        if (workingConditionsMap.find("electronTemperature") != workingConditionsMap.end()){
            electronTemperatureArray = Parse::evalVectorExpress(workingConditionsMap["electronTemperature"]);
            electronTemperature = electronTemperatureArray[0];
        }
        if (workingConditionsMap.find("chamberLength") != workingConditionsMap.end()){
            chamberLength = Parse::str2value(workingConditionsMap["chamberLength"]);
        }
        if (workingConditionsMap.find("chamberRadius") != workingConditionsMap.end()){
            chamberRadius = Parse::str2value(workingConditionsMap["chamberRadius"]);
        }
        if (workingConditionsMap.find("reducedElecField") != workingConditionsMap.end()){
            std::string ENString = workingConditionsMap["reducedElecField"];
            reducedElecFieldArray = Parse::evalVectorExpress(workingConditionsMap["reducedElecField"]);
            reducedElecField = reducedElecFieldArray[0];
            reducedElecFieldSI = reducedElecField*1e-21;
        }

        std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");
        if(eedfType == "boltzmannMC"){
            if (reducedElecFieldArray.empty()){
                Message::error("Error in the configuration of the working conditions. When choosing 'boltzmannMC' eedfType, the 'reducedElecField' condition must be defined.");
            }
            variableCondition = "reducedElecField";
        }        
        else if(eedfType == "prescribedEedf"){
            if (electronTemperatureArray.empty()){
                Message::error("Error in the configuration of the working conditions. When choosing 'prescribedEedf' eedfType, the 'electronTemperature' condition must be defined.");
            }
            variableCondition = "electronTemperature";
        }

        checkPositiveConditions();
    }

    void checkPositiveConditions(){
        if ((gasPressure < 0 && gasPressure != Constant::NON_DEF) || (gasTemperature < 0 && gasTemperature != Constant::NON_DEF) || (wallTemperature < 0 && wallTemperature != Constant::NON_DEF) ||
            (electronDensity < 0 && electronDensity != Constant::NON_DEF) || (chamberLength < 0 && chamberLength != Constant::NON_DEF) || (chamberRadius < 0 && chamberRadius != Constant::NON_DEF)){
            Message::error("Error in the configuration of the working conditions. The working conditions must be all non-negative.");
        }
        for (auto elecTemp: electronTemperatureArray){
            if (elecTemp <= 0){
                Message::error("Error in the configuration of the working conditions. 'electronTemperature' must be positive.");
            }
        }
        for (auto redElecField: reducedElecFieldArray){
            if (redElecField < 0){
                Message::error("Error in the configuration of the working conditions. 'reducedElecField' must be non-negative.");
            }
        }
    }

    void update(std::vector<std::string> propertiesToUpdate, std::vector<double> newValues){
        for (int i = 0; i < propertiesToUpdate.size(); ++i){
            std::string property = propertiesToUpdate[i];
            if (property == "gasPressure"){
                gasPressure = newValues[i];
                gasDensity = gasPressure/(Constant::boltzmann*gasTemperature);
                updatedGasPressure();
                updatedGasDensity();
            }
            else if (property == "electronDensity"){
                electronDensity = newValues[i];
                updatedElectronDensity();
            }
            else if (property == "electronTemperature"){
                electronTemperature = newValues[i];
                updatedElectronTemperature();
            }
            else if (property == "reducedElecField"){
                reducedElecField = newValues[i];
                reducedElecFieldSI = reducedElecField*1e-21;
                updatedReducedElecField();
            }
        }
    }

    double getValue(std::string property){
        double value;
        if (property == "gasPressure"){
            value = gasPressure;
        }
        else if (property == "gasTemperature"){
            value = gasTemperature;
        }
        else if (property == "wallTemperature"){
            value = wallTemperature;
        }
        else if (property == "gasDensity"){
            value = gasDensity;
        }
        else if (property == "electronDensity"){
            value = electronDensity;
        }
        else if (property == "electronTemperature"){
            value = electronTemperature;
        }
        else if (property == "chamberLength"){
            value = chamberLength;
        }
        else if (property == "chamberRadius"){
            value = chamberRadius;
        }
        else if (property == "reducedElecField"){
            value = reducedElecField;
        }
        else if (property == "reducedElecFieldSI"){
            value = reducedElecFieldSI;
        }
        else{
            Message::error(std::string("Trying to access to the value of the working condition '") + property + "', which is not defined in the code.");
        }
        return value;
    }
};

#endif