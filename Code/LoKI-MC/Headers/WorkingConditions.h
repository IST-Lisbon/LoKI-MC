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
    double reducedMagField = Constant::NON_DEF;       // in Hx
    double reducedMagFieldSI = Constant::NON_DEF;     // In T.m^3
    double elecFieldAngle = Constant::NON_DEF;        // in degrees
    double excitationFrequency = Constant::NON_DEF;   // in Hz

    std::vector<double> electronTemperatureArray;
    std::vector<double> reducedElecFieldArray;
    std::vector<double> reducedMagFieldArray;
    std::vector<double> elecFieldAngleArray;
    std::vector<double> excitationFrequencyArray;

    // signals that are equivalent to the matlab events
    boost::signals2::signal<void ()> updatedGasPressure;
    boost::signals2::signal<void ()> updatedGasDensity;  
    boost::signals2::signal<void ()> updatedElectronDensity;
    boost::signals2::signal<void ()> updatedElectronTemperature;
    boost::signals2::signal<void ()> updatedReducedElecField;
    boost::signals2::signal<void ()> updatedReducedMagField;
    boost::signals2::signal<void ()> updatedElecFieldAngle;
    boost::signals2::signal<void ()> updatedExcitationFrequency;

    std::string variableCondition; // can be "reducedElecField", "reducedMagField", "elecFieldAngle", "excitationFrequency", "electronTemperature"
    bool isCylindricallySymmetric;

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
            reducedElecFieldArray = Parse::evalVectorExpress(ENString);
            reducedElecField = reducedElecFieldArray[0];
            reducedElecFieldSI = reducedElecField*1e-21;
        }
        if (workingConditionsMap.find("reducedMagField") != workingConditionsMap.end()){
            reducedMagFieldArray = Parse::evalVectorExpress(workingConditionsMap["reducedMagField"]);
            reducedMagField = reducedMagFieldArray[0];
            reducedMagFieldSI = reducedMagField*1e-27;
        }
        if (workingConditionsMap.find("elecFieldAngle") != workingConditionsMap.end()){
            elecFieldAngleArray = Parse::evalVectorExpress(workingConditionsMap["elecFieldAngle"]);
            elecFieldAngle = elecFieldAngleArray[0];
        }
        if (workingConditionsMap.find("excitationFrequency") != workingConditionsMap.end()){
            excitationFrequencyArray = Parse::evalVectorExpress(workingConditionsMap["excitationFrequency"]);
            excitationFrequency = excitationFrequencyArray[0];
        }

        // check if there are multiple values for magnetic field, reduced field and elec-field angle (BoltzmannMC)
        std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");
        if (eedfType == "boltzmannMC"){
            if (reducedElecFieldArray.size() > 1 && (reducedMagFieldArray.size() > 1 || elecFieldAngleArray.size() > 1 || excitationFrequencyArray.size() > 1) ){
                Message::error("Error in the configuration of the working conditions. When 'reducedElecField' has multiple values, 'reducedMagField', 'elecFieldAngle' and 'excitationFrequency' must be single-valued.");
            }
            else if (reducedMagFieldArray.size() > 1 && (reducedElecFieldArray.size() > 1 || elecFieldAngleArray.size() > 1 || excitationFrequencyArray.size() > 1) ){
                Message::error("Error in the configuration of the working conditions. When 'reducedMagField' has multiple values, 'reducedElecField', 'elecFieldAngle' and 'excitationFrequency' must be single-valued.");
            }
            else if (elecFieldAngleArray.size() > 1 && (reducedMagFieldArray.size() > 1 || reducedElecFieldArray.size() > 1 || excitationFrequencyArray.size() > 1) ){
                Message::error("Error in the configuration of the working conditions. When 'elecFieldAngle' has multiple values, 'reducedMagField', 'reducedElecField' and 'excitationFrequency' must be single-valued.");
            }
            else if (excitationFrequencyArray.size() > 1 && (reducedMagFieldArray.size() > 1 || elecFieldAngleArray.size() > 1 || reducedElecFieldArray.size() > 1) ){
                Message::error("Error in the configuration of the working conditions. When 'excitationFrequency' has multiple values, 'reducedMagField', 'elecFieldAngle' and 'reducedElecField' must be single-valued.");
            }          
            else if (reducedElecFieldArray.size() > 1){
                variableCondition = "reducedElecField";
            }
            else if (reducedMagFieldArray.size() > 1){
                variableCondition = "reducedMagField";
            }
            else if (elecFieldAngleArray.size() > 1){
                variableCondition = "elecFieldAngle";
            }
            else if (excitationFrequencyArray.size() > 1){
                variableCondition = "excitationFrequency";
            }
            if (reducedElecFieldArray.empty() || reducedMagFieldArray.empty() || elecFieldAngleArray.empty() || excitationFrequencyArray.empty()){
                Message::error("Error in the configuration of the working conditions. When choosing 'boltzmannMC' eedfType, the 'reducedElecField', 'reducedMagField', 'elecFieldAngle' and 'excitationFrequency' must be defined.");
            }
            // if there are no multiple values, define 'reducedElecField' as the variable condition
            // to avoid problems in the output
            if (variableCondition.empty()){
                variableCondition = "reducedElecField";
            }
        }
        // (PrescribedEedf)
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
        for (auto redMagField: reducedMagFieldArray){
            if (redMagField < 0){
                Message::error("Error in the configuration of the working conditions. 'reducedMagField' must be non-negative.");
            }
        }
        isCylindricallySymmetric = true;
        for (auto elecAngle: elecFieldAngleArray){
            if (elecAngle < 0 && elecAngle > 180){
                Message::error("Error in the configuration of the working conditions. 'elecFieldAngle' must be in the interval [0,180] degrees.");
            }
            else if (elecAngle != 180){
                isCylindricallySymmetric = false;
            }
        }
        for (auto excitFreq: excitationFrequencyArray){
            if (excitFreq < 0){
                 Message::error("Error in the configuration of the working conditions. 'excitationFrequency' must be non-negative.");
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
            else if (property == "reducedMagField"){
                reducedMagField = newValues[i];
                reducedMagFieldSI = reducedMagField*1e-27;
                updatedReducedMagField();
            }
            else if (property == "elecFieldAngle"){
                elecFieldAngle = newValues[i];
                updatedElecFieldAngle();
            }
            else if (property == "excitationFrequency"){
                excitationFrequency = newValues[i];
                updatedExcitationFrequency();
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
        else if (property == "reducedMagField"){
            value = reducedMagField;
        }
        else if (property == "reducedMagFieldSI"){
            value = reducedMagFieldSI;
        }
        else if (property == "elecFieldAngle"){
            value = elecFieldAngle;
        }
        else if (property == "excitationFrequency"){
            value = excitationFrequency;
        }
        else{
            Message::error(std::string("Trying to access to the value of the working condition '") + property + "', which is not defined in the code.");
        }
        return value;
    }
};

#endif