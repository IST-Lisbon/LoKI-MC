#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Message.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "External/MathParser/parser.h"
#include <map>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/algorithm/string.hpp>
#include <filesystem>

void Parse::setupFile(std::string fileName){
	// 'setupFile' parses the information in the setup file. Note that the indentation of the file is fundamental!

	std::ifstream setupFile("Input/"+fileName);
	std::string setupLine, cleanLine;
	int lineNumber = 0;
	if (setupFile.is_open()){
		while( getline(setupFile,setupLine) ){
			cleanLine = removeComments(setupLine);
			if (cleanLine.empty() || hasJustSpaces(cleanLine)){
				continue;
			}
			else{
				// add the setup field in the static array
				FieldInfo::fieldArray.push_back(new FieldInfo(lineNumber, cleanLine));
			}
			++lineNumber;
		}
	}
	else{
		Message::error(std::string("The setup file '") + fileName + "' could not be opened. Put it inside the 'Input' folder.");
	}
	setupFile.close();
}

std::vector<LXCatEntryStruct> Parse::LXCatFiles(std::vector<std::string> LXCatFileNames){
	// LXCatFiles reads LXCat files, parse their content and returns a structure array 'LXCatEntryArray' with all the information

	std::string line, descriptionLine;
	std::vector<std::string> tokens;
	int beginThreshold, endThreshold;
	std::vector<double> xCross, yCross;
	std::vector<std::vector<double>> rawCrossSection;
	double threshold = 0;
	std::vector<LXCatEntryStruct> LXCatEntryArray;

    // loop over different "LXCatFiles" files that have to be read
    for (auto fileName : LXCatFileNames){
    	std::ifstream file(std::string("Input/") + fileName);
    	if (file.is_open()){
			while( getline(file,line) ){
				if (line.empty()){
					continue;
				}
				else{
					// find a process
					if (line.find("PROCESS:") != std::string::npos){

						std::string processLine = line;
						xCross.clear();
						yCross.clear();
						rawCrossSection.clear();
						threshold = 0;

						// read the threshold (it must be in the line immediately after the process)
						getline(file, line);
						beginThreshold = line.find("E =");
						endThreshold = line.find("eV");
						if (beginThreshold != std::string::npos && endThreshold != std::string::npos){
							threshold = str2value( line.substr(beginThreshold + 3, endThreshold-beginThreshold-3) );
						}
						else if (processLine.find("Elastic") == std::string::npos && processLine.find("Effective") == std::string::npos){
							Message::error(std::string("Error! Could not find a threshold for the following process at file '") + fileName + "':\n" + processLine + "\nThe threshold must be in the line immediately after the definition of the process. For example:\nPARAM.: E = someValue eV\nPlease add it and run the code again.");
						}
				
						// get the line with the description of the LXCat collision (the one that LoKI reads)
						getline(file, line);
						descriptionLine = line;
						
						// get the raw section
						// find the beggining of the cross section
						while(getline(file, line)){
							if (line.find("-----") != std::string::npos){
								break;
							}
						}
						// read until the end of the cross section
						while (getline(file, line) && line.find("-----") == std::string::npos){
							tokens = tokenizeSpaces(line);
							if (tokens.size() != 2){
								Message::error(std::string("Error when reading the cross section of the LXCat collision described in the following line:\n") + descriptionLine +"Check the LXCat file '" + fileName + "'");
							}
							xCross.push_back(str2value(tokens[0]));
							yCross.push_back(str2value(tokens[1]));
						}
						rawCrossSection.push_back(xCross);
						rawCrossSection.push_back(yCross);
						addLXCatEntry(descriptionLine, threshold, rawCrossSection, LXCatEntryArray);
					}
				}
			}
    	}
		else{
			Message::error(std::string("The file '") + fileName + "' could not be opened.");
		}
		file.close();
    }
    return LXCatEntryArray;
}

void Parse::addLXCatEntry(std::string descriptionLine, double threshold, std::vector<std::vector<double>> rawCrossSection, std::vector<LXCatEntryStruct> &LXCatEntryArray){ 
	// addLXCatEntry analyses the information of a particular LXCat entry and adds it to the structure array LXCatEntryArray

	LXCatEntryStruct LXCatEntry;

	// eliminate the spaces from the string
	std::string spacedLine = descriptionLine;
	descriptionLine = removeSpaces(spacedLine);

	// start by finding the description part: collision, type
	int descriptionBegin = descriptionLine.find_first_of('[')+1, descriptionEnd = descriptionLine.find_last_of(']')-1;
	std::string description = descriptionLine.substr(descriptionBegin, descriptionEnd-descriptionBegin+1);

	// find the type of collision
	int collisionEnd = description.find_last_of(",");
	std::string type = description.substr(collisionEnd+1);
	if (type != "Elastic" && type != "Effective" && type != "Attachment" && type != "Ionization" &&
		type != "Excitation" && type != "Vibrational" && type != "Rotational"){
		Message::error(std::string("Error in the parsing of one LXCat file. Invalid type of collision in the following line: \n") + spacedLine + 
					   "\nValid types are: Elastic, Effective, Attachment, Ionization, Excitation, Vibrational and Rotational.");
	}
	LXCatEntry.type = type;

	// get the collision string
	std::string collisionPart = description.substr(0,collisionEnd);

	// parsing the collision part

	std::string collisionLeftPart, collisionRightPart;
	int separatorFound;
	if (collisionPart.find("<->") != std::string::npos){
		LXCatEntry.isReverse = true;
		separatorFound = collisionPart.find("<->");
		collisionLeftPart = collisionPart.substr(0, separatorFound);
		collisionRightPart = collisionPart.substr(separatorFound+3, std::string::npos);
	}
	else if (collisionPart.find("->") != std::string::npos){
		separatorFound = collisionPart.find("->");
		collisionLeftPart = collisionPart.substr(0, separatorFound);
		collisionRightPart = collisionPart.substr(separatorFound+2, std::string::npos);	
	}
	else{
		Message::error(std::string("Error in the parsing of one LXCat file. Check the direction of the collision presented in the following line: \n") + spacedLine);
	}


	double stoiCoeff = 1;
	std::string stateString(""), character;
	bool leftParenthesesFound = false;

	// parsing the left part of the collision: target

	int collisionLeftPartSize = collisionLeftPart.size();
	for (int i = 0; i < collisionLeftPartSize; ++i){
		character = collisionLeftPart[i];
		if ((character == "+" && !leftParenthesesFound) || i == collisionLeftPartSize-1){
			leftParenthesesFound = false;

			if (i == collisionLeftPartSize-1){
				stateString += character;
			}

			// get the stoiCoeff and eliminate the correspondent part
			stoiCoeff = getStoiCoeff(stateString); 

			// parsing electrons orregular reactants (states)
			if (stateString == "e" || stateString == "E"){
				LXCatEntry.reactantElectrons += stoiCoeff;
			}
			else{
				LXCatEntry.targetStoiCoeff.push_back(stoiCoeff);
				LXCatEntry.target.push_back(getRawState(stateString));
			}
			stateString = "";		
		}
		else if (character == "("){
			leftParenthesesFound = true;
			stateString += character;
		}
		else if (character == ")"){
			leftParenthesesFound = false;
			stateString += character;
		}
		else{
			stateString += character;
		}
	}

	if (LXCatEntry.target.empty()){
		Message::error(std::string("Could not find a target in the collision presented in the following LXCat line:\n") + spacedLine + "\n Please check your LXCat files.");
	}

	// parsing the right part of the collision: products

	int collisionRightPartSize = collisionRightPart.size();
	stateString = "";
	for (int i = 0; i < collisionRightPartSize; ++i){
		character = collisionRightPart[i];
		if ((character == "+" && !leftParenthesesFound) || i == collisionRightPartSize-1){
			leftParenthesesFound = false;

			if (i == collisionRightPartSize-1){
				stateString += character;
			}
			// get the stoiCoeff and eliminate the correspondent part
			stoiCoeff = getStoiCoeff(stateString); 

			// parsing electrons or regular products (states)
			if (stateString == "e" || stateString == "E"){
				LXCatEntry.productElectrons += stoiCoeff;
			}
			else{
				LXCatEntry.productStoiCoeff.push_back(stoiCoeff);
				LXCatEntry.productArray.push_back(getRawState(stateString));
			}
			stateString = "";		
		}
		else if (character == "("){
			leftParenthesesFound = true;
			stateString += character;
		}
		else if (character == ")"){
			leftParenthesesFound = false;
			stateString += character;
		}
		else{
			stateString += character;
		}
	}

	removeDuplicatedStates(LXCatEntry.productArray, LXCatEntry.productStoiCoeff);

	if (LXCatEntry.target.empty()){
		Message::error(std::string("Could not find a product in the collision presented in the following LXCat line:\n") + spacedLine + "\n Please check your LXCat files.");
	}

	LXCatEntry.threshold = threshold;
	LXCatEntry.rawCrossSection = rawCrossSection;


	LXCatEntryArray.push_back(LXCatEntry);
}

void Parse::removeDuplicatedStates(std::vector<rawStateStruct>& stateArray, std::vector<double>& stoiCoeff){

	if (stateArray.empty()){
		return;
	}

	std::vector<rawStateStruct> newStateArray;
	std::vector<double> newStoiCoeff;

	int numStates = stateArray.size();
	int numNewStates;
	newStateArray.push_back(stateArray[0]);
	newStoiCoeff.push_back(stoiCoeff[0]);


	for (int i=1; i<numStates; ++i){
		numNewStates = newStateArray.size();
		for (int j=0; j<numNewStates; ++j){
			if (stateArray[i].gasName == newStateArray[j].gasName && stateArray[i].ionCharg == newStateArray[j].ionCharg &&
				stateArray[i].eleLevel == newStateArray[j].eleLevel && stateArray[i].vibLevel == newStateArray[j].vibLevel &&
				stateArray[i].rotLevel == newStateArray[j].rotLevel){
				newStoiCoeff[j] +=  stoiCoeff[i];
				break;
			}
			if (j == (numNewStates-1) ){
				newStateArray.push_back(stateArray[i]);
				newStoiCoeff.push_back(stoiCoeff[i]);
			}
		}
	}


	stateArray = newStateArray;
	stoiCoeff = newStoiCoeff;
}

void Parse::findCatalysts(std::vector<rawStateStruct> &reactantArray, std::vector<double> &reactantStoiCoeff, std::vector<rawStateStruct> &productArray, std::vector<double> &productStoiCoeff, std::vector<rawStateStruct> &catalystArray, std::vector<double> &catalystStoiCoeff){ 

	if (reactantArray.empty() || productArray.empty()){
		return;
	}

	int i = 0, j;
	while(i<reactantArray.size()){
		j = 0;
		while(j<productArray.size()){
			if (reactantArray[i].gasName == productArray[j].gasName && reactantArray[i].ionCharg == productArray[j].ionCharg &&
				reactantArray[i].eleLevel == productArray[j].eleLevel && reactantArray[i].vibLevel == productArray[j].vibLevel &&
				reactantArray[i].rotLevel == productArray[j].rotLevel){
				catalystArray.push_back(reactantArray[i]);
				catalystStoiCoeff.push_back(reactantStoiCoeff[i]);
				if (reactantStoiCoeff[i] < productStoiCoeff[j]){
					productStoiCoeff[j] = productStoiCoeff[j] - reactantStoiCoeff[i];
					reactantArray.erase(reactantArray.begin()+i);
					reactantStoiCoeff.erase(reactantStoiCoeff.begin()+i);
					--i;
					break;
				}
				else if (reactantStoiCoeff[i] > productStoiCoeff[j]){
					catalystStoiCoeff.back() = productStoiCoeff[j];
					reactantStoiCoeff[i] = reactantStoiCoeff[i] - productStoiCoeff[j];
					productArray.erase(productArray.begin()+j);
					productStoiCoeff.erase(productStoiCoeff.begin()+j);
					break;
				}
				else{
					reactantArray.erase(reactantArray.begin()+i);
					reactantStoiCoeff.erase(reactantStoiCoeff.begin()+i);
					productArray.erase(productArray.begin()+j);
					productStoiCoeff.erase(productStoiCoeff.begin()+j);
					--i;
					break;
				}
			}
			++j;
		}
		++i;
	}
}

double Parse::getStoiCoeff(std::string &stateToken){
	// gets the stoichiometric coefficient from a string with the stoiCoeff + stateName.
	// once the coefficient is parsed, the correspondent part is eliminated, leaving only the stateName
	int countCoeffChar = 0;
	std::string coeffString;

	for (auto car : stateToken){
		if (!isalpha(car)){
			++countCoeffChar;
			coeffString.append(1,car);
		}
		else{
			break;
		}
	}
	for (int i=0; i<countCoeffChar; ++i){
		stateToken.erase(0,1);
	}
	if (countCoeffChar == 0){
		return 1;
	}
	else{
		return str2value(coeffString);
	}
}


rawStateStruct Parse::getRawState(std::string stateName){
	// returns a structure with the information parsed from the state name: gasName, eleLevel, vibLevel and rotLevel (and vibRange)
		rawStateStruct rawState;

		std::vector<std::string> tokens = tokenizeCharacters(stateName,(char*)"()");
		if (tokens.size() == 1){
			Message::error(std::string("The states must have the electronic state defined between parentheses!\nPlease check '") + stateName + "'");
		}
		rawState.gasName = tokens[0];
		std::vector<std::string> stateFields = tokenizeCharacters(tokens[1],(char*)",");
		

		switch ( stateFields.size() ){
			case 1 :
				rawState.eleLevel = stateFields[0];
				break;
			case 2 :
				if (stateFields[1].find("v=") != std::string::npos){
					rawState.eleLevel = stateFields[0];
					rawState.vibLevel = tokenizeCharacters(stateFields[1],(char*)"=")[1];
					rawState.vibRange = "v";
				}
				else if (stateFields[1].find("w=") != std::string::npos){
					rawState.eleLevel = stateFields[0];
					rawState.vibLevel = tokenizeCharacters(stateFields[1],(char*)"=")[1];
					rawState.vibRange = "w";
				}
				else{
					rawState.ionCharg = stateFields[0];
					rawState.eleLevel = stateFields[1];
				}
				break;
			case 3 :
				if (stateFields[2].find("J=") != std::string::npos){
					rawState.eleLevel = stateFields[0];
					if (stateFields[1].find("v=") != std::string::npos){
						rawState.vibLevel = tokenizeCharacters(stateFields[1],(char*)"=")[1];
						rawState.vibRange = "v";
					}
					else if (stateFields[1].find("w=") != std::string::npos){
						rawState.vibLevel = tokenizeCharacters(stateFields[1],(char*)"=")[1];
						rawState.vibRange = "w";
					}
					rawState.rotLevel = tokenizeCharacters(stateFields[2],(char*)"=")[1];
				}
				else{
					rawState.ionCharg = stateFields[0];
					rawState.eleLevel = stateFields[1];
					rawState.vibLevel = stateFields[2];
				}
				break;
			case 4 :
				rawState.ionCharg = stateFields[0];
				rawState.eleLevel = stateFields[1];
				if (stateFields[2].find("v=") != std::string::npos){
					rawState.vibLevel = tokenizeCharacters(stateFields[2],(char*)"=")[1];
					rawState.vibRange = "v";
				}
				else if (stateFields[2].find("w=") != std::string::npos){
					rawState.vibLevel = tokenizeCharacters(stateFields[2],(char*)"=")[1];
					rawState.vibRange = "w";
				}
				rawState.rotLevel = tokenizeCharacters(stateFields[3],(char*)"=")[1];
				break;
			default:
				Message::error(std::string("Error! Check the state '") + stateName + "'.");
		}

		return rawState;
}

void Parse::modifyPropertyMap(std::string fileName, std::map<std::string,std::string> &propertyMap){
	// modifies a property map using an input file. E.g. Input/Databases/masses.txt
	std::ifstream file(fileName);
	std::string line, cleanLine;
	std::vector<std::string> tokens;

	if (file.is_open()){

		while( getline(file,line) ){

			cleanLine = removeComments(line);
			tokens = tokenizeSpaces(cleanLine);
			if (tokens.size()==0){
				continue;
			}
			else if (tokens.size()!=2){
				Message::error(std::string("Error in the parsing of the following property file: ") + fileName + ". Check the following line: \n" + line);
			}
			else{
				propertyMap[tokens[0]] = tokens[1];
			}
		}	
	}

	else{
		Message::error(std::string("The file ") + fileName + " could not be opened.");
	}
	file.close();
}

void Parse::modifyPropertyMap(std::string fileName, std::map<std::string,double> &propertyMap){
	// modifies a property map using an input file. E.g. Input/Databases/masses.txt
	std::ifstream file(fileName);
	std::string line, cleanLine;
	std::vector<std::string> tokens;

	if (file.is_open()){

		while( getline(file,line) ){

			cleanLine = removeComments(line);
			tokens = tokenizeSpaces(cleanLine);
			if (tokens.size()==0){
				continue;
			}
			else if (tokens.size()!=2){
				Message::error(std::string("Error in the parsing of the following property file: ") + fileName + ". Check the following line: \n" + line);
			}
			else{
				propertyMap[tokens[0]] = Parse::str2value(tokens[1]);
			}
		}	
	}

	else{
		Message::error(std::string("The file ") + fileName + " could not be opened.");
	}
	file.close();
}

std::map<std::string,std::string> Parse::createPropertyMap(std::string fileName){
	// creates a property map using an input file. E.g. Input/Databases/masses.txt
	std::ifstream file(fileName);
	std::string line, cleanLine;
	std::map<std::string,std::string> propertyMap;
	std::vector<std::string> tokens;

	if (file.is_open()){

		while( getline(file,line) ){

			cleanLine = removeComments(line);
			tokens = tokenizeSpaces(cleanLine);
			if (tokens.size()==0){
				continue;
			}
			else if (tokens.size()!=2){
				Message::error(std::string("Error in the parsing of the following property file: ") + fileName + ". Check the following line: \n" + line);
			}
			else{
				propertyMap[tokens[0]] = tokens[1];
			}
		}	
	}

	else{
		Message::error(std::string("The file ") + fileName + " could not be opened.");
	}
	file.close();
	return propertyMap;
}

std::vector<std::string> Parse::getMapKeys(std::map<std::string,std::map<std::string,double>> const& input_map) {
  std::vector<std::string> retKey;
  for (auto const& element : input_map) {
    retKey.push_back(element.first);
  }
  return retKey;
}

std::vector<std::string> Parse::getMapKeys(std::map<std::string,std::string> const& input_map) {
  std::vector<std::string> retKey;
  for (auto const& element : input_map) {
    retKey.push_back(element.first);
  }
  return retKey;
}

std::vector<std::string> Parse::getMapKeys(std::map<std::string,double> const& input_map) {
  std::vector<std::string> retKey;
  for (auto const& element : input_map) {
    retKey.push_back(element.first);
  }
  return retKey;
}

std::vector<double> Parse::evalVectorExpress(std::string expression){
	std::vector<double> vec;
	std::string parametersString;
	std::vector<std::string> tokens;
	std::vector<double> parameters;
	if (expression.find("linspace") != std::string::npos){
		parametersString = expression.substr(expression.find("(")+1, expression.find(")")-expression.find("(")-1);
		tokens = tokenizeCharacters(parametersString, (char*)",");
		if (tokens.size() != 3){
			Message::error(std::string("Invalid use of 'linspace' in the following expression: ") + expression);
		}
		for (auto& token: tokens){
			parameters.push_back( str2value(token) );
		}
		return linSpace(parameters[0], parameters[1], parameters[2]);
	}
	else if (expression.find("logspace") != std::string::npos){
		parametersString = expression.substr(expression.find("(")+1, expression.find(")")-expression.find("(")-1);
		tokens = tokenizeCharacters(parametersString, (char*)",");
		if (tokens.size() != 3){
			Message::error(std::string("Invalid use of 'logspace' in the following expression: ") + expression);
		}
		for (auto& token: tokens){
			parameters.push_back( str2value(token) );
		}
		return logSpace(parameters[0], parameters[1], parameters[2]);
	}
	else if (count(expression.begin(), expression.end(), ':') == 1){
		tokens = tokenizeCharacters(expression, (char*)":");
		int ini = str2value(tokens[0]);
		int fin = str2value(tokens[1]);
		if (ini > fin){
			Message::error(std::string("Error while using 'evalVectorExpress'. When an expression like 'initial:final' is to be evaluated, 'final' cannot be smaller than 'initial'. Check this expression: ") + expression);
		}
		for (int i = ini; i <= fin; ++i){
			vec.push_back(i);
		}
		return vec;
	}
	else if (count(expression.begin(), expression.end(), ':') == 2){
		tokens = tokenizeCharacters(expression, (char*)":");
		double ini = str2value(tokens[0]);
		double step = str2value(tokens[1]);
		double fin = str2value(tokens[2]);
		if (ini > fin || step <= 0){
			Message::error(std::string("Error while using 'evalVectorExpress'. When an expression like 'initial:step:final' is to be evaluated, 'final' cannot be smaller than 'initial' and 'step' must be larger than 0. Check this expression: ") + expression);
		}
		for (double i = ini; i <= fin; i += step){
			vec.push_back(i);
		}
		return vec;
	}
	else if (expression.find("[") != std::string::npos && expression.find("]") != std::string::npos){
		tokens = tokenizeCharacters(expression, (char*)"[],");
		for (auto token: tokens){
			vec.push_back(str2value(token));
		}
		return vec;
	}
	else{
		return std::vector<double> (1, str2value(expression));
	}
}

std::vector<double> Parse::logSpace(double startExponent, double endExponent, double num){
	// generates num points between decades 10^startExponent and 10^endExponent
	std::vector<double> logSpaced = linSpace(startExponent, endExponent, num);
	for (auto& value: logSpaced){
		value = std::pow(10, value);
	}
	return logSpaced;
}

std::vector<double> Parse::linSpace(double start, double end, double num){
	std::vector<double> linSpaced;

	if (num == 0){
		return linSpaced;
	}
	else if (num == 1){
		linSpaced.push_back(start);
		return linSpaced;
	}

	double delta = (end - start) / (num - 1);

	for(double i = 0; i < num-1; ++i){
		linSpaced.push_back(start + delta * i);
	}
	linSpaced.push_back(end); // I want to ensure that start and end are exactly the same as the input

	return linSpaced;
}


std::vector<int> Parse::evalRange(std::string stringRange){
	std::vector<int> vec;
	if (stringRange.find(":")==std::string::npos){
		vec.push_back(str2value(stringRange));
		return vec;
	}
	std::vector<std::string> tokens = tokenizeCharacters(stringRange, (char*)":");
	int ini = str2value(tokens[0]);
	int fin = str2value(tokens[1]);
	boost::push_back(vec, boost::irange(ini, fin+1));
	return vec;
}

std::vector<double> Parse::str2valueArray(std::vector<std::string> strings){
	std::vector<double> values;
	for (auto str: strings){
		values.push_back(str2value(str));
	}
	return values;
}

template <class type>
double Parse::str2value(type mathExpression){
	// parse a number using an open source library
	static Parser prs;
	char* result;

	if (mathExpIsValid(mathExpression)){
		result = prs.parse( (char*)mathExpression.c_str() );
	}
	else{
		Message::error("Could not parse the mathematical expression: " + std::string(mathExpression) + ". Please check the input files.");
	}

	double value;
    std::sscanf(result, "Ans = %lf", &value);
    if (std::isnan(value)){
    	Message::error("Could not parse the mathematical expression: " + std::string(mathExpression) + ". Please check the input files.");
    }
    return value;
}


template <class type>
bool Parse::mathExpIsValid(type mathExpression){
	// check the parentheses
	int leftParantheses = 0, rightParentheses = 0;
	for (auto charac: mathExpression){
		if (charac == '('){
			++leftParantheses;
		}
		else if (charac == ')'){
			++rightParentheses;
		}
	}
	if (leftParantheses != rightParentheses){
		return false;
	}

	// check if the constants and functions used are valid
	static std::vector<std::string> validConstants = {"E","PI","M_PI"};
	static std::vector<std::string> validFunctions = {"ABS","EXP","SIGN","SQRT","LOG","LOG10","SIN","COS","TAN","ASIN","ACOS","ATAN","FACTORIAL"};
	std::vector<std::string> mathParts = tokenizeCharacters(mathExpression, (char*)"()*+-/ .^0123456789%");
	for (auto part: mathParts){
		bool validPart = false;
		for (auto constant: validConstants){
			if (boost::iequals(constant, part)){
				validPart = true;
				break;
			}
		}
		for (auto function: validFunctions){
			if (boost::iequals(function, part)){
				validPart = true;
				break;
			}
		}
		if (!validPart){
			return false;
		}
	}
	return true;
}

bool Parse::str2bool(std::string expression){
	if (expression == "true" || expression == "True" || expression == "1"){
		return true;
	}
	else{
		return false;
	}
}

bool Parse::isLogical(std::string expression){
	if (expression == "true" || expression == "True" || expression == "1" || expression == "0" || expression == "false" || expression == "False"){
		return true;
	}
	else{
		return false;
	}
}

std::string Parse::removeComments(std::string textLine){
	// returns a string without comments, which are recognized with the '%' character
	size_t idx = textLine.find_first_of('%');
	if (idx != std::string::npos){
		return textLine.substr(0,idx);
	}
	else{
		return textLine;
	}
}

std::string Parse::removeSpaces(std::string spacedLine){
	std::string nonSpacedLine;
	std::vector<std::string> tokens = tokenizeSpaces(spacedLine);
	for (auto token: tokens){
		nonSpacedLine += token;
	}
	return nonSpacedLine;
}

std::vector<std::string> Parse::tokenizeSpaces(std::string textLine){
	
	static char *token;
	static char delimiters[]="   \r\n\v\f\t";
	std::vector<std::string> tokens;
	token = std::strtok( (char*)textLine.c_str(),delimiters);
	while (token != NULL){
		tokens.push_back(std::string(token));
		token = std::strtok(NULL,delimiters);
	}
	return tokens;
}

std::vector<std::string> Parse::tokenizeCharacters(std::string textLine, char* delimiters){
	// tokenize with a desired set of characters
	static char *token;
	std::vector<std::string> tokens;
	token = std::strtok( (char*)textLine.c_str(),delimiters);
	while (token != NULL){
		tokens.push_back(std::string(token));
		token = std::strtok(NULL,delimiters);
	}
	return tokens;
}

bool Parse::hasJustSpaces(std::string textLine){
	bool justSpaces = true;
	for (auto c: textLine){
		if (!std::isspace(c)){
			justSpaces = false;
			break;
		}
	}
	return justSpaces;
}

