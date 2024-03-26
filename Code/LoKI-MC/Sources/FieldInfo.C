#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Message.h"
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>

std::vector<FieldInfo*> FieldInfo::fieldArray; // vector with the pointers of all FieldInfo objects parsed from the setup file

FieldInfo::FieldInfo(int lineNumber1, std::string rawLine1){
	static int lastID = -1;
	++lastID;
	// assign the ID of the general field array
	ID = lastID;
	// assign the raw line from the setup
	rawLine = Parse::removeComments(rawLine1);
	// assign the line number of the setup file
	lineNumber = lineNumber1;
	// evaluate name and value
	evaluateNameAndValue();
	// evaluate level and parent of the field
	evaluateLevelAndParent();
}

void FieldInfo::evaluateNameAndValue(){
	std::vector<std::string> tokens = Parse::tokenizeSpaces(rawLine);
	if (tokens.empty()){
		Message::error("Error! Trying to evaluate a FieldInfo in an empty line.");
	}
	// check if it is an enumeration
	if (tokens[0] == "-"){
		isEnumeration = true;
	}

	if (isEnumeration){
		// case: "- name"
		if (tokens.size() == 2){
			name = tokens[1];
		}
		// case: "- name = value"
		else if (tokens.size() == 4){
			name = tokens[1];
			value = tokens[3];
		}
		else{
			Message::error(std::string("Could not parse line ") + std::to_string(lineNumber) + " of the setup file:\n" + rawLine);
		}
	}
	else{
		int foundTwoPoints = tokens[0].find(":");
		if (foundTwoPoints != std::string::npos){
			// case: "name:"
			if (tokens.size() == 1){
				name = tokens[0].substr(0, foundTwoPoints);
			}
			// case: "name: value"
			else if (tokens.size() == 2){
				name = tokens[0].substr(0, foundTwoPoints);
				value = tokens[1];
			}
			else{
				Message::error(std::string("Could not parse line ") + std::to_string(lineNumber) + " of the setup file:\n" + rawLine);
			}
		}
		else{
			Message::error(std::string("Could not parse line ") + std::to_string(lineNumber) + " of the setup file:\n" + rawLine);
		}
	}
}

void FieldInfo::evaluateLevelAndParent(){
	// count the number of space characters at the left
	numberOfSpaces = 0;
	for (auto charac: rawLine){
		if (isspace(charac)){
			++numberOfSpaces;
		}
		else{
			break;
		}
	}
	// find the parent and evaluate the level
	level = 0;

	for (int searchID = ID-1; searchID != -1; --searchID){
		FieldInfo* possibleParent = fieldArray[searchID];
		// find the first field which has no value assigned and with a smaller spacing
		if (numberOfSpaces > possibleParent->numberOfSpaces && possibleParent->value.empty() && !possibleParent->isEnumeration){
			// assign the parent
			parent = possibleParent;
			// assign the level
			level = possibleParent->level + 1;
			// find the siblings
			siblingArray = parent->childArray;
			// update the siblingArray of the siblings
			for (auto& sibling: siblingArray){
				sibling->siblingArray.push_back(this);
			}
			// update the childArray of the parent
			parent->childArray.push_back(this);
			break;
		}
	}
}

void FieldInfo::printSetupInfo(){
	for (auto& field: fieldArray){
		// evaluate the number of spaces
		std::string lineToPrint("");
		for (int i = 0; i < field->level; ++i){
			lineToPrint += "  ";
		}
		// evaluate the rest
		if (field->isEnumeration){
			if (field->value.empty()){
				lineToPrint += "- " + field->name;
			}
			else{
				lineToPrint += "- " + field->name + " = " + field->value;
			}
		}
		else{
			lineToPrint += field->name + ": " + field->value;
		}
		std::cout<<lineToPrint<<std::endl;
	}
}

void FieldInfo::saveSetupInfo(std::string fileName){
	std::ofstream setupFile(fileName);
	if (setupFile.is_open()){
		for (auto& field: fieldArray){
			// evaluate the number of spaces
			std::string lineToPrint("");
			for (int i = 0; i < field->level; ++i){
				lineToPrint += "  ";
			}
			// evaluate the rest
			if (field->isEnumeration){
				if (field->value.empty()){
					lineToPrint += "- " + field->name;
				}
				else{
					lineToPrint += "- " + field->name + " = " + field->value;
				}
			}
			else{
				lineToPrint += field->name + ": " + field->value;
			}
			setupFile<<lineToPrint<<std::endl;
		}
	}
	else{
		Message::error(std::string("The file '") + fileName + "' could not be created. Please check the correspondent directory.");
	}
}

std::map<std::string,double> FieldInfo::getFieldNumericMap(std::string chainFieldString){
	// getFieldMap obtains a numeric map correspondent to a given chain of field names. For example, if we want a map regarding the fraction of the gases used
	// in the electron kinetics, we should use fractionMap = getFieldMap("electronKinetics.gasProperties.fraction", fieldArray)

	std::map<std::string,double> fieldMap;
	// get the field correspondent to this chain
	FieldInfo* field = getField(chainFieldString);
	std::string fileName;

	if (field){
		// if the map is to be evaluated through the children
		if (field->value.empty()){
			for (auto& childField: field->childArray){
				if (childField->value.empty()){
					fileName = std::string("Input/") + childField->name;
					Parse::modifyPropertyMap(fileName, fieldMap);
				}
				else{
					fieldMap[childField->name] = Parse::str2value(childField->value);
				}
			}
		}
		// if the map is to be evaluated through the field value (file name)
		else{
			fileName = std::string("Input/") + field->value;
			Parse::modifyPropertyMap(fileName, fieldMap);
		}
	}

	return fieldMap;
}

std::map<std::string,double> FieldInfo::getFieldNumericMap(FieldInfo* field){

	std::map<std::string,double> fieldMap;
	std::string fileName;

	// if the map is to be evaluated through the children
	if (field->value.empty()){
		for (auto& childField: field->childArray){
			if (childField->value.empty()){
				fileName = std::string("Input/") + childField->name;
				Parse::modifyPropertyMap(fileName, fieldMap);
			}
			else{
				fieldMap[childField->name] = Parse::str2value(childField->value);
			}
		}
	}
	// if the map is to be evaluated through the field value (file name)
	else{
		fileName = std::string("Input/") + field->value;
		Parse::modifyPropertyMap(fileName, fieldMap);
	}

	return fieldMap;
}

std::map<std::string,std::string> FieldInfo::getFieldMap(std::string chainFieldString){
	// getFieldMap obtains a map correspondent to a given chain of field names. For example, if we want a map regarding the fraction of the gases used
	// in the electron kinetics, we should use fractionMap = getFieldMap("electronKinetics.gasProperties.fraction", fieldArray)

	std::map<std::string,std::string> fieldMap;
	// get the field correspondent to this chain
	FieldInfo* field = getField(chainFieldString);
	std::string fileName;

	if (field){
		// if the map is to be evaluated through the children
		if (field->value.empty()){
			for (auto& childField: field->childArray){
				if (childField->value.empty()){
					fileName = std::string("Input/") + childField->name;
					Parse::modifyPropertyMap(fileName, fieldMap);
				}
				else{
					fieldMap[childField->name] = childField->value;
				}
			}
		}
		// if the map is to be evaluated through the field value (file name)
		else{
			fileName = std::string("Input/") + field->value;
			Parse::modifyPropertyMap(fileName, fieldMap);
		}
	}

	return fieldMap;
}

std::map<std::string,std::string> FieldInfo::getFieldMap(FieldInfo* field){

	std::map<std::string,std::string> fieldMap;
	std::string fileName;

	// if the map is to be evaluated through the children
	if (field->value.empty()){
		for (auto& childField: field->childArray){
			if (childField->value.empty()){
				fileName = std::string("Input/") + childField->name;
				Parse::modifyPropertyMap(fileName, fieldMap);
			}
			else{
				fieldMap[childField->name] = childField->value;
			}
		}
	}
	// if the map is to be evaluated through the field value (file name)
	else{
		fileName = std::string("Input/") + field->value;
		Parse::modifyPropertyMap(fileName, fieldMap);
	}

	return fieldMap;
}

std::vector<std::string> FieldInfo::getFieldChildNames(std::string chainFieldString){
	// getFieldChildNames obtains the children of the field corresponding to the last name of the chain

	std::vector<std::string> childNames;
	// get the field correspondent to this chain
	FieldInfo* field = getField(chainFieldString);

	if (field){
		if (field->value.empty()){
			for (auto& childField: field->childArray){
				childNames.push_back(childField->name);
			}
		}
		else{
			childNames.push_back(field->value);
		}
	}

	return childNames;
}

void FieldInfo::changeFieldValue(std::string chainFieldString, std::string newValue){
	// changeFieldValue changes the value of the field corresponding to the last name of the chain

	FieldInfo* field = getField(chainFieldString);

	if (field){
		field->value = newValue;
	}
}

std::string FieldInfo::getFieldValue(std::string chainFieldString){
	// getFieldValue obtains the value of the field corresponding to the last name of the chain

	FieldInfo* field = getField(chainFieldString);

	if (!field){
		return std::string("");
	}
	else{
		return field->value;
	}
}

double FieldInfo::getFieldNumericValue(std::string chainFieldString){
	// getFieldValue obtains the numerical value of the field corresponding to the last name of the chain

	FieldInfo* field = getField(chainFieldString);

	if (!field){
		return 0;
	}
	else{
		return Parse::str2value(field->value);
	}
}

FieldInfo* FieldInfo::getField(std::string chainFieldString){
	// getField obtains the field corresponding to the last name of the chain

	std::vector<std::string> chainFieldNames = Parse::tokenizeCharacters(chainFieldString, (char*)".");

	FieldInfo* currentField = getLevel0Field(chainFieldNames[0]);
	FieldInfo* nextField;

	for (int i = 1; i < chainFieldNames.size() && currentField; ++i){
		nextField = currentField->getChild(chainFieldNames[i]);
		currentField = nextField;
	}
	return currentField;
}

FieldInfo* FieldInfo::getLevel0Field(std::string fieldName){
	FieldInfo* foundField = NULL;
	for (auto& field: fieldArray){
		if (field->name == fieldName && field->level == 0){
			foundField = field;
		}
	}
	return foundField;
}

FieldInfo* FieldInfo::getChild(std::string childName){
	FieldInfo* childPointer = NULL;
	for (auto& child: childArray){
		if (child->name == childName){
			childPointer = child;
			break;
		}
	}
	return childPointer;
}