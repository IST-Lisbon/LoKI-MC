#ifndef __FieldInfo__
#define __FieldInfo__

#include <string>
#include <vector>
#include <map>

class FieldInfo{
public:
	static std::vector<FieldInfo*> fieldArray; // vector with the pointers of all FieldInfo objects parsed from the setup file

	int ID; 						     // position in the general field array
	int lineNumber;						 // number of the correspondent line in the setup file
	std::string rawLine;				 // line of the setup file
	std::string name;	
	std::string value;
	int numberOfSpaces = 0;				 // number of spaces at the LHS of the field
	int level = -1;						 // field level: 0,1,2,... depending on the spacing
	FieldInfo* parent = NULL;			 
	std::vector<FieldInfo*> siblingArray;
	std::vector<FieldInfo*> childArray;
	bool isEnumeration = false;			 // indicates if it is an enumeration or not. E.g.: "- O2(X) = 5"

	FieldInfo(int lineNumber1, std::string rawLine1);
	void evaluateNameAndValue();
	void evaluateLevelAndParent();

	static void printSetupInfo();
	static void saveSetupInfo(std::string fileName);

	static std::map<std::string,double> getFieldNumericMap(std::string chainFieldString);
	static std::map<std::string,double> getFieldNumericMap(FieldInfo* field);
	static std::map<std::string,std::string> getFieldMap(std::string chainFieldString);
	static std::map<std::string,std::string> getFieldMap(FieldInfo* field);
	static std::vector<std::string> getFieldChildNames(std::string chainFieldString);
	static void changeFieldValue(std::string chainFieldString, std::string newValue);
	static double getFieldNumericValue(std::string chainFieldString);
	static std::string getFieldValue(std::string chainFieldString);
	static FieldInfo* getField(std::string chainFieldString);
	static FieldInfo* getLevel0Field(std::string fieldName);

	FieldInfo* getChild(std::string childName);
};

#endif