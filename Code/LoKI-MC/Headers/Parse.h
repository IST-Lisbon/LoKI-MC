#ifndef __Parse__
#define __Parse__

#include "LoKI-MC/Headers/Constant.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include <string>
#include <vector>
#include <map>

namespace Parse{

	struct LXCatEntryStruct;		//
	struct rawStateStruct;			//

	void setupFile(std::string);

	std::vector<LXCatEntryStruct> LXCatFiles(std::vector<std::string>);
	void addLXCatEntry(std::string, double, std::vector<std::vector<double>>, std::vector<LXCatEntryStruct>&);

	void removeDuplicatedStates(std::vector<rawStateStruct>&, std::vector<double>&);
	void findCatalysts(std::vector<rawStateStruct>&, std::vector<double>&, std::vector<rawStateStruct>&, std::vector<double>&, std::vector<rawStateStruct>&, std::vector<double>&);
	double getStoiCoeff(std::string&);
	rawStateStruct getRawState(std::string);

	std::vector<std::string> readFileStrings(std::string fileName);
	
	void modifyPropertyMap(std::string, std::map<std::string,std::string>&);
	void modifyPropertyMap(std::string, std::map<std::string,double>&);	
	std::map<std::string,std::string> createPropertyMap(std::string);
	std::vector<std::string> getMapKeys(std::map<std::string,std::string> const&);
	std::vector<std::string> getMapKeys(std::map<std::string,double> const&);
	std::vector<std::string> getMapKeys(std::map<std::string,std::map<std::string,double>> const&);

	std::vector<double> evalVectorExpress(std::string);
	std::vector<double> logSpace(double, double, double);
	std::vector<double> linSpace(double, double, double);
	std::vector<int> evalRange(std::string);
	std::vector<double> str2valueArray(std::vector<std::string>);
	template <class type>
	double str2value(type);
	template <class type>
	bool mathExpIsValid(type);
	bool str2bool(std::string);
	bool isLogical(std::string);

	std::string removeComments(std::string);
	std::string removeSpaces(std::string);
	std::vector<std::string> tokenizeSpaces(std::string);
	std::vector<std::string> tokenizeCharacters(std::string,char*);
	bool hasJustSpaces(std::string);

	struct LXCatEntryStruct{
		std::string type;

		int reactantElectrons = 0;
		int productElectrons = 0;
		std::vector<rawStateStruct> target;
		std::vector<double> targetStoiCoeff;
		std::vector<rawStateStruct> productArray;
		std::vector<double> productStoiCoeff;

		bool isReverse = false;

		double threshold = 0;
		std::vector<std::vector<double>> rawIntegralCrossSection;
		std::vector<std::vector<double>> rawMomTransfCrossSection;
	};

	struct rawStateStruct{
		std::string gasName;
		std::string ionCharg;
		std::string eleLevel;
		std::string vibLevel;
		std::string rotLevel;
		std::string vibRange;
	};	
};

#endif