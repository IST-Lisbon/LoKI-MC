#include "LoKI-MC/Headers/Parse.h"
#include "LoKI-MC/Headers/Setup.h"
#include "LoKI-MC/Headers/BoltzmannMC.h"
#include "LoKI-MC/Headers/PrescribedEedf.h"
#include "LoKI-MC/Headers/FieldInfo.h"
#include "LoKI-MC/Headers/Message.h"
#include "External/eigen-3.4.0/Eigen/Core"
#include "External/eigen-3.4.0/Eigen/Dense"
#include <iostream>
#include <string>
#include <vector>

template <class ElectronKineticsType> 
void LoKISimulation(std::string setupFileName){
	// ----- CREATING SETUP OBJECT -----

	Setup<ElectronKineticsType> setup(setupFileName);

	// ----- INITIALIZING SIMULATION -----

	setup.initializeSimulation();

	// ----- MAIN BODY OF THE SIMULATION -----

	// loop over the different jobs specified in the setup
	while (setup.currentJobID < setup.numberOfJobs){

		// --- run a particular job 
		if (setup.enableElectronKinetics){
			setup.electronKinetics->solve();
		}

		// --- set up next job
		setup.nextJob();
	}

	// ----- FINISHING SIMULATION -----

	setup.finishSimulation();	
}


int main(int argc, char* argv[]){

	// clear the error log file
	int ret = std::remove("errorLog.txt");

	// ----- GET THE SETUP FILE NAME -----

	std::string setupFileName;
	int numberOfCores = 1; // number of cpu cores (1 by default)
	if (argc == 2){
	    setupFileName.assign(argv[1]);
	}
	else if (argc == 3){
		setupFileName.assign(argv[1]);
		numberOfCores = std::atoi(argv[2]);
	}
	else{
		//Message::error("The name of the setup file should be put when the executable file is called.");
		std::cout<<"Insert the name of the setup file and the number of threads in the following form:\nSETUP_FILE  NUM_THREADS\n";
		char setupName[100];
		int ret = std::scanf("%s %d", setupName, &numberOfCores);
		setupFileName.assign(setupName);
	}

	omp_set_num_threads(numberOfCores);
	Eigen::setNbThreads(numberOfCores);

	// ----- PARSING THE SETUP INFO -----

	Parse::setupFile(setupFileName);

	std::string eedfType = FieldInfo::getFieldValue("electronKinetics.eedfType");

	// (this has to be repeated for each electron kinetics type)

	if (eedfType == "boltzmannMC"){
		LoKISimulation<BoltzmannMC>(setupFileName);
	}

	else if (eedfType == "prescribedEedf"){
		LoKISimulation<PrescribedEedf>(setupFileName);
	}

	else{
		Message::error("Please choose a valid 'electronKinetics->eedfType' in the setup file: 'boltzmannMC' or 'prescribedEedf'.");
	}


	return 0;
}


