#!/bin/bash

if [ "$1" = "brewAll" ]
then
	# ************* installation commands using Homebrew package manager *****************
	# install Homebrew 
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
	# install gcc compiler
	brew install gcc@13
	# install the gsl libraries
	brew install gsl
	# install the boost libraries
	brew install boost
	# install gnuplot
	brew install gnuplot
	# install cmake
	brew install cmake
elif [ "$1" = "portAll" ]
then
	# ************* installation commands using Macports package manager *****************
	# install gcc compiler
	sudo port install gcc13
	# install the gsl libraries
	sudo port install gsl
	# install the boost libraries
	sudo port install boost
	# install gnuplot
	sudo port install gnuplot +qt
	# install cmake
	sudo port install cmake
elif [ "$1" != "brewOnlyCode" ] && [ "$1" != "portOnlyCode" ]
then
	echo "The option '$1' is not valid. Choose 'brewAll', 'portAll', 'brewOnlyCode', 'portOnlyCode'."
	exit 0
fi

# create/clean the 'Code/build' directory 
cd ../Code
sudo rm -r build
mkdir build
# go to the 'Code/build' directory
cd build
if [ "$1" = "brewOnlyCode" ] || [ "$1" = "brewAll" ]
then
	# run cmake. If other compiler version than gc13 is installed, this must be changed. AppleClang is not supported 
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 ..
elif [ "$1" = "portOnlyCode" ] || [ "$1" = "portAll" ]
then
	cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-mp-13 -DCMAKE_CXX_COMPILER=g++-mp-13 ..
fi
# run make
make
# move the executable to the 'Code' folder
sudo mv lokimc ..
# go to the 'Code' folder
cd ..