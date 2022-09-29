#!/bin/bash
if [ "$1" = "all" ]
then
	# update
	apt-get update
	# install gcc and g++ compilers and make
	apt install build-essential
	# assuring that g++9 compiler is installed
	# Details in: https://linuxize.com/post/how-to-install-gcc-compiler-on-ubuntu-18-04/
	apt install software-properties-common
	add-apt-repository ppa:ubuntu-toolchain-r/test
	apt -y install gcc-9 g++-9
	update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9
	update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9
	# install gsl library
	apt-get install libgsl-dev
	# install boost libraries
	apt-get install libboost-all-dev
	# install gnuplot
	apt-get install gnuplot
	# install cmake
	apt-get install cmake
elif [ "$1" != "onlyCode" ]
then
	echo "The option '$1' is not valid. Choose 'all' or 'onlyCode'."
	exit 0
fi
# create/clean the 'Code/build' directory 
cd ../Code
rm -r build
mkdir build
# go to the 'Code/build' directory
cd build
# run cmake
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..
# run make
make
# move the executable to the 'Code' folder
mv lokimc ..
# go to the 'Code' folder
cd ..