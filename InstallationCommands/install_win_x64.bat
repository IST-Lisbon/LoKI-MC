@echo off
if "%1" == "all" (
	set codeInstDir=%cd%
	echo ************************************************************************
	echo ************  Create the directory C:/DEV ************
	mkdir C:\DEV
	cd C:\DEV
	echo ************************************************************************
	echo ************  Get the latest vcpkg version and install it at C:/DEV ************
	set GIT_TRACE_PACKET=1
	set GIT_TRACE=1
	set GIT_CURL_VERBOSE=1
	git clone https://github.com/microsoft/vcpkg.git
	cd vcpkg
	call .\bootstrap-vcpkg.bat
	call .\vcpkg integrate install
	echo ************************************************************************
	echo ************  Install gsl library using vcpkg ************
	call .\vcpkg install gsl gsl:x64-windows
	cd ..
	echo ************************************************************************
	echo ************  Download the boost library **************
	powershell -Command "(New-Object Net.WebClient).DownloadFile('https://sourceforge.net/projects/boost/files/boost/1.77.0/boost_1_77_0.zip/download', 'C:\DEV\boost_1_77_0.zip')"
	echo ************************************************************************
	echo ************  Extract the boost zip file at C:/DEV************
	powershell -Command Expand-Archive -Force boost_1_77_0.zip .
	echo ************************************************************************
	echo ************ Return to the directory Code\InstalllationCommands ************
	cd %codeInstDir%
) else (
	if "%1" neq "onlyCode" (
		echo The option '%1' is not valid. Valid options are 'all' or 'onlyCode'.
		exit /b
	)
)
echo ************************************************************************
echo ************ Compile lokimc.exe using cmake ************
mkdir ..\Code\build
cd ..\Code\build
del /F /Q /S *
cmake -DBOOST_ROOT=C:/DEV/boost_1_77_0 -DGSL_ROOT=C:/DEV/vcpkg/packages/gsl_x64-windows ..
cmake --build . --config Release
echo ************************************************************************
echo ************ Move the executable and dll files to the Code folder ************
move /y  Release\lokimc.exe ..
move /y  Release\gsl.dll ..
move /y  Release\gslcblas.dll ..
echo ************************************************************************
echo ************ JOB DONE! ************
echo ************ You can try to run the code: .\lokimc.exe default_setup.in [number_of_threads]
cd ..