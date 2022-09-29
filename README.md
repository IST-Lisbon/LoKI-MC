<p align="center">----------------------------- README FILE FOR THE LoKI-MC SIMULATION TOOL -----------------------------<br>
<align="center">-----------------------------  10 steps to get acquainted with the tool  -----------------------------<br>
<align="center">(Version LoKI-MC_v1.0)</p>

1. What's LoKI-MC distribution license ?   
   The LisbOn KInetics Monte Carlo (LoKI-MC) is an open-source tool, licensed under the GNU general public license.  
   LoKI-MC is freely available for users to perform electron kinetics calculations, and for expert researchers who are invited to continue testing the tool and/or to contribute for its development and improvement.
   Note that the tool makes use of other external codes, present at "Code/External/eigen-3.4.0" and "Code/External/MathParser", which are distributed under the licenses Mozilla Public License version 2.0 and Apache License version 2.0, respectively.

2. How to contact the developers ?   
   You are much welcome to report any problem or bug you may find using the code, or to send any feedback with suggestions for further developments.
   After downloading LoKI-MC, and especially if you intend to interact with us, you are invited to send a short message   
   to: loki@tecnico.ulisboa.pt   
   with subject: <b>LoKI-MC</b>   
   just giving your <b>name</b> and <b>affiliation</b>.   

3. What's LoKI-MC ?   
   LoKI-MC is a Monte Carlo simulation tool that solves the electron kinetics for non-magnetized low-temperature plasmas excited by homogeneous DC electric fields in any complex gas mixture.

4. What's the programming language of LoKI-MC ?   
   LoKI-MC is developed in C++, with a highly-performant object-oriented structure and adopting an ontology that privileges the separation between tool and data.

5. What are the input data of LoKI-MC ?
   On input, the code requires the working conditions (e.g. the applied reduced electric field, the gas pressure and temperature), the gas-mixture composition, the distributions of populations for the levels of the atomic/molecular gases considered, and the relevant sets of electron-scattering cross sections obtained from the open-access website LXCat (http://www.lxcat.net/).

6. What are the output results of LoKI-MC ?   
   On output, it yields the electron energy distribution function, the electron velocity distribution function, the electron swarm parameters, the collision rate-coefficients, and the electron power absorbed from the electric field and transferred to the different collisional channels.

7. How to find your way in the code ?   
   After pulling the files in the repository, the LoKI-MC folder contains   
   A) Subfolder "Documentation", with important documentation files - PLEASE READ THEM BEFORE USING THE CODE !!!
   B) Subfolder "CompiledExecutables", with compiled executables of the code for different operating systems. Please do not use them without reading the user manual in the documentation.
   C) Subfolder "InstallationCommands", with installation scripts for different operating systems. The utilization of these scripts is described in the 'Quick-start guide' and the appendix 'Installation guide' of the user manual.   
   D) Subfolder "Code" containing    
   &ensp;(a) A subfolder "Input", containing the input files required for the simulations, organized as follows   
   &ensp;&ensp;i. Default configuration file 'default_lokimc_setup.in' 
   &ensp;&ensp;ii. A subfolder "Databases" with '\*.txt' files, containing different properties (masses, energies of levels, atomic/molecular constants, ...) for the gases used in the simulations.   
   &ensp;&ensp;iii. Several subfolders "Helium", "Nitrogen", ... with '\*.txt' files, containing the electron-scattering cross sections for the different gases used in the simulations, usually obtained from the open-access website LXCat (http://www.lxcat.net/).
   &ensp;(b) A subfolder "LoKI-MC", containing the code files of LoKI-MC, organized as follows
   &ensp;&ensp;i. A subfolder "Headers" with all ".h" files. Along them, you can find, for example: in "Constant.h" the physical constants used by code; in "GasPropertyFunctions.h" and "StatePropertyFunctions.h" the functions used to assign the properties of the gases and states, respectively.
   &ensp;&ensp;ii. A subfolder "Sources" with all ".C" files. The main function is located at "lokimc.C".
   &ensp;&ensp;iii. A subfolder "GnuplotFunctions" with Gnuplot functions used for the graphical interface.
   &ensp;(c ) A subfolder "External", containing the external codes used by LoKI-MC (see also step 1 for their distribution licenses).           
   &ensp;(d) A subfolder "Output", where LoKI-MC will write the output files resulting from the simulations. The subfolder "swarm_O2" contains an example of the output provided by running the default configuration file.

8. How to run LoKI-MC ?   
   Before running the code, the user should read the manual present in the documentation folder.
   The LoKI-MC installation is detailed in the "Quick-start guide" and the appendix "Installation guide".

   LoKI-MC runs upon calling the command './lokimc setupFile numberOfThreads'.   
   The end-user interacts with the code by specifying a particular "setupFile" for the simulation and by defining the number of threads that will be used for the calculations (ideally, equal to the number of CPU cores).      
   The setup files should be located in [repository folder]/LoKI-MC/Code/Input/ with a '.in' extension   
   (this is just a recommendation in order to keep the input folder organized; the setup files are just plain text files).   

   After the installation, the distribution of LoKI-MC includes a couple of default configuration files to help you make a first run of the code, following the sequence of steps below:   
   A) Open the terminal   
   B) Navigate to the "Code" folder of your local copy of the repository:   
   &ensp;&ensp;>> cd [repository folder]/LoKI-MC/Code/   
   C) Execute the following command:   
   &ensp;&ensp;>> './lokimc default_setup.in numberOfThreads'  
   D) When the calculations are finished, the graphical user interface (GUI) should show the solution(s) for the default setup file.

9. How to reference the code and the cross sections?   
   LoKI-MC is the result of the efforts of the Portuguese group N-PRiME, that decided to share the outcome of its research with the members of the Low-Temperature Plasmas community.

   When using LoKI-MC in your work, please give proper credits to the main developers, by adding the following citation:     
   [] Dias T C et al The LisbOn KInetics Monte Carlo solver 2022 Comput. Phys. Commun. 108554, doi: https://doi.org/10.1016/j.cpc.2022.108554

   Additionally, do not forget to reference properly the LXCat databases used for the cross sections, which can be found in the beginning of each LXCat file.

10. Acknowledgments   
   This work was funded by Portuguese FCT - Fundação para a Ciência e a Tecnologia, under Projects UIDB/50010/2020 and UIDP/50010/2020 and Grant PD/BD/150414/2019 (PD-F APPLAuSE).
