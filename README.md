<p align="center">----------------------------- README FILE FOR THE LoKI-MC SIMULATION TOOL -----------------------------<br>
<align="center">-----------------------------  10 steps to get acquainted with the tool  -----------------------------<br>
<align="center">(Version LoKI-MC_v1.1.0)</p>

1. What's LoKI-MC distribution license ?<br>
   The LisbOn KInetics Monte Carlo (LoKI-MC) is an open-source tool, licensed under the GNU general public license.<br>
   LoKI-MC is freely available for users to perform electron kinetics calculations, and for expert researchers who are invited to continue testing the tool and/or to contribute for its development and improvement.<br>
   Note that the tool makes use of other external codes, present at "Code/External/eigen-3.4.0" and "Code/External/MathParser", which are distributed under the licenses Mozilla Public License version 2.0 and Apache License version 2.0, respectively.

2. How to contact the developers ?<br> 
   You are much welcome to report any problem or bug you may find using the code, or to send any feedback with suggestions for further developments.
   After downloading LoKI-MC, and especially if you intend to interact with us, you are invited to send a short message   
   to: loki@tecnico.ulisboa.pt   
   with subject: <b>LoKI-MC</b>   
   just giving your <b>name</b> and <b>affiliation</b>.   

3. What's LoKI-MC ?<br>
   LoKI-MC is a Monte Carlo simulation tool that solves the electron kinetics for low-temperature plasmas excited by AC/DC electric and DC magnetic fields in any complex gas mixture.<br>
   One very interesting feature is the inclusion of general anisotropic scattering for any collision type.

4. What's the programming language of LoKI-MC ?<br>
   LoKI-MC is developed in C++, with a highly-performant object-oriented structure and adopting an ontology that privileges the separation between tool and data.

5. What are the input data of LoKI-MC ?<br>
   On input, the code requires the working conditions (e.g. the reduced values of the electric field and magnetic field, the AC frequency, the gas pressure and temperature), the gas-mixture composition, the distributions of populations for the levels of the atomic/molecular gases considered, and the relevant sets of electron-scattering cross sections obtained from the open-access website LXCat (http://www.lxcat.net/).

6. What are the output results of LoKI-MC ?<br>
   On output, it yields the electron energy distribution function, the electron velocity distribution function, the electron swarm parameters, the collision rate-coefficients, and the electron power absorbed from the electric field and transferred to the different collisional channels.

7. How to find your way in the code ?<br>
   After pulling the files in the repository, the LoKI-MC folder contains<br>
   A) Subfolder "Documentation", with important documentation files - PLEASE READ THEM BEFORE USING THE CODE !!!<br>
   B) Subfolder "CompiledExecutables", with compiled executables of the code for different operating systems. Please do not use them without reading the user manual in the documentation. <br>
   C) Subfolder "InstallationCommands", with installation scripts for different operating systems. The utilization of these scripts is described in the 'Quick-start guide' and the appendix 'Installation guide' of the user manual.<br>
   D) Subfolder "Code" containing <br>
   &ensp;(a) A subfolder "Input", containing the input files required for the simulations, organized as follows <br>
   &ensp;&ensp;i. Default configuration file 'default_lokimc_setup.in' <br>
   &ensp;&ensp;ii. A subfolder "Databases" with '\*.txt' files, containing different properties (masses, energies of levels, atomic/molecular constants, ...) for the gases used in the simulations.<br>
   &ensp;&ensp;iii. Several subfolders "Helium", "Nitrogen", ... with '\*.txt' files, containing the electron-scattering cross sections for the different gases used in the simulations, usually obtained from the open-access website LXCat (http://www.lxcat.net/). <br>
   &ensp;(b) A subfolder "LoKI-MC", containing the code files of LoKI-MC, organized as follows <br>
   &ensp;&ensp;i. A subfolder "Headers" with all ".h" files. Along them, you can find, for example: in "Constant.h" the physical constants used by code; in "GasPropertyFunctions.h" and "StatePropertyFunctions.h" the functions used to assign the properties of the gases and states, respectively. <br>
   &ensp;&ensp;ii. A subfolder "Sources" with all ".C" files. The main function is located at "lokimc.C". <br>
   &ensp;&ensp;iii. A subfolder "GnuplotFunctions" with Gnuplot functions used for the graphical interface. <br>
   &ensp;(c ) A subfolder "External", containing the external codes used by LoKI-MC (see also step 1 for their distribution licenses). <br>        
   &ensp;(d) A subfolder "Output", where LoKI-MC will write the output files resulting from the simulations. The subfolder "swarm_O2_short" contains an example of the output provided by running the default configuration file. <br>

8. How to run LoKI-MC ? <br>
   Before running the code, the user should read the manual present in the documentation folder. <br>
   The LoKI-MC installation is detailed in the "Quick-start guide" and the appendix "Installation guide". <br>

   LoKI-MC runs upon calling the command './lokimc setupFile numberOfThreads'. <br>
   The end-user interacts with the code by specifying a particular "setupFile" for the simulation and by defining the number of threads that will be used for the calculations (ideally, equal to the number of CPU cores). <br>
   The setup files should be located in [repository folder]/LoKI-MC/Code/Input/ with a '.in' extension <br>
   (this is just a recommendation in order to keep the input folder organized; the setup files are just plain text files). <br>

   After the installation, the distribution of LoKI-MC includes a couple of default configuration files to help you make a first run of the code, following the sequence of steps below: <br>
   A) Open the terminal <br>
   B) Navigate to the "Code" folder of your local copy of the repository: <br>
   &ensp;&ensp;>> cd [repository folder]/LoKI-MC/Code/ <br>
   C) Execute the following command: <br>
   &ensp;&ensp;>> './lokimc default_lokimc_setup.in numberOfThreads' <br>
   D) When the calculations are finished, the graphical user interface (GUI) should show the solution(s) for the default setup file. <br>

9. How to reference the code and the cross sections? <br>
   LoKI-MC is the result of the efforts of the Portuguese group N-PRiME, that decided to share the outcome of its research with the members of the Low-Temperature Plasmas community. 

   When using LoKI-MC in your work, please give proper credits to the main developers, by adding the following citations: <br>
   [] Dias T C, Pintassilgo C D and Guerra V 2023 Plasma Sources Sci. Technol. 32 095003 <br>
   [] Dias T C, Tejero-del-Caz A, Alves L L and Guerra V 2023 Comput. Phys. Commun. 282 108554 <br>

   Additionally, do not forget to reference properly the LXCat databases used for the cross sections, which can be found in the beginning of each LXCat file.

10. Acknowledgments <br>
   This work was funded by Portuguese FCT - Fundação para a Ciência e a Tecnologia, under Projects UIDB/50010/2020 and UIDP/50010/2020 and Grant PD/BD/150414/2019 (PD-F APPLAuSE). 
