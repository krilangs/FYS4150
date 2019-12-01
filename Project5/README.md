FYS4150 - Project 5

Report: Contains the report files for the project.

Data: Contains all the data files created by the C++-programs

Figures: Contains figures produced by the Python code `plot_SolarSystem` by reading in the data from the text-files in the Data folder.

Codes: Conatins all the codes for this project. All the code was done on a Windows 10 computer with 8GB RAM. The C++ codes are done Qt Creator, and are compiled and run with the following terminal command: `c++ -floop-parallelize-all -o main.exe main.cpp classes.cpp utils.cpp SolarSystem.cpp`.
- plot_SolarSystem.py: Creates all the figures in the Figures folder. This was done in Python 3.7 with Spyder 3.3.
- main.cpp: This is the main program for running the functions and simulating the Solar System tasks.
- classes.cpp: Contains virtual class Obj for n-body systems. This program simulates the positions of the bodies in the system using the forward Euler and velocity Verlet algorithms, and writes the results to files in the Data folder.
- SolarSystem.cpp: Contains the class SolarSystem, which is a subclass of Obj, and the class SolarRelativistic, which is a subclass of SolarSystem.
- utils.cpp: Contains supporting functions.
- classes.h: Header file for the class Obj and subclasses. Also contains the struct MassObject
- planetsLib.h: Contains the StellarObjectsLibrary, which has the masses, initial positions and initial velocities of the Sun and the planets in the Solar System.
