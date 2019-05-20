******************************************************************************************************************
Author:	Al Timofeyev
Date:	May 3, 2019
Desc:	This is how to compile and execute the code.
******************************************************************************************************************

**********************************************************************************
NOTES TO PROFESSOR (if any) ARE AT THE VERY BOTTOM OF THIS README!!
**********************************************************************************


**********************************************************************************
---------------------------- ENVIRONMENT USED TO CODE ----------------------------
Windows 10
CLion version 2019.1
cygwin version 3.0.4
cygwin GDB version 8.1.1
gcc version 7.4.0 
g++ version 7.4.0

******
NOTE:
1)	CLion generated a CMakeLists.txt file included with the source code.
	cmake_minimum_required(VERSION 3.13)
2)	The program was written in C++.
******
**********************************************************************************


**********************************************************************************
--------------------------- SETUP CONFIGURATION FILES ----------------------------
---- Structure of Genetic Algorithm configuration file
<Parameter Header>		---- First line is headers only (DO NOT CHANGE HEADER LINE).
<list of parameters>		---- Second line is a list of parameters.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
Dims,Pop,Gen,CR,MProb,MRange,MPrec,ER,Selection,CRPs		---- First line is the parameter headers.
30,200,100,0.8,0.005,0.1,5,0.2,2,1				---- Second line is a list of parameters.

-- What Each Parameter Means:
Dims	=	The number of dimensions.
Pop	=	The population size.
Gen	=	The maximum number of generations.
CR	=	The Crossover Probability (between [0, 1]).
MProb	=	The Mutation Probability (between [0, 1]).
MRange	=	The Mutation Range.
MPrec	=	The Mutation Precision (between [1,5]).
ER	=	The Elitism Rate (between [0,1]).
Selection=	The Selection Type (between [1,2]).
CRPs	=	The number of crossover points (between [1,2]).

----------------------------------------------------------------------------------
---- Structure of Differential Evolution Algorithm configuration file
<Parameter Header>		---- First line is headers only (DO NOT CHANGE HEADER LINE).
<list of parameters>		---- Second line is a list of parameters.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
Dims,Pop,Gen,CR,F,Lambda,Strategy		---- First line is the parameter headers.
30,200,100,0.8,0.6,1.6,1			---- Second line is a list of parameters.

-- What Each Parameter Means:
Dims	=	The number of dimensions.
Pop	=	The population size.
Gen	=	The maximum number of generations.
CR	=	The Crossover Probability between [0, 1].
F	=	A scaling factor between (0, 1.2].
Lambda	=	A scaling factor between (0,2].
Strategy=	One of the 10 implemented strategies.


-------------------------- SETUP INITIALIZATION FILES ----------------------------
---- Structure of configuration file
<list of dimensions>				---- First line only
<list of function IDs and their bounds>		---- All other lines.
Use only a comma (,) delimiter, no spaces between values.

-- Example:
1,-500,40		---- First line is for Benchmark Function 1, with -500/40 min/max bounds.
5,-32,100		---- Second line is for Benchmark Function 5, with -32/100 min/max bounds.
8,0,pi,0		---- Third line is Benchmark Function 8, with 0/pi min/max bounds.
			     \--> PLEASE LOOK AT NOTE 3 IN THIS SECTION.
******
NOTE:
1)	Depending on which IDE you are running, conig files should be either
	in the same folder as source code or in build folder.
2)	Initialization files can be passed as command line parameters or use the default
	configuration file (just alter it).
3)	Please note the extra zero (0) value after pi (in above example). On any line
	that contains the value pi, please include an extra value, like zero. This is
	for conversion purposes, otherwise the program will not run.
******
**********************************************************************************


**********************************************************************************
------------------------------ COMPILE AND EXECUTE -------------------------------
---- To compile for an IDE project.
To Compile:
You could use CMake to compile CMakeLists.txt file that's included with source code.

To Execute:
run main.cpp


---- I'm assuming it could also be compiled and run from command line:
To Compile:
g++ -o main main.cpp

To Execute:
./main				---- Default ini.txt is used as initialization file.
./main inifrog.txt		---- inifrog.txt is initialization file example.
./main inifigFile2.txt		---- inifigFile2.txt is initialization file example.
./main blabla.txt		---- blabla.txt is initialization file example.
**********************************************************************************


**********************************************************************************
------------------------------ NOTES TO PROFESSOR --------------------------------
1)	You can have only ONE configuration file but mutlitple initialization files.

2)	The initialization file used to start the program can be passed as a command
	line parameter or you can change the default iniFile.txt file.

3)	Some of the functions in the DifferentialEvolution.h file have not been defined.

4)	Even though the 2-crossover point is defined for Genetic Algorithm, it does
	not work. The code just stops working in the reduce function when it tries
	to add to the new population the old individuals from the previous generation.
**********************************************************************************