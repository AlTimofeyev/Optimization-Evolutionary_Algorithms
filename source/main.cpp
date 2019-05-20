#include <iostream>
#include <fstream>
#include "GeneticAlgorithm.h"
#include "DifferentialEvolution.h"
#include "utilities.h"

using namespace std;

void test_GA(string configFilename, string iniFilename);
void test_DEA(string configFilename, string iniFilename);

int main(int argc, char **argv)
{
    // Default configuration and initialization filenames.
    string gaConfigFilename = "GA-Config.txt";
    string eaConfigFilename = "DEA-Config.txt";
    string iniFilename = "iniFile.txt";

    // If an initialization filename was provided, reassign the variable to user input.
    if (argc == 2)
        iniFilename = argv[1];

    // Test the Genetic Algorithm.
    //test_GA(gaConfigFilename, iniFilename);

    // Test the Differential Evolution
    test_DEA(eaConfigFilename, iniFilename);

    return 0;
}

// ******************************************************************************************************
// ******************************************************************************************************

void test_GA(string configFilename, string iniFilename)
{
    // Open the config file.
    ifstream configFile;
    configFile.open(configFilename);
    if(configFile.fail())
    {
        cout << "Failed to open file: " << configFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory\n";
        cout << "or does not exist.\n";
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Get the configuration values from file.
    string line;
    while(configFile.good())
    {
        getline(configFile, line);
    }
    vector<double> tokens = parseStringDbl(line, ",");

    // Close the config file.
    configFile.close();

    // If there aren't enough values.
    if(tokens.size() < 10)
    {
        cout << "Not Enough Values in " << configFilename << endl;
        cout << "There Should Be 10 Total Values.\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Create Genetic Algorithm object and configure it.
    GeneticAlgorithm geneticAlgorithm(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5], tokens[6], tokens[7], tokens[8], tokens[9]);

    // Open the ini file.
    ifstream iniFile;
    iniFile.open(iniFilename);
    if(iniFile.fail())
    {
        cout << "Failed to open file: " << iniFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory,\n";
        cout << "does not exist, or was not provided as\n";
        cout << "a command line argument." << endl;
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Run Genetic Algorithm on all values in initialization file.
    while(iniFile.good())
    {
        getline(iniFile, line);
        tokens = parseStringDbl(line, ",");
        vector<double> bounds = {tokens[1], tokens[2]};

        // Initialize the Genetic Algorithm.
        geneticAlgorithm.setPopulationParams(tokens[0], bounds);

        // Run the Genetic Algorithm.
        geneticAlgorithm.runGeneticAlgorithm();

        // Print the results.
        geneticAlgorithm.printGAResults();
    }

    geneticAlgorithm.analyzeGAResults();
    geneticAlgorithm.printGAAnalysis();

    // Save the results.
    geneticAlgorithm.saveGAResults(iniFilename);
    geneticAlgorithm.saveGAAnalysis();

    // Close the ini file.
    iniFile.close();
}

// ******************************************************************************************************
// ******************************************************************************************************

void test_DEA(string configFilename, string iniFilename)
{
    // Open the config file.
    ifstream configFile;
    configFile.open(configFilename);
    if(configFile.fail())
    {
        cout << "Failed to open file: " << configFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory\n";
        cout << "or does not exist.\n";;
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Get the configuration values from file.
    string line;
    while(configFile.good())
    {
        getline(configFile, line);
    }
    vector<double> tokens = parseStringDbl(line, ",");

    // Close the config file.
    configFile.close();

    // If there aren't enough values.
    if(tokens.size() < 7)
    {
        cout << "Not Enough Values in " << configFilename << endl;
        cout << "There Should Be 7 Total Values.\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Create Differential Evolution Algorithm object and configure it.
    DifferentialEvolution deAlgorithm(tokens[0], tokens[1], tokens[2], tokens[3], tokens[4], tokens[5], tokens[6]);

    //cout << tokens[0] << "\t" << tokens[1] << "\t" << tokens[2] << "\t" << tokens[3]<< "\t" << tokens[4]<< "\t" << tokens[5]<< "\t" << tokens[6] << endl;

    // Open the ini file.
    ifstream iniFile;
    iniFile.open(iniFilename);
    if(iniFile.fail())
    {
        cout << "Failed to open file: " << iniFilename << endl;
        cout << "-----------------------------------------------\n";
        cout << "File is either not in the right directory,\n";
        cout << "does not exist, or was not provided as\n";
        cout << "a command line argument." << endl;
        cout << "-----------------------------------------------\n";
        cout << "Accepted File Formats: .txt" << endl;
        cout << "-----------------------------------------------\n";
        cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
        exit(1);
    }

    // Run Differential Evolution Algorithm on all values in initialization file.
    while(iniFile.good())
    {
        getline(iniFile, line);
        tokens = parseStringDbl(line, ",");

        // Run the Differential Evolution Algorithm.
        deAlgorithm.runDifferentialEvolution(tokens[0], tokens[1], tokens[2]);
    }

    deAlgorithm.analyzeDEAResults();
    deAlgorithm.printDEAAnalysis();

    // Save the results.
    deAlgorithm.saveDEAResults();
    deAlgorithm.saveDEAAnalysis();

    // Close the ini file.
    iniFile.close();
}