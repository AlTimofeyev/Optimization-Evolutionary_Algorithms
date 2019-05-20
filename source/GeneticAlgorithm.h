/**
 * @file GeneticAlgorithm.h
 * @author  Al Timofeyev
 * @date    April 25, 2019
 * @brief   Implementation of the Genetic Algorithm.
 */

#ifndef EVOLUTIONARYALGORITHMS_GENETICALGORITHM_H
#define EVOLUTIONARYALGORITHMS_GENETICALGORITHM_H

#include <fstream>
#include <chrono>
#include "EA_Utilities.h"

using namespace std;


/**
 * @brief Holds all the user defined variables.
 * Genetic Algorithm Configuration Structure, where all user defined
 * variables that are used to configure the Genetic Algorithm are stored.
 */
struct GA_Config
{
    int dimensions;     /**< Number of dimensions per individual in population. */
    int solutions;      /**< Population size. */
    int generations;    /**< Maximum number of generations. */
    double cr;          /**< Crossover Probability. */
    double mutProb;     /**< Mutation Probability. */
    double mutRange;    /**< Mutation Range. */
    double mutPrec;     /**< Mutation Precision. */
    double er;          /**< Elitism Rate. */
    int eliteIndex;     /**< Elitism ending Index in the population (0 - eliteIndex). */
    int selectionID;    /**< The selection type to use for the Genetic Algorithm. */
    int crPoints;       /**< The number of crossover points. */
};

/**
 * @brief Holds all the population information.
 * Genetic Algorithm Population Structure, holds all the data related
 * to the population of the Genetic Algorithm.
 */
struct GA_Population
{
    int functionID;                 /**< The ID determines which benchmark function to call. */
    vector<double> bounds;          /**< Holds the (min,max) bounds of the values for each individual in the population. */
    int functionCounter = 0;        /**< The function counter keeps track of how many times the benchmark function was called. */
    vector<double> fitness;         /**< The fitness for each vector in the population matrix. */
    vector<vector<double>> pop;     /**< The population matrix. */
    vector<double> bestGenFitness;  /**< Keeps track of best fitness from each generation. */
    double totalFitness;            /**< The summed up total of the fitness values. */
    double executionTime = -1.0;    /**< Time(ms) it took to run the Genetic Algorithm on this population. */
};

/**
 * @brief Genetic Algorithm Analysis
 * Genetic Algorithm Analysis Structure, to keep track of the analysis
 * performed on each population in the population list.
 */
struct GAAnalysis
{
    string header = "Function ID,Average Fitness,Standard Deviation,Range(min),Range(max),Median,Time(ms),Function Calls\n"; /**< Header used when saving the data.*/
    vector<int> functionIDs;                /**< List of function IDs.*/
    vector<double> avgFunctionFitness;      /**< List of the average fitness per FunctionData structure.*/
    vector<double> standardDeviation;       /**< List of standard fitness deviations.*/
    vector<vector<double>> ranges;          /**< List of ranges for each fitness result in resultsOfFunctions.*/
    vector<double> medianFunctionFitness;   /**< List of the Median fitness from each FunctionData structure.*/
    vector<double> executionTimes;          /**< List of execution times in ms for all functions.*/
    vector<int> functionCalls;              /**< List of the amount of times a function was called. */
};


class GeneticAlgorithm
{
public:
    // ---------------------- CONSTRUCTORS ----------------------
    GeneticAlgorithm(int dim, int sol, int gen, double cr, double mProb, double mR, double mPrec, double er, int selectionID, int crPoints);

    // ------------------------- METHODS ------------------------
    void setPopulationParams(int functionID, vector<double> bounds);    /**< Sets up the parameters for the population. */
    double runGeneticAlgorithm();                                       /**< Runs the Genetic Algorithm with set parameters. */
    void analyzeGAResults();                                            /**< Analyzes the results of the Genetic Algorithm. */

    void printGAResults();                                              /**< Prints the Results of the Genetic Algorithm. */
    void printGAAnalysis();                                             /**< Prints the Analysis of the Genetic Algorithm. */

    void saveGAResults(string iniFilename);                             /**< Saves all Genetic Algorithm Results to file. */
    void saveGAAnalysis();                                              /**< Saves the Analysis of the Genetic Algorithm to file. */

private:
    // ------------------------ VARIABLES -----------------------
    GA_Config gaConfig;
    vector<GA_Population> popList;
    GA_Population population;
    GAAnalysis gaAnalysis;

    // ------------------------- METHODS ------------------------
    void generateRandPopulation();                              /**< Generates the initial population for Genetic Algorithm. */
    void evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter); /**< Calculates fitness of all solutions in population. */

    vector<int> select(int selectionID, vector<double> &popFitness, mt19937 &randGenerator);    /**< Selects two parents from the population. */
    int t_select(vector<double> &popFitness, int numOfInid, mt19937 &randGenerator);            /**< Uses Tournament Selection to select one parent from population. */
    int rw_select(vector<double> &popFitness, mt19937 &randGenerator);                          /**< Uses Roulette Wheel Selection to select one parent from population. */

    vector<vector<double>> crossover(int crPoints, vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator);  /**< Returns two children that are a crossover of the parents. */
    vector<vector<double>> crossover1(vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator);               /**< 1 Crossover - Returns two children that are a crossover of the parents. */
    vector<vector<double>> crossover2(vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator);               /**< 2 Crossovers - Returns two children that are a crossover of the parents. */

    void mutate(vector<double> &child, double mProb, double mRange, double mPrec, vector<double> bounds, mt19937 &randGenerator);   /**< Mutates the child. */

    void reduce(GA_Population &pop, vector<vector<double>> &newPop, vector<double> &newFit, int eliteIndex); /**< Combines new population with the old one. */

    void normalizeFitness(int functionID, vector<double> &popFitness);    /**< Normalizes all fitness values for a population. */
    void recordBestFitness(GA_Population &pop);           /**< Saves the best fitness value from the population. */
};


#endif //EVOLUTIONARYALGORITHMS_GENETICALGORITHM_H
