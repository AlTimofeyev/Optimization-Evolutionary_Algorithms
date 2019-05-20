/**
 * @file DifferentialEvolution.h
 * @author  Al Timofeyev
 * @date    May 3, 2019
 * @brief   Implementation of the Differential Evolution Algorithm.
 */

#ifndef EVOLUTIONARYALGORITHMS_DIFFERENTIALEVOLUTION_H
#define EVOLUTIONARYALGORITHMS_DIFFERENTIALEVOLUTION_H

#include <fstream>
#include <chrono>
#include "DEA_Strategies.h"
#include "EA_Utilities.h"

/**
 * @brief Holds all the user defined variables.
 * Differential Evolution Algorithm Configuration Structure, where all user defined
 * variables that are used to configure the Differential Evolution Algorithm are stored.
 */
struct DEA_Config
{
    int dimensions;     /**< Number of dimensions per individual in population. */
    int NP;             /**< Population size. */
    int generations;    /**< Maximum number of generations. */
    double cr;          /**< Crossover Probability. */
    double f;           /**< A scaling factor. */
    double lambda;      /**< A scaling factor. */
    int strategy;       /**< The Differential Evolution strategy to use for mutation/crossover. */
};

/**
 * @brief Holds all the population information.
 * Differential Evolution Algorithm Population Structure, holds all the
 * data related to the population of the Differential Evolution Algorithm.
 */
struct DEA_Population
{
    int functionID;                 /**< The ID determines which benchmark function to call. */
    vector<double> bounds;          /**< Holds the (min,max) bounds of the values for each individual in the population. */
    int functionCounter = 0;        /**< The function counter keeps track of how many times the benchmark function was called. */
    vector<double> fitness;         /**< The fitness for each vector in the population matrix. */
    vector<vector<double>> pop;     /**< The population matrix. */
    vector<double> bestGenFitness;  /**< Keeps track of best fitness from each generation. */
    double executionTime = -1.0;    /**< Time(ms) it took to run the Genetic Algorithm on this population. */
};

/**
 * @brief Differential Evolution Algorithm Analysis
 * Differential Evolution Algorithm Analysis Structure, to keep track
 * of the analysis performed on each population in the population list.
 */
struct DEAAnalysis
{
    string header = "Function ID,Average Fitness,Standard Deviation,Range(min),Range(max),Median,Time(ms),Function Calls,Strategy\n"; /**< Header used when saving the data.*/
    vector<int> functionIDs;                /**< List of function IDs.*/
    vector<double> avgFunctionFitness;      /**< List of the average fitness per FunctionData structure.*/
    vector<double> standardDeviation;       /**< List of standard fitness deviations.*/
    vector<vector<double>> ranges;          /**< List of ranges for each fitness result in resultsOfFunctions.*/
    vector<double> medianFunctionFitness;   /**< List of the Median fitness from each FunctionData structure.*/
    vector<double> executionTimes;          /**< List of execution times in ms for all functions.*/
    vector<int> functionCalls;              /**< List of the amount of times a function was called. */
};


class DifferentialEvolution
{
public:
    // ---------------------- CONSTRUCTORS ----------------------
    DifferentialEvolution(int dim, int sol, int gen, double cr, double f, double lambda, int strategy);

    // ------------------------- METHODS ------------------------
    double runDifferentialEvolution(int functionID, double minBound, double maxBound);  /**< Runs the Differential Evolution Algorithm with set parameters. */
    void analyzeDEAResults();                                                           /**< Analyzes the results of the Differential Evolution Algorithm. */

    void printDEAResults();                                                             /**< Prints the Results of the Differential Evolution Algorithm. */
    void printDEAAnalysis();                                                            /**< Prints the Analysis of the Differential Evolution Algorithm. */

    void saveDEAResults();                                                              /**< Saves all Differential Evolution Algorithm Results to file. */
    void saveDEAAnalysis();                                                             /**< Saves the Analysis of the Differential Evolution Algorithm to file. */

private:
    // ------------------------ VARIABLES -----------------------
    DEA_Config deConfig;
    vector<DEA_Population> popList;
    DEAAnalysis deAnalysis;

    // ------------------------- METHODS ------------------------
    void generateRandDEAPopulation(DEA_Population &population);                                                             /**< Generates the initial population for Differential Evolution Algorithm. */
    void evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter);    /**< Calculates fitness of all solutions in population. */
    void evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter);           /**< Calculate the fitness of an individual solution of the population. */

    vector<double> mutateAndCrossover(const double &cr, const double &f, const double &lambda, vector<double> &x, vector<vector<double>> indiv, mt19937 &randGenerator);    /**< Mutate and Crossover produces one new individual. */
    void select(int functionID, vector<double> &x, double &fitness, vector<double> &newSol, int &functionCounter);                                                          /**< Selects a new individual to be part of the next generation. */

    void saveBestFitness(DEA_Population &pop);  /**< Saves the best fitness from the population. */
};


#endif //EVOLUTIONARYALGORITHMS_DIFFERENTIALEVOLUTION_H
