/**
 * @file EA_Utilities.h
 * @author  Al Timofeyev
 * @date    April 27, 2019
 * @brief   This utilities file is used by the Evolutionary Algorithms,
 *          and to create matrices using the Mersenne Twister.
 */

#ifndef EVOLUTIONARYALGORITHMS_EA_UTILITIES_H
#define EVOLUTIONARYALGORITHMS_EA_UTILITIES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "BenchmarkFunctions.h"

using namespace std;


/** Prints all the possible Function IDs to the screen. */
void printAllFunctionIDs();
/** Prints all the possible GA Selection IDs to the screen. */
void printAllGASelectionIDs();

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Creates a matrix with the given min/max bound for the given number of rows/columns.*/
vector<vector<double>> createGAMatrix(int rows, int columns, double minBound, double maxBound);
vector<vector<double>> createDEAMatrix(int rows, int columns, double minBound, double maxBound);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Calculates the fitness of a single vector.*/
double calculateFitnessOfVector(vector<double> &vect, int functionID);
/** Calculates the fitness of all vectors in matrix.*/
vector<double> calculateFitnessOfMatrix(vector<vector<double>> matrix, int functionID);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Calculates the average value of a vector of doubles.*/
double calculateAverage(vector<double> vect);

/** Calculates the standard deviation value of a vector of doubles.*/
double calculateStandardDeviation(vector<double> vect);

// ------------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------
/** Special Quicksort implementation for fitness/matrices.*/
void quicksort(vector<double> &fitnessList, vector<vector<double>> &matrix, int L, int R);
/** Swap function for the Quicksort.*/
void swap(vector<double> &fitnessList, vector<vector<double>> &matrix, int x, int y);

// ------------------------------------------------------------------------------------------------------
/** Normal Quicksort implementation for vector arrays.*/
void quicksort(vector<double> &vec, int L, int R);
void swap(vector<double> &v, int x, int y);

#endif //EVOLUTIONARYALGORITHMS_EA_UTILITIES_H
