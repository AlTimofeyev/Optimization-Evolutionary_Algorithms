/**
 * @file DEA_Strategies.h
 * @author  Al Timofeyev
 * @date    May 2, 2019
 * @brief   A library of strategies used in the Differential Evolution Algorithm.
 *
 * The general notation used for these strategies is DE/x/y/z: where DE stands for Differential
 * Evolution algorithm, x represents a string denoting the vector to be perturbed, y is the
 * number of difference vectors considered for perturbation of x, and z is the type of crossover
 * being used. Two types of crossovers: exp (exponential) and  bin (binomial).
 */

#ifndef EVOLUTIONARYALGORITHMS_DEA_STRATEGIES_H
#define EVOLUTIONARYALGORITHMS_DEA_STRATEGIES_H


#include <vector>
#include <random>

using namespace std;


/** Strategy 1: DE/best/1/exp */
vector<double> de_Strategy1(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator);

/** Strategy 2: DE/rand/1/exp */
vector<double> de_Strategy2(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator);

/** Strategy 3: DE/rand-to-best/1/exp */
vector<double> de_Strategy3(const double &cr, const double &F, const double &lambda, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, mt19937 &randGenerator);

/** Strategy 4: DE/best/2/exp */
vector<double> de_Strategy4(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, mt19937 &randGenerator);

/** Strategy 5: DE/rand/2/exp */
vector<double> de_Strategy5(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, vector<double> xRand5, mt19937 &randGenerator);

/** Strategy 6: DE/best/1/bin */
vector<double> de_Strategy6(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator);

/** Strategy 7: DE/rand/1/bin */
vector<double> de_Strategy7(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator);

/** Strategy 8: DE/rand-to-best/1/bin */
vector<double> de_Strategy8(const double &cr, const double &F, const double &lambda, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, mt19937 &randGenerator);

/** Strategy 9: DE/best/2/bin */
vector<double> de_Strategy9(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, mt19937 &randGenerator);

/** Strategy 10: DE/rand/2/bin */
vector<double> de_Strategy10(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, vector<double> xRand5, mt19937 &randGenerator);


#endif //EVOLUTIONARYALGORITHMS_DEA_STRATEGIES_H
