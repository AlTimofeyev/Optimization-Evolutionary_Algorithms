/**
 * @file DEA_Strategies.cpp
 * @author  Al Timofeyev
 * @date    May 2, 2019
 * @brief   A library of strategies used in the Differential Evolution Algorithm.
 *
 * The general notation used for these strategies is DE/x/y/z: where DE stands for Differential
 * Evolution algorithm, x represents a string denoting the vector to be perturbed, y is the
 * number of difference vectors considered for perturbation of x, and z is the type of crossover
 * being used. Two types of crossovers: exp (exponential) and  bin (binomial).
 */

#include "DEA_Strategies.h"

/**
 * @brief Strategy 1: DE/best/1/exp.
 *
 * @note xBest != xRand2 != xRand3 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy1(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xBest[i] + F * (xRand2[i] - xRand3[i]));

    // Generate a random crossover probability.
    double randCR = crDis(randGenerator);

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    int i = 0;
    while(randCR < cr && i < newSolutionU.size())
    {
        // Do crossover.
        newSolutionU[i] = trialV[i];

        // Update the random cr and i index.
        randCR = crDis(randGenerator);
        i++;
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 2: DE/rand/1/exp.
 *
 * @note xRand1 != xRand2 != xRand3 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy2(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xRand1[i] + F * (xRand2[i] - xRand3[i]));

    // Generate a random crossover probability.
    double randCR = crDis(randGenerator);

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    int i = 0;
    while(randCR < cr && i < newSolutionU.size())
    {
        // Do crossover.
        newSolutionU[i] = trialV[i];

        // Update the random cr and i index.
        randCR = crDis(randGenerator);
        i++;
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 3: DE/rand-to-best/1/exp.
 *
 * @note xBest != xRand1 != xRand2 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param lambda A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy3(const double &cr, const double &F, const double &lambda, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(x[i] + lambda * (xBest[i] - x[i]) + F * (xRand1[i] - xRand2[i]));

    // Generate a random crossover probability.
    double randCR = crDis(randGenerator);

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    int i = 0;
    while(randCR < cr && i < newSolutionU.size())
    {
        // Do crossover.
        newSolutionU[i] = trialV[i];

        // Update the random cr and i index.
        randCR = crDis(randGenerator);
        i++;
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 4: DE/best/2/exp.
 *
 * @note xBest != xRand1 != xRand2 != xRand3 != xRand4 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param xRand4 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy4(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xBest[i] + F * (xRand1[i] + xRand2[i] - xRand3[i] - xRand4[i]));

    // Generate a random crossover probability.
    double randCR = crDis(randGenerator);

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    int i = 0;
    while(randCR < cr && i < newSolutionU.size())
    {
        // Do crossover.
        newSolutionU[i] = trialV[i];

        // Update the random cr and i index.
        randCR = crDis(randGenerator);
        i++;
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 5: DE/rand/2/exp.
 *
 * @note xRand1 != xRand2 != xRand3 != xRand4 != xRand5 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param xRand4 A randomly chosen solution from the population.
 * @param xRand5 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy5(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, vector<double> xRand5, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xRand5[i] + F * (xRand1[i] + xRand2[i] - xRand3[i] - xRand4[i]));

    // Generate a random crossover probability.
    double randCR = crDis(randGenerator);

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    int i = 0;
    while(randCR < cr && i < newSolutionU.size())
    {
        // Do crossover.
        newSolutionU[i] = trialV[i];

        // Update the random cr and i index.
        randCR = crDis(randGenerator);
        i++;
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 6: DE/best/1/bin.
 *
 * @note xBest != xRand2 != xRand3 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy6(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> dimDis(0, x.size()-1);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xBest[i] + F * (xRand2[i] - xRand3[i]));

    // Generate a random index withing range of x-dimensions.
    int randIndex = dimDis(randGenerator);

    // Declare a random crossover probability variable.
    double randCR;

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    for(int i = 0; i < trialV.size(); i++)
    {
        // Generate a random crossover probability.
        randCR = crDis(randGenerator);

        // Crossover when random cr is within the cr range.
        if(randCR < cr || i == randIndex)
            newSolutionU[i] = trialV[i];
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 7: DE/rand/1/bin.
 *
 * @note xRand1 != xRand2 != xRand3 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy7(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> dimDis(0, x.size()-1);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xRand1[i] + F * (xRand2[i] - xRand3[i]));

    // Generate a random index withing range of x-dimensions.
    int randIndex = dimDis(randGenerator);

    // Declare a random crossover probability variable.
    double randCR;

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    for(int i = 0; i < trialV.size(); i++)
    {
        // Generate a random crossover probability.
        randCR = crDis(randGenerator);

        // Crossover when random cr is within the cr range.
        if(randCR < cr || i == randIndex)
            newSolutionU[i] = trialV[i];
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 8: DE/rand-to-best/1/bin.
 *
 * @note xBest != xRand1 != xRand2 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param lambda A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy8(const double &cr, const double &F, const double &lambda, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> dimDis(0, x.size()-1);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(x[i] + lambda * (xBest[i] - x[i]) + F * (xRand1[i] - xRand2[i]));

    // Generate a random index withing range of x-dimensions.
    int randIndex = dimDis(randGenerator);

    // Declare a random crossover probability variable.
    double randCR;

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    for(int i = 0; i < trialV.size(); i++)
    {
        // Generate a random crossover probability.
        randCR = crDis(randGenerator);

        // Crossover when random cr is within the cr range.
        if(randCR < cr || i == randIndex)
            newSolutionU[i] = trialV[i];
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 9: DE/best/2/bin.
 *
 * @note xBest != xRand1 != xRand2 != xRand3 != xRand4 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xBest The best solution of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param xRand4 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy9(const double &cr, const double &F, vector<double> x, vector<double> xBest, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> dimDis(0, x.size()-1);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xBest[i] + F * (xRand1[i] + xRand2[i] - xRand3[i] - xRand4[i]));

    // Generate a random index withing range of x-dimensions.
    int randIndex = dimDis(randGenerator);

    // Declare a random crossover probability variable.
    double randCR;

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    for(int i = 0; i < trialV.size(); i++)
    {
        // Generate a random crossover probability.
        randCR = crDis(randGenerator);

        // Crossover when random cr is within the cr range.
        if(randCR < cr || i == randIndex)
            newSolutionU[i] = trialV[i];
    }

    // Return the new solution.
    return newSolutionU;
}

// --------------------------------------------------------------------------------------------------------------------------------------
// **************************************************************************************************************************************
// --------------------------------------------------------------------------------------------------------------------------------------

/**
 * @brief Strategy 10: DE/rand/2/bin.
 *
 * @note xRand1 != xRand2 != xRand3 != xRand4 != xRand5 (where "!=" means "not equal to").
 *
 * @param cr The crossover probability.
 * @param F A scaling factor.
 * @param x The current solution (individual) of the population.
 * @param xRand1 A randomly chosen solution from the population.
 * @param xRand2 A randomly chosen solution from the population.
 * @param xRand3 A randomly chosen solution from the population.
 * @param xRand4 A randomly chosen solution from the population.
 * @param xRand5 A randomly chosen solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new solution (individual).
 */
vector<double> de_Strategy10(const double &cr, const double &F, vector<double> x, vector<double> xRand1, vector<double> xRand2, vector<double> xRand3, vector<double> xRand4, vector<double> xRand5, mt19937 &randGenerator)
{
    // Create a Crossover Distribution for Mersenne Twister.
    uniform_real_distribution<double> crDis(0.0, 1.0);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> dimDis(0, x.size()-1);

    // Declare a trial vector and new solution vector.
    vector<double> trialV;
    vector<double> newSolutionU;

    // Mutation.
    for(int i = 0; i < x.size(); i++)
        trialV.push_back(xRand5[i] + F * (xRand1[i] + xRand2[i] - xRand3[i] - xRand4[i]));

    // Generate a random index withing range of x-dimensions.
    int randIndex = dimDis(randGenerator);

    // Declare a random crossover probability variable.
    double randCR;

    // Crossover - Make trial vector into a new solution.
    newSolutionU = x;
    for(int i = 0; i < trialV.size(); i++)
    {
        // Generate a random crossover probability.
        randCR = crDis(randGenerator);

        // Crossover when random cr is within the cr range.
        if(randCR < cr || i == randIndex)
            newSolutionU[i] = trialV[i];
    }

    // Return the new solution.
    return newSolutionU;
}