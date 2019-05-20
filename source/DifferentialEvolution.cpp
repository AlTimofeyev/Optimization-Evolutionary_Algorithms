/**
 * @file DifferentialEvolution.cpp
 * @author  Al Timofeyev
 * @date    May 3, 2019
 * @brief   Implementation of the Differential Evolution Algorithm.
 */

#include "DifferentialEvolution.h"


// **********************************************************************************************************
// ************************************************ PUBLIC **************************************************
// **********************************************************************************************************
// ---------------------- CONSTRUCTORS ----------------------
/**
 * @brief The only constructor available to initialize the Differential Evolution Algorithm.
 *
 * @note No Default (zero-param) Constructor exists.
 *
 * @param dim The number of dimensions each individual in the population has.
 * @param sol The population size.
 * @param gen The max number of generations possible.
 * @param cr The crossover probability.
 * @param f A scaling factor.
 * @param lambda A scaling factor.
 * @param strategy The Differential Evolution strategy to use.
 */
DifferentialEvolution::DifferentialEvolution(int dim, int sol, int gen, double cr, double f, double lambda, int strategy)
{
    // ---------- PERFORM PARAMETER CHECK FIRST ----------
    // Check the population size.
    if(sol < 6)
    {
        cout << "\n***************************************************************\n";
        cout << "*****            Population Size Is Too Small            ******\n";
        cout << "***** Population size must be greater than or equal to 6 ******\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    if(cr < 0 || cr > 1)
    {
        cout << "\n***************************************************************\n";
        cout << "********    Crossover (cr) Value Is Too Big/Small     *********\n";
        cout << "******** Crossover (cr) Value Must Be In Range [0, 1] *********\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    if(f <= 0 || f > 1.2)
    {
        cout << "\n***************************************************************\n";
        cout << "*********    Scaling Factor F Is Too Big/Small       **********\n";
        cout << "********* Scaling Factor F Must Be In Range (0, 1.2] **********\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    if(lambda <= 0 || lambda > 2)
    {
        cout << "\n***************************************************************\n";
        cout << "*******    Scaling Factor Lambda Is Too Big/Small     *********\n";
        cout << "******* Scaling Factor Lambda Must Be In Range (0, 2] *********\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    if(strategy < 1 || lambda > 10)
    {
        cout << "\n***************************************************************\n";
        cout << "*********    The DE Strategy Is Too Big/Small      ************\n";
        cout << "********* The DE Strategy Must Be In Range [1, 10] ************\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    // If all parameter checks pass successfully, then configure the Differential Evolution Algorithm.
    deConfig.dimensions = dim;
    deConfig.NP = sol;
    deConfig.generations = gen;
    deConfig.cr = cr;
    deConfig.f = f;
    deConfig.lambda = lambda;
    deConfig.strategy = strategy;
}

// ------------------------- METHODS ------------------------
/**
 * Runs the Differential Evolution Algorithm with a set of parameters.
 *
 * @param functionID The ID of the benchmark function to use.
 * @param minBound, maxBound The minimum and maximum bounds of the population.
 *
 * @return Returns the best result of the Differential Evolution Algorithm.
 */
double DifferentialEvolution::runDifferentialEvolution(int functionID, double minBound, double maxBound)
{
    // Create a Mersenne Twister pseudo-random number generator.
    mt19937 randGenerator(time(NULL));

    // Record the start time of Differential Evolution Algorithm.
    auto startTime = chrono::high_resolution_clock::now();

    // Create population with given parameters and evaluate and sort it.
    DEA_Population population;
    population.functionID = functionID;
    population.bounds.push_back(minBound);
    population.bounds.push_back(maxBound);
    generateRandDEAPopulation(population);
    evaluatePopulation(functionID, population.pop, population.fitness, population.functionCounter);
    quicksort(population.fitness, population.pop, 0, population.fitness.size()-1);

    // Create a Dimensional Index Distribution for Mersenne Twister.
    uniform_int_distribution<int> indivDis(1, population.pop.size()-1);

    // Start the Differential Evolution Algorithm.
    int generation = 0;
    while(generation < deConfig.generations)
    {
        // For every individual in the population.
        for(int i = 0; i < population.pop.size(); i++)
        {
            // --------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------
            // Construct a list of the best and 5 other random individuals.
            vector<vector<double>> someIndiv;
            vector<int> someIndivChosenIndex;

            // Save the best individual first.
            someIndiv.push_back(population.pop[0]);
            someIndivChosenIndex.push_back(0);

            // Choose 5 more random individuals from the population.
            int tempIndex;
            bool found = false;
            while(someIndiv.size() < 6) // Random individuals must be unique.
            {
                // Generate a random index from the population.
                tempIndex = indivDis(randGenerator);

                for(int j = 0; j < someIndivChosenIndex.size(); j++) // If the index was chosen previously, set found flag.
                {
                    if(tempIndex == someIndivChosenIndex[j])
                    {
                        found = true;
                        break;
                    }
                }

                // If the index was not previously chosen, then add it tho chosen list.
                if(found == false)
                {
                    someIndiv.push_back(population.pop[tempIndex]);
                    someIndivChosenIndex.push_back(tempIndex);
                }

                // Reset the flag to false if it was set to true.
                found = false;
            }
            // --------------------------------------------------------------------------------------------------------
            // --------------------------------------------------------------------------------------------------------

            // Mutate and do Crossover.
            vector<double> newSolution = mutateAndCrossover(deConfig.cr, deConfig.f, deConfig.lambda, population.pop[i], someIndiv, randGenerator);

            // Select an individual to be part of the next generation.
            select(functionID, population.pop[i], population.fitness[i], newSolution, population.functionCounter);
        }

        // Sort the new population;
        quicksort(population.fitness, population.pop, 0, population.fitness.size()-1);

        // Save the best fitness to list.
        saveBestFitness(population);

        // Increment the generation.
        generation++;
    }

    // Record the end time of Genetic Algorithm.
    auto endTime = chrono::high_resolution_clock::now();

    // Calculate elapsed time in milliseconds it took to execute the benchmark function.
    auto elapsedTime = endTime - startTime;
    double elapsedTimeMS = chrono::duration_cast<chrono::milliseconds>(elapsedTime).count();

    // Save elapsed time to the population.
    population.executionTime = elapsedTimeMS;

    // Add the population to the population list.
    popList.push_back(population);

    // Return the best fitness of the Differential Evolution Algorithm.
    int bestFitnessSize = population.bestGenFitness.size();
    return population.bestGenFitness[bestFitnessSize-1];
}

/**
 * @brief Analyzes the results of the Differential Evolution Algorithm.
 */
void DifferentialEvolution::analyzeDEAResults()
{
    if(popList.size() == 0)
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze  ***********\n";
        cout << "********** Please Run The Differential Evolution Algorithm First ***********";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new analysis object.
    deAnalysis = DEAAnalysis();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].bestGenFitness.size();

        // Save the function ID.
        deAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].bestGenFitness);
        deAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].bestGenFitness);
        deAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].bestGenFitness[fitnessSize-1]);
        range.push_back(popList[i].bestGenFitness[0]);
        deAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        deAnalysis.medianFunctionFitness.push_back(popList[i].bestGenFitness[fitnessSize / 2]);

        // Save the execution time of data.
        deAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        deAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Prints the Results of the Differential Evolution Algorithm.
 */
void DifferentialEvolution::printDEAResults()
{
    cout << "****************************************************************************\n";
    cout << "**** Printing Results of Differential Evolution on Current Population ******\n";
    cout << "----------------------------------------------------------------------------\n";

    // If the popList is empty, then Differential Evolution Algorithm has not been run yet.
    if(popList.size() == 0)
    {
        cout << "********** NO RESULTS FOR THIS POPULATION\n";
        cout << "********** PLEASE RUN THE DIFFERENTIAL EVOLUTION\n";
        cout << "----------------------------------------------------------------------------\n\n";
        return;
    }

    for(int i = 0; i < popList.size(); i++)
    {

        cout << "Function ID: " << popList[i].functionID << endl;
        cout << "Generation\t\tBest Fitness of Generation\n";
        cout.precision(12);
        for (int i = 0; i < popList[i].bestGenFitness.size(); i++)
            cout << i << "\t\t\t\t" << popList[i].bestGenFitness[i] << endl;

        cout << "\nElapsed Time (ms) : " << popList[i].executionTime << endl;
        cout << "----------------------------------------------------------------------------\n\n";
    }
}

/**
 * @brief Prints the Analysis of the Differential Evolution Algorithm.
 */
void DifferentialEvolution::printDEAAnalysis()
{
    cout << "\n\n********************************************************\n";
    cout << "************** Printing Analysis Results ***************\n";
    cout << "--------------------------------------------------------\n";

    cout << "Function ID\t\tAverage Fitness\t\t\tStandard Deviation\t\t\tRange(min)\t\t\tRange(max)\t\t\t\tMedian\t\t\t\tTime(ms)\t\t\tFunction Calls\n";
    cout.precision(12);
    for(int row = 0; row < deAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << deAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(deAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << deAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(deAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << deAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(deAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << deAnalysis.ranges[row][0] << "\t\t\t";
        if(deAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << deAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(deAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << deAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << deAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << deAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Saves all Differential Evolution Algorithm Results to file.
 */
void DifferentialEvolution::saveDEAResults()
{
    // If the popList is empty, exit the function.
    if(popList.size() == 0)
    {
        cout << "\n******************************************************\n";
        cout << "THERE IS NO DIFFERENTIAL EVOLUTION ALGORITHM DATA TO SAVE";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "DifferentialEvolution-Results-Strategy" + to_string(deConfig.strategy) + ".csv";

    // Initialize the number of rows (generations.).
    int rows = deConfig.generations;

    // Create the file to where the matrix is saved.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    string header = "Generations,";
    for(int pIndex = 0; pIndex < popList.size(); pIndex++)
    {
        header += "f" + to_string(popList[pIndex].functionID) + "(x)";
        if(pIndex == popList.size()-1)
            header += "\n";
        else
            header += ",";
    }
    outputFile << header;

    // Save the data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the generation.
        line += to_string(row) + ",";

        // Save the best fitness from generation <row> of each population.
        for(int pIndex = 0; pIndex < popList.size(); pIndex++)
        {
            line += to_string(popList[pIndex].bestGenFitness[row]);

            if(pIndex == popList.size()-1)
                line += "\n";
            else
                line += ",";
        }

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Save the averages if they exist.
    if(deAnalysis.functionIDs.size() > 0)
    {
        line = "Average,";
        for(int i = 0; i < deAnalysis.avgFunctionFitness.size(); i++)
        {
            line += to_string(deAnalysis.avgFunctionFitness[i]);

            if(i == deAnalysis.avgFunctionFitness.size()-1)
                line += "\n";
            else
                line += ",";
        }
        outputFile << line;
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Saves the Analysis of the Differential Evolution Algorithm to file.
 */
void DifferentialEvolution::saveDEAAnalysis()
{
    if(deAnalysis.functionIDs.size() == 0)
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Save *********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = deAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create filename based on strategy used.
    string filename = "DifferentialEvolution-Analysis-Strategy" + to_string(deConfig.strategy) + ".csv";


    // Create the file to where the matrix is saved.
    ofstream outputFile;
    outputFile.open (filename);

    // Save the header line first.
    outputFile << deAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(deAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(deAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(deAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(deAnalysis.ranges[row][0]) + ",";
        line += to_string(deAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(deAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(deAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(deAnalysis.functionCalls[row]) + ",";

        // Save the strategy used.
        line += to_string(deConfig.strategy) + "\n";

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}



// **********************************************************************************************************
// ************************************************ PRIVATE *************************************************
// **********************************************************************************************************
// ------------------------- METHODS ------------------------
/**
 * @brief Generates the initial population for Differential Evolution Algorithm.
 *
 * @param population The population.
 */
void DifferentialEvolution::generateRandDEAPopulation(DEA_Population &population)
{
    population.pop = createDEAMatrix(deConfig.NP, deConfig.dimensions, population.bounds[0], population.bounds[1]);
    population.fitness.resize(population.pop.size());
}

/**
 * @brief Calculates fitness of all solutions in population.
 *
 * @note Makes function call to EA_Utilities.h --> calculateFitnessOfVector().
 *
 * @param functionID The ID of the benchmark function to use.
 * @param pop The matrix population of the Differential Evolution Algorithm.
 * @param fitness The fitness vector for each solution from the population.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void DifferentialEvolution::evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter)
{
    for(int i = 0; i < pop.size(); i++)
    {
        fitness[i] = calculateFitnessOfVector(pop[i], functionID);
        functionCounter++;
    }
}

/**
 * @brief Calculate the fitness of an individual solution of the population.
 *
 * @note Makes function call to EA_Utilities.h --> calculateFitnessOfVector().
 *
 * @param functionID The ID of the benchmark function to use.
 * @param indiv The individual of the population.
 * @param fitness The fitness variable for the individual.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void DifferentialEvolution::evaluateIndividual(const int &functionID, vector<double> &indiv, double &fitness, int &functionCounter)
{
    fitness = calculateFitnessOfVector(indiv, functionID);
    functionCounter++;
}

/**
 * @brief Mutate and Crossover produces one new individual.
 *
 * @note Makes function call to DEA_Strategies.h.
 * @note indiv list of solution has the best solution at index 0 and 5 random solutions at indices 1 - 5.
 *
 * @param cr The crossover probability.
 * @param f A scaling factor.
 * @param lambda A scaling factor.
 * @param x The initial solution of the population.
 * @param indiv A list of the best solution and 5 random solution from the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A new potential individual of the population.
 */
vector<double> DifferentialEvolution::mutateAndCrossover(const double &cr, const double &f, const double &lambda, vector<double> &x, vector<vector<double>> indiv, mt19937 &randGenerator)
{
    // Make sure there are 6 individuals in the indiv vector list.
    if(indiv.size() < 6)
    {
        cout << "\n***************************************************************\n";
        cout << "***** The indiv List in MutateAndCrossover()Is Too Small ******\n";
        cout << "*****    The Individual List must have 6 Individuals     ******\n";
        cout << "*****             Please Fix The Source Code             ******\n";
        cout << "***************************************************************\n";
        cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
        cout << "***************************************************************\n\n";
        exit(1);
    }

    switch(deConfig.strategy)
    {
        case 1:
            return de_Strategy1(cr, f, x, indiv[0], indiv[2], indiv[3], randGenerator);
        case 2:
            return de_Strategy2(cr, f, x, indiv[1], indiv[2], indiv[3], randGenerator);
        case 3:
            return de_Strategy3(cr, f, lambda, x, indiv[0], indiv[1], indiv[2], randGenerator);
        case 4:
            return de_Strategy4(cr, f, x, indiv[0], indiv[1], indiv[2], indiv[3], indiv[4], randGenerator);
        case 5:
            return de_Strategy5(cr, f, x, indiv[1], indiv[2], indiv[3], indiv[4], indiv[5], randGenerator);
        case 6:
            return de_Strategy6(cr, f, x, indiv[0], indiv[2], indiv[3], randGenerator);
        case 7:
            return de_Strategy7(cr, f, x, indiv[1], indiv[2], indiv[3], randGenerator);
        case 8:
            return de_Strategy8(cr, f, lambda, x, indiv[0], indiv[1], indiv[2], randGenerator);
        case 9:
            return de_Strategy9(cr, f, x, indiv[0], indiv[1], indiv[2], indiv[3], indiv[4], randGenerator);
        case 10:
            return de_Strategy10(cr, f, x, indiv[1], indiv[2], indiv[3], indiv[4], indiv[5], randGenerator);
        default:
            cout << "\n***************************************************************\n";
            cout << "****************** Strategy " << deConfig.strategy << " Does Not Exist ******************\n";
            cout << "***************************************************************\n";
            cout << "---------------- TERMINATING PROGRAM EXECUTION ----------------\n";
            cout << "***************************************************************\n\n";
            exit(1);
    }
}

/**
 * @brief Selects a new individual to be part of the next generation.
 *
 * @param functionID The ID of the benchmark function to use.
 * @param x The original individual (solution) of the population.
 * @param fitness The fitness of the original individual of the population.
 * @param newSol A new potential individual of the population.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void DifferentialEvolution::select(int functionID, vector<double> &x, double &fitness, vector<double> &newSol, int &functionCounter)
{
    // Evaluate the new solution and get it's fitness value.
    double newFitness;
    evaluateIndividual(functionID, newSol, newFitness, functionCounter);

    // If the new solution is better than the original solution,
    // then replace original solution with the new solution.
    if(newFitness <= fitness)
    {
        x = newSol;
        fitness = newFitness;
    }
}

/**
 * @brief Saves the best fitness from the population.
 *
 * @note The best fitness is at the top (index 0) of the population.
 *
 * @param pop
 */
void DifferentialEvolution::saveBestFitness(DEA_Population &pop)
{
    pop.bestGenFitness.push_back(pop.fitness[0]); // Save the best fitness (at index 0) to best fitness list.
}