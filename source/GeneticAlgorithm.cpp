/**
 * @file GeneticAlgorithm.cpp
 * @author  Al Timofeyev
 * @date    April 25, 2019
 * @brief   Implementation of the Genetic Algorithm.
 */

#include "GeneticAlgorithm.h"


// **********************************************************************************************************
// ************************************************ PUBLIC **************************************************
// **********************************************************************************************************
// ---------------------- CONSTRUCTORS ----------------------
/**
 * @brief The only constructor available to initialize the Genetic Algorithm.
 *
 * @note No Default (zero-param) Constructor exists.
 *
 * @param dim The number of dimensions each individual in the population has.
 * @param sol The population size.
 * @param gen The max number of generations possible.
 * @param cr The crossover probability.
 * @param mProb The mutation probability.
 * @param mR The mutation range.
 * @param mPrec The mutation precision.
 * @param er The elitism rate.
 * @param selectionID The selection type to use for the Genetic Algorithm.
 * @param crPoints The number of crossover points.
 */
GeneticAlgorithm::GeneticAlgorithm(int dim, int sol, int gen, double cr, double mProb, double mR, double mPrec, double er, int selectionID, int crPoints)
{
    gaConfig.dimensions = dim;
    gaConfig.solutions = sol;
    gaConfig.generations = gen;
    gaConfig.cr = cr;
    gaConfig.mutProb = mProb;
    gaConfig.mutRange = mR;
    gaConfig.mutPrec = mPrec;
    gaConfig.er = er;
    gaConfig.eliteIndex = er * sol;
    gaConfig.selectionID = selectionID;
    gaConfig.crPoints = crPoints;
}

// ------------------------- METHODS ------------------------
/**
 * @brief Sets up the population parameters.
 *
 * @note The bounds vector list must have min bound at index 0 and max bound at index 1.
 *
 * @param functionID The ID of which benchmark function to use (IDs: 1 - 18).
 * @param bounds A pair of min/max boundaries for the individuals in the population.
 */
void GeneticAlgorithm::setPopulationParams(int functionID, vector<double> bounds)
{
    population = GA_Population();
    population.functionID = functionID;
    population.bounds = bounds;
}

/**
 * @brief Runs the Genetic Algorithm with set parameters.
 *
 * @return Returns the best result of the Genetic Algorithm.
 */
double GeneticAlgorithm::runGeneticAlgorithm()
{
    // Create a Mersenne Twister pseudo-random number generator.
    mt19937 randGenerator(time(NULL));

    // Record the start time of Genetic Algorithm.
    auto startTime = chrono::high_resolution_clock::now();

    // Generate the initial population and evaluate it.
    generateRandPopulation();
    evaluatePopulation(population.functionID, population.pop, population.fitness, population.functionCounter);

    // Starting generation to 1.
    int generation = 1;

    // Start the Genetic Algorithm.
    while(generation <= gaConfig.generations)
    {
        // Declare a new population and it's fitness vector.
        vector<vector<double>> newPopulation;
        vector<double> newFitness;
        newFitness.resize(population.fitness.size());

        // Normalize the fitness values.
        normalizeFitness(population.functionID, population.fitness);

        for(int i = 0; i < gaConfig.solutions; i += 2)
        {
            // Select parents.
            vector<int> parentIndex = select(gaConfig.selectionID, population.fitness, randGenerator);

            // Perform Crossover.
            vector<vector<double>> children = crossover(gaConfig.crPoints, population.pop[parentIndex[0]], population.pop[parentIndex[1]], gaConfig.cr, randGenerator);

            // Mutate and add children to population.
            for(int cIndex = 0; cIndex < children.size(); cIndex++)
            {
                // Mutate child.
                mutate(children[cIndex], gaConfig.mutProb, gaConfig.mutRange, gaConfig.mutPrec, population.bounds, randGenerator);

                // Add child to new population.
                newPopulation.push_back(children[cIndex]);
            }
        }

        // Evaluate the new population.
        evaluatePopulation(population.functionID, newPopulation, newFitness, population.functionCounter);

        // Generate the new population by combining new population with the old one.
        reduce(population, newPopulation, newFitness, gaConfig.eliteIndex);

        // Save the best solution.
        recordBestFitness(population);

        // Update generation counter.
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

    // Return the best fitness result.
    int bestGenFitSize = population.bestGenFitness.size();
    return population.bestGenFitness[bestGenFitSize-1];
}

/**
 * @brief Analyzes the results of the Genetic Algorithm.
 */
void GeneticAlgorithm::analyzeGAResults()
{
    if(popList.size() == 0)
    {
        cout << "****************************************************************************\n";
        cout << "********** Analysis Could NOT Be Completed - No Data To Analyze ************\n";
        cout << "**********         Please Run The Genetic Algorithm First       ************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Create new analysis object.
    gaAnalysis = GAAnalysis();

    // Perform analysis on all populations stored in popList.
    for(int i = 0; i < popList.size(); i++)
    {

        int fitnessSize = popList[i].bestGenFitness.size();

        // Save the function ID.
        gaAnalysis.functionIDs.push_back(popList[i].functionID);

        // Save the average fitness of data.
        double averageFitness = calculateAverage(popList[i].bestGenFitness);
        gaAnalysis.avgFunctionFitness.push_back(averageFitness);

        // Save the standard deviation fitness of data
        double stdDeviationFitness = calculateStandardDeviation(popList[i].bestGenFitness);
        gaAnalysis.standardDeviation.push_back(stdDeviationFitness);

        // Save the fitness ranges.
        vector<double> range;
        range.push_back(popList[i].bestGenFitness[fitnessSize-1]);
        range.push_back(popList[i].bestGenFitness[0]);
        gaAnalysis.ranges.push_back(range);

        // Save the median fitness of data.
        gaAnalysis.medianFunctionFitness.push_back(popList[i].bestGenFitness[fitnessSize / 2]);

        // Save the execution time of data.
        gaAnalysis.executionTimes.push_back(popList[i].executionTime);

        // Save the function counter.
        gaAnalysis.functionCalls.push_back(popList[i].functionCounter);
    }
}

/**
 * @brief Prints the Results of the Genetic Algorithm.
 */
void GeneticAlgorithm::printGAResults()
{
    cout << "****************************************************************************\n";
    cout << "******* Printing Results of Genetic Algorithm on Current Population ********\n";
    cout << "----------------------------------------------------------------------------\n";

    // If there are no best fitness values from at least 1 generation,
    // then Genetic Algorithm has not been run yet.
    if(population.bestGenFitness.size() == 0)
    {
        cout << "********** NO RESULTS FOR THIS POPULATION\n";
        cout << "********** PLEASE RUN THE GENETIC ALGORITHM\n";
        cout << "----------------------------------------------------------------------------\n\n";
        return;
    }

    cout << "Function ID: " << population.functionID << endl;
    cout << "Generation\t\tBest Fitness of Generation\n";
    cout.precision(12);
    for(int i = 0; i < population.bestGenFitness.size(); i++)
        cout << i << "\t\t\t\t" << population.bestGenFitness[i] << endl;

    cout << "\nElapsed Time (ms) : " << population.executionTime << endl;
    cout << "----------------------------------------------------------------------------\n\n";
}

/**
 * @brief Prints the Analysis of the Genetic Algorithm.
 */
void GeneticAlgorithm::printGAAnalysis()
{
    cout << "\n\n********************************************************\n";
    cout << "************** Printing Analysis Results ***************\n";
    cout << "--------------------------------------------------------\n";

    cout << "Function ID\t\tAverage Fitness\t\t\tStandard Deviation\t\t\tRange(min)\t\t\tRange(max)\t\t\t\tMedian\t\t\t\tTime(ms)\t\t\tFunction Calls\n";
    cout.precision(12);
    for(int row = 0; row < gaAnalysis.functionIDs.size(); row++)
    {
        // Print function ID.
        cout << gaAnalysis.functionIDs[row] << "\t\t\t\t";

        // Print average fitness.
        if(gaAnalysis.avgFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << gaAnalysis.avgFunctionFitness[row] << "\t\t\t";

        // Print the standard deviation.
        if(gaAnalysis.standardDeviation[row] >= 0.0)
            cout << " ";
        cout << gaAnalysis.standardDeviation[row] << "\t\t\t";

        // Print the range.
        if(gaAnalysis.ranges[row][0] >= 0.0)
            cout << " ";
        cout << gaAnalysis.ranges[row][0] << "\t\t\t";
        if(gaAnalysis.ranges[row][1] >= 0.0)
            cout << " ";
        cout << gaAnalysis.ranges[row][1] << "\t\t\t";

        // Print the median.
        if(gaAnalysis.medianFunctionFitness[row] >= 0.0)
            cout << " ";
        cout << gaAnalysis.medianFunctionFitness[row] << "\t\t\t";

        // Print the Time in milliseconds.
        cout << gaAnalysis.executionTimes[row] << "\t\t\t";

        // Print the number of function calls.
        cout << gaAnalysis.functionCalls[row] << "\n";
    }

    cout << "********************************************************\n\n";
}

/**
 * @brief Saves all Genetic Algorithm Results to file.
 *
 * @param iniFilename   The name of the initialization file that was used to set up
 *                      the population parameters of the Genetic Algorithm.
 */
void GeneticAlgorithm::saveGAResults(string iniFilename)
{
    // If the popList is empty, exit the function.
    if(popList.size() == 0)
    {
        cout << "\n******************************************************\n";
        cout << "THERE IS NO GENETIC ALGORITHM DATA TO SAVE";
        cout << "\n******************************************************\n\n";
        return;
    }

    // Setup the output filename.
    string filename = "GeneticAlgorithm-Results.csv";

    // Initialize the number of rows (generations.).
    int rows = gaConfig.generations;

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
    if(gaAnalysis.functionIDs.size() > 0)
    {
        line = "Average,";
        for(int i = 0; i < gaAnalysis.avgFunctionFitness.size(); i++)
        {
            line += to_string(gaAnalysis.avgFunctionFitness[i]);

            if(i == gaAnalysis.avgFunctionFitness.size()-1)
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
 * @brief Saves the Analysis of the Genetic Algorithm to file.
 */
void GeneticAlgorithm::saveGAAnalysis()
{
    if(gaAnalysis.functionIDs.size() == 0)
    {
        cout << "****************************************************************************\n";
        cout << "******************** There Is No Analysis Data To Save *********************";
        cout << "\n****************************************************************************\n\n";
        return;
    }

    // Rows.
    int rows = gaAnalysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create the file to where the matrix is saved.
    ofstream outputFile;
    outputFile.open ("GeneticAlgorithm-Analysis.csv");

    // Save the header line first.
    outputFile << gaAnalysis.header;

    // Save data to file.
    string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += to_string(gaAnalysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += to_string(gaAnalysis.avgFunctionFitness[row]) + ",";

        // Save the standard deviation.
        line += to_string(gaAnalysis.standardDeviation[row]) + ",";

        // Save the range.
        line += to_string(gaAnalysis.ranges[row][0]) + ",";
        line += to_string(gaAnalysis.ranges[row][1]) + ",";

        // Save the median.
        line += to_string(gaAnalysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += to_string(gaAnalysis.executionTimes[row]) + ",";

        // Save the Function Calls Counter.
        line += to_string(gaAnalysis.functionCalls[row]) + "\n";

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
/**
 * @brief Generates a random initial population for the Genetic Algorithm.
 *
 * @note Makes function call to EA_Utilities.h --> createMatrix().
 */
void GeneticAlgorithm::generateRandPopulation()
{
    population.pop = createGAMatrix(gaConfig.solutions, gaConfig.dimensions, population.bounds[0], population.bounds[1]);
    population.fitness.resize(population.pop.size());
}

/**
 * @brief Calculates fitness of all solutions in population.
 *
 * @note Makes function call to EA_Utilities.h --> calculateFitnessOfVector().
 *
 * @param functionID The ID of the benchmark function to use.
 * @param pop The matrix population of the Genetic Algorithm.
 * @param fitness The fitness vector for each solution from the population.
 * @param functionCounter A counter to keep track of how many times fitness function was called.
 */
void GeneticAlgorithm::evaluatePopulation(int functionID, vector<vector<double>> &pop, vector<double> &fitness, int &functionCounter)
{
    for(int i = 0; i < pop.size(); i++)
    {
        fitness[i] = calculateFitnessOfVector(pop[i], functionID);
        functionCounter++;
    }
}

/**
 * @brief Selects two parents from the population.
 *
 * @note Makes function call to EA_Utilities.h --> printAllGASelectionIDs().
 *
 * @param selectionID The ID of the selection type to use.
 * @param popFitness The fitness list of the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return A list of parent indices.
 */
vector<int> GeneticAlgorithm::select(int selectionID, vector<double> &popFitness, mt19937 &randGenerator)
{
    // Set the tournament size (Only for Tournament Selection).
    int tournamentSize = popFitness.size()/2;

    // Start the selection process.
    vector<int> parentIndex;
    for(int i = 0; i < 2; i++)
    {
        switch(selectionID)
        {
            // Tournament Selection.
            case 1:
                parentIndex.push_back(t_select(popFitness, tournamentSize, randGenerator));
                break;
            // Roulette Wheel Selection.
            case 2:
                parentIndex.push_back(rw_select(popFitness, randGenerator));
                break;
            default:
                cout << "\n***********************************************\n";
                cout << "Invalid Selection ID was used in GA-Config.txt";
                printAllGASelectionIDs();
                cout << "******** TERMINATING PROGRAM EXECUTION ********\n\n";
                exit(1);
        }
    }

    // Return the parent indices.
    return parentIndex;
}

/**
 * @brief Uses Tournament Selection to select one parent from population.
 *
 * @note All fitness values have been normalized to between 0 to 1.
 *
 * @param popFitness The fitness list of the population.
 * @param numOfInid The number of individuals to select for tournament (max = size of population).
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return The index of the parent in the population.
 */
int GeneticAlgorithm::t_select(vector<double> &popFitness, int numOfInid, mt19937 &randGenerator)
{
    // Create a Mersenne Twister pseudo-random number generator.
    uniform_int_distribution<int> dis(0, popFitness.size()-1);

    // Select the specified number of individuals for the tournament.
    vector<double> tourFitList;
    vector<int> tourIndexList;
    for(int i = 0, index; i < numOfInid; i++)
    {
        index = dis(randGenerator);                 // Generate random index in range.
        tourIndexList.push_back(index);             // Save the index generated.
        tourFitList.push_back(popFitness[index]);   // Save the fitness at that index.
    }

    // Set the initial values for the parent (best) fitness and its index.
    int parentIndex = tourIndexList[0];
    double parentFitness = tourFitList[0];

    // Select best (minimum) fitness from tournament.
    for(int i = 1; i < tourFitList.size(); i++)
    {
        // If the next fitness in the tournament is better than the current one.
        // We're looking at normalized fitness values, not actual fitness values.
        if(tourFitList[i] > parentFitness)
        {
            parentIndex = tourIndexList[i]; // Update the parent Index.
            parentFitness = tourFitList[i]; // Update the parent Fitness.
        }
    }

    // Return the parent index of the parent.
    return parentIndex;
}

/**
 * @brief Uses Roulette Wheel Selection to select one parent from population.
 *
 * @note All fitness values have been normalized to between 0 to 1.
 *
 * @param popFitness The fitness list of the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return The index of the parent in the population.
 */
int GeneticAlgorithm::rw_select(vector<double> &popFitness, mt19937 &randGenerator)
{
    // Sum up the fitness values.
    double fitnessSum = 0;
    for(int i = 0; i < popFitness.size(); i++)
        fitnessSum += popFitness[i];

    // Calculate the probability of each fitness value.
    vector<double> probabilities;
    double probability;
    for(int i = 0; i < popFitness.size(); i++)
    {
        probability = popFitness[i] / fitnessSum;
        probabilities.push_back(probability);
    }

    // Create a Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Generate a random probability.
    double randProbability = dis(randGenerator);

    // Start Roulette Wheel Selection.
    double probabilitySum = 0;
    int parentIndex = 0;
    for(int i = 0; i < probabilities.size(); i = ((i+1) % probabilities.size()))
    {
        // Add probability to probability sum.
        probabilitySum += probabilities[i];

        // The first probability found that's greater than the random probability
        // gets chosen as the parent.
        if(probabilitySum > randProbability)
        {
            parentIndex = i;    // Update parent index.
            break;
        }
    }

    // Return the parent index.
    return parentIndex;
}

/**
 * @brief Creates two children using the two parent individuals.
 *
 * @param crPoints The number of crossover points.
 * @param parent1 The first parent.
 * @param parent2 The second parent.
 * @param cr The crossover probability.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return Returns two children that are a crossover of the parents.
 */
vector<vector<double>> GeneticAlgorithm::crossover(int crPoints, vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator)
{
    switch(crPoints)
    {
        case 1:
            return crossover1(parent1, parent2, cr, randGenerator);
        case 2:
            return crossover2(parent1, parent2, cr, randGenerator);
        default:
            cout << "\n*****************************************************\n";
            cout << "******* Invalid Crossover Point Was Selected ******** \n";
            cout << "*********** TERMINATING PROGRAM EXECUTION ***********\n";
            cout << "*****************************************************\n\n";
            exit(1);
    }

}

/**
 * @brief Returns two children that are a crossover of the parents using 1 crossover point.
 *
 * @param parent1 The first parent.
 * @param parent2 The second parent.
 * @param cr The crossover probability.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return Returns two children that are a crossover of the parents.
 */
vector<vector<double>> GeneticAlgorithm::crossover1(vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator)
{
    // Declare children vectors.
    vector<vector<double>> children;
    vector<double> child1;
    vector<double> child2;

    // Create a Mersenne Twister pseudo-random number generator for random CR.
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Create a Mersenne Twister pseudo-random number generator for crossover index.
    uniform_int_distribution<int> crossoverDis(0, parent1.size()-1);

    // Generate a random crossover probability.
    double randCR = dis(randGenerator);

    // Declare a crossover index.
    int crIndex;

    // If the random crossover is less than the initial crossover probability.
    if(randCR < cr)
    {
        // Assign a crossover index.
        crIndex = crossoverDis(randGenerator);

        // Add the first half of each parent two it's respective child.
        for(int i = 0; i < crIndex; i++)
        {
            child1.push_back(parent1[i]);
            child2.push_back(parent2[i]);
        }
        // Add the second half of each parent to the opposite child.
        for(int i = crIndex; i < parent1.size(); i++)
        {
            child1.push_back(parent2[i]);
            child2.push_back(parent1[i]);
        }
    }
        // Else the children are the same as the parents.
    else
    {
        child1 = parent1;
        child2 = parent2;
    }

    // Add the two children the the list.
    children.push_back(child1);
    children.push_back(child2);

    // Return the two children.
    return children;
}

/**
 * @brief Returns two or six children that are a crossover of the parents using 2 crossover points.
 *
 * @param parent1 The first parent.
 * @param parent2 The second parent.
 * @param cr The crossover probability.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 *
 * @return Returns two or six children that are a crossover of the parents.
 */
vector<vector<double>> GeneticAlgorithm::crossover2(vector<double> parent1, vector<double> parent2, double cr, mt19937 &randGenerator)
{
    // Declare children vectors.
    vector<vector<double>> children;
    vector<double> child1;
    vector<double> child2;
    vector<double> child3;
    vector<double> child4;
    vector<double> child5;
    vector<double> child6;

    // Create a Mersenne Twister pseudo-random number generator for random CR.
    uniform_real_distribution<double> dis(0.0, 1.0);

    // Create a Mersenne Twister pseudo-random number generator for crossover index.
    uniform_int_distribution<int> crossoverDis(0, parent1.size()/2);
    uniform_int_distribution<int> crossoverDis2((parent1.size()/2) + 1, parent1.size()-1);

    // Generate a random crossover probability.
    double randCR = dis(randGenerator);

    // Declare a crossover indices.
    int crIndex;
    int crIndex2;

    // If the random crossover is less than the initial crossover probability.
    if(randCR < cr)
    {
        // Assign crossover indices.
        crIndex = crossoverDis(randGenerator);
        crIndex2 = crossoverDis2(randGenerator);

        // Add the first half of each parent two it's respective child. (head)
        for(int i = 0; i < crIndex; i++)
        {
            child1.push_back(parent1[i]);
            child2.push_back(parent2[i]);
            child3.push_back(parent2[i]);

            child4.push_back(parent2[i]);
            child5.push_back(parent1[i]);
            child6.push_back(parent1[i]);
        }
        // Add the second half of each parent to the opposite child. (middle)
        for(int i = crIndex; i < crIndex2; i++)
        {
            child1.push_back(parent2[i]);
            child2.push_back(parent1[i]);
            child3.push_back(parent2[i]);

            child4.push_back(parent1[i]);
            child5.push_back(parent2[i]);
            child6.push_back(parent1[i]);
        }
        // Add the third half of each parent to the opposite child. (tail)
        for(int i = crIndex2; i < parent1.size(); i++)
        {
            child1.push_back(parent2[i]);
            child2.push_back(parent2[i]);
            child3.push_back(parent1[i]);

            child4.push_back(parent1[i]);
            child5.push_back(parent1[i]);
            child6.push_back(parent2[i]);
        }
    }
        // Else the children are the same as the parents.
    else
    {
        child1 = parent1;
        child2 = parent2;
    }

    // Add the two children the the list.
    children.push_back(child1);
    children.push_back(child2);

    // If the 6th child is not empty, add the other four children.
    if(child6.size() > 0)
    {
        children.push_back(child3);
        children.push_back(child4);
        children.push_back(child5);
        children.push_back(child6);
    }

    // Return the two children.
    return children;
}

/**
 * @brief Mutates the child.
 *
 * @param child The child vector to be mutated.
 * @param mProb Mutation Probability.
 * @param mRange Mutation Range.
 * @param mPrec Mutation Precision.
 * @param bounds The min/max bounds of the population.
 * @param randGenerator A Mersenne Twister pseudo-random number generator.
 */
void GeneticAlgorithm::mutate(vector<double> &child, double mProb, double mRange, double mPrec, vector<double> bounds, mt19937 &randGenerator)
{
    // Create a Mersenne Twister pseudo-random number generator for random mutation probability.
    uniform_real_distribution<double> randMProb(0.0, 1.0);

    // Create a Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis(-1.0, 1.0);

    // Create a Mersenne Twister pseudo-random number generator.
    uniform_real_distribution<double> dis2(0.0, 1.0);

    // Declare variables needed for mutation loop.
    double randMutProb;
    double randTemp1;
    double randTemp2;

    // Randomly mutate some parts of the child.
    for(int i = 0; i < child.size(); i++)
    {
        // Generate a random mutation probability.
        randMutProb = randMProb(randGenerator);

        // If the random mutation probability is less than initial
        // mutation probability, then mutate the child's value at index i.
        if(randMutProb < mProb)
        {
            randTemp1 = dis(randGenerator);
            randTemp2 = dis2(randGenerator);
            child[i] += randTemp1 * (bounds[1] - bounds[0]) * mRange * pow(2, (-1 * randTemp2 * mPrec));
        }
    }
}

/**
 * @brief Combines new population with the old one.
 *
 * @note Makes function call to EA_Utilities.h --> quicksort().
 *
 * @param pop The GA_Population population structure.
 * @param newPop The new population.
 * @param newFit The new set of fitness values for the new population.
 * @param eliteIndex The ending index of how much of the old population to rertain.
 */
void GeneticAlgorithm::reduce(GA_Population &pop, vector<vector<double>> &newPop, vector<double> &newFit, int eliteIndex)
{
    // Re-evaluate old population.
    evaluatePopulation(pop.functionID, pop.pop, pop.fitness, pop.functionCounter);

    // Sort the old populations based on fitness.
    quicksort(pop.fitness, pop.pop, 0, pop.fitness.size()-1);

    // Add the best <eliteIndex> solutions from the old population to the new population.
    // DON'T REPLACE ANYTHING YET!!
    for(int i = 0; i < eliteIndex; i++)
    {
        newPop.push_back(pop.pop[i]);
        newFit.push_back(pop.fitness[i]);
    }

    // Evaluate the new population, and sort it.
    //evaluatePopulation(pop.functionID, newPop, newFit, pop.functionCounter);
    quicksort(newFit, newPop, 0, newFit.size()-1);

    // Resize the new population, leaving only the most optimal values.
    newPop.resize(pop.pop.size());
    newFit.resize(pop.fitness.size());

    // Set the new population data.
    pop.pop = newPop;
    pop.fitness = newFit;
}

/**
 * @brief Normalizes all fitness values for a population.
 *
 * @param functionID The ID of the benchmark function to normalize to.
 * @param popFitness A vector of fitness values from the population.
 */
void GeneticAlgorithm::normalizeFitness(int functionID, vector<double> &popFitness)
{
    // Declare a shift value.
    // Shift value is optimal solution for a function.
    // Basically, shifting all the population fitness values
    // by the shift value will make sure that the optimal solution
    // would then be located at 0. After that, we can normalize.
    double shiftValue;

    switch(functionID)
    {
        case 6:
            // Shift value is optimal solution for this function.
            shiftValue = -1.4915 * (gaConfig.solutions - 1);

            // Normalize the fitness values of the population using shiftValue.
            for (int i = 0; i < popFitness.size(); i++)
            {
                // First, shift the population fitness value towards 0.
                popFitness[i] -= shiftValue;

                // Then normalize it.
                if (popFitness[i] >= 0)
                    popFitness[i] = 1 / (1 + popFitness[i]);
                else
                    popFitness[i] = 1 / (1 + fabs(popFitness[i]));
            }
            break;
        case 8:
            // Shift value is optimal solution for this function.
            shiftValue = -7.54276 - 2.91867 * (gaConfig.solutions - 3);

            // Normalize the fitness values of the population using shiftValue.
            for (int i = 0; i < popFitness.size(); i++)
            {
                // First, shift the population fitness value towards 0.
                popFitness[i] -= shiftValue;

                // Then normalize it.
                if (popFitness[i] >= 0)
                    popFitness[i] = 1 / (1 + popFitness[i]);
                else
                    popFitness[i] = 1 / (1 + fabs(popFitness[i]));
            }
            break;
        case 13:
            // Shift value is optimal solution for this function.
            shiftValue = 0.966 * gaConfig.solutions;

            // Normalize the fitness values of the population using shiftValue.
            for (int i = 0; i < popFitness.size(); i++)
            {
                // First, shift the population fitness value towards 0.
                popFitness[i] -= shiftValue;

                // Then normalize it.
                if (popFitness[i] >= 0)
                    popFitness[i] = 1 / (1 + popFitness[i]);
                else
                    popFitness[i] = 1 / (1 + fabs(popFitness[i]));
            }
            break;
        case 14:
            // Shift value is optimal solution for this function.
            shiftValue = 1 - gaConfig.solutions;

            // Normalize the fitness values of the population using shiftValue.
            for (int i = 0; i < popFitness.size(); i++)
            {
                // First, shift the population fitness value towards 0.
                popFitness[i] -= shiftValue;

                // Then normalize it.
                if (popFitness[i] >= 0)
                    popFitness[i] = 1 / (1 + popFitness[i]);
                else
                    popFitness[i] = 1 / (1 + fabs(popFitness[i]));
            }
            break;
        default:
            // Normalize the fitness values of the population.
            for (int i = 0; i < popFitness.size(); i++)
            {
                if (popFitness[i] >= 0)
                    popFitness[i] = 1 / (1 + popFitness[i]);
                else
                    popFitness[i] = 1 / (1 + fabs(popFitness[i]));
            }
            break;
    }
}

/**
 * @brief Saves the best fitness value from the population into the GA_Population population structure.
 *
 * @note Makes function call to EA_Utilities.h --> quicksort().
 *
 * @param pop The GA_Population population structure.
 */
void GeneticAlgorithm::recordBestFitness(GA_Population &pop)
{
    // Sort a copy of the fitness list so as not to change the original.
    vector<double> tempFitnessHolder = pop.fitness;
    quicksort(tempFitnessHolder, 0, tempFitnessHolder.size()-1);

    // Save the best fitness (at index 0) to best fitness list.
    pop.bestGenFitness.push_back(tempFitnessHolder[0]);
}