cmake_minimum_required(VERSION 3.13)
project(EvolutionaryAlgorithms)

set(CMAKE_CXX_STANDARD 14)

add_executable(EvolutionaryAlgorithms main.cpp BenchmarkFunctions.cpp BenchmarkFunctions.h utilities.cpp utilities.h GeneticAlgorithm.cpp GeneticAlgorithm.h EA_Utilities.cpp EA_Utilities.h DEA_Strategies.cpp DEA_Strategies.h DifferentialEvolution.cpp DifferentialEvolution.h)