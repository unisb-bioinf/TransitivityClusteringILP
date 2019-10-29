#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "DenseMatrixReader.h"
#include "DenseMatrix.h"
#include "macros.h"
#define IL_STD
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#ifndef ClusteringILP_ClusteringILP2_H
#define ClusteringILP_ClusteringILP2_H


class GT2_EXPORT ClusteringILP2{
private:

    double threshold;

    GeneTrail::DenseMatrix similarity_matrix;
    GeneTrail::DenseMatrix connectivity_matrix;

    IloEnv env;
    IloModel model;
    IloCplex cplex;

    // Variables
    IloNumVarArray y_ij;
    IloNumVarArray z_ij;	
    // Solution
    IloNumArray sol_y_ij;
    IloNumArray sol_z_ij;

    // Save all constraints in order to remove them later
    IloExtractableArray constraints_yij;
    IloExtractableArray constraints_zij;
    std::map<std::string, IloConstraint> edge_constraints;


    bool sol_requested{false}; //needed to ensure that solutions are only generated once
    const std::string name_delim_cplex{"____"};

    auto getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int;

    auto doubleEpsilonEqual(double a, double b) -> bool;

public:

    explicit ClusteringILP2(GeneTrail::DenseMatrix & similarity_matrix, double threshold);

    auto initializeModel() -> bool;
    auto updateVariableConstraints(size_t i, size_t j) -> bool;
    auto updateCycleConstraints(size_t i, size_t j, size_t k) -> bool;
    auto updateCycleConstraints(size_t i, size_t j) -> bool;
    auto extendModel(const std::vector<std::tuple<size_t, size_t, double>>& edges, size_t begin, size_t end) -> bool;
    //auto buildModel() -> bool;
    auto solveModel(size_t i) -> bool;
    GeneTrail::DenseMatrix getSolution();
};


#endif //ClusteringILP_ClusteringILP2_H
