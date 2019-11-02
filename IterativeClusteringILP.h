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

#ifndef ITERATIVE_CLUSTERING_ILP_H
#define ITERATIVE_CLUSTERING_ILP_H


class GT2_EXPORT IterativeClusteringILP {
private:

    double threshold;

    GeneTrail::DenseMatrix similarity_matrix;
    GeneTrail::DenseMatrix connectivity_matrix;

    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloExpr obj_expr;
    IloObjective obj;

    // Variables
    IloNumVarArray y_ij;
    IloNumVarArray z_ij;	
    // Solution
    IloNumArray sol_y_ij;
    IloNumArray sol_z_ij;

    // Check if variable is already set
    std::vector<bool> isVariableUsed;

    // Save all constraints in order to remove them later
    IloExtractableArray constraints_yij;
    IloExtractableArray constraints_zij;
    std::map<std::string, IloConstraint> edge_constraints;

    const std::string name_delim_cplex{"____"};

    auto getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int;

public:

    explicit IterativeClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, double threshold);

    auto initializeModel() -> bool;
    auto addVariableConstraints(size_t i, size_t j) -> bool;
    auto updateCycleConstraints(size_t i, size_t j, size_t k) -> bool;
    auto updateCycleConstraints(size_t i, size_t j) -> bool;
    auto updateObjectiveFunction(size_t i, size_t j) -> bool;
    auto extendModel(const std::vector<std::tuple<size_t, size_t, double>>& edges, size_t begin, size_t end) -> bool;
    
    auto solveModel() -> bool;
    
    GeneTrail::DenseMatrix getSolution();
};


#endif //ITERATIVE_CLUSTERING_ILP_H
