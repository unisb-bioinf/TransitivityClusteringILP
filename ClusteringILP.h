//
// Created by klenhof on 23.10.19.
//

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <DenseMatrixReader.h>
#include <DenseMatrix.h>
#define IL_STD
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#ifndef ClusteringILP_ClusteringILP_H
#define ClusteringILP_ClusteringILP_H


class ClusteringILP{

    bool valid_input;
    double threshold;

public:
    enum Metric{similarity, distance};
private:
    ClusteringILP::Metric sim_or_dist;//specify whether similarity or distance measure was used to create the 'similarity_matrix'
    GeneTrail::DenseMatrix similarity_matrix;//can be a similarity or distance matrix depending on setting
    GeneTrail::DenseMatrix connectivity_matrix;

    IloEnv env;
    IloModel model;
    IloCplex cplex;

    IloNumVarArray y_ij;
    IloNumVarArray z_ij;
    IloNumArray sol_y_ij;
    IloNumArray sol_z_ij;

    bool sol_requested{false}; //needed to ensure that solutions are only generated once
    const std::string name_delim_cplex{"____"};

    auto getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int;

    auto doubleEpsilonEqual(double a, double b) -> bool;

public:

    ClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, double threshold, Metric sim_or_dist);
    ~ClusteringILP();
    //Forbidden until decided otherwise
    ClusteringILP(ClusteringILP & ilp) = delete; //copy ctor
    ClusteringILP(ClusteringILP && ilp) = delete; //move ctor
    ClusteringILP & operator=(ClusteringILP && ilp) = delete; //move operator

    auto buildModel() -> bool;
    auto solveModel() -> void;
    auto getSolution(GeneTrail::DenseMatrix & output_mtx) -> void;//output matrix must have the same dimensions as similarity matrix
};


#endif //ClusteringILP_ClusteringILP_H