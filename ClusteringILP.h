/** 
 * Copyright (C) 2020 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *                    Kerstin Lenhof <klenhof@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "DenseMatrixReader.h"
#include "DenseMatrix.h"
#include "macros.h"
#define IL_STD
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#ifndef ClusteringILP_ClusteringILP_H
#define ClusteringILP_ClusteringILP_H


class GT2_EXPORT ClusteringILP{

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

    // Variables
    IloNumVarArray y_ij;
    IloNumVarArray z_ij;	
    // Solution
    IloNumArray sol_y_ij;
    IloNumArray sol_z_ij;

    bool sol_requested{false}; //needed to ensure that solutions are only generated once
    const std::string name_delim_cplex{"____"};

    auto getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int;

    auto doubleEpsilonEqual(double a, double b) -> bool;

public:

    ClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, double threshold, Metric sim_or_dist);
    ~ClusteringILP(){}
    //Forbidden until decided otherwise
    ClusteringILP(ClusteringILP & ilp) = delete; //copy ctor
    ClusteringILP(ClusteringILP && ilp) = delete; //move ctor
    ClusteringILP & operator=(ClusteringILP && ilp) = delete; //move operator

    auto buildModel() -> bool;
    auto solveModel() -> void;
    auto getSolution() -> GeneTrail::DenseMatrix;//output matrix must have the same dimensions as similarity matrix
};


#endif //ClusteringILP_ClusteringILP_H
