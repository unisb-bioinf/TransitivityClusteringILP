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

#include<string>
#include <time.h>
#include "ClusteringILP2.h"

ClusteringILP2::ClusteringILP2(GeneTrail::DenseMatrix & similarity_matrix, double threshold)
: threshold(threshold),
similarity_matrix(similarity_matrix), 
connectivity_matrix(similarity_matrix),
model(env), 
cplex(env), 
y_ij(env), 
z_ij(env),
sol_y_ij(env), 
sol_z_ij(env),
constraints_yij(env),
constraints_zij(env)
{}

auto ClusteringILP2::getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int{
    if(column < row) std::swap(row, column);
    if(row == column || row >= number_of_elements -1 || column >= number_of_elements){
        std::cout << row << " - " << column << " - " << number_of_elements << std::endl;
        throw std::out_of_range ("Index out of range: Only upper triangle matrix without diagonal is used");
    }

    return (number_of_elements -2)*row + column -1 - (row*(row-1))/2;
}

auto ClusteringILP2::doubleEpsilonEqual(double a, double b) -> bool{

    double epsilon(0.1);

    return abs(a-b) < epsilon;
}

auto ClusteringILP2::initializeModel() -> bool{
    std::cout << "INFO: Initializing model - Defining variables ..." << std::endl;
    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){
        for(unsigned int j = 0; j< connectivity_matrix.rows(); j++){
            connectivity_matrix(i,j) = 0.0;
        }
    }

    //Defining y_ij and z_ij using only a 1D array structure

    unsigned int number_of_elements_whole_matrix(getIndexUpperTriangle(connectivity_matrix.rows()-2, connectivity_matrix.rows()-1, connectivity_matrix.rows())+1);

    y_ij = IloNumVarArray(env, number_of_elements_whole_matrix);
    z_ij = IloNumVarArray(env, number_of_elements_whole_matrix);

    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));

            std::string variable_namey{"y" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};
            std::string variable_namez{"z" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};

            y_ij[mat_index] = IloBoolVar(env, variable_namey.c_str());
            z_ij[mat_index] = IloBoolVar(env, variable_namez.c_str());
        }
    }

    std::cout << "INFO: Initializing model - Number of variables : " << 2 * number_of_elements_whole_matrix << std::endl;
    std::cout << "INFO: Initializing model - Generating constraints on decision variables ..." << std::endl;

    //Constraints for y_ij and z_ij
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));
            IloExpr expr(env);
            expr += y_ij[mat_index];
            IloConstraint cons = expr <= round(connectivity_matrix(i,j));
            model.add(cons);
            constraints_yij.add(cons);

            IloExpr exprz(env);
            exprz += z_ij[mat_index] + y_ij[mat_index];
            IloConstraint consz = exprz == round(connectivity_matrix(i,j));
            model.add(consz);
            constraints_zij.add(consz);
        }
    }
    
    std::cout << "INFO: Initializing model - Defining objective function ..." << std::endl;

    //Objective function
    IloExpr obj(env);

    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            obj += y_ij[mat_index_ij] * similarity_matrix(i,j);
        }
    }

    model.add(IloMinimize(env, obj));
    
    //
    cplex.extract(model);

    //Running on testosterone?
    cplex.setParam(IloCplex::Threads, 32); // max number of threads
    cplex.setParam(IloCplex::ParallelMode, 1); // deterministic calculations enforced
    cplex.setParam(IloCplex::Param::TimeLimit, 60 * 10); // Set time limit
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);

    return true;
}

auto ClusteringILP2::updateVariableConstraints(size_t i, size_t j) -> bool {
    unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));
    model.remove(constraints_yij[mat_index]);
    IloExpr expr(env);
    expr += y_ij[mat_index];
    IloConstraint cons = expr <= round(connectivity_matrix(i,j));
    model.add(cons);
    constraints_yij[mat_index] = cons;
    
    model.remove(constraints_zij[mat_index]);
    IloExpr exprz(env);
    exprz += z_ij[mat_index] + y_ij[mat_index];
    IloConstraint consz = exprz == round(connectivity_matrix(i,j));
    model.add(consz);
    constraints_zij[mat_index] = consz;

    return true;
}

auto ClusteringILP2::updateCycleConstraints(size_t i, size_t j, size_t k) -> bool {
    // Make sure they are sorted i < j < k
    if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(j,k), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){//is that correct?
        unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
        unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
        unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

        // Remove the old constraint that is replaced
	std::string key;
        if(mat_index_ik < mat_index_jk) {
            key = std::to_string(mat_index_ik) + "-" + std::to_string(mat_index_jk);
        } else {
            key = std::to_string(mat_index_jk) + "-" + std::to_string(mat_index_ik);
        }
        IloConstraint to_del = edge_constraints[key];
        model.remove(to_del);

        // Add the new constraints
        // All three constraints needed because there exists a triangle
        IloExpr expr_three1(env);
        expr_three1 = z_ij[mat_index_ij] + z_ij[mat_index_jk] - z_ij[mat_index_ik];
        model.add(expr_three1 <= 1);
        expr_three1.end();

        IloExpr expr_three2(env);
        expr_three2 = z_ij[mat_index_ij] + z_ij[mat_index_ik] - z_ij[mat_index_jk];
        model.add(expr_three2 <=1);
        expr_three2.end();

        IloExpr expr_three3(env);
        expr_three3 = z_ij[mat_index_jk] + z_ij[mat_index_ik] - z_ij[mat_index_ij];
        model.add(expr_three3 <=1);
        expr_three3.end();

	return false;
    }
    else  {
	bool keepExtending = true; 
        if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(j,k), 1.0)){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
    
            IloExpr expr_one(env);
            expr_one = z_ij[mat_index_ij] + z_ij[mat_index_jk];
            IloConstraint cons = expr_one <= 1;
            model.add(cons);
            std::string key;
            if(mat_index_ij < mat_index_jk) {
                key = std::to_string(mat_index_ij) + "-" + std::to_string(mat_index_jk);
            } else {
                key = std::to_string(mat_index_jk) + "-" + std::to_string(mat_index_ij);
            }
            edge_constraints[key] = cons;
            keepExtending = false;
        }
        if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));
    
            IloExpr expr_one2(env);
            expr_one2 = z_ij[mat_index_ij] + z_ij[mat_index_ik];
            IloConstraint cons = expr_one2 <= 1;
            model.add(cons);
            std::string key;
            if(mat_index_ij < mat_index_ik) {
                key = std::to_string(mat_index_ij) + "-" + std::to_string(mat_index_ik);
            } else {
                key = std::to_string(mat_index_ik) + "-" + std::to_string(mat_index_ij);
            }   
            edge_constraints[key] = cons;
            keepExtending = false;
        }
	return keepExtending;
    }
    return true;
    // We don't need this anymore
        /*else if(doubleEpsilonEqual(connectivity_matrix(j,k), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){
            unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
            unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

            IloExpr expr_one3(env);
            expr_one3 = z_ij[mat_index_jk] + z_ij[mat_index_ik];
            model.add(expr_one3 <=1);
            expr_one3.end();
        }*/
}

auto ClusteringILP2::updateCycleConstraints(size_t i, size_t j) -> bool {
    bool keepExtending = true;
    for(unsigned int k = 0; k< connectivity_matrix.rows(); ++k){
        if(k == i || k == j) continue;
        keepExtending = updateCycleConstraints(i, j ,k) ? keepExtending : false;
    }
    return keepExtending;
}

auto ClusteringILP2::extendModel(const std::vector<std::tuple<size_t, size_t, double>>& edges, size_t begin, size_t end) -> bool {
    bool keepExtending = true;
    for(size_t i = begin; i<end; ++i) {
        connectivity_matrix(std::get<0>(edges[i]), std::get<1>(edges[i])) = 1;
        connectivity_matrix(std::get<1>(edges[i]), std::get<0>(edges[i])) = 1;
        updateVariableConstraints(std::get<0>(edges[i]), std::get<1>(edges[i]));
        keepExtending = updateCycleConstraints(std::get<0>(edges[i]), std::get<1>(edges[i])) ? keepExtending : false;
    }

    return keepExtending;
}

auto ClusteringILP2::solveModel(size_t i) -> bool {
    
    bool feasible{false};

    try{
        feasible = cplex.solve();
    if(cplex.getCplexStatus() == IloCplex::CplexStatus::AbortTimeLim) {
        feasible = false;
        std::cout << "INFO: Computation aborted due to time limit." << std::endl;
    }
    if(feasible) {
            std::cout << "INFO: ILP successfully solved." << std::endl;
    }

    }
    catch(IloException & e){
        std::cerr << "ERROR: CPLEX raised an exception!" << std::endl;
        std::cerr << e << std::endl;
        e.end();

    }
    catch(std::exception & e){
        std::cerr << "ERROR: Exception occurred!" << std::endl;
        std::cerr << e.what() << std::endl;
    }

    std::string s = "out_" + std::to_string(i) + ".lp";
    cplex.exportModel(s.c_str());

    /*if(!feasible){
        std::cout << "ILP was not feasible, additional information is written to Not_feasible.txt" << std::endl;
        std::ofstream output("Not_feasible.txt");
        output << "ILP was not feasible" << std::endl;
        output << "Status\t" << cplex.getStatus() << std::endl;
        output << "Solver Status\t" << cplex.getCplexStatus() << std::endl;
        output.close();
    }
    else{
        std::cout << "ILP was feasible" << std::endl;
        std::ofstream sol_info("Feasible.txt");
        sol_info << "ILP was feasible" << std::endl;
        sol_info << "Status\t" << cplex.getStatus() << std::endl;
        sol_info << "Solver Status\t" << cplex.getCplexStatus() << std::endl;
        sol_info.close();
        //Solution as is by cplex
        cplex.writeSolution("Solution.sol");

    }*/

    return feasible;
}

auto ClusteringILP2::getSolution() -> GeneTrail::DenseMatrix
{
    GeneTrail::DenseMatrix output_mtx(connectivity_matrix.rows(), connectivity_matrix.rows());
    cplex.getValues(sol_y_ij, y_ij);
    cplex.getValues(sol_z_ij, z_ij);

    //Put zeroes everywhere, diagonal cannot be touched in next iteration
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = 0; j< connectivity_matrix.cols(); j++){
            output_mtx(i,j) = 0.0;
        }
	output_mtx(i,i) = 1.0;
    }

    for(unsigned int i= 0; i<connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j<connectivity_matrix.rows(); j++){

            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));

            output_mtx(i,j) = connectivity_matrix(i,j) - sol_y_ij[mat_index_ij];
            output_mtx(j,i) = connectivity_matrix(j,i) - sol_y_ij[mat_index_ij];
        }
    }
    return std::move(output_mtx);
}

