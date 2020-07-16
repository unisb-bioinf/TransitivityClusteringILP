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
#include "IterativeClusteringILP.h"

IterativeClusteringILP::IterativeClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, size_t cplex_threads, size_t cplex_time_limit)
:
similarity_matrix(similarity_matrix), 
connectivity_matrix(similarity_matrix),
output_matrix(similarity_matrix),
cplex_threads(cplex_threads),
cplex_time_limit(cplex_time_limit),
model(env), 
cplex(env), 
obj_expr(env),
obj(env),
y_ij(env), 
z_ij(env),
sol_y_ij(env), 
sol_z_ij(env),
hasSolution(false),
constraints_yij(env),
constraints_zij(env)
{
    for(unsigned int i = 0; i< output_matrix.rows(); i++){
        for(unsigned int j = 0; j< output_matrix.cols(); j++){
            output_matrix(i,j) = 0.0;
        }
        output_matrix(i,i) = 1.0;
    }
}

auto IterativeClusteringILP::getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int{
    if(column < row) std::swap(row, column);
    if(row == column || row >= number_of_elements -1 || column >= number_of_elements){
        std::cout << row << " - " << column << " - " << number_of_elements << std::endl;
        throw std::out_of_range ("Index out of range: Only upper triangle matrix without diagonal is used");
    }

    return (number_of_elements -2)*row + column -1 - (row*(row-1))/2;
}

auto IterativeClusteringILP::initializeModel() -> bool{
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
    constraints_yij = IloExtractableArray(env, number_of_elements_whole_matrix);
    constraints_zij = IloExtractableArray(env, number_of_elements_whole_matrix);
    isVariableUsed = std::vector<bool>(number_of_elements_whole_matrix, false);
    isFixed = std::vector<bool>(number_of_elements_whole_matrix, false);

    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));

            std::string variable_namey{"y" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};
            std::string variable_namez{"z" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};

            y_ij[mat_index] = IloBoolVar(env, variable_namey.c_str());
            z_ij[mat_index] = IloBoolVar(env, variable_namez.c_str());
        }
    }

    std::cout << "INFO: Initializing model - Potential number of variables : " << 2 * number_of_elements_whole_matrix << std::endl;
    std::cout << "INFO: Initializing model - Generating constraints on decision variables ..." << std::endl;

    //Constraints for y_ij and z_ij
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));
            // Only add used variables
	    if(!isVariableUsed[mat_index]) continue;

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
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            // Only add used variables
            if(!isVariableUsed[mat_index_ij]) continue;
            obj_expr += y_ij[mat_index_ij] * similarity_matrix(i,j);
        }
    }
 
    obj = IloMinimize(env, obj_expr);
    model.add(obj);
    
    //
    cplex.extract(model);

    //Running on testosterone?
    cplex.setParam(IloCplex::Threads, cplex_threads); // max number of threads
    cplex.setParam(IloCplex::ParallelMode, 1); // deterministic calculations enforced
    if(cplex_time_limit > 0) cplex.setParam(IloCplex::Param::TimeLimit, 60 * cplex_time_limit); // Set time limit
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);

    return true;
}

auto IterativeClusteringILP::addVariableConstraints(size_t i, size_t j) -> bool {
    unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));
    if(isVariableUsed[mat_index]) {
        return true;
    }

    // Add constraints
    IloExpr expr(env);
    expr += y_ij[mat_index];
    IloConstraint cons = expr <= round(connectivity_matrix(i,j));
    model.add(cons);
    constraints_yij[mat_index] = cons;    

    IloExpr exprz(env);
    exprz += z_ij[mat_index] + y_ij[mat_index];
    IloConstraint consz = exprz == round(connectivity_matrix(i,j));
    model.add(consz);
    constraints_zij[mat_index] = consz;
    
    isVariableUsed[mat_index] = true;
    return true;
}

//Might actually not be needed anymore
/*auto ClusteringILP3::updateVariableConstraints(size_t i, size_t j) -> bool {
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
}*/

auto IterativeClusteringILP::updateCycleConstraints(size_t i, size_t j, size_t k) -> bool {
    unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
    unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
    unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

    // Make sure they are sorted i < j < k
    if(isVariableUsed[mat_index_ij] && isVariableUsed[mat_index_jk] && isVariableUsed[mat_index_ik]){
        // Remove the old constraint that is replaced
        // TODO: Make this pretty
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
        // TODO: Refactor this !!!
        bool keepExtending = true; 
        if(isVariableUsed[mat_index_ij] && isVariableUsed[mat_index_jk]){
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
        // No else, we need to check both
        if(isVariableUsed[mat_index_ij] && isVariableUsed[mat_index_ik]){
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
}

auto IterativeClusteringILP::updateCycleConstraints(size_t i, size_t j) -> bool {
    bool keepExtending = true;
    for(unsigned int k = 0; k< connectivity_matrix.rows(); ++k){
        if(k == i || k == j) continue;
        keepExtending = updateCycleConstraints(i, j ,k) ? keepExtending : false;
    }
    return keepExtending;
}

auto IterativeClusteringILP::updateObjectiveFunction(size_t i, size_t j) -> bool {
    unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
    // Update expression and objective function
    obj_expr += y_ij[mat_index_ij] * similarity_matrix(i,j);
    obj.setExpr(obj_expr);
    return true;
}

auto IterativeClusteringILP::fixPreviousSolution() -> void {
    for(unsigned int i= 0; i<connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j<connectivity_matrix.rows(); j++){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            if(!isVariableUsed[mat_index_ij]) continue;
            if(isFixed[mat_index_ij]) continue;
            isFixed[mat_index_ij] = output_matrix(i,j);
            if(isFixed[mat_index_ij]) {
                IloExpr expr(env);
                expr += z_ij[mat_index_ij];
                IloConstraint cons = expr == isFixed[mat_index_ij];
                model.add(cons);
            }
        }
    }
}

auto IterativeClusteringILP::extendModel(const std::vector<std::tuple<size_t, size_t, double>>& edges, size_t begin, size_t end, bool fixSolution) -> bool {
    bool keepExtending = true;
    if(fixSolution) {
	fixPreviousSolution();
    }
    for(size_t i = begin; i<end; ++i) {
        connectivity_matrix(std::get<0>(edges[i]), std::get<1>(edges[i])) = 1;
        connectivity_matrix(std::get<1>(edges[i]), std::get<0>(edges[i])) = 1;
        addVariableConstraints(std::get<1>(edges[i]), std::get<0>(edges[i]));
        //updateVariableConstraints(std::get<0>(edges[i]), std::get<1>(edges[i]));
        keepExtending = updateCycleConstraints(std::get<0>(edges[i]), std::get<1>(edges[i])) ? keepExtending : false;
        updateObjectiveFunction(std::get<0>(edges[i]), std::get<1>(edges[i]));
    }

    return keepExtending;
}

auto IterativeClusteringILP::solveModel() -> bool {
    
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

    return feasible;
}

auto IterativeClusteringILP::getSolution() -> GeneTrail::DenseMatrix
{
    for(unsigned int i= 0; i<connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j<connectivity_matrix.rows(); j++){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            if(!isVariableUsed[mat_index_ij]) continue;
            output_matrix(i,j) = connectivity_matrix(i,j) - cplex.getValue(y_ij[mat_index_ij]);
            output_matrix(j,i) = connectivity_matrix(j,i) - cplex.getValue(y_ij[mat_index_ij]);
        }
    }
    return output_matrix;
}

