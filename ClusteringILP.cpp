//
// Created by klenhof on 25.10.19.
//


#include "ClusteringILP.h"

ClusteringILP::ClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, double threshold, Metric sim_or_dist): valid_input(true), threshold(threshold), sim_or_dist(sim_or_dist),
                                                                                                                similarity_matrix(similarity_matrix), connectivity_matrix(similarity_matrix),
                                                                                                                model(env), cplex(env), y_ij(env), z_ij(env),
                                                                                                                sol_y_ij(env), sol_z_ij(env){

    if(similarity_matrix.rows() != similarity_matrix.cols()){
        valid_input = false;
        std::cout << "Please use a quadratic, symmetric similarity/distance matrix as input" << std::endl;
    }

    if(sim_or_dist == similarity){

        for(unsigned int i = 0; i<similarity_matrix.rows(); i++){

            for(unsigned int j = 0; j<similarity_matrix.cols(); j++){

                if(i == j){
                    connectivity_matrix(i,j) = 0.0; //the diagonal does not have to be included
                    continue;
                }
                if(similarity_matrix(i,j)>= threshold){
                    connectivity_matrix(i,j) = 1.0;
                }
                else{
                    connectivity_matrix(i,j) = 0.0;
                }
            }
        }
    }
    else if(sim_or_dist == distance){

        for(unsigned int i = 0; i<similarity_matrix.rows(); i++){

            for(unsigned int j = 0; j<similarity_matrix.cols(); j++){

                if(i==j){
                    connectivity_matrix(i,j) = 0.0; //diagonal not included
                    continue;
                }
                if(similarity_matrix(i,j)<= threshold){
                    connectivity_matrix(i,j) = 1.0;
                }
                else{
                    connectivity_matrix(i,j) = 0.0;
                }
            }
        }

    }
    else{
        valid_input = false;
        std::cout << "Please specify whether you use a similarity or distance measure for clustering. " << std::endl;
    }

}

auto ClusteringILP::getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int{

    if(row == column || column> row|| row>= number_of_elements -1 || column>=number_of_elements){
        throw std::out_of_range ("Index out of range: Only upper triangle matrix without diagonal is used");
    }

    return (number_of_elements -2)*row + column -1 - (row*(row-1))/2;
}

auto ClusteringILP::doubleEpsilonEqual(double a, double b){

    double epsilon(0.1);

    return abs(a-b) < epsilon;
}


auto ClusteringILP::buildModel() -> void{

    if(!valid_input){
        std::cout << "Please check your input parameters for validity. Calculation not possible." << std::endl;
        return;
    }

    //Defining y_ij and z_ij using only a 1D array structure

    unsigned int number_of_elements(getIndexUpperTriangle(connectivity_matrix.rows()-2, connectivity_matrix.rows()-1)+1);

    y_ij = IloNumVarArray(env,number_of_elements);
    z_ij = IloNumVarArray(env, number_of_elements);

    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){

        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));

            std::string variable_namey{"y" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};
            std::string variable_namenz{"z" + name_delim_cplex + connectivity_matrix.rowName(i) + name_delim_cplex + connectivity_matrix.colName(j)};

            y_ij[mat_index] = IloBoolVar(env, variable_namey.c_str());
            z_ij[mat_index] = IloBoolVar(env, variable_namez.c_str());
        }
    }


    //Constraints for y_ij and z_ij
    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){

        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));

            IloExpr expr(env);
            expr += y_ij[mat_index];
            model.add(y_ij <= round(connectivity_matrix(i,j)));
            expr.end();


            IloExpr exprz(env);
            exprz += z_ij[mat_index] + y_ij[mat_index];
            model.add(exprz == round(connectivity_matrix(i,j)));
            exprz.end();
        }
    }


    //Further constraints for z_ij

    for(unsigned int i = 0; i< connectivity_matrix.rows()-1; i++){

        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            for(unsigned int k = j+1; k< connectivity_matrix.rows(); k++){


                if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(j,k), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k))){//is that correct?

                    unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
                    unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
                    unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

                    //All three constraints needed because there exists a triangle
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
                }
                else{
                    un
                }
                IloExpr expr_
            }
        }
    }
}

auto ClusteringILP::solveModel() {


}

auto ClusteringILP::getSolution(GeneTrail::DenseMatrix &output_mtx) {


}

