//
// Created by klenhof on 25.10.19.
//

#include <time.h>
#include "ClusteringILP.h"

ClusteringILP::ClusteringILP(GeneTrail::DenseMatrix & similarity_matrix, double threshold, Metric sim_or_dist): 
valid_input(true), 
threshold(threshold), 
sim_or_dist(sim_or_dist),                                                                                                             similarity_matrix(similarity_matrix), 
connectivity_matrix(similarity_matrix),                                                                                                                model(env), 
cplex(env), 
y_ij(env), 
z_ij(env),
sol_y_ij(env), 
sol_z_ij(env){

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

    int null = 0;
    for(unsigned int i = 0; i<connectivity_matrix.rows(); i++){
        int row = 0;
	for(unsigned int j = i+1; j<connectivity_matrix.cols(); j++){
		row += connectivity_matrix(i,j);
        }
	if(row == 0) {
		null += 1;
	}
    }

    std::cout << "Null: " << null << std::endl;

}

auto ClusteringILP::getIndexUpperTriangle(unsigned int row, unsigned int column, unsigned int number_of_elements) -> unsigned int{
    if(column <= row || row >= number_of_elements -1 || column >= number_of_elements){
        throw std::out_of_range ("Index out of range: Only upper triangle matrix without diagonal is used");
    }

    return (number_of_elements -2)*row + column -1 - (row*(row-1))/2;
}

auto ClusteringILP::doubleEpsilonEqual(double a, double b) -> bool{

    double epsilon(0.1);

    return abs(a-b) < epsilon;
}


auto ClusteringILP::buildModel() -> bool{

    if(!valid_input){
        std::cout << "Please check your input parameters for validity. Calculation not possible." << std::endl;
        return false;
    }

    std::cout << "INFO: Building model - Defining variables ..." << std::endl;
    //Defining y_ij and z_ij using only a 1D array structure

    unsigned int number_of_elements_whole_matrix(getIndexUpperTriangle(connectivity_matrix.rows()-2, connectivity_matrix.rows()-1, connectivity_matrix.rows())+1);

    y_ij = IloNumVarArray(env,number_of_elements_whole_matrix);
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

    std::cout << "INFO: Building model - Number of variables : " << 2 * number_of_elements_whole_matrix << std::endl;

    std::cout << "INFO: Building model - Generating constraints on decision variables ..." << std::endl;
    //Constraints for y_ij and z_ij
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){

        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            unsigned int mat_index(getIndexUpperTriangle(i,j,connectivity_matrix.rows()));

            IloExpr expr(env);
            expr += y_ij[mat_index];
            model.add(y_ij[mat_index] <= round(connectivity_matrix(i,j)));
            expr.end();


            IloExpr exprz(env);
            exprz += z_ij[mat_index] + y_ij[mat_index];
            model.add(exprz == round(connectivity_matrix(i,j)));
            exprz.end();
        }
    }


    std::cout << "INFO: Building model - Generation triangle constraints ..." << std::endl;  
    //Further constraints for z_ij
    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){

        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            for(unsigned int k = j+1; k< connectivity_matrix.rows(); k++){


                if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(j,k), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){//is that correct?

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
                else if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(j,k), 1.0)){
                    unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
                    unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));

                    //Only one expression needed that ensures that only one edge is taken

                    IloExpr expr_one(env);
                    expr_one = z_ij[mat_index_ij] + z_ij[mat_index_jk];
                    model.add(expr_one <=1);
                    expr_one.end();
                }
                else if(doubleEpsilonEqual(connectivity_matrix(i,j), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){
                    unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
                    unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

                    IloExpr expr_one2(env);
                    expr_one2 = z_ij[mat_index_ij] + z_ij[mat_index_ik];
                    model.add(expr_one2 <=1);
                    expr_one2.end();
                }
                else if(doubleEpsilonEqual(connectivity_matrix(j,k), 1.0) && doubleEpsilonEqual(connectivity_matrix(i,k), 1.0)){
                    unsigned int mat_index_jk(getIndexUpperTriangle(j,k, connectivity_matrix.rows()));
                    unsigned int mat_index_ik(getIndexUpperTriangle(i,k, connectivity_matrix.rows()));

                    IloExpr expr_one3(env);
                    expr_one3 = z_ij[mat_index_jk] + z_ij[mat_index_ik];
                    model.add(expr_one3 <=1);
                    expr_one3.end();
                }
            }
        }
    }

    
    std::cout << "INFO: Building model - Setting objective function ..." << std::endl;

    //Objective function
    IloExpr obj(env);

    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j< connectivity_matrix.rows(); j++){

            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));

            if(sim_or_dist == similarity){
                obj += y_ij[mat_index_ij] * similarity_matrix(i,j);
            }
            else{
                obj += y_ij[mat_index_ij] * (1.0/(1.0 + similarity_matrix(i,j)));
            }
        }
    }

    model.add(IloMinimize(env, obj));

    return true;
}

auto ClusteringILP::solveModel() -> void {

    cplex.extract(model);

    //Running on testosterone?
    cplex.setParam(IloCplex::Threads, 32); //max number of threads
    cplex.setParam(IloCplex::ParallelMode, 1); //deterministic calculations enforced
    cplex.setParam(IloCplex::Param::TimeLimit, 60 * 10);
   
    bool feasible{false};

    try{
        feasible = cplex.solve();
	if(cplex.getCplexStatus() == IloCplex::CplexStatus::AbortTimeLim) {
		feasible = false;
		std::cout << "Computation aborted due to time limit." << std::endl;
	}
	if(feasible) {
        	std::cout << "ILP successfully solved" << std::endl;
	}

    }
    catch(IloException & e){
        std::cout << "CPLEX raised an exception" << std::endl;
        std::cout << e << std::endl;
        e.end();

    }
    catch(std::exception & e){
        std::cout << "Exception occurred" << std::endl;
        std::cout << e.what() << std::endl;
    }

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
}

auto ClusteringILP::getSolution() -> GeneTrail::DenseMatrix
{
    GeneTrail::DenseMatrix output_mtx(connectivity_matrix.rows(), connectivity_matrix.rows());

    for(unsigned int i = 0; i< connectivity_matrix.rows(); i++){
        for(unsigned int j = 0; j< connectivity_matrix.cols(); j++){
            output_mtx(i,j) = 0.0;
        }
        output_mtx(i,i) = 1.0;
    }

    for(unsigned int i= 0; i<connectivity_matrix.rows(); i++){
        for(unsigned int j = i+1; j<connectivity_matrix.rows(); j++){
            unsigned int mat_index_ij(getIndexUpperTriangle(i,j, connectivity_matrix.rows()));
            output_mtx(i,j) = connectivity_matrix(i,j) - cplex.getValue(y_ij[mat_index_ij]);
            output_mtx(j,i) = connectivity_matrix(j,i) - cplex.getValue(y_ij[mat_index_ij]);
        }
    }
    return std::move(output_mtx);
}

