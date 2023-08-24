This comparison has not been automated; as it involves setting some preprocessor directives and recompiling BCM3. The steps used to make the comparison are described below.

1) Set the prediction loop in bcminf to skip over some samples to reduce the overall computation time
bcminf/main.cpp:
		//for (size_t i = total_samples / 2; i < total_samples; i++) {
		for (size_t i = total_samples / 2; i < total_samples; i += 50) {
        
        
2) To switch to CVODEs default dense solver
odecommon/LinearAlgebraSelector.h:
set CVODE_USE_EIGEN_SOLVER to 0


3) To switch to Eigen's default solver
odecommon/sunlinsol_dense_eigen.cpp:
//EIGSOL(S).compute_optimized(EIGMAT(A));
EIGSOL(S).compute(EIGMAT(A));


4) To disable code generation, 
cellpop/Cell.cpp, line 13, set Cell::use_generated_code to 0
(This is only possible when CVODE_USE_EIGEN_SOLVER is also set to 0, since without code generation we don't have a procedure to calculate the Jacobian and have to resort to CVODEs difference quotient Jacobian)


5) To disable cell parallelization, run with "-k 1"