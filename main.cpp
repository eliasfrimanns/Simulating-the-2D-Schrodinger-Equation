#define ARMA_USE_SUPERLU
#include "initialize.hpp"     // also includes: armadillo, cmath, cstdio, generate_AB.hpp, iomanip, iostream, string, tuple


int main(int argc, char* argv[]){

  if (argc != 13){
    throw std::runtime_error("Wrong number of input parameters. Check README for proper use.");
  }

  double h            = atof(argv[1 ]);
  double dt           = atof(argv[2 ]);
  double T            = atof(argv[3 ]);
  double xc           = atof(argv[4 ]);
  double sigma_x      = atof(argv[5 ]);
  double p_x          = atof(argv[6 ]);
  double yc           = atof(argv[7 ]);
  double sigma_y      = atof(argv[8 ]);
  double p_y          = atof(argv[9 ]);
  double v_0          = atof(argv[10]);
  int number_of_slits = atoi(argv[11]);


  int M     = 1 / h;
  int M_2   = M - 2;
  int t_max = T / dt;


  arma::mat V(M_2,M_2);

  arma::sp_cx_mat A, B;

  arma::cx_mat quantum_state(M_2*M_2,t_max);

  arma::cx_vec quantum_state_iteration; // this is the state of the previous/current iteration
  arma::cx_vec b; // the b-vector in A*x_{i+1} = b = B*x_i

  std::string custom_slits = argv[12];
  std::string yes = "y"; std::string no = "n";

  double slit_width, slit_center, wall_width, wall_center;

  if (custom_slits == yes)
  {
    std::cout << "Enter slit width, slit center, wall width and wall center: ";
    std::cin  >> slit_width >> slit_center >> wall_width >> wall_center;
    std::cout << std::endl;

    if (number_of_slits >= 0){
      initialize_n_slit_potential(V, number_of_slits, v_0, slit_width, wall_width, slit_center, wall_center);
    }

  } else if (custom_slits == no){

      if (number_of_slits >= 0){
        initialize_n_slit_potential(V, number_of_slits, v_0);
      }
  } else {
    throw std::runtime_error("You must type y or n in the last variable space. Check README for proper use.");
  }


  std::tie(A,B) = generate_AB_matrices(V,M,h,dt);

  initialize_quantum_state(quantum_state, xc, yc, sigma_x, sigma_y, p_x, p_y);

  quantum_state_iteration = quantum_state.col(0);

  arma::spsolve_factoriser SF; // initializing the LU-factorizer that will be used to solve the wave-equation
  bool factorization_status = SF.factorise(A); // LU-factorizes A and reuses that one factorization, also checks the state of the factorization

  if(factorization_status == false) { throw std::runtime_error("LU-factorization of A failed."); }

  for (int iteration = 1; iteration < t_max; iteration++){

    b = B*quantum_state_iteration;
    SF.solve(quantum_state_iteration,b); // solves Ax = b, where x = quantum_state_iteration is unknown

    check_if_normalized(quantum_state_iteration, iteration); // checks if x is normalized within a tolerance, from initialize.hpp
    quantum_state.col(iteration) = quantum_state_iteration;
  }
  std::string filename;
  if (number_of_slits >= 0){
    filename = "quantum_state_vec_" + std::to_string(number_of_slits) + "_slit(s).bin";
  } else {
    filename = "quantum_state_vec_no_potential.bin";
  }
  quantum_state.save(filename);
  return 0;
}
