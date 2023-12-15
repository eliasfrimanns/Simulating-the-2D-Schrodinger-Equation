#include <cmath>
#include <cstdio>
#include "generate_AB.hpp"
#include <string>

void normalize(arma::cx_vec& vector){
  arma::cx_vec vector_conj = arma::conj(vector);
  arma::vec     abs_vector = arma::abs(vector % vector_conj);

  double sum = arma::accu(abs_vector);
  arma::cx_double sum_sqrt = 1 / std::sqrt(sum);

  vector *= sum_sqrt;
}

void check_if_normalized(arma::cx_vec& vector, const int iteration){
  arma::cx_vec vector_conj = arma::conj(vector);
  arma::vec abs_vector = arma::abs(vector % vector_conj);

  double sum = arma::accu(abs_vector);

  double tolerance = 1e-11;

  if (sum < 1. - tolerance || sum > 1. + tolerance){
    printf("Problem occured at iteration %i, got a sum of %1.18f.\n", iteration, sum);
    throw std::runtime_error("The vector is not normalized");
  }
}

void initialize_n_slit_potential(arma::mat& V, const int number_of_slits, const double v_0, const double slit_width=0.05, const double wall_width=0.02, const double slit_center=0.5, const double wall_center=0.5){

  int M_2 = V.n_cols;
  int M   = M_2 + 2;

  int slit_center_index = M*(1-slit_center) - 1; // to account for the problem with indexing in generate_AB_matrices()
  int wall_center_index = wall_center * M - 1;

  int slit_width_index  = slit_width  * M;
  int wall_width_index  = wall_width  * M;

  int wall_start_index = wall_center_index - wall_width_index/2;
  int wall_end_index   = wall_start_index  + wall_width_index;

  int iteration_start = slit_center_index  -  slit_width_index * number_of_slits + slit_width_index/2;
  int iteration_end   = iteration_start + 2 * slit_width_index * number_of_slits;

  if (wall_width_index == 0) {
    std::cout << "The wall width is smaller than steplength, wall width is rounded up to 1 steplength." << "\n";
    wall_width_index = 1;
  }

  if (slit_width_index == 0) {
    std::cout << "The slit width is smaller than steplength, slit width is rounded up to 1 steplength." << "\n";
    slit_width_index = 1;
  }

  // add wall
  V.cols(wall_start_index, wall_end_index).fill(v_0);

  // add slits

  if (number_of_slits == 0) { return; } // stops the function here since there are no slits

  if (iteration_start < 0 || iteration_end > M_2) {
    throw std::runtime_error("Too many slits relative to the centering of the slits on the y-axis.");
  }

  V.cols(wall_center*M - wall_width_index/2, wall_center*M + wall_width_index/2).fill(v_0);

  for (int i = iteration_start; i < iteration_end; i += 2*slit_width_index){
    V.rows(i,i+slit_width_index-1).fill(0.);
  }
}

void initialize_quantum_state(arma::cx_mat& quantum_state, const double xc, const double yc, const double sigma_x, const double sigma_y, const double p_x, const double p_y){

  arma::cx_vec initial_quantum_state( arma::size( quantum_state.col(0) ) );

  int length_of_initial_quantum_state = initial_quantum_state.n_elem;
  int sidelength_of_matrix = std::sqrt(length_of_initial_quantum_state);

  arma::cx_double exponent;
  double x, y, x_subtract_xc, y_subtract_yc;
  int i,j,k;

  for (j=0;j<sidelength_of_matrix;j++){
    for (i=0;i<sidelength_of_matrix;i++){

      x = ((double)i) / sidelength_of_matrix;
      y = ((double)j) /sidelength_of_matrix;

      x_subtract_xc = (x - xc) / sigma_x;
      y_subtract_yc = (y - yc) / sigma_y;

      exponent = arma::cx_double( (-x_subtract_xc*x_subtract_xc - y_subtract_yc*y_subtract_yc) * 0.5 , x*p_x + y*p_y );

      k = ij_to_k(i,j,sidelength_of_matrix);
      initial_quantum_state(k) = std::exp(exponent);
    }
  }

  normalize(initial_quantum_state);
  quantum_state.col(0) = initial_quantum_state;
}
