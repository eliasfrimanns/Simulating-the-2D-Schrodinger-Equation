#include <armadillo>
#include <iomanip>
#include <tuple>

int ij_to_k(const int i, const int j, const int sidelength_of_matrix){

  if ((i >= sidelength_of_matrix) && (j>= sidelength_of_matrix))
  {
      throw std::runtime_error("i and j cannot be larger or equal to the sidelength of the matrix");

  } else if (i >= sidelength_of_matrix){
      throw std::runtime_error("i cannot be larger or equal to the sidelength of the matrix");

  } else if (j >= sidelength_of_matrix){
      throw std::runtime_error("j cannot be larger or equal to the sidelength of the matrix");
  }

  int k = i + j * sidelength_of_matrix;

  return k;
}

std::tuple<arma::sp_cx_mat,arma::sp_cx_mat> generate_AB_matrices(const arma::mat V, const int M, const double h, const double dt)
{
  int M_2 = M - 2;
  int V_sidelength = V.n_cols;

  if (V_sidelength != M_2){
    throw std::runtime_error("V must have sidelenths of size (M-2)");

  }else if (V_sidelength != V.n_rows){
    throw std::runtime_error("V must be a square matrix");
  }

  int matrix_sidelength = M_2 * M_2;

  arma::sp_cx_mat A(matrix_sidelength, matrix_sidelength);
  arma::sp_cx_mat B(matrix_sidelength, matrix_sidelength);

  arma::cx_double r(0, dt/(2*h*h));

  // the values
  arma::cx_vec r_diag1(matrix_sidelength-1,   arma::fill::value(r));
  arma::cx_vec r_diag2(matrix_sidelength-M_2, arma::fill::value(r));

  for (int idx = M_2-1; idx < matrix_sidelength; idx += M_2){
    r_diag1.at(idx) = 0;
  }
  arma::cx_double i_dt_half(0.,dt*0.5);

  int k;
  // the values along the main diagonal
  for (int j = 0; j < M_2; j++){
    for (int i = 0; i < M_2; i++){

      k = ij_to_k(i,j,M_2);

      A(k,k) = 1. + 4.*r + i_dt_half * V(j,i); // for some reason, the indexing needs to be reversed
      B(k,k) = 1. - 4.*r - i_dt_half * V(j,i); // for the potential to be implemented properly.
    }
  }
  // filling the diagonals of A and B with their
  A.diag(1)   = -r_diag1; A.diag(-1)   = -r_diag1;
  A.diag(M_2) = -r_diag2; A.diag(-M_2) = -r_diag2;

  B.diag(1)   = r_diag1; B.diag(-1)   = r_diag1;
  B.diag(M_2) = r_diag2; B.diag(-M_2) = r_diag2;

  return std::make_tuple(A,B);
}
