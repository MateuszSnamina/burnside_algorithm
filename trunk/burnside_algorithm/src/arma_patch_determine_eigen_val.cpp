#include<armadillo>
#include<vector>
#include<string>
#include<stdexcept>
#include<arma_patch.hpp>

namespace armaPatch {

    const char* NotEigenVectorError::what() const noexcept {
        return "The tested vector is not an eigen vector.";
    };

    arma::cx_double determine_eigen_val(const arma::cx_mat & M, const arma::cx_vec & v, double threshold) {
        if (M.n_rows != M.n_cols) {
            std::string str = "The given matrix is not a square one. "
                    "(Note: The eigenproblem may be considered only for the square matrices.)";
            throw std::invalid_argument(str);
        }
        if (M.n_cols != v.n_rows) {
            std::string str = "The given tested vector and the given matrix sizes are incompatible. "
                    "(Note: The tested eigenvector and the matrix must have corresponding sizes).";
            throw std::invalid_argument(str);
        }
        const int dim = M.n_rows;
        const arma::cx_vec res = M * v;
        const arma::cx_vec lambdas = res / v;
        arma::cx_double lambda_numerator;
        int lambda_denominator = 0;
        for (int i = 0; i < dim; i++)
            if (std::abs(v(i)) > threshold) {
                lambda_numerator += lambdas(i);
                lambda_denominator++;
            }
        const arma::cx_double lambda = lambda_numerator / double(lambda_denominator);
        if (arma::norm(res - lambda * v) > threshold) {
            throw NotEigenVectorError();
        }
        return lambda;
    }

} // end of namespace armaPatch 