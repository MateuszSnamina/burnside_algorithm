#include<armadillo>
#include<vector>
#include<arma_patch.hpp>

namespace armaPatch {

    Decomposition _eig_gen(const arma::cx_mat & M) {
        const int dim = M.n_rows;
        const arma::cx_mat unity_matrix = arma::eye<arma::cx_mat>(dim, dim);
        std::vector<arma::cx_mat> Ms = {M, M + 1 * unity_matrix, M + arma::cx_double(0, 1) * unity_matrix};
        Decomposition decomposition = common_svn(Ms);
        //std::vector<arma::cx_vec> basis = decomposition.get_basis();
        return decomposition;
    }

    Decomposition eig_gen(const arma::cx_mat & M) {
        Decomposition decomposition = _eig_gen(M);
        const std::vector<arma::cx_vec> basis = decomposition.get_basis();
        // If basis_vec is not an eigenvector than
        // determine_eigen_val function throws an exception.
        for (const arma::cx_vec & basis_vec : basis)
            determine_eigen_val(M, basis_vec);
        return decomposition;
    }

    void eig_gen(arma::cx_vec & eigval, arma::cx_mat & eigvec, const arma::cx_mat & M) {
        Decomposition decomposition = _eig_gen(M);
        std::vector<arma::cx_vec> basis = decomposition.get_basis();
        const unsigned dim = basis.size();
        eigvec.resize(dim, dim);
        eigval.resize(dim);
        for (unsigned i = 0; i < dim; i++) {
            eigvec.col(i) = basis[i] / basis[i](0);
            eigval(i) = determine_eigen_val(M, basis[i]);
        }
    }

} // end of namespace armaPatch 