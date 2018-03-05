#include<armadillo>
#include<vector>
#include<arma_patch.hpp>

namespace armaPatch {

    Decomposition common_eig_gen(std::vector<arma::cx_mat> Ms) {
        std::vector<Decomposition> decompositions;
        for (arma::cx_mat M : Ms) {
            const Decomposition decomposition = eig_gen(M);
            decompositions.push_back(decomposition);
        }
        return Decomposition::common_decomposition(decompositions);
    }

} // end of namespace armaPatch 