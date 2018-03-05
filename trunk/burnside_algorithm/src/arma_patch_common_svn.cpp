#include<armadillo>
#include<vector>
#include<arma_patch.hpp>

namespace armaPatch {

    Decomposition common_svn(std::vector<arma::cx_mat> Ms) {
        std::vector<Decomposition> decompositions;
        for (arma::cx_mat M : Ms) {
            Decomposition M_dec = Decomposition::decomposition_from_rightEigVecs(M); // to tworzy rozklad zgodny z prawymi wektorami wlasnymi. 	
            decompositions.push_back(M_dec);
        }
        return Decomposition::common_decomposition(decompositions);
    }


} // end of namespace armaPatch 