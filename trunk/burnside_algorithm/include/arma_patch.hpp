#ifndef ARMA_PATCH_HPP
#define ARMA_PATCH_HPP

#include<armadillo>
#include<vector>
#include<exception>

namespace armaPatch {

    // --------- Klasa reprezentujaca: rozklad jednosci --------------------------------------------
    // ------------------------------  [ lub matematycznie rownowaznie ] ---------------------------
    // ------------------------------  rozkald przestrzeni wektorowej na sume prosta ---------------

    class Decomposition {
    public:
        Decomposition();
        Decomposition(unsigned dim);
        static Decomposition common_decomposition(Decomposition decomposition1, Decomposition decomposition2);
        static Decomposition common_decomposition(const std::vector<Decomposition> & decompositions);
        static Decomposition decomposition_from_eigval_and_eigvec(const arma::vec & eig_vals, const arma::cx_mat & eig_vecs, double threshold = 1e-5);
        static Decomposition decomposition_from_rightEigVecs(const arma::cx_mat & M, double threshold = 1e-5);
        static Decomposition decomposition_from_eigVecs_of_symMat(const arma::mat & M, const char* method = "std", double threshold = 1e-5);
        std::vector<arma::cx_vec> get_basis() const;
        void print() const;

        unsigned get_dim() const {
            return _m_dim;
        };
    private:
        std::vector<arma::cx_mat> _m_projections;
        unsigned _m_dim;
    };

    // -------------------------- sprawdzanie wartosci wlasnej ---------------------------------
    // Gdy v to wektor wlasny M to zwracana jest wartosc wlasna,
    // Gdy v to nie wektor wlasny M to rzucany jest wyjatek notEigenVectorError.
    arma::cx_double determine_eigen_val(const arma::cx_mat & M, const arma::cx_vec & v, double threshold = 1e-5);

    struct NotEigenVectorError : std::exception {
        const char* what() const noexcept;
    };

    // -------------------------- robust eig_gen ------------------------------------------------
    // Funkcja buduje rozklad i nie sprawdza jego poprawnosci:
    Decomposition _eig_gen(const arma::cx_mat & M);
    // Funkcja buduje rozklad i sprawdza jego poprawnosc:
    Decomposition eig_gen(const arma::cx_mat & M);
    // Funkcja buduje rozklad i sprawdza jego poprawnosc, API armadillo:
    void eig_gen(arma::cx_vec & eigval, arma::cx_mat & eigvec, const arma::cx_mat & M);

    // -------------------------- common_eig_gen ------------------------------------------------
    // Funkcja robioca jednoczesna diagonalizacje macierzy symetrycznych
    Decomposition common_eig_gen(std::vector<arma::cx_mat> Ms);

    // -------------------------- common_eig_SVN ------------------------------------------------
    // Funkcja znajdujaca wspolna baze prawych wektorow singulanrych 
    Decomposition common_svn(std::vector<arma::cx_mat> Ms);

} // end of namespace armaPatch 

#endif