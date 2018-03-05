#include<armadillo>
#include<cmath>
#include<vector>
#include<list>
#include<string>
#include<iostream>
#include<stdexcept>

#include<arma_patch.hpp>

/*
 *  Konstruktor budujacy rozklad identycznosci dla zerowymiarowej przestrzeni.
 */
armaPatch::Decomposition::Decomposition() : _m_dim(0) {
};

/*
 *  Konstruktor budujacy trywialny rozklad identycznosci (tj. id=id).
 */
armaPatch::Decomposition::Decomposition(unsigned dim) : _m_dim(dim) {
    _m_projections.push_back(arma::diagmat(arma::cx_vec(dim, arma::fill::ones)));
}

/*
 * Funkcja budujaca rozklad identycznosci prosta bedacy ,,przecieciem'' dwoch zadanych rozkladow
 */
armaPatch::Decomposition armaPatch::Decomposition::common_decomposition(Decomposition decomposition1, Decomposition decomposition2) {
    if (decomposition1._m_dim != decomposition2._m_dim) {
        std::string str = "The two given decomposition are not defined in spaces of the same dimension. "
                "(Note: The two given decompositions to combine must have the same dimension.)";
                throw std::invalid_argument(str);
    }
    armaPatch::Decomposition finalDecomposition;
    finalDecomposition._m_dim = decomposition1.get_dim();
    for (arma::cx_mat projection1 : decomposition1._m_projections)
        for (arma::cx_mat projection2 : decomposition2._m_projections) {
            arma::cx_mat projection1projection2 = projection1 * projection2;
            if (std::abs(trace(projection1projection2)) > 1e-6)
                finalDecomposition._m_projections.push_back(projection1projection2); // nowy operator dodajemy gdy jest on czyms innym operatorem zerowym.
        }
    return finalDecomposition;
}

/*
 * Funkcja budujaca rozklad identycznosci bedacy ,,przecieciem'' wielu zadanych rozkladow
 */
armaPatch::Decomposition armaPatch::Decomposition::common_decomposition(const std::vector<Decomposition> & decompositions) {
    armaPatch::Decomposition finalDecomposition(decompositions[0].get_dim());
    for (const Decomposition & decomposition : decompositions)
        finalDecomposition = Decomposition::common_decomposition(finalDecomposition, decomposition);
    return finalDecomposition;
}

/*
 * Funkcja budujaca rozklad identycznosci 
 * rownowazy
 * rozkladowi na sume prosta podprzestrzeni
 * zdefiniowanych przez wektory z eigVecs o tych samych wartosciach eigVals.
 */
armaPatch::Decomposition armaPatch::Decomposition::decomposition_from_eigval_and_eigvec(const arma::vec & eig_vals, const arma::cx_mat & eig_vecs, double threshold) {
    armaPatch::Decomposition decomposition;
    decomposition._m_dim = eig_vecs.n_rows;
    std::list<unsigned> elements;
    for (unsigned i = 0; i < eig_vecs.n_cols; i++) elements.push_back(i);
    while (!elements.empty()) {
        arma::mat projectionDiag = arma::mat(decomposition._m_dim, decomposition._m_dim, arma::fill::zeros);
        double eigVal = eig_vals(*(elements.begin()));
        for (std::list<unsigned>::iterator it = elements.begin(); it != elements.end();)
            if (std::abs(eig_vals(*it) - eigVal) < 0.0001) {
                projectionDiag(*it, *it) = 1;
                it = elements.erase(it);
            } else
                it++;
        arma::cx_mat projection = eig_vecs * projectionDiag * inv(eig_vecs);
        decomposition._m_projections.push_back(projection);
    }
    return decomposition;
}

/*
 *  Funkcja budujaca rozklad identycznosci zwiazana z prawymi wektorami singularnymi:
 */
armaPatch::Decomposition armaPatch::Decomposition::decomposition_from_rightEigVecs(const arma::cx_mat & M, double threshold) {
    if (M.n_rows != M.n_cols) {
        std::string str = "M powinna byc macierza kwadratowa.";
        throw std::invalid_argument(str);
    }
    arma::cx_mat U, V;
    arma::vec sing_vals;
    arma::svd(U, sing_vals, V, M);
    return decomposition_from_eigval_and_eigvec(sing_vals, V, threshold);
}

/*
 *  Funkcja budujaca rozklad identycznosci zwiazana z wektorami wlasnymi macierzy symetrycznej:
 */
armaPatch::Decomposition armaPatch::Decomposition::decomposition_from_eigVecs_of_symMat(const arma::mat & M, const char* method, double threshold) {
    arma::vec eigVals;
    arma::mat eigVecs;
    arma::eig_sym(eigVals, eigVecs, M, method);
    arma::cx_mat cx_eigVecs(eigVecs, arma::mat(eigVecs.n_rows, eigVecs.n_cols, arma::fill::zeros));
    return decomposition_from_eigval_and_eigvec(eigVals, cx_eigVecs, threshold);
}

/*
 *  Funkcja zwracajaca wektory, na ktore rzutuja operatory rzutowe wystepujace w rozkladzie identycznosci
 */
std::vector<arma::cx_vec> armaPatch::Decomposition::get_basis() const {
    std::vector<arma::cx_vec> basis;
    for (arma::cx_mat projection : _m_projections) {
        arma::cx_vec eigVals;
        arma::cx_mat beta;
        arma::eig_gen(eigVals, beta, projection);
        for (unsigned i = 0; i < _m_dim; i++)
            if (std::abs(eigVals(i) - 1.0) < 1e-6) basis.push_back(beta.col(i));
            else if (std::abs(eigVals(i)) > 1e-6) {
                std::string str = "zmaleziono operator rzutowy o wartosci wlasnej."
                        + std::to_string(std::real(eigVals(i))) + " " + std::to_string(std::imag(eigVals(i))) + "i.";
                throw std::logic_error(str);
            }
    }
    return basis;
}

void armaPatch::Decomposition::print() const {
    std::cout << std::string(100, '#') << std::endl;
    for (arma::cx_mat projection : _m_projections) {
        std::cout << std::string(40, '-') << std::endl;
        projection.print("projection");
        arma::cx_vec eigVals;
        arma::cx_mat beta;
        arma::eig_gen(eigVals, beta, projection);
        for (unsigned i = 0; i < _m_dim; i++)
            if (std::abs(eigVals(i) - 1.0) < 0.0001) beta.col(i).print("wektor bazowy");
            else if (std::abs(eigVals(i)) > 0.0001) std::cout << "COS NIE TAK JEST Z TYM OPERATOREM RZUTOWYM, jego wartosc wlasna wynosi:" << eigVals(i) << std::endl;
    }
    std::cout << std::string(100, '#') << std::endl;
}