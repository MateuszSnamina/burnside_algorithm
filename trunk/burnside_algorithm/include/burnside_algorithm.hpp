#ifndef BURNSIDE_ALGORITHM_HPP
#define BURNSIDE_ALGORITHM_HPP

#include<armadillo>
#include<vector>

std::vector<std::vector<unsigned>> determine_conj_classes(const arma::Mat<unsigned> & multiplication);
// funkcja ta buduje macierze (M_R)_{ST} o elementach c_{RST}:
std::vector<arma::mat> bulid_Ms(const arma::Mat<unsigned> & multiplication, std::vector<std::vector<unsigned> > conj_classes);
std::vector<arma::cx_vec> build_character_table(const arma::Mat<unsigned> & multiplication, const std::vector<std::vector<unsigned>> &conj_classes);
std::vector<arma::cx_vec> build_character_table(const arma::Mat<unsigned> & multiplication);

#endif