#include<cassert>
#include<armadillo>
#include<cmath>
#include<vector>
#include<iostream>
#include<iomanip>
#include<memory>

#include<read_file_tools.hpp>
#include<arma_patch.hpp>
#include<burnside_algorithm.hpp>

void print_character_table(const std::vector<arma::cx_vec> & character_table, unsigned precision = 2) {
    auto f = std::cout.flags();
    std::cout << "Character table:" << std::endl;
    std::cout << "+";
    for (unsigned i = 0; i < character_table.size(); i++)
        std::cout << std::string(2 * precision + 6 + 1, '-') << '+';
    std::cout << std::endl;
    std::cout << "|";
    for (unsigned i = 0; i < character_table.size(); i++)
        std::cout << std::setw(2 * precision + 6 + 1) << "[" + std::to_string(i) + "]" << '|';
    std::cout << std::endl;
    std::cout << "+";
    for (unsigned i = 0; i < character_table.size(); i++)
        std::cout << std::string(2 * precision + 6 + 1, '-') << '+';
    std::cout << std::endl;
    std::cout << std::fixed << std::right << std::showpos << std::setprecision(precision);
    for (const arma::cx_vec& chis_row : character_table) {
        std::cout << "|";
        assert(chis_row.n_rows == character_table.size());
        for (auto chi : chis_row) {
            std::cout << std::setw(precision + 3) << real(chi);
            if (std::abs(imag(chi)) > pow(10, -int(precision)))
                std::cout << std::setw(precision + 3) << imag(chi) << "i" << "|";
            else
                std::cout << std::setw(precision + 4) << ' ' << "|";
        }
        std::cout << std::endl;
    }
    std::cout << "+";
    for (unsigned i = 0; i < character_table.size(); i++)
        std::cout << std::string(2 * precision + 6 + 1, '-') << '+';
    std::cout << std::endl;
    std::cout.flags(f);
}

void print_conj_classes(const std::vector<std::vector<unsigned> > & conjClasses) {
    auto f = std::cout.flags();
    std::cout << "Conjugacy classes:" << std::endl;
    std::cout << std::right;
    for (unsigned i = 0; i < conjClasses.size(); i++) {
        std::cout << std::setw(5) << "[" + std::to_string(i) + "]" << " = {";
        for (std::vector<unsigned>::const_iterator it = conjClasses[i].cbegin(); it != conjClasses[i].cend() - 1; it++) {
            std::cout << std::setw(3) << *it << ',';
        }
        std::cout << std::setw(3) << std::right << conjClasses[i].back();
        std::cout << "}" << std::endl;
    }
    std::cout.flags(f);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "The program determine the character table for given group." << std::endl;
        std::cerr << "Input:  the group multiplication table." << std::endl;
        std::cerr << "Output: the group character table." << std::endl;
        std::cerr << "" << std::endl;
        std::cerr << "Synopis:" << std::endl;
        std::cerr << argv[0] << " path_to_matrix_file_with_multiplication_table" << std::endl;                
        return 1;
    }
    // Czytanie pliku zawierajacego tabele mnozenia grupowego:
    std::string fileName = argv[1];
    std::shared_ptr<arma::Mat<unsigned>> multiplication_ptr(file_to_mat<unsigned>(fileName));
    std::cout << std::endl;
    multiplication_ptr->print("Group multiplication table:");
    std::cout << std::endl;
    // Znajdujemy klasy sprzezonosci:
    std::vector<std::vector<unsigned> > conj_classes = determine_conj_classes(*multiplication_ptr);
    print_conj_classes(conj_classes);
    std::cout << std::endl;
    // Znajdujemy tabele charakterow:
    std::vector<arma::cx_vec> character_table = build_character_table(*multiplication_ptr);
    // Wypisywanie wynikow:
    print_character_table(character_table);
}


/*
int main(){

        arma::cx_mat M;

        const double t = 3.4641;
        M << 0.0 << 0.0 << 0.0 << 2.0 << 0.0 << 0.0 << 0.0 << 0.0 << arma::endr
          << 0.0 << 0.0 << 2.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << arma::endr
          << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 << 4.0 << 0.0 << 0.0 << arma::endr
          << 0.0 << 0.0 << 0.0 << 0.0 << 4.0 << 0.0 << 0.0 << 0.0 << arma::endr
          << 0.0 << 2.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 <<   t << arma::endr
          << 2.0 << 0.0 << 0.0 << 0.0 << 0.0 << 0.0 <<   t << 0.0 << arma::endr
          << 0.0 << 0.0 << 0.0 <<   t << 0.0 << 0.0 << 0.0 << 0.0 << arma::endr
          << 0.0 << 0.0 <<   t << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;

        arma::cx_vec eigval;
        arma::cx_mat eigvec;

        armaPatch::eig_gen(eigval, eigvec, M);

        eigval.print("eigval");
        eigvec.print("eigvec");
	
}
 */