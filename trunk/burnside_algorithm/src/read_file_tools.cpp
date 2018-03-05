#include<armadillo>
#include<string>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<stdexcept>

const unsigned my_width = 50;

const std::string color_red("\033[31m");
const std::string color_green("\033[32m");
const std::string color_blue("\033[34m");
const std::string color_end("\033[0m");

template<typename T>
arma::Mat<T>* file_to_mat(std::string file_name) {
    auto f = std::cout.flags();
    std::cout << std::setw(my_width + color_blue.size() + color_end.size())
            << std::setfill('.') << std::left << "Otwieranie pliku: "
            << color_blue + file_name + color_end;
    std::ifstream fcin(file_name.data());
    if (!fcin) {
        std::cerr << color_red + " -> Problem z otwieraniem pliku." + color_end << std::endl;
        std::string str = "Blad podczas otwierania pliku.";
        throw std::runtime_error(str);
    }
    std::cout << color_green << " -> Plik poprawnie otwarty." << color_end << std::endl;
    std::cout << std::setw(my_width) << std::setfill('.') << std::left << "Czytamy rozmiar macierzy ...";
    unsigned N, M;
    fcin >> N >> M;
    if (!fcin) {
        std::cerr << color_red + " -> Wystopil blad podczas czytania macierzy!" + color_end;
        if (fcin.eof()) std::cerr << " [EoF bit set]";
        if (fcin.fail()) std::cerr << " [fail bit set]";
        if (fcin.bad()) std::cerr << " [bad bit set]";
        std::cerr << std::endl;
        std::string str = "Blad podczas czytania ilosci wierszy i kolumn.";
        throw std::runtime_error(str);
    }
    std::cout << color_green << " -> Rozmiar to: " << N << "x" << M << "." << color_end << std::endl;
    arma::Mat<T>* result = new arma::Mat<T>(N, M);
    std::cout << std::setw(my_width) << std::setfill('.') << std::left << "Czytamy macierz ...";
    for (unsigned i = 0; i < N; i++)
        for (unsigned j = 0; j < M; j++) {
            T temp;
            fcin >> temp;
            result->at(i, j) = temp;
        }
    if (!fcin) {
        std::cerr << color_red + " -> Wystopil blad podczas czytania macierzy!" + color_end;
        if (fcin.eof()) std::cerr << " [EoF bit set]";
        if (fcin.fail()) std::cerr << " [fail bit set]";
        if (fcin.bad()) std::cerr << " [bad bit set]";
        std::cerr << std::endl;
        std::string str = "Blad podczas czytania macierzy.";
        throw std::runtime_error(str);
    }
    std::cout << color_green + " -> Macierz przeczytana bez BledowStrumienia." << color_end << std::endl;
    fcin.close();
    std::cout << std::setfill(' ');
    std::cout.flags(f);
    return result;
}

template arma::Mat<double>* file_to_mat<double>(std::string file_name);
template arma::Mat<int>* file_to_mat<int>(std::string file_name);
template arma::Mat<unsigned>* file_to_mat<unsigned>(std::string file_name);