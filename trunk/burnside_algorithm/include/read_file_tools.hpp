#ifndef READ_FILE_TOOLS_HPP
#define READ_FILE_TOOLS_HPP

#include<armadillo>
#include<string>

template<typename T>
arma::Mat<T>* file_to_mat(std::string file_name);

#endif