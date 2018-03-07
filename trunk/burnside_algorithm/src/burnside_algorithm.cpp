#include<armadillo>
#include<vector>
#include<list>

#include<arma_patch.hpp>
#include<burnside_algorithm.hpp>

// The adopted convention for group's theory related code:
// (i) Conjugacy classes are indexed by unsigned variables,
// the letters R,S,T,I,... are used.
// (ii) Groups elements are indexed by unsigned variables,
// the letters r,s,t,i,j, ... are used

/*
 * The functions that determines all conjugacy classes.
 */
std::vector<std::vector<unsigned> > determine_conj_classes(const arma::Mat<unsigned> & multiplication) {
    std::vector<std::vector<unsigned> > conj_classes;
    double size_of_group = multiplication.n_rows;
    std::vector<unsigned> Inv(size_of_group);
    for (unsigned i = 0; i < size_of_group; i++)
        for (unsigned j = 0; j < size_of_group; j++)
            if (multiplication(i, j) == 0) {
                Inv[i] = j;
                break;
            }
    std::list<unsigned> elements;
    for (unsigned i = 0; i < size_of_group; i++) elements.push_back(i);
    while (!elements.empty()) {
        std::vector<unsigned> conj_class;
        unsigned classRepresentant = *(elements.begin());
        for (std::list<unsigned>::iterator it = elements.begin(); it != elements.end();) {
            bool areConjugate = false;
            for (unsigned j = 0; j < size_of_group; j++)
                if (multiplication(multiplication(j, *it), Inv[j]) == classRepresentant) {
                    areConjugate = true;
                    break;
                }
            if (areConjugate) {
                conj_class.push_back(*it);
                it = elements.erase(it);
            } else it++;
        }
        conj_classes.push_back(conj_class);
    }
    return conj_classes;
}

// *********************************************************************************************
// ************************************  Burnside's algorithm  *********************************
// *********************************************************************************************

/*
 * The functions that determines (M_R)_{ST} matrices with c_{RST} matrix elements.
 */

std::vector<arma::mat> bulid_Ms(const arma::Mat<unsigned> & multiplication, std::vector<std::vector<unsigned> > conj_classes) {
    unsigned number_of_conj_classes = conj_classes.size();
    std::vector<unsigned> sizes_of_conj_classes;
    for (std::vector<unsigned> conj_class : conj_classes) sizes_of_conj_classes.push_back(conj_class.size());
    std::vector<arma::mat> Ms(number_of_conj_classes, arma::mat(number_of_conj_classes, number_of_conj_classes, arma::fill::zeros));
    for (unsigned R = 0; R < number_of_conj_classes; R++)
        for (unsigned S = 0; S < number_of_conj_classes; S++)
            for (unsigned r : conj_classes[R])
                for (unsigned s : conj_classes[S])
                    for (unsigned T = 0; T < number_of_conj_classes; T++)
                        for (unsigned t : conj_classes[T])
                            if (multiplication(r, s) == t) Ms[R](S, T) += 1;
    for (unsigned R = 0; R < number_of_conj_classes; R++)
        for (unsigned S = 0; S < number_of_conj_classes; S++)
            for (unsigned T = 0; T < number_of_conj_classes; T++)
                Ms[R](S, T) /= sqrt(sizes_of_conj_classes[S]) * sqrt(sizes_of_conj_classes[T]);
    return Ms;
}

/* 
 * The functions that determines character_table:
 * 
 * The arguments:
 *
 * multiplication:
 * the group multiplication table.
 * conventions and constrains:
 * the left index (row index) corresponds to the left operand in the group multiplication.
 * the right index (coll index) corresponds to the right operand in the group multiplication.
 * The neutral group element is indexed by the index equals to 0.
 *
 * conj_classes:
 * all the conjugacy classes for the given group.
 * A conjugacy class is represented by a vector of indices
 * corresponding to the groups elements belonging to the conjugacy class.
 * 
 * Returns:
 * A vector of vector of the character table
 * The inner vectors are the rows of the table.
 * Convention: there is one entry per one conjugacy class
 * (not one entry for one group element).
 */
std::vector<arma::cx_vec> build_character_table(const arma::Mat<unsigned> & multiplication, const std::vector<std::vector<unsigned> > & conj_classes) {
    // Tu bedzie zapisywany wynik tej tunkcji.
    std::vector<arma::cx_vec> character_table;
    // Pomocnicza wielkosci: rzad grupy, liosc klas sprzezonosci oraz ilosci elelentow w kolejnych klasach sprzezonosci
    const double size_of_group = multiplication.n_rows;
    const unsigned number_of_conj_classes = conj_classes.size();
    std::vector<unsigned> sizes_of_conj_classes;
    for (std::vector<unsigned> conjClas : conj_classes) sizes_of_conj_classes.push_back(conjClas.size());
    // Budujemy macierze M wystepujace w algorytmie Burnsidea i (wspolnie) je diagonalizujemy
    std::vector<arma::mat> Ms = bulid_Ms(multiplication, conj_classes);
    // Techniczne przerobienia max -> cx_max:
    std::vector<arma::cx_mat> cxMs;
    const arma::mat zero_mat(number_of_conj_classes, number_of_conj_classes, arma::fill::zeros);
    for (arma::mat M : Ms)
        cxMs.push_back(arma::cx_mat(M, zero_mat));
    // Analizyjemy wspolne wektory wlasne i wyciagamy z nich informacje o charakterach (jeden wektro wlasny - jeden wiersz w tab charakterow)
    armaPatch::Decomposition decomposition = armaPatch::common_eig_gen(cxMs);
    for (arma::cx_vec chis_row : decomposition.get_basis()) {
        // Najpierw z wektora bazowego robimy wektor, ktory bedzie proporcjonalny do wiersza w tab charakterow
        for (unsigned i = 0; i < number_of_conj_classes; i++) chis_row(i) /= sqrt(sizes_of_conj_classes[i]);
        // Normowanie z norma jak dla charakterow
        arma::cx_double group_norm = 0;
        for (unsigned I = 0; I < number_of_conj_classes; I++) group_norm += (sizes_of_conj_classes[I]) * std::real(conj(chis_row(I)) * chis_row(I));
        group_norm = sqrt(group_norm);
        chis_row /= (group_norm / sqrt(size_of_group));
        // Robimy tak, by charakter elementu neutralnego byl dodatni
        arma::cx_double z = chis_row(0) / abs(chis_row(0));
        chis_row /= z;
        // Gotowy wreszczie wiersz tabeli charakterow dodajemy do pojemnika na wyniki:
        character_table.push_back(chis_row);
    }
    // Zwracamy wynik:
    return character_table;
}

std::vector<arma::cx_vec> build_character_table(const arma::Mat<unsigned> & multiplication) {
    // Argumenty funkcji:
    //
    // mat &multiplication to tabela mnozenia w grupie; 
    // Konwencje:
    // lewy indeks(indeks wiersza)  odpowiada elementowi stojącemu po lewej  stronie znaku dzialania grupowego.
    // prawy indeks(indeks kolumny) odpowiada elementowi stojącemu po prawej stronie znaku dzialania grupowego.
    // Elementem grupy odpowiadającym indeksowi równemu zero jest element neutralny.
    std::vector<std::vector<unsigned> > conj_classes = determine_conj_classes(multiplication);
    return build_character_table(multiplication, conj_classes);
}