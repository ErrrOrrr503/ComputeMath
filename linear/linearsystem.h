#include "matrix.h"

template <typename num_t>
class linearsystem {
    private:
        int _verboseness = 1;
        size_t _print_width = 4;
        size_t _N;
        matrix<num_t> _matrix;
        std::vector<num_t> _colon;
        std::vector<num_t> _ans;
    public:
        linearsystem () = delete;
        linearsystem (size_t N);
        ~linearsystem ();
        void solve_gauss ();
        void solve_jacobi (size_t iterations, num_t precision);
        void solve_seidel (size_t iterations, num_t precision);
        void solve_fastest_descend (size_t iterations, num_t precision);
        bool check_jacobi ();
        void print (std::ostream& out);
        void scan (std::istream& in);
        void print ();
        void scan ();
        void print_ans (std::ostream& out);
        void print_ans ();
        void set_verboseness (int verboseness);
};

template <typename num_t>
linearsystem<num_t>::linearsystem (size_t N) : _N(N), _matrix (matrix<num_t> (N, N)), _colon (std::vector<num_t> (N)), _ans (std::vector<num_t> (N))
{

}

template <typename num_t>
linearsystem<num_t>::~linearsystem ()
{
    
}

template <typename num_t>
void linearsystem<num_t>::print (std::ostream& out)
{
    for (size_t i = 0; i < _N; i++) {
        for (size_t j = 0; j < _N; j++) {
            out.width (_print_width);
            out.setf (std::ios::right);
            out << _matrix[i][j] << "   ";
        }
        out << "   |   " << _colon[i];
        out << std::endl;
    }
}

template <typename num_t>
void linearsystem<num_t>::print ()
{
    print (std::cout);
}

template <typename num_t>
void linearsystem<num_t>::scan (std::istream& in)
{
    for (size_t i = 0; i < _N; i++) {
        for (size_t j = 0; j < _N; j++)
            in >> _matrix[i][j];
        in >> _colon[i];
    }
}

template <typename num_t>
void linearsystem<num_t>::scan ()
{
    scan (std::cin);
}

template <typename num_t>
void linearsystem<num_t>::solve_gauss ()
{
    matrix<num_t> solving_matrix (_matrix);
    std::vector<num_t> solving_colon (_colon);
    solving_matrix.conv_to_upper_triangle_with_colon (solving_colon);
std::cout << "debug: upt = " << std::endl;
solving_matrix.print (); 
    //now all we have to do is solve upper-triangle system
    if (!solving_matrix.float_cmp (solving_matrix[_N - 1][_N - 1], 0))
        throw std::runtime_error ("can't resolve incompatible system or system with multiple solutions");
    for (ssize_t i = _N - 1; i >= 0; i--) {
        for (ssize_t j = _N - 1; j > i; j--)
            solving_colon[i] -= solving_matrix[i][j] * _ans[j];
        _ans[i] = solving_colon[i] / solving_matrix[i][i];
    }
}

template <typename num_t>
bool linearsystem<num_t>::check_jacobi ()
{
    for (size_t i = 0; i < _N; i++) {
        num_t sum_non_diag = 0;
        for (size_t j = 0; j < _N; j++)
            if (j != i)
                sum_non_diag += fabs (_matrix[i][j]);
        if (sum_non_diag > fabs (_matrix[i][i]))
            return false;
    }
    return true;
}

template <typename num_t>
void linearsystem<num_t>::solve_jacobi (size_t iterations, num_t precision)
{
    std::vector<num_t> x_even (_N), x_uneven (_N);
    for (size_t i = 0; i < _N; i++) {
        x_uneven[i] = 0;
        x_even[i] = precision + precision;
    }
    size_t i = 0;
    num_t metrics = 0;
    std::vector<num_t> *xi = nullptr, *xi_1 = nullptr;
    for (i = 0; (metrics = _matrix.vector_metrics_compare (x_uneven, x_even)) > precision && i < iterations; i++) {
        if (i % 2) {
            xi = &x_even;
            xi_1 = &x_uneven;
        }
        else {
            xi = &x_uneven;
            xi_1 = &x_even;
        }
        for (size_t j = 0; j < xi->size (); j++) {
            (*xi_1)[j] = _colon[j];
            for (size_t k = 0; k < xi->size (); k++) {
                if (k != j)
                    (*xi_1)[j] -= _matrix[j][k] * (*xi)[k];
            }
            (*xi_1)[j] /= _matrix[j][j];
        }
std::cout << "x" << i + 1 << " = ";
for (size_t i = 0; i < (*xi_1).size (); i++)
    std::cout << (*xi_1)[i] << " ";
std::cout << std::endl;
    }
std::cout << std::endl <<"iterations = " << i << std::endl << "precision = " << metrics << std::endl;
    _ans = *xi_1;
}

template <typename num_t>
void linearsystem<num_t>::solve_seidel (size_t iterations, num_t precision)
{
    std::vector<num_t> x_even (_N), x_uneven (_N);
    for (size_t i = 0; i < _N; i++) {
        x_uneven[i] = 0;
        x_even[i] = precision + precision;
    }
    size_t i = 0;
    num_t metrics = 0;
    std::vector<num_t> *xi = nullptr, *xi_1 = nullptr;
    for (i = 0; (metrics = _matrix.vector_metrics_compare (x_uneven, x_even)) > precision && i < iterations; i++) {
        if (i % 2) {
            xi = &x_even;
            xi_1 = &x_uneven;
        }
        else {
            xi = &x_uneven;
            xi_1 = &x_even;
        }
        for (size_t j = 0; j < xi->size (); j++) {
            (*xi_1)[j] = _colon[j];
            for (size_t k = 0; k < j; k++)
                (*xi_1)[j] -= _matrix[j][k] * (*xi_1)[k];
            for (size_t k = j + 1; k < xi->size (); k++)
                (*xi_1)[j] -= _matrix[j][k] * (*xi)[k];    
            (*xi_1)[j] /= _matrix[j][j];
        }
std::cout << "x" << i + 1 << " = ";
for (size_t i = 0; i < (*xi_1).size (); i++)
    std::cout << (*xi_1)[i] << " ";
std::cout << std::endl;
    }
std::cout << std::endl <<"iterations = " << i << std::endl << "precision = " << metrics << std::endl;
    _ans = *xi_1;
}

template<typename num_t>
void linearsystem<num_t>::solve_fastest_descend (size_t iterations, num_t precision)
{
    std::vector<num_t> x_even (_N), x_uneven (_N);
    for (size_t i = 0; i < _N; i++) {
        x_uneven[i] = 0;
        x_even[i] = precision + precision;
    }
    bool watch_precision = 1;
    if (precision == 0)
        watch_precision = 0;
    size_t i = 0;
    num_t metrics = 0;
    std::vector<num_t> *xi = nullptr, *xi_1 = nullptr;
    for (i = 0; (((metrics = _matrix.vector_metrics_compare (x_uneven, x_even)) > precision) || !watch_precision) && i < iterations; i++) {
        if (i % 2) {
            xi = &x_even;
            xi_1 = &x_uneven;
        }
        else {
            xi = &x_uneven;
            xi_1 = &x_even;
        }
        std::vector<num_t> ri (_matrix.mul_colon (*xi));
        for (size_t j = 0; j < ri.size (); j++)
            ri[j] -= _colon[j];
        num_t scalar_mul_ri_ri = 0, scalar_mul_Ari_ri;
        std::vector<num_t> Ari (_matrix.mul_colon (ri));
        for (size_t j = 0; j < ri.size (); j++) {
            scalar_mul_ri_ri += ri[j] * ri[j];
            scalar_mul_Ari_ri += Ari[j] * ri[j];
        }
        num_t ti = scalar_mul_ri_ri / scalar_mul_Ari_ri;
        for (size_t j = 0; j < xi->size (); j++)
            (*xi_1)[j] = (*xi)[j] - ti * ri[j];

        if (_verboseness == 2) {
            std::cout << "x" << i + 1 << " = ";
            for (size_t i = 0; i < (*xi_1).size (); i++)
                std::cout << (*xi_1)[i] << " ";
            std::cout << std::endl;
        }
    }
    if (_verboseness == 1)
        std::cout << std::endl <<"iterations = " << i << std::endl << "precision = " << metrics << std::endl;
    _ans = *xi_1;
}

template <typename num_t>
void linearsystem<num_t>::print_ans (std::ostream& out)
{
    for (size_t i = 0; i < _N; i++)
        out << "x" << i << " = " << _ans[i] << std::endl;
}

template <typename num_t>
void linearsystem<num_t>::print_ans ()
{
    print_ans (std::cout);
}

template <typename num_t>
void linearsystem<num_t>::set_verboseness (int verboseness)
{
    _verboseness = verboseness;
}