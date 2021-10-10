#include "matrix.h"

template <typename num_t>
class linearsystem {
    private:
        matrix<num_t> _matrix;
        std::vector<num_t> _colon;
        std::vector<num_t> _ans;
        size_t _N;
        size_t _print_width = 4;
        num_t _prec = 1e-8;
    public:
        linearsystem () = delete;
        linearsystem (size_t N);
        ~linearsystem ();
        void solve_gauss ();
        int float_cmp (num_t fl1, num_t fl2);
        void print (std::ostream& out);
        void scan (std::istream& in);
        void print ();
        void scan ();
        void print_ans (std::ostream& out);
        void print_ans ();
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
int linearsystem<num_t>::float_cmp (num_t fl1, num_t fl2)
{
    if (fl1 - fl2 > _prec)
        return 1;
    if (fl1 - fl2 < -_prec)
        return -1;
    return 0;
}

template <typename num_t>
void linearsystem<num_t>::solve_gauss ()
{
    matrix<num_t> solving_matrix (_matrix);
    std::vector<num_t> solving_colon (_colon);
    solving_matrix.conv_to_upper_triangle_with_colon (solving_colon);
    //now all we have to do is solve upper-triangle system
    if (!float_cmp (solving_matrix[_N - 1][_N - 1], 0))
        throw std::runtime_error ("can't resolve incompatible system or system with multiple solutions");
    for (ssize_t i = _N - 1; i >= 0; i--) {
        for (ssize_t j = _N - 1; j > i; j--)
            solving_colon[i] -= solving_matrix[i][j] * _ans[j];
        _ans[i] = solving_colon[i] / solving_matrix[i][i];
    }
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