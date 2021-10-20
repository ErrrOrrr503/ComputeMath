#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstring>
#include <cmath>

template <typename num_t>
class matrix {
    private:
        num_t *_arr;
        size_t _lines;
        size_t _colons;
        size_t _print_width = 4;
        num_t _prec = 1e-10;

        num_t _det (matrix& m);
        num_t _relative_distance_to_one (num_t a);
    public:
        matrix () = delete;
        matrix (size_t lines, size_t colons);
        matrix (matrix<num_t>& m);
        ~matrix ();

        int float_cmp (num_t fl1, num_t fl2);
        num_t vector_metrics_compare (const std::vector<num_t>& v1, const std::vector<num_t>& v2);

        size_t lines ();
        size_t colons ();
        const num_t* data();
        void print (std::ostream& out);
        void scan (std::istream& in);
        void print ();
        void scan ();
        void swap_colons (size_t first, size_t second);
        void change_colon (size_t colon_num, const std::vector<num_t>& new_colon);
        void conv_to_upper_triangle ();
        std::vector<num_t> mul_colon (const std::vector<num_t> colon);
        std::vector<num_t>& conv_to_upper_triangle_with_colon (std::vector<num_t> &colon);
        num_t det ();
        int is_symmetric ();
        num_t calculate_machinery_precision ();
        bool is_3_diag ();

        // doorka inducted successfully, for 'operator[][]'
        num_t& elem (size_t i, size_t j);
        class matrix_line {
            private:
                matrix &_matrix;
                num_t _line;
            public:
                matrix_line (matrix &matrix, num_t line) : _matrix(matrix), _line(line) {}
                num_t& operator[] (size_t j) {
                    return _matrix.elem (_line, j);
                }
        };
        matrix_line operator[] (size_t i)
        {
            return matrix_line (*this, i);
        }
};

template <typename num_t>
matrix<num_t>::matrix (size_t lines, size_t colons)
{
    if (!(lines && colons))
        throw std::runtime_error ("Matrix dimensions must be grater 0!");
    _arr = new num_t[lines * colons];
    _lines = lines;
    _colons = colons;
    _prec = 5 * calculate_machinery_precision ();
}

template <typename num_t>
matrix<num_t>::matrix (matrix<num_t>& m)
{
    _arr = new num_t[m.lines () * m.colons ()];
    _lines = m.lines ();
    _colons = m.colons ();
    std::memcpy (_arr, m.data (), _lines * _colons * sizeof (num_t));
}

template <typename num_t>
matrix<num_t>::~matrix ()
{
    delete[] _arr;
}

template <typename num_t>
int matrix<num_t>::float_cmp (num_t fl1, num_t fl2)
{
    if (fl1 - fl2 > _prec)
        return 1;
    if (fl1 - fl2 < -_prec)
        return -1;
    return 0;
}

template <typename num_t>
void matrix<num_t>::print (std::ostream& out)
{
    for (size_t i = 0; i < _lines; i++) {
        for (size_t j = 0; j < _colons; j++) {
            out.width (_print_width);
            out.setf (std::ios::right);
            out << _arr[i * _colons + j] << "   ";
        }
        out <<std::endl;
    }
}

template <typename num_t>
void matrix<num_t>::scan (std::istream& in)
{
    for (size_t i = 0; i < _lines; i++) {
        for (size_t j = 0; j < _colons; j++)
            in >> _arr[i * _colons + j];
    }
}

template <typename num_t>
void matrix<num_t>::print ()
{
    print (std::cout);
}

template <typename num_t>
void matrix<num_t>::scan ()
{
    scan (std::cin);
}

template <typename num_t>
num_t& matrix<num_t>::elem (size_t i, size_t j) {
    if (i >= _lines || j >= _colons)
        throw std::runtime_error ("Invalid index!");
    return _arr[i * _colons + j];
}

template <typename num_t>
void matrix<num_t>::swap_colons (size_t first, size_t second)
{
    if (first >= _colons || second >= _colons)
        throw std::runtime_error ("Invalid colons for swapping");
    num_t tmp;
    for (size_t i = 0; i < _lines; i++) {
        tmp = _arr[i][first];
        _arr[i][first] = _arr[i][second];
        _arr[i][second] = tmp;
    }
}
template <typename num_t>
void matrix<num_t>::change_colon (size_t colon_num, const std::vector<num_t>& new_colon)
{
    if (new_colon.size () != _lines)
        throw std::runtime_error ("Invalid new_colon size");
    if (colon_num >= _colons)
        throw std::runtime_error ("Invalid colon_num");
    for (size_t i = 0; i < _lines; i++) {
        _arr[i][colon_num] = new_colon;
    }
}

template <typename num_t>
num_t matrix<num_t>::det ()
{
    if (_lines != _colons)
        throw std::runtime_error ("Matrix should be square to calculate determinant!");
    return _det (*this);
}

template <typename num_t>
const num_t* matrix<num_t>::data ()
{
    return _arr;
}

template <typename num_t>
size_t matrix<num_t>::lines ()
{
    return _lines;
}

template <typename num_t>
size_t matrix<num_t>::colons ()
{
    return _colons;
}

template <typename num_t>
num_t matrix<num_t>::_det (matrix& m)
{
    if (m.lines () == 2)
        return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    num_t det = 0;
    for (size_t i = 0; i < m.lines (); i++) {

        //gen minor
        matrix<num_t> mi (m.lines () - 1, m.colons () - 1);
        for (size_t k = 0; k < i; k++) {
            for (size_t l = 0; l < m.colons () - 1; l++)
                mi[k][l] = m[k][l+1];
        }
        for (size_t k = i + 1; k < m.lines (); k++) {
            for (size_t l = 0; l < m.colons () - 1; l++)
                mi[k - 1][l] = m[k][l+1];
        }

        if (i % 2)
            det -= m[i][0] * _det (mi);
        else
            det += m[i][0] * _det (mi);
    }
    return det;
}

template <typename num_t>
void matrix<num_t>::conv_to_upper_triangle ()
{
    std::vector<num_t> tmp (_lines);
    conv_to_upper_triangle_with_colon (tmp);
}

template <typename num_t>
num_t matrix<num_t>::_relative_distance_to_one (num_t a)
{
    if (a > 1)
        return 1 - 1/a;
    if (a < -1)
        return 1 + 1/a;
    if (a > 0)
        return 1 - a;
    return 1 + a;
}

template <typename num_t>
num_t matrix<num_t>::vector_metrics_compare (const std::vector<num_t>& v1, const std::vector<num_t>& v2)
{
    if (v1.size () != v2.size ())
        throw std::runtime_error ("Vectors for comparing should be equally dimensioned!");
    num_t metrics = 0;
    for (size_t i = 0; i < v1.size (); i++) {
        num_t tmp = fabs (v1[i] - v2[i]);
        if (tmp > metrics)
            metrics = tmp;
    }
    return metrics;
}

template <typename num_t>
std::vector<num_t>& matrix<num_t>::conv_to_upper_triangle_with_colon (std::vector<num_t> &colon)
{
    for (size_t i = 0; i < std::min (_lines, _colons); i++) {
        // locate most favorable line for division by M[i][i]
        num_t dist = _relative_distance_to_one (_arr[i * _colons + i]);
        size_t new_i = i;
        for (size_t j = i + 1; j < _lines; j++) {
            num_t tmp_dist = _relative_distance_to_one (_arr[j * _colons + i]);
            if (tmp_dist < dist) {
                dist = tmp_dist;
                new_i = j;
            }
        }
        // swap [i] and most favorable lines
        if (new_i != i) {
            num_t tmp;
            for (size_t k = 0; k < _colons; k++) {
                tmp = _arr[i * _colons + k];
                _arr[i * _colons + k] = _arr[new_i * _colons + k];
                _arr[new_i * _colons + k] = tmp;
            }
            tmp = colon[i];
            colon[i] = colon[new_i];
            colon[new_i] = tmp;
        }
        // convert
        for (size_t j = i + 1; j < _lines; j++) {
            num_t mul = _arr[j * _colons + i] / _arr[i * _colons + i];
            for (size_t k = 0; k < _colons; k++)
                _arr[j * _colons + k] = _arr[j * _colons + k] - _arr[i * _colons + k] * mul;
            colon[j] = colon[j] - colon[i] * mul;
        }
    }
    return colon;
}

template<typename num_t>
std::vector<num_t> matrix<num_t>::mul_colon (const std::vector<num_t> colon)
{
    if (_colons != colon.size ())
        throw std::runtime_error ("wrong vector dimension for multiplication");
    std::vector<num_t> ans (_lines);
    for (size_t i = 0; i < _lines; i++) {
        ans[i] = 0;
        for (size_t j = 0; j < _colons; j++)
            ans[i] += _arr[i * _colons + j] * colon[j];
    }
    return ans;
}

template <typename num_t>
int matrix<num_t>::is_symmetric ()
{
    if (_lines != _colons)
        return 0;
    if (_lines == 1)
        return 1;
    int sym_antisym = 1;
    if (_arr[0 * _colons + 1] == -_arr[1 * _colons + 0])
        sym_antisym = -1;
    for (size_t i = 0; i < _lines - 1; i++) {
        for (size_t j = i + 1; j < _colons; j++) {
            if (_arr[i * _colons + j] != sym_antisym * _arr[j * _colons + i])
                return 0;
        }
    }
    return sym_antisym;
}

template <typename num_t>
num_t matrix<num_t>::calculate_machinery_precision ()
{
    num_t a = 1;
    num_t b = 1;
    while (a + b != a)
        b /= 10;
    return b;
}

template <typename num_t>
bool matrix<num_t>::is_3_diag ()
{
    if (_lines != _colons)
        return false;
    if (_lines <= 2)
        return true;
    for (size_t i = 2; i < _lines; i++)
        for (size_t j = 0; j < i - 1; j++)
            if (_arr[i * _colons + j] || _arr[j * _colons + i])
                return false;
    return true;
}