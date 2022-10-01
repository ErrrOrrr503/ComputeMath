#include <iostream>
#include <stdexcept>
#include <vector>
#include <cstring>
#include <cmath>
#ifdef OMP
#include <omp.h>
#endif

template <typename num_t>
class matrix {
    private:
        num_t *_arr;
        size_t _lines;
        size_t _colons;
        size_t _print_width = 4;
        num_t _prec = 1e-10;
        int _verboseness = 1;

        num_t _det (matrix& m);
        num_t _calculate_minor (matrix& m, size_t i, size_t j);
        num_t _relative_distance_to_one (num_t a);
    public:
        bool _arr_is_custom = 0;
        matrix () = delete;
        matrix (size_t lines, size_t colons);
        matrix (size_t lines, size_t colons, num_t prec);
        matrix (size_t lines, size_t colons, num_t *custom_arr);
        matrix (const matrix<num_t>& m);
        ~matrix ();

        int float_cmp (num_t fl1, num_t fl2);
        num_t vector_metrics_compare (const std::vector<num_t>& v1, const std::vector<num_t>& v2);
        num_t matrix_metrics_1 ();

        size_t lines () const;
        size_t colons () const;
        const num_t* data();
        void print (std::ostream& out) const ;
        void scan (std::istream& in);
        void random_generate ();
        void print () const;
        void scan ();
        void swap_colons (size_t first, size_t second);
        void change_colon (size_t colon_num, const std::vector<num_t>& new_colon);
        void conv_to_upper_triangle ();
        std::vector<num_t> mul_colon (const std::vector<num_t> &colon);
        std::vector<num_t>& conv_to_upper_triangle_with_colon (std::vector<num_t> &colon);
        matrix<num_t> transpose () const;
        num_t det ();
        matrix<num_t> inverse_matrix ();
        num_t calculate_conditionality_number ();
        int is_symmetric ();
        num_t calculate_machinery_precision ();
        bool is_3_diag ();
        void set_verboseness (int verboseness);
        num_t* matrix_c_array_WARNING_rw ();
        const num_t* matrix_c_array () const;

        //friend matrix<num_t> operator* (const matrix<num_t>& M1, const matrix<num_t>& M2);

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

        friend matrix<num_t> operator+ (const matrix<num_t>& M1, const matrix<num_t>& M2)
        {
            return add (M1, M2);
        }
        friend matrix<num_t> operator- (const matrix<num_t>& M1, const matrix<num_t>& M2)
        {
            return sub (M1, M2);
        }
        friend matrix<num_t> operator* (const matrix<num_t>& M1, const matrix<num_t>& M2)
        {
            return line_mul (M1, M2);
        }
};

template <typename num_t>
matrix<num_t> add (const matrix<num_t>& M1, const matrix<num_t>& M2)
{
    if (M1.colons () != M2.colons () || M1.lines () != M2.lines ())
        throw std::runtime_error ("Summed matrixes should have same sizes!");
    size_t colons = M1.colons ();
    size_t lines = M1.lines ();
    matrix<num_t> M = matrix<num_t> (M1.lines (), M1.colons ());
    num_t *p_M = M.matrix_c_array_WARNING_rw ();
    const num_t *p_M1 = M1.matrix_c_array ();
    const num_t *p_M2 = M2.matrix_c_array ();
    #ifdef OMP
    #pragma omp parallel for simd shared (p_M1, p_M2, p_M, colons, lines)
    #endif
    for (size_t i = 0; i < lines; i++) {
        size_t temp = i * colons;
        num_t *p_M_i = p_M + temp;
        const num_t *p_M1_i = p_M1 + temp;
        const num_t *p_M2_i = p_M2 + temp;
        for (size_t j = 0; j < colons; j++) {
            p_M_i[j] = p_M1_i[j] + p_M2_i[j];
        }
    }
}

template <typename num_t>
matrix<num_t> sub (const matrix<num_t>& M1, const matrix<num_t>& M2)
{
    if (M1.colons () != M2.colons () || M1.lines () != M2.lines ())
        throw std::runtime_error ("Summed matrixes should have same sizes!");
    size_t colons = M1.colons ();
    size_t lines = M1.lines ();
    matrix<num_t> M = matrix<num_t> (M1.lines (), M1.colons ());
    num_t *p_M = M.matrix_c_array_WARNING_rw ();
    const num_t *p_M1 = M1.matrix_c_array ();
    const num_t *p_M2 = M2.matrix_c_array ();
    #ifdef OMP
    #pragma omp parallel for simd shared (p_M1, p_M2, p_M, colons, lines)
    #endif
    for (size_t i = 0; i < lines; i++) {
        size_t temp = i * colons;
        num_t *p_M_i = p_M + temp;
        const num_t *p_M1_i = p_M1 + temp;
        const num_t *p_M2_i = p_M2 + temp;
        for (size_t j = 0; j < colons; j++) {
            p_M_i[j] = p_M1_i[j] - p_M2_i[j];
        }
    }
}

template <typename num_t>
matrix<num_t> line_mul (const matrix<num_t>& M1, const matrix<num_t>& M2)
{
    if (M1.colons () != M2.lines () || M1.lines () != M2.colons ())
        throw std::runtime_error ("Multiplying matrixes should have corresponding sizes!");
    matrix<num_t> T = M2.transpose ();
    const num_t *p_T = T.matrix_c_array ();
    matrix<num_t> M = matrix<num_t> (M1.lines (), M2.colons ());
    num_t *p_M = M.matrix_c_array_WARNING_rw ();
    const num_t *p_M1 = M1.matrix_c_array ();
    size_t lines = M.lines ();
    size_t colons = M.colons ();
    size_t reduced_colons = M1.colons ();
    #ifdef OMP
    #pragma omp parallel for simd shared (p_M1, p_T, p_M, colons, lines, reduced_colons)
    #endif
    for (size_t i = 0; i < lines; i++) {
        num_t *p_M_i = p_M + i * colons;
        const num_t *p_M1_i = p_M1 + i * reduced_colons;
        for (size_t j = 0; j < colons; j++) {
            num_t *p_M_ij = p_M_i + j;
            const num_t *p_T_j = p_T + j * reduced_colons;
            for (size_t k = 0; k < reduced_colons; k++)
                *p_M_ij += p_M1_i[k] * p_T_j[k];
        }
    }
    return M;
}

template <typename num_t>
matrix<num_t> plain_mul (const matrix<num_t>& M1, const matrix<num_t>& M2)
{
    if (M1.colons () != M2.lines () || M1.lines () != M2.colons ())
        throw std::runtime_error ("Multiplying matrixes should have corresponding sizes!");
    const num_t *p_M2 = M2.matrix_c_array ();
    matrix<num_t> M = matrix<num_t> (M1.lines (), M2.colons ());
    num_t *p_M = M.matrix_c_array_WARNING_rw ();
    const num_t *p_M1 = M1.matrix_c_array ();
    size_t lines = M.lines ();
    size_t colons = M.colons ();
    size_t reduced_colons = M1.colons ();
    #ifdef OMP
    #pragma omp parallel for simd shared (p_M1, p_M2, p_M, colons, lines, reduced_colons)
    #endif
    for (size_t i = 0; i < lines; i++) {
        num_t *p_M_i = p_M + i * colons;
        const num_t *p_M1_i = p_M1 + i * reduced_colons;
        for (size_t j = 0; j < colons; j++) {
            num_t *p_M_ij = p_M_i + j;
            const num_t *p_M2_j = p_M2 + j;
            for (size_t k = 0; k < reduced_colons; k++)
                *p_M_ij += p_M1_i[k] * p_M2_j[k * reduced_colons];
        }
    }
    return M;
}

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
matrix<num_t>::matrix (size_t lines, size_t colons, num_t prec) : _lines (lines), _colons (colons), _prec (prec)
{
    _arr = new num_t[lines * colons];
}

template <typename num_t>
matrix<num_t>::matrix (size_t lines, size_t colons, num_t *custom_arr)
{
    if (!(lines && colons))
        throw std::runtime_error ("Matrix dimensions must be grater 0!");
    _arr = custom_arr;
    _arr_is_custom = 1;
    _lines = lines;
    _colons = colons;
    _prec = 5 * calculate_machinery_precision ();
}

template <typename num_t>
matrix<num_t>::matrix (const matrix<num_t>& m)
{
    _lines = m.lines ();
    _colons = m.colons ();
    if (!m._arr_is_custom)
        _arr = new num_t[m.lines () * m.colons ()];
    std::memcpy (_arr, m.matrix_c_array (), _lines * _colons * sizeof (num_t));
}

template <typename num_t>
matrix<num_t>::~matrix ()
{
    if (!_arr_is_custom)
        delete[] _arr;
}

template <typename num_t>
num_t* matrix<num_t>::matrix_c_array_WARNING_rw ()
{
    return _arr;
}

template <typename num_t>
const num_t* matrix<num_t>::matrix_c_array () const
{
    return _arr;
}

template <typename num_t>
void matrix<num_t>::random_generate ()
{
    #ifdef OMP
    #pragma omp parallel for shared (_arr)
    #endif
    for (size_t i = 0; i < _lines; i++) {
        num_t *_arr_i = _arr + i * _colons;
        for (size_t j = 0; j < _colons; j++) {
            _arr_i[j] = rand();
        }
    }
}

template <typename num_t>
void matrix<num_t>::set_verboseness (int verboseness)
{
    _verboseness = verboseness;
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
void matrix<num_t>::print (std::ostream& out) const
{
    if (_verboseness > 0) {
        for (size_t i = 0; i < _lines; i++) {
            for (size_t j = 0; j < _colons; j++) {
                out.width (_print_width);
                out.setf (std::ios::right);
                out << _arr[i * _colons + j] << "   ";
            }
            out <<std::endl;
        }
    }
    if (_verboseness == 0) {
        for (size_t i = 0; i < _lines; i++) {
            for (size_t j = 0; j < _colons; j++) {
                out.width (12);
                out.setf (std::ios::right);
                out << _arr[i * _colons + j] << " ";
            }
            out <<std::endl;
        }
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
void matrix<num_t>::print () const
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
    if (_arr == NULL)
        throw std::runtime_error ("trying to index uninit matrix");
    if (i >= _lines || j >= _colons)
        throw std::runtime_error ("Invalid index at elem method: i = " + std::to_string (i) + " j = " + std::to_string (j)+ " lines = " + std::to_string (_lines)+ " colons = " + std::to_string (_colons));
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
matrix<num_t> matrix<num_t>::transpose () const
{
    matrix<num_t> T (_colons, _lines, _prec);
    num_t *p_T = T.matrix_c_array_WARNING_rw ();
    #ifdef OMP
    #pragma omp parallel for simd shared (_arr, p_T, _colons, _lines)
    #endif
    for (size_t i = 0; i < _lines; i++) {
        num_t *p_T_i = p_T + i;
        const num_t *_arr_i = _arr + i * _colons;
        for (size_t j = 0; j < _colons; j++) {
            p_T_i[j * _lines] = _arr_i[j];
        }
    }
    return T;
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
size_t matrix<num_t>::lines () const
{
    return _lines;
}

template <typename num_t>
size_t matrix<num_t>::colons () const
{
    return _colons;
}

template <typename num_t>
num_t matrix<num_t>::_calculate_minor (matrix& m, size_t i, size_t j)
{
    matrix<num_t> mij(m.lines() - 1, m.colons() - 1);
    for (size_t k = 0; k < i; k++) {
        for (size_t l = 0; l < j; l++)
            mij[k][l] = m[k][l];
        for (size_t l = j; l < m.colons () - 1; l++)
            mij[k][l] = m[k][l + 1];
    }
    for (size_t k = i; k < m.lines() - 1; k++) {
        for (size_t l = 0; l < j; l++)
            mij[k][l] = m[k + 1][l];
        for (size_t l = j; l < m.colons () - 1; l++)
            mij[k][l] = m[k + 1][l + 1];
    }
    if ((i + j) % 2)
        return -_det (mij);
    return _det (mij);
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
std::vector<num_t> matrix<num_t>::mul_colon (const std::vector<num_t> &colon)
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

template <typename num_t>
matrix<num_t> matrix<num_t>::inverse_matrix ()
{
    if (_lines != _colons)
        throw std::runtime_error ("Matrix should be square to calculate invrse matrix!");
    matrix<num_t> inverse (_lines, _colons);
    num_t D = det ();
    for (size_t k = 0; k < _lines; k++)
        for (size_t l = 0; l < _colons; l++)
            inverse[k][l] = _calculate_minor (*this, k, l) / D;
    return inverse;
}

template <typename num_t>
num_t matrix<num_t>::matrix_metrics_1 ()
{
    num_t metric = 0, temp = 0;
    for (size_t i = 0; i < _lines; i++) {
        temp = 0;
        for (size_t j = 0; j < _colons; j++)
            temp += _arr[i * _colons + j];
        if (temp > metric)
            metric = temp;
    }
    return metric;
}

template <typename num_t>
num_t matrix<num_t>::calculate_conditionality_number ()
{
    matrix<num_t> inverse (inverse_matrix ());
    return matrix_metrics_1 () * inverse.matrix_metrics_1 ();
}