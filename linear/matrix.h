#include <iostream>
#include <vector>

template <typename num_t>
class matrix {
    private:
        num_t *_arr;
        size_t _lines;
        size_t _colons;
        size_t _print_width = 4;
        num_t _det (matrix& m);
    public:
        matrix () = delete;
        matrix (size_t lines, size_t colons);
        ~matrix ();

        size_t lines ();
        size_t colons ();
        void print (std::ostream& out);
        void scan (std::istream& in);
        void swap_colons (size_t first, size_t second);
        void change_colon (size_t colon_num, const std::vector<num_t>& new_colon);
        num_t det ();

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
}

template <typename num_t>
matrix<num_t>::~matrix ()
{
    delete[] _arr;
}

template <typename num_t>
void matrix<num_t>::print (std::ostream& out)
{
    //out.width (_print_width);
    for (size_t i = 0; i < _lines; i++) {
        for (size_t j = 0; j < _colons; j++)
            out << _arr[i * _colons + j] << "   ";
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