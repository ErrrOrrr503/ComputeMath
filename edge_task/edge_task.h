#include "../linear/linearsystem.h"

template <typename num_t>
class difeq_2_order {
    // y`` + p(x)y` + q(x)y = f(x)
    private:
        num_t (*_px) (num_t);
        num_t (*_qx) (num_t);
        num_t (*_fx) (num_t);
        size_t _net_size;
        num_t _h;
        num_t _start = 0;
        num_t _end = 0;
        num_t _y_start = 0;
        num_t _y_end = 0;
        num_t _a = 0;
        num_t _b = 0;
        int _edge_is_set = 0;
        int _verboseness = 1;
        size_t _print_width = 6;
        std::vector<num_t> _ans;
    public:
        difeq_2_order (num_t (*px) (num_t), num_t (*qx) (num_t), num_t (*fx) (num_t), num_t start, num_t end, size_t net_size);
        difeq_2_order () = delete;
        ~difeq_2_order ();
        void set_edge_1_kind (num_t y_start, num_t y_end);
        // (y' + ay) (start) = y_start
        // (y' + by) (end) = y_end
        void set_edge_23_kind (num_t y_start, num_t y_end, num_t a, num_t b);
        void set_verboseness (int verboseness);
        void solve_fin_dif_1_kind ();
        void solve_fin_dif_23_kind ();
        void print_ans ();
};

template <typename num_t>
difeq_2_order<num_t>::difeq_2_order (num_t (*px) (num_t), num_t (*qx) (num_t), num_t (*fx) (num_t), num_t start, num_t end, size_t net_size) : \
    _px (px), _qx (qx), _fx (fx), \
    _net_size (net_size), \
    _start (start), _end (end), \
    _ans (std::vector<num_t> (net_size))
{
    if (end <= start || net_size <= 0)
        throw std::runtime_error ("end must be greater than start and net size should be greater than 0");
    _h = (end - start) / (net_size - 1);
}

template <typename num_t>
difeq_2_order<num_t>::~difeq_2_order ()
{

}

template <typename num_t>
void difeq_2_order<num_t>::set_edge_1_kind (num_t y_start, num_t y_end)
{
    _y_end = y_end;
    _y_start = y_start;
    _edge_is_set = 1;
}

template <typename num_t>
void difeq_2_order<num_t>::set_edge_23_kind (num_t y_start, num_t y_end, num_t a, num_t b)
{
    _y_end = y_end;
    _y_start = y_start;
    _a = a;
    _b = b;
    _edge_is_set = 23;
}

template <typename num_t>
void difeq_2_order<num_t>::set_verboseness (int verboseness)
{
    _verboseness = verboseness;
}

template <typename num_t>
void difeq_2_order<num_t>::solve_fin_dif_1_kind ()
{
    if (_edge_is_set != 1)
        throw std::runtime_error ("edge task 1 kind must be set");
    matrix<num_t> M (_net_size, _net_size);
    std::vector<num_t> colon (_net_size);
    // calculate matrix
    num_t _1_h2 = 1 / _h / _h;
    num_t __2_h2 = -2 / _h / _h;
    // edge condition
    M[0][0] = 1;
    colon[0] = _y_start;
    // n-2 eqv in center
    for (size_t i = 1; i < _net_size - 1; i++) {
        num_t xi = _start + (_end - _start) * (num_t) i / (num_t) (_net_size - 1);
        M[i][i - 1] = _1_h2 - _px (xi) / 2 / _h;
        M[i][i] = __2_h2 + _qx (xi);
        M[i][i + 1] = _1_h2 + _px (xi) / 2 / _h;
        colon[i] = _fx (xi);
    }
    // edge condition
    M[_net_size - 1][_net_size - 1] = 1;
    colon[_net_size - 1] = _y_end;
    linearsystem<num_t> linsys (_net_size, M, colon);
    if (_verboseness) {
        linsys.print ();
    }
    linsys.solve_run_through ();
    _ans = linsys.ans ();
}

template <typename num_t>
void difeq_2_order<num_t>::solve_fin_dif_23_kind ()
{
    if (_edge_is_set != 23)
        throw std::runtime_error ("edge task 2 or 3 kind must be set");

/*
    std:: cout << "_start: " << _start << std::endl;
    std:: cout << "_end: " << _end << std::endl;
    std:: cout << "_y_start: " << _y_start << std::endl;
    std:: cout << "_y_end: " << _y_end << std::endl;
    std:: cout << "_a: " << _a << std::endl;
    std:: cout << "_b: " << _b << std::endl;
*/

    matrix<num_t> M (_net_size, _net_size);
    std::vector<num_t> colon (_net_size);
    // calculate matrix
    num_t _1_h2 = 1 / _h / _h;
    num_t __2_h2 = -2 / _h / _h;
    // edge condition
    M[0][0] = __2_h2 + _qx (_start) + 2 * _h * _a * (_1_h2 - _px (_start) / 2 / _h);
    M[0][1] = _1_h2 + _1_h2;
    colon[0] = _fx (_start) + 2 * _h * _y_start * (_1_h2 - _px (_start) / 2 / _h);
    // n-2 eqv in center
    for (size_t i = 1; i < _net_size - 1; i++) {
        num_t xi = _start + (_end - _start) * (num_t) i / (num_t) (_net_size - 1);
        M[i][i - 1] = _1_h2 - _px (xi) / 2 / _h;
        M[i][i] = __2_h2 + _qx (xi);
        M[i][i + 1] = _1_h2 + _px (xi) / 2 / _h;
        colon[i] = _fx (xi);
    }
    // edge condition
    M[_net_size - 1][_net_size - 1] = __2_h2 + _qx (_end) - 2 * _h * _b * (_1_h2 + _px (_end) / 2 / _h);
    M[_net_size - 1][_net_size - 2] = _1_h2 + _1_h2;
    colon[_net_size - 1] = _fx (_end) - 2 * _h * _y_end * (_1_h2 + _px (_end) / 2 / _h);
    linearsystem<num_t> linsys (_net_size, M, colon);
    if (_verboseness) {
        linsys.print ();
    }
    linsys.solve_run_through ();
    _ans = linsys.ans ();
}

template <typename num_t>
void difeq_2_order<num_t>::print_ans ()
{
    if (_verboseness == 0) {
        for (size_t i = 0; i < _net_size; i++) {
            std::cout.precision (10);
            std::cout << _ans[i] << std::endl;
        }
    }
}