#include <iostream>
#include <vector>
#include <cmath>

template <typename num_t>
struct root {
    num_t val = 0;
    num_t froot = 0;
    num_t multiplicity = -1;
    root (num_t val_, num_t froot_, num_t multiplicity_) :
        val (val_), froot (froot_), multiplicity (multiplicity_)
    {
    }
};

template <typename num_t>
class nonlinearexpr
{
    private:
        num_t (*_expr) (num_t x) = nullptr;
        num_t (*_differential_expr) (num_t x) = nullptr;
        std::vector<root<num_t>> _roots;
        std::ostream& _out = std::cout;
        num_t _start;
        num_t _end;
        size_t _verboseness = 0;
        num_t _calculate_xi_1_from_xi_newton (num_t xi);
        num_t _calculate_xi_1_from_xi_secant (num_t xi);
        num_t _calculate_xi_1_from_xi_chord (num_t xi);
        void _solve_iterated (num_t (nonlinearexpr<num_t>::*_calculate_xi_1_from_xi) (num_t), num_t strt_x, size_t iterations, num_t precision);
        void _solve_divtwo (num_t start, num_t end, size_t depth, num_t precision);
    public:
        nonlinearexpr () = delete;
        nonlinearexpr (num_t (*expression) (num_t x), num_t start, num_t end);
        ~nonlinearexpr ();
        void print_roots ();
        void set_out (std::ostream& out);
        void set_verboseness (size_t verboseness);
        void solve_divtwo (size_t depth, num_t precision);
        void solve_newton (size_t iterations, num_t precision);
        void solve_secant (size_t iterations, num_t precision);
        void solve_chord (size_t iterations, num_t precision);
        void set_differential_expr (num_t (*differential_expr) (num_t x));
};

template <typename num_t>
nonlinearexpr<num_t>::nonlinearexpr (num_t (*expression) (num_t x), num_t start, num_t end) : _expr (expression), _start(start), _end(end)
{

}

template <typename num_t>
nonlinearexpr<num_t>::~nonlinearexpr ()
{

}

template <typename num_t>
void nonlinearexpr<num_t>::set_out (std::ostream& out)
{
    _out = out;
}

template <typename num_t>
void nonlinearexpr<num_t>::set_verboseness (size_t verboseness)
{
    _verboseness = verboseness;
}

template <typename num_t>
void nonlinearexpr<num_t>::set_differential_expr (num_t (*differential_expr) (num_t x))
{
    _differential_expr = differential_expr;
}

template <typename num_t>
void nonlinearexpr<num_t>::print_roots ()
{
    if (!_roots.size ()) {
        _out << "no roots" << std::endl;
        return;
    }
    for (size_t i = 0; i < _roots.size(); i++) {
        _out << "x" << i << " = " << _roots[i].val;
        if (_roots[i].multiplicity > 0)
            _out << " multiplicity = " << _roots[i].multiplicity << std::endl;
        else
            _out << " multiplicity undefined" << std::endl;
    }
}

template <typename num_t>
void nonlinearexpr<num_t>::solve_divtwo (size_t depth, num_t precision)
{
    _roots.clear ();
    _solve_divtwo (_start, _end, depth, precision);
}

template <typename num_t>
void nonlinearexpr<num_t>::solve_newton (size_t iterations, num_t precision)
{
    if (_differential_expr == nullptr)
        throw std::runtime_error ("differential expr is not set!");
    _roots.clear ();
    _solve_iterated (&nonlinearexpr<num_t>::_calculate_xi_1_from_xi_newton, (_end + _start) / 2, iterations, precision);
}

template <typename num_t>
void nonlinearexpr<num_t>::solve_secant (size_t iterations, num_t precision)
{
    _roots.clear ();
    _calculate_xi_1_from_xi_secant ((_end + _start) / 2 + 10 * precision); //first_run to init i-1;
    _solve_iterated (&nonlinearexpr<num_t>::_calculate_xi_1_from_xi_secant, (_end + _start) / 2 - 10 * precision, iterations, precision);
    _roots[0].multiplicity = -1;
}

template <typename num_t>
void nonlinearexpr<num_t>::solve_chord (size_t iterations, num_t precision)
{
    _roots.clear ();
    _calculate_xi_1_from_xi_chord ((_end + _start) / 2 + 10 * precision); //first_run to init i-1;
    _solve_iterated (&nonlinearexpr<num_t>::_calculate_xi_1_from_xi_chord, (_end + _start) / 2 - 10 * precision, iterations, precision);
    _roots[0].multiplicity = -1;
}

template <typename num_t>
void nonlinearexpr<num_t>::_solve_divtwo (num_t start, num_t end, size_t depth, num_t precision)
{
    // assume that there are no more than 1 root on 20x precision interval.
    num_t fstart = _expr (start);
    num_t fend = _expr (end);
    if (fstart <= precision && fstart >= -precision) {
        if (_roots.size ()) {
            const root<num_t> &last = _roots[_roots.size() - 1];
            if (last.val - start > 20 * precision || last.val - start < -20 * precision)
                _roots.push_back (root<num_t> (start, fstart, -1));
        }
        else
            _roots.push_back (root<num_t> (start, fstart, -1));
    }

    if (end - start > precision + precision) {
        if (fstart * fend < 0) {
            _solve_divtwo (start, (start + end) / 2, depth, precision);
            _solve_divtwo ((start + end) / 2, end, depth, precision);
        }
        else if (depth) {
            _solve_divtwo (start, (start + end) / 2, depth - 1, precision);
            _solve_divtwo ((start + end) / 2, end, depth - 1, precision);
        }
    }

    if (fend <= precision && fend >= -precision) {
        if (_roots.size ()) {
            const root<num_t> &last = _roots[_roots.size() - 1];
            if (last.val - end > 20 * precision || last.val - end < -20 * precision)
                _roots.push_back (root<num_t> (end, fend, -1));
        }
        else
            _roots.push_back (root<num_t> (end, fend, -1));
    }
}

template <typename num_t>
void nonlinearexpr<num_t>::_solve_iterated (num_t (nonlinearexpr<num_t>::*_calculate_xi_1_from_xi) (num_t), num_t strt_x, size_t iterations, num_t precision)
{
    num_t x_even = strt_x, x_uneven = strt_x + precision + precision, x__1 = 0;
    num_t multiplicity = -1;
    bool watch_precision = 1;
    if (precision == 0)
        watch_precision = 0;
    size_t i = 0;
    num_t *xi = nullptr, *xi_1 = nullptr;
    for (i = 0; ((x_even - x_uneven > precision || x_even - x_uneven < -precision) || !watch_precision) && i < iterations; i++) {
        if (i % 2) {
            xi = &x_even;
            xi_1 = &x_uneven;
        }
        else {
            xi = &x_uneven;
            xi_1 = &x_even;
        }
        
        *xi_1 = (this->*_calculate_xi_1_from_xi) (*xi);
        if (i >= 1)
            multiplicity = 1 / (1 - (*xi_1 - *xi) / (*xi - x__1));

        if (_verboseness == 2) {
            _out << "[dbg]: x" << i + 1 << " = " << *xi_1 << " multiplicity = " << multiplicity << std::endl;
        }
        x__1 = *xi;
    }
    if (_verboseness == 1)
        _out << std::endl <<"iterations = " << i << std::endl << "precision = " << _expr (*xi_1) << std::endl;
    _roots.push_back (root (*xi_1, _expr (*xi_1), multiplicity));
}

template <typename num_t>
num_t nonlinearexpr<num_t>::_calculate_xi_1_from_xi_newton (num_t xi)
{
    return xi - _expr (xi) / _differential_expr (xi);
}

template <typename num_t>
num_t nonlinearexpr<num_t>::_calculate_xi_1_from_xi_chord (num_t xi)
{
    static num_t x0;
    static bool first_run = 0;
    if (!first_run) {
        x0 = xi;
        first_run = 1;
        return 0;
    }
    return xi - _expr (xi) * (xi - x0) / (_expr (xi) - _expr (x0));
}

template <typename num_t>
num_t nonlinearexpr<num_t>::_calculate_xi_1_from_xi_secant (num_t xi)
{
    static num_t xi_1;
    static bool first_run = 0;
    if (!first_run) {
        xi_1 = xi;
        first_run = 1;
        return 0;
    }
    num_t res = xi - _expr (xi) * (xi - xi_1) / (_expr (xi) - _expr (xi_1));
    xi_1 = xi;
    return res;
}