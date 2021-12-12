#include "nonlinearexpr.h"
using namespace std;

typedef double num_t;

num_t myfunc (num_t x)
{
    return x * x * x;
}

num_t myfunc_differential (num_t x)
{
    return 3 * x * x;
}

int main ()
{
    nonlinearexpr<num_t> nonlinexpr (myfunc, -1, 10);
    //nonlinearexpr<num_t> nonlinexpr (sin, 1.6, 3.2);
    nonlinexpr.solve_divtwo (0, 1e-12);
    cout << "divtwo:" << endl;
    nonlinexpr.print_roots ();
    nonlinexpr.set_verboseness (1);
    nonlinexpr.set_differential_expr (myfunc_differential);
    //nonlinexpr.set_differential_expr (cos);
    nonlinexpr.solve_newton (100, 1e-14);
    cout << "newton:" << endl;
    nonlinexpr.print_roots ();
    nonlinexpr.solve_secant (1000, 1e-14);
    cout << "secant:" << endl;
    nonlinexpr.print_roots ();
    nonlinexpr.solve_chord (1000, 1e-14);
    cout << "chord:" << endl;
    nonlinexpr.print_roots ();
    return 0;
}