#include "integrate.h"
using namespace std;

typedef double flt;

flt fin_left (flt x)
{
    return cos (x) / sqrt (x);
}

flt fin_helper_left (flt x)
{
    return 1 * sqrt (x) - sqrt (x * x * x * x * x) / 5;
}

flt fin_right (flt x)
{
    return cos (x) / sqrt (-x);
}

flt fin_helper_right (flt x)
{
    return -1 * sqrt (-x) + sqrt (-x * x * x * x * x) / 5;
}

flt inf_left (flt x)
{
    return exp(x);
}

flt inf_helper_left (flt x)
{
    return exp(x);
}

flt inf_right (flt x)
{
    return exp(-x);
}

flt inf_helper_right (flt x)
{
    return -exp(-x);
}

flt F (flt x)
{
    return  3 * x * x;
}

int main ()
{
    //proper 
    cout.precision (15);
    netfunc<flt> f (F, 0, 1, 10);
    integrator<flt> I (f);
    I.integrate_rectangles ();
    cout << "Integral rectangles = " << I.ans << endl;
    I.integrate_trapezoids ();
    cout << "Integral trapezoids = " << I.ans << endl;
    if (I.integrate_simpson ())
        cout << I.get_err () << endl;
    else 
        cout << "Integral simpson = " << I.ans << endl;
    I.enable_runge (1e-10);
    I.integrate_rectangles ();
    cout << "Integral rectangles runge = " << I.ans << " runge = " << I.runge << " iterations = " << I.iterations << endl;
    I.integrate_trapezoids ();
    cout << "Integral trapezoids runge = " << I.ans << " runge = " << I.runge << " iterations = " << I.iterations << endl;
    if (I.integrate_simpson ())
        cout << I.get_err () << endl;
    else
        cout << "Integral simpson runge = " << I.ans << " runge = " << I.runge << " iterations = " << I.iterations << endl;
    I.assume_pure_net ();
    I.integrate_rectangles ();
    cout << "Integral rectangles pure net = " << I.ans << endl;
    if (I.integrate_simpson ())
        cout << I.get_err () << endl;
    else
        cout << "Integral simpson pure net = " << I.ans << endl;

    //improper
    //finite left
    netfunc<flt> net_fin_left (fin_left, fin_helper_left, 0, 1, improper_left);
    integrator<flt> Im_fin_left (net_fin_left);
    if (Im_fin_left.integrate_improper (1e-3))
        cout << Im_fin_left.get_err () << endl;
    else
        cout << "Integral improper finite left = " << Im_fin_left.ans << " edge = " << Im_fin_left.improper_edge << endl;
    
    //finite right
    netfunc<flt> net_fin_right (fin_right, fin_helper_right, -1, 0, improper_right);
    integrator<flt> Im_fin_right (net_fin_right);
    if (Im_fin_right.integrate_improper (1e-3))
        cout << Im_fin_right.get_err () << endl;
    else
        cout << "Integral improper finite right = " << Im_fin_right.ans << " edge = " << Im_fin_right.improper_edge << endl;

    //infinite left
    netfunc<flt> net_inf_left (inf_left, inf_helper_left, NAN, 0, improper_left);
    integrator<flt> Im_inf_left (net_inf_left);
    if (Im_inf_left.integrate_improper (1e-5))
        cout << Im_inf_left.get_err () << endl;
    else
        cout << "Integral improper infinite left = " << Im_inf_left.ans << " edge = " << Im_inf_left.improper_edge << endl;

    //infinite right
    netfunc<flt> net_inf_right (inf_right, inf_helper_right, 0, NAN, improper_right);
    integrator<flt> Im_inf_right (net_inf_right);
    if (Im_inf_right.integrate_improper (1e-5))
        cout << Im_inf_right.get_err () << endl;
    else
        cout << "Integral improper infinite right = " << Im_inf_right.ans << " edge = " << Im_inf_right.improper_edge << endl;
    return 0;
}