#include "linearsystem.h"

/* TODO:
prec = 0;
verboseness;
check for det > 0 && symmetry
get rid of copypaste in iteration and variation methods
*/

int main ()
{
    size_t dim;
    std::cin >> dim;
    linearsystem<double> syslin (dim);
    syslin.scan ();
    std::cout << "syslin = " << std::endl;
    syslin.print ();
    syslin.solve_gauss ();
    syslin.print_ans ();
    std::cout << std::endl << "jacobi: " <<std::endl;
    if (!syslin.check_jacobi ())
        std::cout << "Hadamard condition is not met! Discrepancy is inevitable!" << std::endl;
    else
        std::cout << "Hadamard condition is met! Discrepancy is highly unlikely." << std::endl;
    syslin.solve_jacobi (100, 1e-6);
    syslin.print_ans ();
    std::cout << std::endl << "Seidel: " <<std::endl;
    syslin.solve_seidel (100, 1e-6);
    syslin.print_ans ();
    std::cout << std::endl << "fastest descend: " <<std::endl;
    syslin.solve_fastest_descend (10000, 0);
    syslin.print_ans ();
    std::cout << std::endl << "least residuals: " <<std::endl;
    syslin.solve_least_residuals (10000, 0);
    syslin.print_ans ();
    return 0;
}