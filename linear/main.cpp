#include "linearsystem.h"

typedef long double floating;

/*
todo: неваязка
прогонка
*/

int main ()
{
    size_t dim;
    std::cin >> dim;
    linearsystem<floating> syslin (dim);
    floating prec_machinery = syslin.calculate_machinery_precision ();
    std::cout.setf (std::ios::scientific);
    std::cout << "Machinery precision = " << prec_machinery << std::endl;
    std::cout.unsetf (std::ios::scientific);
    syslin.set_verboseness (1);
    syslin.scan ();
    std::cout << "syslin = " << std::endl;
    syslin.print ();
    //gauss
    std::cout << std::endl << "Gauss" << std:: endl;
    syslin.solve_gauss ();
    syslin.print_ans ();
    //iteration
    if (!syslin.check_jacobi ())
        std::cout << std::endl << "Hadamard condition is not met! Discrepancy is likely!" << std::endl;
    else
        std::cout << std::endl << "Hadamard condition is met! Discrepancy is highly unlikely." << std::endl;
    std::cout << "jacobi: " <<std::endl;
    syslin.solve_jacobi (1000, 5 * prec_machinery);
    syslin.print_ans ();
    std::cout << std::endl << "Seidel: " <<std::endl;
    syslin.solve_seidel (1000, 5 * prec_machinery);
    syslin.print_ans ();
    //variation
    if (!syslin.check_variation_methods ())
        std::cout << std::endl << "Enough condition for variation methods is not met! Discrepancy is likely!" << std::endl;
    else
        std::cout << std::endl << "Enough condition for variation methods is met! Discrepancy is highly unlikely." << std::endl;
    std::cout << "fastest descend: " <<std::endl;
    syslin.solve_fastest_descend (10000, 20 * prec_machinery);
    syslin.print_ans ();
    std::cout << std::endl << "least residuals: " <<std::endl;
    syslin.solve_least_residuals (10000, 20 * prec_machinery);
    syslin.print_ans ();
    std::cout << std::endl << "Run-through: " <<std::endl;
    if (!syslin.check_run_through ())
        std::cout << "matrix is not 3-diag, aborting" << std::endl;
    else {
        syslin.solve_run_through ();
        syslin.print_ans ();
    }
    return 0;
}