#include "linearsystem.h"

int main ()
{
    linearsystem<float> syslin (3);
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
    return 0;
}