#include "linearsystem.h"

int main ()
{
    linearsystem<float> syslin (3);
    syslin.scan ();
    std::cout << "syslin = " << std::endl;
    syslin.print ();
    syslin.solve_gauss ();
    syslin.print_ans ();
    return 0;
}