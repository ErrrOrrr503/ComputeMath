#include "../linear/linearsystem.h"

typedef long double floating;

int main ()
{
    size_t dim;
    std::cin >> dim;
    linearsystem<floating> syslin (dim);

    syslin.set_verboseness (0);

    syslin.scan ();

    syslin.solve_run_through ();
    syslin.print_ans ();
    return 0;
}