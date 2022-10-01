#include <iostream>

#define OMP
#include "matrix.h"

typedef double num_t;

int main (int argc, char *argv[])
{
    bool human_readable = 0;
    bool time = 0;
    if (argc > 1) {
        argv[1][0] == 'h' ? human_readable = 1 : human_readable = 0;
        argv[1][0] == 't' ? time = 1 : time = 0;
    }
    if (argc > 2) {
        argv[2][0] == 'h' ? human_readable = 1 : human_readable = human_readable;
        argv[2][0] == 't' ? time = 1 : time = time;
    }
    if (human_readable)
        std::cout << "Enter M1 lines colons (colons and lines for M2)" << std::endl;
    size_t lines, colons;
    std::cin >> lines >> colons;
    matrix<num_t> M1(lines, colons);
    matrix<num_t> M2(colons, lines);
    if (human_readable)
        std::cout << "Enter M1" << std::endl;
    M1.scan ();
    if (human_readable)
        std::cout << "Enter M2" << std::endl;
    M2.scan ();
    double t_start = omp_get_wtime ();
    matrix<num_t> M = line_mul (M1, M2);
    double time_mul = omp_get_wtime () - t_start;
    if (time) {
        if (human_readable)
            std::cout << "Consumed: " << time_mul << std::endl;
        else
            std::cout << time_mul << std::endl;
    }
    if (human_readable)
        std::cout << "M1 * M2 = " << std::endl;
    if (!time || human_readable)
        M.print ();
    return 0;
}