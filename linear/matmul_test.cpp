#include <iostream>

#define OMP
#include "matrix.h"

typedef double num_t;

int main (int argc, char *argv[])
{
    size_t lines = atoll (argv[1]), colons= atoll (argv[1]);
    int num_threads = atoi (argv[2]);
    matrix<num_t> M1(lines, colons, num_threads);
    matrix<num_t> M2(colons, lines, num_threads);
    M1.random_generate ();
    M2.random_generate ();

    double t_start = omp_get_wtime ();
    matrix<num_t> M = plain_mul (M1, M2);
    double time_mul = omp_get_wtime () - t_start;
    std::cout << time_mul << std::endl;

    t_start = omp_get_wtime ();
    matrix<num_t> N = line_mul (M1, M2);
    time_mul = omp_get_wtime () - t_start;
    std::cout << time_mul << std::endl;

    return 0;
}