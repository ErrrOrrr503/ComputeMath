#include <iostream>

#define OMP
#define OMP_TARGET
#include "matrix.h"

#ifndef OMP
#include <omp.h>
#endif

typedef float num_t;

int main (int argc, char *argv[])
{
    size_t lines = atoll (argv[1]), colons= atoll (argv[2]);
    int num_threads = atoi (argv[3]);
    int algo = 0;
    if (argc == 5)
        algo = atoi (argv[4]);
    #ifdef OMP
    matrix<num_t> M1(lines, colons, num_threads);
    matrix<num_t> M2(colons, lines, num_threads);
    #else
    matrix<num_t> M1(lines, colons);
    matrix<num_t> M2(colons, lines);
    #endif
    M1.random_generate ();
    M2.random_generate ();
    double t_start, time_mul;

    if (!algo || (algo & 1)) {
        t_start = omp_get_wtime ();
        matrix<num_t> M (plain_mul (M1, M2));
        time_mul = omp_get_wtime () - t_start;
        std::cout << time_mul << std::endl;
    }

    if (!algo || (algo & 2)) {
        t_start = omp_get_wtime ();
        matrix<num_t> N (line_mul (M1, M2));
        time_mul = omp_get_wtime () - t_start;
        std::cout << time_mul << std::endl;
    }

    if (!algo || (algo & 4)) {
        t_start = omp_get_wtime ();
        matrix<num_t> O (strassen_mul (M1, M2));
        time_mul = omp_get_wtime () - t_start;
        std::cout << time_mul << std::endl;
    }

    if (!algo || (algo & 8)) {
        t_start = omp_get_wtime ();
        matrix<num_t> O (omp_target_mul (M1, M2));
        time_mul = omp_get_wtime () - t_start;
        std::cout << time_mul << std::endl;
    }

    return 0;
}