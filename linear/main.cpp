#include "matrix.h"

int main ()
{
    matrix<float> A (3, 3);
    A.scan (std::cin);
    std::cout << "A = " << std::endl;
    A.print (std::cout);
    std::cout << "A[1][1] = " << A[1][1] << std::endl;
    float det = A.det ();
    std::cout << "|A| = " << det << std::endl;
}