#include "matrix.h"

int main ()
{
    matrix<float> A (3, 4);
    A.scan (std::cin);
    std::cout << "A = " << std::endl;
    A.print (std::cout);
    //float det = A.det ();
    //std::cout << "|A| = " << det << std::endl;
    A.conv_to_upper_triangle ();
    std::cout << "A.up_t = " << std::endl;
    A.print (std::cout);
}