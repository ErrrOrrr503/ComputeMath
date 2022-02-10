#include "edge_task.h"

typedef long double num_t;

num_t px (num_t x)
{
    return -4;
}

num_t qx (num_t x)
{
    return 8;
}

num_t fx (num_t x)
{
    return exp (x) * (2 * sin (x) - cos (x));
}

int main (int argc, char *argv[])
{
    if (argc < 7) {
        std::cout << "usage: " << argv[0] << " <type[1 or 23]> <start> <y_start> <a, 23> <end> <y_end> <b, 23> <net_size>" << std::endl;
        return 0;
    }
    if (atoi (argv[1]) == 1) {
        if (argc != 7) {
            std::cout << "usage: " << argv[0] << " <type[1 or 23]> <start> <y_start> <a, 23> <end> <y_end> <b, 23> <net_size>" << std::endl;
            return 0;
        }
        difeq_2_order<num_t> d_2_o (px, qx, fx, (num_t) atof (argv[2]), (num_t) atof (argv[4]), (size_t) atoll (argv[6]));
        d_2_o.set_edge_1_kind ((num_t) atof (argv[3]), (num_t) atof (argv[5]));
        d_2_o.set_verboseness (0);
        d_2_o.solve_fin_dif_1_kind ();
        d_2_o.print_ans ();
    }
    if (atoi (argv[1]) == 23) {
        if (argc != 9) {
            std::cout << "usage: " << argv[0] << " <type[1 or 23]> <start> <y_start> <a, 23> <end> <y_end> <b, 23> <net_size>" << std::endl;
            return 0;
        }
        difeq_2_order<num_t> d_2_o (px, qx, fx, (num_t) atof (argv[2]), (num_t) atof (argv[5]), (size_t) atoll (argv[8]));
        d_2_o.set_edge_23_kind ((num_t) atof (argv[3]), (num_t) atof (argv[6]), (num_t) atof (argv[4]), (num_t) atof (argv[7]));
        d_2_o.set_verboseness (0);
        d_2_o.solve_fin_dif_23_kind ();
        d_2_o.print_ans ();
    }
    if (atoi (argv[1]) != 23 && atoi (argv[1]) != 1) {
        std::cout << "usage: " << argv[0] << " <type[1 or 23]> <start> <y_start> <a, 23> <end> <y_end> <b, 23> <net_size>" << std::endl;
        return 0;
    }
    return 0;
}