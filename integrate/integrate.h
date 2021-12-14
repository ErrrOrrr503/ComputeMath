#include <iostream>
#include <vector>
#include <cmath>
#include <string>

enum left_right {
    improper_right,
    improper_left
};

template <typename num_t>
class netfunc {
public:
    num_t (*ref_func) (num_t);
    num_t (*eps_func) (num_t);
    size_t netsize;
    std::vector<num_t> x_data;
    std::vector<num_t> f_data;
    num_t h;
    num_t x_max;
    num_t x_min;
    bool has_reference;
    left_right end;
    bool is_improper;

    netfunc () = delete;
    netfunc (num_t (*reference_func) (num_t), num_t min, num_t max, size_t netsz) : 
        ref_func (reference_func), 
        eps_func (NULL),
        x_max (max), x_min (min), 
        has_reference (1),
        is_improper (0)
    {
        build_net (netsz);
    }
    netfunc (num_t (*reference_func) (num_t), num_t (*epsilon_func) (num_t), num_t min, num_t max, left_right feature_end) : 
        ref_func (reference_func), 
        eps_func (epsilon_func), 
        x_max (max), x_min (min), 
        has_reference (1), 
        end (feature_end),
        is_improper (1)
    {
        if ((isnanl (x_max) && end != improper_right) || (isnanl (x_min) && end != improper_left))
            throw std::runtime_error ("proper end must correspond to nan");
    }
    void build_net (size_t netsz)
    {
        if (netsz < 2)
            throw std::runtime_error ("number of points in netsize must be greater than 1");
        netsize = netsz;
        h = (x_max - x_min) / (netsz - 1);
        x_data.clear ();
        f_data.clear ();
        for (size_t i = 0; i < netsz; i++) {
            x_data.push_back (x_min + (num_t) (i) / (num_t) (netsz - 1) * (x_max - x_min));
            f_data.push_back (ref_func (x_data[i]));
        }
    }
    ~netfunc ()
    {
    }
};

template <typename num_t>
class integrator {
private:
    netfunc<num_t> func;
    bool use_runge = 0;
    size_t max_iterations = 10;
    num_t precision = 0;
    std::string errstr;

public:
    num_t ans = 0;
    num_t runge = 0;
    size_t iterations = 0;
    num_t improper_edge = 0;

    integrator () = delete;
    integrator (netfunc<num_t> f) : func (f)
    {
    }
    ~integrator ()
    {
    }

    void assume_pure_net ()
    {
        func.has_reference = 0;
        use_runge = 0;
    }

    const std::string & get_err ()
    {
        return errstr;
    }

    int enable_runge (num_t precision)
    {
        if (func.has_reference) {
            use_runge = 1;
            this->precision = precision;
            return 0;
        }
        else {
            errstr = "func is pure net, runge needs more data";
            return 1;
        }
    }

    int integrate_rectangles ()
    {
        if (func.is_improper) {
            errstr = "improper must use improper methods";
            return 1;
        }
        num_t I = 0;
        if (func.has_reference) {
            for (size_t i = 0; i < func.netsize - 1; i++)
                I += func.ref_func ((func.x_data[i] + func.x_data[i + 1]) / 2);
            I *= func.h;
            if (use_runge) {
                num_t I_2h = I;
                I = 0;
                size_t old_netsize = func.netsize;
                func.build_net (func.netsize + func.netsize - 1); // x2 number of h-s
                for (size_t i = 0; i < func.netsize - 1; i++)
                    I += func.ref_func ((func.x_data[i] + func.x_data[i + 1]) / 2);
                I *= func.h;
                for (size_t i = 0; i < max_iterations && precision < fabsl ((I - I_2h)); i++) {
                    I_2h = I;
                    I = 0;
                    func.build_net (func.netsize + func.netsize - 1);
                    for (size_t i = 0; i < func.netsize - 1; i++)
                        I += func.ref_func ((func.x_data[i] + func.x_data[i + 1]) / 2);
                    I *= func.h;
                    iterations = i + 2;
                }
                func.build_net (old_netsize);
                runge = (I - I_2h) / 3;
                I += runge;
            }
        }
        else {
            for (size_t i = 0; i < func.netsize - 1; i += 2)
                I += (func.f_data[i + 1]);
            I *= 2 * func.h;
        }
        ans = I;
        return 0;
    }

    int integrate_trapezoids ()
    {
        if (func.is_improper) {
            errstr = "improper must use improper methods";
            return 1;
        }
        num_t I = 0;
        for (size_t i = 0; i < func.netsize - 1; i++)
            I += (func.f_data[i] + func.f_data[i + 1]);
        I *= func.h / 2;
        if (func.has_reference && use_runge) {
            num_t I_2h = I;
            I = 0;
            size_t old_netsize = func.netsize;
            func.build_net (func.netsize + func.netsize - 1); // x2 number of h-s
            for (size_t i = 0; i < func.netsize - 1; i++)
                I += (func.f_data[i] + func.f_data[i + 1]);
            I *= func.h / 2;
            for (size_t i = 0; i < max_iterations && precision < fabsl ((I - I_2h)); i++) {
                I_2h = I;
                I = 0;
                func.build_net (func.netsize + func.netsize - 1);
                for (size_t i = 0; i < func.netsize - 1; i++)
                    I += (func.f_data[i] + func.f_data[i + 1]);
                I *= func.h / 2;
                iterations = i + 2;
            }
            func.build_net (old_netsize);
            runge = (I - I_2h) / 3;
            I += runge;
        }
        ans = I;
        return 0;
    }

    int integrate_simpson ()
    {
        if (func.is_improper) {
            errstr = "improper must use improper methods";
            return 1;
        }
        num_t I = 0;
        if (func.has_reference) {
            for (size_t i = 0; i < func.netsize - 1; i++)
                I += (func.f_data[i] + func.f_data[i + 1] + 4 * func.ref_func((func.x_data[i] + func.x_data[i + 1]) / 2));
            I *= func.h / 6;
            if (use_runge) {
                num_t I_2h = I;
                I = 0;
                size_t old_netsize = func.netsize;
                func.build_net (func.netsize + func.netsize - 1); // x2 number of h-s
                for (size_t i = 0; i < func.netsize - 1; i++)
                    I += (func.f_data[i] + func.f_data[i + 1] + 4 * func.ref_func((func.x_data[i] + func.x_data[i + 1]) / 2));
                I *= func.h / 6;
                for (size_t i = 0; i < max_iterations && precision < fabsl (I - I_2h); i++) {
                    I_2h = I;
                    I = 0;
                    func.build_net (func.netsize + func.netsize - 1);
                    for (size_t i = 0; i < func.netsize - 1; i++)
                        I += func.h * (func.f_data[i] + func.f_data[i + 1] + 4 * func.ref_func((func.x_data[i] + func.x_data[i + 1]) / 2));
                    I /= 6;
                    iterations = i + 2;
                }
                func.build_net (old_netsize);
                runge = (I - I_2h) / 15;
                I += runge;
            }
        }
        else {
            if (func.netsize < 3) {
                errstr = "simpson pure net requires at least 3 knots";
                return 1;
            }
            for (size_t i = 0; i < func.netsize - 2; i += 2)
                I += (func.f_data[i] + func.f_data[i + 2] + 4 * func.f_data[i + 1]);
            I *= func.h / 3;
        }
        ans = I;
        return 0;
    }

    int integrate_improper (num_t precision)
    {
        if (!func.is_improper) {
            errstr = "proper must use proper methods";
            return 1;
        }
        num_t I1 = 0, I2 = 0;
        if (func.end == improper_right) {
            if (!isnanl (func.x_max)) {
                num_t edge = func.x_min;
                do {
                    edge = (func.x_max + edge) / 2;
                    I1 = func.eps_func (func.x_max) - func.eps_func (edge);
                } while (fabsl (I1) > precision / 2);
                improper_edge = edge;
            }
            else {
                num_t edge = 1;
                do {
                    edge *= 2;
                    I1 = func.eps_func (1e8) - func.eps_func (edge);
                } while (fabsl (I1) > precision / 2);
                improper_edge = edge;
            }
            netfunc<num_t> proper_f (func.ref_func, func.x_min, improper_edge, 100);
            integrator<num_t> I2_integrator (proper_f);
            I2_integrator.enable_runge (precision / 2);
            I2_integrator.integrate_rectangles ();
            I2 = I2_integrator.ans;
        }
        else {
            if (!isnanl (func.x_min)) {
                num_t edge = func.x_max;
                do {
                    edge = (func.x_min + edge) / 2;
                    I1 = func.eps_func (edge) - func.eps_func (func.x_min);
                } while (fabsl (I1) > precision / 2);
                improper_edge = edge;
            }
            else {
                num_t edge = -1;
                do {
                    edge *= 2;
                    I1 = func.eps_func (edge) - func.eps_func (-1e8);
                } while (fabsl (I1) > precision / 2);
                improper_edge = edge;
            }
            netfunc<num_t> proper_f (func.ref_func, improper_edge, func.x_max, 100);
            integrator<num_t> I2_integrator (proper_f);
            I2_integrator.enable_runge (precision / 2);
            I2_integrator.integrate_rectangles ();
            I2 = I2_integrator.ans;
        }
        ans = I1 + I2;
        return 0;
    }
};