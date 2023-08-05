//
// Created by mblomquist on 8/3/23.
//

#ifndef LDG_POISSON_MULTIGRID_HPP
#define LDG_POISSON_MULTIGRID_HPP

#include "uvector.hpp"
#include "multiloop.hpp"
#include "smatrix.hpp"
#include "legendre.hpp"
#include "math_tools.hpp"

template <int P, int N>
smatrix<double, ipow(P,N)> transform(algoim::uvector<algoim::uvector<double, N>, 2> rect_s,
                                                   algoim::uvector<algoim::uvector<double, N>, 2> rect_d)
{
    GaussQuad quad;

    smatrix<double, ipow(P,N)> C;

    algoim::uvector<double, ipow(P,N)> basis_s, basis_d;
    algoim::uvector<double, N> pos_d, pos_s;

    double weights;

    for (algoim::MultiLoop<N> i(0, P); ~i; ++i)
    {
        weights = 1.;

        for (int dim = 0; dim < N; ++dim) {
            weights *= quad.w(P,i(dim));
            pos_d(dim) = quad.x(P, i(dim));
        }

        for (int dim = 0; dim < N; ++dim) {
            pos_s(dim) = (pos_d(dim) * (rect_s(1)(dim) - rect_s(0)(dim)) + rect_s(0)(dim) - rect_d(0)(dim)) / (rect_d(1)(dim) - rect_d(0)(dim));
        }

        compute_basis_at_point<P,N>(basis_s, pos_s);
        compute_basis_at_point<P,N>(basis_d, pos_d);

        C += weights * outer_prod<ipow(P,N)>(basis_d, basis_s);
    }
    return C;
}





template<int P, int N>
void test_transform()
{
    algoim::uvector<algoim::uvector<double, N>, 2> source, dest;
    smatrix<double, ipow(P,N)> C_sd;

    source(0) = -1.;
    source(1) = 1.;

    dest(0) = 1.0;
    dest(1) = 3.0;

    C_sd = transform<P,N>(source, dest);

    constexpr int num_pts = 10;

    algoim::uvector<algoim::uvector<double, N>, num_pts> rnd_pts, x_s, x_d;

    for (int i = 0; i < num_pts; ++i) {
        for (int dim = 0; dim < N; ++dim) {
            rnd_pts(i)(dim) = static_cast<double>(std::rand()) / RAND_MAX;
        }
    }

    for (int i = 0; i < num_pts; ++i) {
        for (int dim = 0; dim < N; ++dim) {
            x_s(i)(dim) = rnd_pts(i)(dim) * (source(1)(dim) - source(0)(dim)) + source(0)(dim);
            x_d(i)(dim) = rnd_pts(i)(dim) * (dest(1)(dim)   - dest(0)(dim))   + dest(0)(dim);
        }
    }

    algoim::uvector<double, ipow(P,N)> u, Cu;
    algoim::uvector<double, num_pts> result0, result1;

    for (int i = 0; i < ipow(P,N); ++i) {
        u(i) = static_cast<double>(std::rand()) / RAND_MAX;
    }

    Cu = matvec<ipow(P,N)>(C_sd, u);

    std::cout << "\n--- Random Points ---" << std::endl;
    for (int i = 0; i < num_pts; ++i) {
        std::cout << "pt " << i << ": " << rnd_pts(i) << ", " << x_s(i) << ", " << x_d(i) << std::endl;
    }

    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "source rectangle: " << source << std::endl;
    std::cout << "destination rectangle: " << dest << std::endl;
    for (int i = 0; i < num_pts; ++i) {
        result0(i) = compute_basis_coefficients_at_point<P,N>(u,x_s(i));
        result1(i) = compute_basis_coefficients_at_point<P,N>(Cu,x_d(i));

        std::cout << "pt " << i << ": " << std::abs(result0(i) - result1(i)) << std::endl;
    }
}

#endif //LDG_POISSON_MULTIGRID_HPP
