//
// Created by mblomquist on 7/7/23.
//
#pragma once

#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include <cmath>
#include <vector>
#include <functional>

template<int P, int N>
int fold(const algoim::uvector<int, N> &indx_t)
{
    int indx = 0;

    for (int dim = 0; dim < N; ++dim)
    {
        indx += indx_t(dim) * ipow(P, dim);
    }

    return indx;
}

template<int P, int N>
algoim::uvector<int, N> unfold(int indx)
{
    algoim::uvector<int, N> indx_t = 0;

    indx_t(N-1) = indx / ipow(P,N-1);

    for (int dim = N-2; dim > 0; --dim) {
        int temp = indx;
        for (int i = dim+1; i < N; ++i) {
            temp -= indx_t(i) * ipow(P,i);
        }
        indx_t(dim) = temp / ipow(P,dim);
    }

    indx_t(0) = indx % P;

    return indx_t;
}

template<int P>
void evaluate_legendre_at_point(const double x, algoim::uvector<double, P> &px)
{
    px(0) = 1.;

    if constexpr (P >= 1)
        px(1) = x;

    for (int i = 1; i < P-1; ++i) {
        px(i+1) = ((2.*i+1.)*x*px(i) - i*px(i-1))/(i+1.);
    }
}

template<int P>
void evaluate_shifted_legendre(const double x, algoim::uvector<double, P> &px)
{
    evaluate_legendre_at_point<P>((2.*x-1.), px);

    // normalize
    for (int i = 0; i < P; ++i) {
        px(i) *= std::sqrt(2.*i+1.);
    }
}

template<int P, int N>
void compute_basis_at_point(algoim::uvector<double, ipow(P,N)> &basis_Nd,
                            algoim::uvector<double, N> &pt)
{
    algoim::uvector<algoim::uvector<double,P>, N> basis_1d = 0.;

    for (int dim = 0; dim < N; ++dim) {
        evaluate_shifted_legendre<P>(pt(dim), basis_1d(dim));
    }

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        double prod = 1.;
        for (int dim = 0; dim < N; ++dim) {
            prod *= basis_1d(dim)(i(dim));
        }
        basis_Nd(fold<P, N>(i())) = prod;
    }

}

template<int P, int N>
void l2_projection_on_reference_element(const std::function<double(const algoim::uvector<double,N>&)> func,
                                 algoim::uvector<double, ipow(P,N)> &p_fun)
{
    constexpr int Q = P;

    // create a stacked variable of the basis at multiple quadrature points
    algoim::uvector<double, ipow(P,N)> basis;

    algoim::uvector<double, N> pos;
    double weights;

    for (algoim::MultiLoop<N> i(0, Q); ~i; ++i)
    {
        weights = 1.;

        for (int dim = 0; dim < N; ++dim) {
            weights *= GaussQuad::w(Q,i(dim));
            pos(dim) = GaussQuad::x(Q,i(dim));
        }

        compute_basis_at_point<P,N>(basis, pos);

        p_fun += basis * weights * func(pos);

    }

}

template<int P, int N>
algoim::uvector<double, ipow(P,N)> evaluate_basis_as_point(algoim::uvector<double, N> pt)
{
    algoim::uvector<algoim::uvector<double,P>, N> basis = 0.;
    algoim::uvector<double, ipow(P,N)> multi_basis;

    for (int dim = 0; dim < N; ++dim) {
        evaluate_shifted_legendre<P>(pt(dim), basis(dim));
    }

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        double prod = 1.;
        for (int dim = 0; dim < N; ++dim) {
            prod *= basis(dim)(i(dim));
        }
        multi_basis(fold<P,N>(i()))= prod;
    }

    return multi_basis;
}

template<int P, int N>
double compute_basis_coefficients_at_point(algoim::uvector<double, ipow(P, N)> coeff,
                                           algoim::uvector<double, N> pt)
{
    algoim::uvector<double, ipow(P,N)> basis = evaluate_basis_as_point<P,N>(pt);

    double sum = 0.;

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        sum += coeff(fold<P, N>(i())) * basis(fold<P,N>(i()));
    }

    return sum;
}

