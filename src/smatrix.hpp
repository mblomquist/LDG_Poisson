//
// Created by mblomquist on 7/7/23.
//

#ifndef DG_UTILS_SMATRIX_HPP
#define DG_UTILS_SMATRIX_HPP



#include "uvector.hpp"

template<typename T, int M, int N = M>
class smatrix
{
    double data_[M*N];

public:

    smatrix()
    {
        *this = 0.;
    }

    T& operator() (int i, int j) {return data_[i*N+j];}
    const T& operator() (int i, int j) const {return data_[i*N+j];}

    smatrix& operator= (const double x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] = x;
        }

        return *this;
    }

    smatrix& operator= (const smatrix& x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] = x.data_[i];
        }

        return *this;
    }

    smatrix& operator* (const double x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] *= x;
        }

        return *this;
    }



};

template<int M, int N, int O>
smatrix<double, M, O> matmat(const smatrix<double, M, N> &A, const smatrix<double, N, O> &B)
{
    smatrix<double, M, O> C;

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < O; ++j) {
             for (int k = 0; k < N; ++k) {
                 C(i,j) += A(i,k) * B(k,j);
             }
        }
    }

    return C;
}

template<int M, int N>
algoim::uvector<double,M> matvec(const smatrix<double, M, N> &A, const algoim::uvector<double,N> &x)
{
    algoim::uvector<double,M> b;

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            b(i) += A(i,j) * x(j);
        }
    }

    return b;
}

template<int M, int N, int O>
smatrix<double, M*O, N*O> kron_A_eyeO(const smatrix<double, M, N> &A)
{
    smatrix<double, M*O, N*O> K;

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < O; ++k) {
                K(i*O+k,j*O+k) = A(i,j);
            }
        }
    }

    return K;
}

#endif //DG_UTILS_SMATRIX_HPP