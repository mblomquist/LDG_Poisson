#pragma once

#include "../../eigen/Eigen/Dense"

#include <iomanip>

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

    smatrix<T, M, N> operator+ (const smatrix& x) const
    {
        smatrix<T, M, N> result;
        for (int i = 0; i < M*N; ++i) {
            result.data_[i] = data_[i] + x.data_[i];
        }

        return result;
    }

    smatrix<T, M, N> operator- (const smatrix& x) const
    {
        smatrix<T, M, N> result;
        for (int i = 0; i < M*N; ++i) {
            result.data_[i] = data_[i] - x.data_[i];
        }

        return result;
    }

    smatrix& operator+= (const smatrix& x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] += x.data_[i];
        }

        return *this;
    }

    smatrix& operator-= (const smatrix& x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] -= x.data_[i];
        }

        return *this;
    }

    smatrix<double, M, N> operator* (const double x) const
    {
        smatrix<double, M, N> result;
        for (int i = 0; i < M*N; ++i) {
                result.data_[i] = data_[i] * x;
        }

        return result;
    }

    smatrix<double, M, N> operator/ (const double x) const
    {
        smatrix<double, M, N> result;
        for (int i = 0; i < M*N; ++i) {
            result.data_[i] = data_[i] / x;
        }

        return result;
    }

    smatrix& operator*= (const double x)
    {
        for (int i = 0; i < M*N; ++i) {
            data_[i] *= x;
        }

        return *this;
    }

    smatrix<double, N, M> transpose()
    {
        smatrix<double, N, M> At;

        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                At(j,i) = data_[i*N+j];
            }
        }
        return At;
    }

    void print()
    {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                std::cout << ((std::abs(data_[i*N+j]) > 1.e-12) ? data_[i*N+j] : 0.0) << ((j == N-1) ? ";" : ", ");
            }
            std::cout << std::endl;
        }
    }

};

template<int M, int N>
smatrix<double, M, N> operator* (const double x, const smatrix<double, M, N> &A)
{
    smatrix<double, M, N> result;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            result(i,j) = x * A(i,j);
        }
    }
    return result;
}

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

template<int N>
smatrix<double, N> outer_prod(const algoim::uvector<double,N> &a, const algoim::uvector<double,N> &b)
{
    smatrix<double, N> prod;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            prod(i,j) = a(i) * b(j);
        }
    }

    return prod;
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

template<int M>
void point_Gauss_Seidel(smatrix<double, M> &A,
                        algoim::uvector<double, M> &x,
                        algoim::uvector<double, M> &b,
                        int max_itrs = 20,
                        double tol = 1.0e-12,
                        double omega = 1.)
{
    if (omega == 1.)
        std::cout << "\n --- Gauss Seidel Iteration --- \n" << std::endl;
    else
        std::cout << "\n --- SOR Iteration --- \n" << std::endl;

    algoim::uvector<double, M> res;

    for (int itr = 0; itr < max_itrs; ++itr) {
        for (int i = 0; i < M; ++i) {
            double sum = 0.;

            for (int j = 0; j < M; ++j) {
                sum += A(i,j) * x(j);
            }

            x(i) = omega * (b(i) - sum + A(i,i) * x(i)) / A(i,i) + (1 - omega) * x(i);
        }

        res = b - matvec(A,x);

        if (norm(res) < tol)
        {
            std::cout << "\nIterations converged at iteration " << itr << " with a residual of " << norm(res) << std::endl;
            return;
        }
        else
            std::cout << "Iteration " << itr << ", norm(res): " << norm(res) << std::endl;
    }
}

template<int M>
void SOR(smatrix<double, M> &A,
         algoim::uvector<double, M> &x,
         algoim::uvector<double, M> &b,
         double omega,
         int max_itrs = 20,
         double tol = 1.0e-12)
{
    point_Gauss_Seidel(A, x, b, tol, max_itrs, omega);
}


template<int M>
smatrix<double, M> pseudo_inverse_with_Eigen(const smatrix<double, M> &A_in)
{

    smatrix<double, M> u, ut, B;
    algoim::uvector<double, M> sigma;

    // create a copy of A_in for eigen
    Eigen::MatrixXd A(M,M);

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            A(i,j) = A_in(i,j);
        }
    }

    // use eigen to compute the svd of A
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd S = svd.singularValues();
    Eigen::MatrixXd V = svd.matrixV();

    // copy eigen svd values to smatrix
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            u(i,j) = U(i,j);
            ut(i,j) = U(j,i);
        }
        sigma(i) = S(i);
    }

    // compute the pseudoinverse of A_in
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (abs(sigma(i)) > max(sigma)*1.e-13)
                ut(i,j) /= sigma(i);
            else
                ut(i,j) *= 0.;
        }
    }

    B = matmat(u, ut);

    return B;
}
