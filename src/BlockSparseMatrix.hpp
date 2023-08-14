#pragma once

#include <cmath>
#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>

template<int P, int N>
using elem_vec = std::unordered_map<int,algoim::uvector<double,ipow(P,N)>>;

// Create struct to be used to hash a pair of integers
struct int_pair
{
    int i, j;

    bool operator==(const int_pair& other) const noexcept
    {
        return i == other.i && j == other.j;
    }
};

// overload std lib to hash an integer pair
namespace std
{
    template<>
    struct hash<int_pair>
    {
        std::size_t operator()(const int_pair& x) const noexcept
        {
            static_assert(sizeof(int_pair)==8);
            return std::hash<int64_t>{}(int64_t(x.i) << 32 | int64_t(x.j));
        }
    };
}

// simple BlockSparseMatrix implementation for operators defined with std::unordered_map
template<typename T>
struct BlockSparseMatrix
{
    // maps (i,j) to A_ij
    std::unordered_map<int_pair, T> A;

    // maps i to the non-zero entries
    std::unordered_map<int, std::unordered_set<int>> row;

    // maps i to the non-zero entries
    std::unordered_map<int, std::unordered_set<int>> col;

    T& operator()(int i, int j)
    {
        row[i].insert(j);
        col[j].insert(i);
        return A[{i,j}];
    }

    const T& operator()(int i, int j) const
    {
        return A.at[{i,j}];
    }
};

template<int P, int N>
void block_Gauss_Seidel(BlockSparseMatrix<smatrix<double, ipow(P, N)>> &A,
                        elem_vec<P,N> &x,
                        elem_vec<P,N> &b,
                        int num_elements,
                        int n_itr = 2)
{
    for (int itr = 0; itr < n_itr; ++itr)
    {
        for (int i = 0; i < num_elements; ++i)
        {
            algoim::uvector<double, ipow(P,N)> sum;

            for (auto j : A.row[i])
            {
                if (i != j)
                    sum += matvec(A(i,j),x[j]);
            }

            smatrix<double, ipow(P,N)> Ap = pseudo_inverse_with_Eigen(A(i, i));
            algoim::uvector<double, ipow(P,N)> result = b[i] - sum;
            x[i] = matvec(Ap,result);
        }
    }
}

template<int P, int N>
elem_vec<P, N> BlockSparse_matvec(BlockSparseMatrix<smatrix<double, ipow(P,N)>> &A,
                                  elem_vec<P,N> &x,
                                  int num_rows_b,
                                  bool A_is_transpose = false)
{
    elem_vec<P,N> b;

    for (int i = 0; i < num_rows_b; ++i) {
        if (!A_is_transpose) {
            for (auto j : A.row[i]) {
                b[i] = matvec(A(i,j), x[j]);
            }
        } else {
            for (auto j : A.col[i]) {
                b[i] = matvec(A(j,i).transpose(),x[j]);
            }
        }
    }

    return b;
}

/*
 * Computes the C = A * B where A and B are BlockSparseMatrices.
 */
template<int P, int N>
BlockSparseMatrix<smatrix<double, ipow(P,N)>> BlockSparse_matmat(BlockSparseMatrix<smatrix<double, ipow(P,N)>> &A,
                                                                               BlockSparseMatrix<smatrix<double, ipow(P,N)>> &B,
                                                                               int num_rows_C,
                                                                               bool A_is_transpose = false)
{
    BlockSparseMatrix<smatrix<double, ipow(P,N)>> C;

    for (int i = 0; i < num_rows_C; ++i) {
        if (!A_is_transpose)
        {
            for (auto k : A.row[i]) {
                for (auto j : B.row[i]) {
                    C(i,j) += matmat(A(i,k), B(k, j));
                }
            }
        } else {
            for (auto k : A.col[i]) {
                for (auto j : B.row[k]) {
                    C(i,j) += matmat(A(k,i).transpose(),B(k,j));
                }
            }
        }
    }

    return C;
}