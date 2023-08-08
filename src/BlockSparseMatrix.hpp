//
// Created by mblomquist on 8/8/23.
//

#ifndef LDG_POISSON_BLOCKSPARSEMATRIX_HPP
#define LDG_POISSON_BLOCKSPARSEMATRIX_HPP

#include <cmath>
#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>

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
void block_Gauss_Seidel(BlockSparseMatrix<smatrix<double, ipow(P,N)>> &A,
                        std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> x,
                        std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> b)
{

}

#endif //LDG_POISSON_BLOCKSPARSEMATRIX_HPP
