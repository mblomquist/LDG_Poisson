//
// Created by mblomquist on 7/26/23.
//

#ifndef LDG_POISSON_TREEGRID_HPP
#define LDG_POISSON_TREEGRID_HPP

#include "uvector.hpp"

template<int N>
class treeGrid
{
    struct Node
    {
        int type = 0;
        algoim::uvector<double, N> xmin, xmax;
    };
};

#endif //LDG_POISSON_TREEGRID_HPP
