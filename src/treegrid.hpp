//
// Created by mblomquist on 7/26/23.
//

#ifndef LDG_POISSON_TREEGRID_HPP
#define LDG_POISSON_TREEGRID_HPP

#include "uvector.hpp"
#include "multiloop.hpp"
#include "math_tools.hpp"

template<int N>
class treeGrid
{
    struct Node
    {
        int type = 1;
        algoim::uvector<double, N> xmin, xmax;
        algoim::uvector<Node, ipow(N,2)> children = nullptr;
    };

    int min_level = 0;
    int max_level = 0;

    algoim::uvector<double, N> domain_min, domain_max;


public:


    void split_node(Node n)
    {
        n.type = 0;

        for (algoim::MultiLoop<N> i(0,2); ~i; ++i)
        {

        }
    }
};

#endif //LDG_POISSON_TREEGRID_HPP
