//
// Created by mblomquist on 7/10/23.
//

#ifndef DG_UTILS_MATH_TOOLS_HPP
#define DG_UTILS_MATH_TOOLS_HPP




constexpr int ipow(int base, int exponent)
{
    if (exponent > 0)
        return ipow(base, exponent - 1) * base;
    else
        return 1;
}


#endif //DG_UTILS_MATH_TOOLS_HPP