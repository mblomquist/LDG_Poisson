//
// Created by mblomquist on 7/7/23.
//

#ifndef DG_UTILS_LEGENDRE_HPP
#define DG_UTILS_LEGENDRE_HPP



#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include <cmath>
#include <vector>
#include <functional>

template<int P, int N>
int unfold(const algoim::uvector<int, N> indx_t)
{
    int indx_f = 0;

    for (int i = 0; i < N; ++i)
    {
        indx_f += indx_t(i) * ipow(P,i);
    }

    return indx_f;
}


template<int P>
void eval_legendre_at_point(const double x, double *px)
{
    px[0] = 1.;

    if constexpr (P >= 1)
        px[1] = x;

    for (int i = 1; i < P-1; ++i) {
        px[i+1] = ((2.*i+1.)*x*px[i] - i*px[i-1])/(i+1.);
    }
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

template<int P>
void eval_shifted_legendre_at_point(const double x, double *px)
{
    eval_legendre_at_point<P>((2.*x-1.), px);

    // normalize
    for (int i = 0; i < P; ++i) {
        px[i] *= std::sqrt(2.*i+1.);
    }
}

template<int P, int N>
double eval_multidimensional_legendre_at_point(const algoim::uvector<double, N> x, const algoim::uvector<int, N> indx)
{
    double p_multi[P][N];
    double p_1d[P];

    for (int i = 0; i < N; ++i) {

        eval_shifted_legendre_at_point<P>(x(i), p_1d);

        for (int j = 0; j < P; ++j) {
            p_multi[j][i] = p_1d[j];
        }
    }

    double prod = 1.;
    for (int i = 0; i < N; ++i) {
        prod *= p_multi[indx(i)][i];
    }

    return prod;
}


template<int P, int N>
void compute_elemental_mass_matrix(int print_me = 0)
{
    smatrix<double, ipow(P,N)> Me;

    GaussQuad quad;

    int Q = ceil((2.*P+1.)/2.);

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i) {
        for (algoim::MultiLoop<N> j(0,P); ~j; ++j) {

            double intg = 0.;

            for (algoim::MultiLoop<N> k(0, Q); ~k; ++k) {

                double w = 1.;
                algoim::uvector<double, N> x;

                for (int dim = 0; dim < N; ++dim)
                {
                    w *= quad.w(Q, k(dim));
                    x(dim) = quad.x(Q, k(dim));
                }

                double f = eval_multidimensional_legendre_at_point<P, N>(x, i());
                double g = eval_multidimensional_legendre_at_point<P, N>(x, j());

                intg += w * f * g;

            }

            Me(unfold<P,N>(i()), unfold<P,N>(j())) = intg;

        }
    }

    // Print Me
    if (print_me == 1){
        for (int i = 0; i < ipow(P,N); ++i) {
            for (int j = 0; j < ipow(P,N); ++j) {

                if (std::abs(Me(i,j)) < 1.0e-12)
                    std::cout << "0 ";
                else
                    std::cout << Me(i,j) << " ";
            }
            std::cout << "; " << std::endl;
        }
    }

}


template<int P, int N>
void project_cf_on_ref_element(const std::function<double(const algoim::uvector<double,N>&)> func,
                               algoim::uvector<double, ipow(P,N)> &p_fun,
                               algoim::uvector<double, N> element_mins = 0.,
                               algoim::uvector<double, N> dx = 1.)
{

    GaussQuad quad;

    int Q = ceil((2.*P+1.)/2.);

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i) {

        double intg = 0.;

        for (algoim::MultiLoop<N> k(0, Q); ~k; ++k) {

            double w = 1.;
            algoim::uvector<double, N> x;

            for (int dim = 0; dim < N; ++dim)
            {
                w *= quad.w(Q, k(dim));
                x(dim) = quad.x(Q, k(dim));
            }

            double f = eval_multidimensional_legendre_at_point<P, N>(x, i());

            algoim::uvector<double, N> global_x;
            for (int dim = 0; dim < N; ++dim) {
                global_x(dim) = x(dim) * dx(dim) + element_mins(dim);
            }

            double g = func(global_x);

            intg += w * f * g;

        }

        p_fun(unfold<P,N>(i())) = intg;

    }
}



template<int P, int N>
void project_cf(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &projected_func,
                const std::function<double(const algoim::uvector<double,N>&)> func,
                uniformGrid<N> &grid)
{
    algoim::uvector<int, N> k_per_dim = grid.get_elements_per_dim();
    algoim::uvector<double, N> k_min, k_dx;
    algoim::uvector<double, ipow(P,N)> local_p_fun;

    for (algoim::MultiLoop<N> i(0, k_per_dim); ~i; ++i)
    {
        k_min = grid.get_nodal_pos(i());
        k_dx = grid.get_dx();

        project_cf_on_ref_element<P, N>(func, local_p_fun, k_min, k_dx);

        projected_func[grid.get_element_id(i())] = local_p_fun;
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
        basis_Nd(unfold<P,N>(i())) = prod;
    }

}

template<int P, int N>
void project_func_on_ref_element(const std::function<double(const algoim::uvector<double,N>&)> func,
                                 algoim::uvector<double, ipow(P,N)> &p_fun)
{
    GaussQuad quad;
    constexpr int Q = int((2.*P+1.)/2.)+1;

    // create a stacked variable of the basis at multiple quadrature points
    algoim::uvector<algoim::uvector<double, ipow(P,N)>, ipow(Q,N)> fat_basis;

    algoim::uvector<double, N> pos;
    double weights;

    for (algoim::MultiLoop<N> i(0, Q); ~i; ++i)
    {
        weights = 1.;

        for (int dim = 0; dim < N; ++dim) {
            weights *= quad.w(Q,i(dim));
            pos(dim) = quad.x(Q,i(dim));
        }

        compute_basis_at_point<P,N>(fat_basis(unfold<Q,N>(i())), pos);

        p_fun += fat_basis(unfold<Q,N>(i())) * weights * func(pos);

//        std::cout << i() << " weight: " << weights << " func: " << func(pos) << std::endl;
//        std::cout << i() << " fat_basis: \n   " << fat_basis(unfold<Q,N>(i())) << std::endl;
//        std::cout << i() << "          : \n   " << fat_basis(unfold<Q,N>(i())) * weights * func(pos) << std::endl;
//        std::cout << i() << " p_fun: \n   " << p_fun << std::endl;
//        std::cout << "\n\n\n --- \n\n\n" << std::endl;

    }

}

template<int P, int N>
void project_func(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &projected_func,
                  const std::function<double(const algoim::uvector<double,N>&)> func,
                  uniformGrid<N> &grid)
{
    for (algoim::MultiLoop<N> i(0, grid.get_elements_per_dim()); ~i; ++i)
    {
        algoim::uvector<double, N> shift, scale;

        shift = grid.get_nodal_pos(i());
        scale = grid.get_dx();

        std::function<double(algoim::uvector<double, N> x)> xfunc = [shift, scale, func](algoim::uvector<double, N> x) {
            return func(shift + x * scale);
        };

        project_func_on_ref_element<P,N>(xfunc, projected_func[grid.get_element_id(i())]);
    }
}

template<int P, int N>
double evaluate_basis_at_point(algoim::uvector<double, ipow(P,N)> coeff,
                               algoim::uvector<double, N> pt)
{
    algoim::uvector<algoim::uvector<double,P>, N> basis = 0.;

    for (int dim = 0; dim < N; ++dim) {
        evaluate_shifted_legendre<P>(pt(dim), basis(dim));
    }

    double sum = 0.;

    for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
    {
        double prod = 1.;
        for (int dim = 0; dim < N; ++dim) {
            prod *= basis(dim)(i(dim));
        }
        sum += coeff(unfold<P,N>(i())) * prod;
    }

    return sum;
}

template<int P, int N>
void evaluate_basis_on_uniform_grid(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &basis_coeffs,
                                    uniformGrid<N> oldGrid,
                                    algoim::uvector<int, N> pts_per_dim,
                                    const std::string& filename,
                                    const std::string& data_name)
{

    // create a new grid for evaluation / printing
    uniformGrid<N> evalGrid(pts_per_dim-1, oldGrid.get_xmin(), oldGrid.get_xmax());

    // create a vector to store the evaluation data
    std::vector<double> point_data;

    point_data.resize(evalGrid.get_total_nodes(), 0.);

    algoim::uvector<double, N> eval_position, ref_pos;


    int basis_indx;

    // loop through all the points in the new grid to evaluate the basis
    for (algoim::MultiLoop<N> i(0,pts_per_dim); ~i; ++i)
    {
        eval_position = evalGrid.get_nodal_pos(i());

        basis_indx = oldGrid.get_element_from_pos(eval_position);

        ref_pos = oldGrid.map_pos_to_ref_element(eval_position);

        point_data[evalGrid.get_node_id(i())] = evaluate_basis_at_point<P,N>(basis_coeffs[basis_indx], ref_pos);

    }

    evalGrid.print_grid_point_data_to_vtk(filename, data_name, point_data);

}


#endif //DG_UTILS_LEGENDRE_HPP