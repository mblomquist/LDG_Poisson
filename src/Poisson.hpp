//
// Created by mblomquist on 7/19/23.
//

#ifndef LDG_POISSON_POISSON_HPP
#define LDG_POISSON_POISSON_HPP

#include <cmath>
#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "legendre.hpp"
#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include "smatrix.hpp"

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


template<int P, int N>
class PoissonSolver
{
    uniformGrid<N> grid;

    GaussQuad quad;

    algoim::uvector<algoim::uvector<int, 2>, N> boundary_conditions;

    double c1 = 1.;
    double c2 = 1.-c1;

    double tau_i = 1;
    double tau_d = 1;

    algoim::uvector<smatrix<double, ipow(P,N)>,N> D;

    std::unordered_map<int_pair,smatrix<double, ipow(P,N)>> L[N], T[N], G[N];

    std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> rhs;

public:

    PoissonSolver()
    {
        boundary_conditions = 0;
        grid.set_periodic(true);
    }

    void set_domain(algoim::uvector<double, N> domain_min_,
                    algoim::uvector<double, N> domain_max_)
    {
        grid.set_domain_min(domain_min_);
        grid.set_domain_max(domain_max_);
    }

    void get_domain()
    {
        std::cout << "Domain min: " << grid.get_xmin() << std::endl;
        std::cout << "Domain max: " << grid.get_xmax() << std::endl;
    }

    void set_elements_per_dim(algoim::uvector<int, N> elements_per_dim_)
    {
        grid.set_elements_per_dim(elements_per_dim_);
    }

    void set_boundary_conditions( algoim::uvector<algoim::uvector<int, 2>, N> bcs)
    {
        for (int dim = 0; dim < N; ++dim) {
            boundary_conditions(dim)(0) = bcs(dim)(0);
            boundary_conditions(dim)(1) = bcs(dim)(1);

            if (bcs(dim) == 0)
            {
                grid.periodic_domain(dim, true);
            }
        }
    }

    void l2_projection(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &projected_func,
                       const std::function<double(const algoim::uvector<double,N>&)> func)
    {
        for (algoim::MultiLoop<N> i(0, grid.get_elements_per_dim()); ~i; ++i)
        {
            algoim::uvector<double, N> shift, scale;

            shift = grid.get_nodal_pos(i());
            scale = grid.get_dx();

            std::function<double(algoim::uvector<double, N> x)> xfunc = [shift, scale, func](algoim::uvector<double, N> x) {
                return func(shift + x * scale);
            };

            l2_projection_on_reference_element<P,N>(xfunc, projected_func[grid.get_element_id(i())]);
        }
    }

    void project_rhs(const std::function<double(const algoim::uvector<double,N>&)> func)
    {
        l2_projection(rhs, func);
    }

    void evaluate_basis_on_uniform_grid(std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> &basis_coeffs,
                                        algoim::uvector<int, N> pts_per_dim,
                                        const std::string& filename,
                                        const std::string& data_name)
    {

        // create a new grid for evaluation / printing
        uniformGrid<N> evalGrid(pts_per_dim-1, grid.get_xmin(), grid.get_xmax());

        // create a vector to store the evaluation data
        std::vector<double> point_data;

        point_data.resize(evalGrid.get_total_nodes(), 0.);

        algoim::uvector<double, N> eval_position, ref_pos;


        int basis_indx;

        // loop through all the points in the new grid to evaluate the basis
        for (algoim::MultiLoop<N> i(0,pts_per_dim); ~i; ++i)
        {
            eval_position = evalGrid.get_nodal_pos(i());

            basis_indx = grid.get_element_from_pos(eval_position);

            ref_pos = grid.map_pos_to_ref_element(eval_position);

            point_data[evalGrid.get_node_id(i())] = compute_basis_coefficients_at_point<P, N>(basis_coeffs[basis_indx],
                                                                                                      ref_pos);

        }

        evalGrid.print_grid_point_data_to_vtk(filename, data_name, point_data);

    }

    void construct_broken_gradient_operator()
    {
        smatrix<double, P> D_1d;

        for (int i = 0; i < P; ++i) {
            for (int j = i+1; j < P; ++j) {
                if ((j-i) % 2 != 0){
                    D_1d(i,j) = 2 * std::sqrt(2*i+1) * std::sqrt(2*j+1);
                } else {
                    D_1d(i,j) = 0.;
                }
            }
        }

        for (int dim = 0; dim < N; ++dim) {
            for (algoim::MultiLoop<N> i(0,P); ~i; ++i)
            {
                for (int j = 0; j < P; ++j) {
                    algoim::uvector<int, N> indx = i();
                    indx(dim) = j;
                    D(dim)(fold<P, N>(i()), fold<P, N>(indx)) =
                            D_1d(i(dim),j);
                }
            }
        }
    }

    void print_Dmat()
    {
        for (int dim = 0; dim < N; ++dim) {
            std::cout << "\nD_" << dim << ":" << std::endl;
            D(dim).print();
        }
    }

    void mult_D(std::function<double(algoim::uvector<double, N> x)> func,
                std::function<double(algoim::uvector<double, N> x)> soln_dx,
                std::function<double(algoim::uvector<double, N> x)> soln_dy,
                std::function<double(algoim::uvector<double, N> x)> soln_dz)
    {
        algoim::uvector<double, ipow(P,N)> func_c;
        l2_projection_on_reference_element<P,N>(func, func_c);

        algoim::uvector<double, ipow(P,N)> soln_cx, soln_cy, soln_cz;
        l2_projection_on_reference_element<P,N>(soln_dx, soln_cx);
        l2_projection_on_reference_element<P,N>(soln_dy, soln_cy);
        l2_projection_on_reference_element<P,N>(soln_dz, soln_cz);

        algoim::uvector<double, ipow(P,N)> func_dx, func_dy, func_dz;

        func_dx = matvec(D(0), func_c);
        func_dy = matvec(D(1), func_c);
        func_dz = matvec(D(2), func_c);

        std::cout << "Difference dx: " << algoim::sqrnorm(func_dx - soln_cx) << std::endl;
        std::cout << "Difference dy: " << algoim::sqrnorm(func_dy - soln_cy) << std::endl;
        std::cout << "Difference dz: " << algoim::sqrnorm(func_dz - soln_cz) << std::endl;
    }

    void print_grid(const std::string& filename)
    {
        grid.print_grid_to_vtk(filename);
    }

    void print_rhs_on_uniform_grid(algoim::uvector<int, N> pts_per_dim,
                                   const std::string& filename)
    {
        evaluate_basis_on_uniform_grid(rhs, pts_per_dim, filename, "rhs_data");
    }

    void compute_lifting_operator_on_ref_face(int dim,
                                              smatrix<double, ipow(P,N)> &A_ii,
                                              smatrix<double, ipow(P,N)> &A_ij,
                                              smatrix<double, ipow(P,N)> &A_ji,
                                              smatrix<double, ipow(P,N)> &A_jj)
    {
        A_ii = A_ij = A_ji = A_jj = 0.;

        constexpr int Q = int((2*P+1)/2)+1;

        algoim::uvector<double, ipow(P,N)> eval_i, eval_j;

        algoim::uvector<double, N> eval_pos_i, eval_pos_j;
        algoim::uvector<double, N-1> pos_Dmo;

        eval_pos_i(dim) = 1.;
        eval_pos_j(dim) = 0.;

        for (algoim::MultiLoop<N-1> i(0, Q); ~i; ++i) {
            double weight = 1.;

            for (int j = 0; j < N - 1; ++j) {
                pos_Dmo(j) = quad.x(Q, i(j));
                weight *= quad.w(Q, i(j));
            }

            int t = 0;
            for (int j = 0; j < N; ++j) {
                if (j != dim) {
                    eval_pos_i(j) = pos_Dmo(t);
                    eval_pos_j(j) = pos_Dmo(t);
                    ++t;
                }
            }

            eval_i = evaluate_basis_as_point<P, N>(eval_pos_i);
            eval_j = evaluate_basis_as_point<P, N>(eval_pos_j);

            A_ii += outer_prod(eval_i, eval_i) * weight;
            A_ij += outer_prod(eval_i, eval_j) * weight;
            A_ji += outer_prod(eval_j, eval_i) * weight;
            A_jj += outer_prod(eval_j, eval_j) * weight;

        }
    }

    void construct_lifting_operator_periodic_grid()
    {

        smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj;

        algoim::uvector<int, N> elements_per_dim = grid.get_elements_per_dim();

        for (int dim = 0; dim < N; ++dim) {

            double dv = 1;
            for (int d = 0; d < N; ++d) {
                if (d != dim)
                    dv *= grid.get_dx(d);
            }

            for (algoim::MultiLoop<N> i(0,elements_per_dim); ~i; ++i)
            {
                algoim::uvector<int, N> element_i, element_j;

                element_i(dim) = (i(dim) - 1 == -1) ? elements_per_dim(dim) - 1 : i(dim) - 1;
                element_j(dim) = i(dim);

                for (int d = 0; d < N; ++d) {
                    if (d != dim){
                        element_i(d) = i(d);
                        element_j(d) = i(d);
                    }
                }

                compute_lifting_operator_on_ref_face(dim, A_ii, A_ij, A_ji, A_jj);

                L[dim][{grid.get_element_id(element_i), grid.get_element_id(element_i)}] += (c1-1.) * A_ii * dv;
                L[dim][{grid.get_element_id(element_i), grid.get_element_id(element_j)}] += c2 * A_ij * dv;
                L[dim][{grid.get_element_id(element_j), grid.get_element_id(element_i)}] += -c1 * A_ji * dv;
                L[dim][{grid.get_element_id(element_j), grid.get_element_id(element_j)}] += (1.-c2) * A_jj * dv;

            }
        }
    }

    void construct_penalty_operator_uniform_grid()
    {

        smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj;

        algoim::uvector<int, N> elements_per_dim = grid.get_elements_per_dim();

        for (int dim = 0; dim < N; ++dim) {

            double dv = 1;
            for (int d = 0; d < N; ++d) {
                if (d != dim)
                    dv *= grid.get_dx(d);
            }

            int starting_face = (grid.periodic_domain()) ? 0 : 1;

            for (algoim::MultiLoop<N> i(starting_face,elements_per_dim); ~i; ++i)
            {
                algoim::uvector<int, N> element_i, element_j;

                element_i(dim) = (i(dim) - 1 == -1) ? elements_per_dim(dim) - 1 : i(dim) - 1;
                element_j(dim) = i(dim);

                for (int d = 0; d < N; ++d) {
                    if (d != dim){
                        element_i(d) = i(d);
                        element_j(d) = i(d);
                    }
                }

                compute_lifting_operator_on_ref_face(dim, A_ii, A_ij, A_ji, A_jj);

                T[dim][{grid.get_element_id(element_i), grid.get_element_id(element_i)}] = tau_i * A_ii * dv;
                T[dim][{grid.get_element_id(element_i), grid.get_element_id(element_j)}] = tau_i * A_ij * dv;
                T[dim][{grid.get_element_id(element_j), grid.get_element_id(element_i)}] = tau_i * A_ji * dv;
                T[dim][{grid.get_element_id(element_j), grid.get_element_id(element_j)}] = tau_i * A_jj * dv;

            }

            // add in dirichlet bc's
        }
    }

    void construct_gradient_operator()
    {
        construct_broken_gradient_operator();
        construct_lifting_operator_periodic_grid();

        // loop through ever dimension
        for (int dim = 0; dim < N; ++dim) {
            // loop through every element
            for (algoim::MultiLoop<N> i(0,grid.get_elements_per_dim()); ~i; ++i)
            {
                G[dim][{grid.get_element_id(i()), grid.get_element_id(i())}] += L[dim][{grid.get_element_id(i()), grid.get_element_id(i())}] + D(dim) * grid.get_dx(dim);

                for (int d = 0; d < N; ++d) {
                    for (int dir = 0; dir < 2; ++dir) {
                        if (i(dim) - 1 > -1 || grid.is_periodic(dim) == true){
                            algoim::uvector<int, N> nbr = i();
                            if (i(dim)-1 == -1)
                                nbr(dim) = grid.get_elements_per_dim()(dim)-1;
                            else
                                --nbr(dim);
                            G[dim][{grid.get_element_id(i()),grid.get_element_id(nbr)}] += L[dim][{grid.get_element_id(i()),grid.get_element_id(nbr)}];
                        }

                        if (i(dim)+1 < grid.get_elements_per_dim()(dim) || grid.is_periodic(dim) == true){
                            algoim::uvector<int, N> nbr = i();
                            if (i(dim)+1 == grid.get_elements_per_dim()(dim))
                                nbr(dim) = 0;
                            else
                                ++nbr(dim);

                            G[dim][{grid.get_element_id(i()),grid.get_element_id(nbr)}] += L[dim][{grid.get_element_id(i()),grid.get_element_id(nbr)}];
                        }
                    }
                }
            }
        }
    }

    void print_lifting_operator()
    {
        {
            for (int dim = 0; dim < N; ++dim) {
                for (int i = 0; i < grid.get_total_elements(); ++i) {
                    for (int j = 0; j < grid.get_total_elements(); ++j) {
                        std::cout << "\nL_" << dim << "(" << i << "," << j << "):" << std::endl;
                        L[dim][{i, j}].print();
                    }
                }
            }
        }
    }

    void print_gradient_operator()
    {
        for (int dim = 0; dim < N; ++dim) {
            for (int i = 0; i < grid.get_total_elements(); ++i) {
                for (int j = 0; j < grid.get_total_elements(); ++j) {
                    std::cout << "\nG_" << dim << "(" << i << "," << j << "):" << std::endl;
                    G[dim][{i, j}].print();
                }
            }
        }
    }
};




#endif //LDG_POISSON_POISSON_HPP
