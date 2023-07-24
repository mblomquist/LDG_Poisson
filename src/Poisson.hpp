//
// Created by mblomquist on 7/19/23.
//

#ifndef LDG_POISSON_POISSON_HPP
#define LDG_POISSON_POISSON_HPP

#include <cmath>
#include <vector>
#include <functional>

#include "legendre.hpp"
#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include "smatrix.hpp"

template<int P, int N>
class PoissonSolver
{
    uniformGrid<N> grid;

    algoim::uvector<algoim::uvector<int, 2>, N> boundary_conditions;

    algoim::uvector<smatrix<double, ipow(P,N)>,N> D, L, G;

    std::unordered_map<int, algoim::uvector<double, ipow(P,N)>> rhs;

public:

    PoissonSolver()
    {
        boundary_conditions = 0;
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

            point_data[evalGrid.get_node_id(i())] = evaluate_basis_at_point<P,N>(basis_coeffs[basis_indx], ref_pos);

        }

        evalGrid.print_grid_point_data_to_vtk(filename, data_name, point_data);

    }

    void compute_D()
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
        std::cout << "\n\nPrinting D\n\n" << std::endl;
        for (int dim = 0; dim < N; ++dim) {
            std::cout << "Dim " << dim << std::endl;
            for (int i = 0; i < ipow(P, N); ++i) {
                for (int j = 0; j < ipow(P, N); ++j) {
                    std::cout << D(dim)(i,j) << " ";
                }
                std::cout << ";" << std::endl;
            }
            std::cout << std::endl;
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
};

#endif //LDG_POISSON_POISSON_HPP
