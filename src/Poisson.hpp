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
#include <iomanip>

#include "legendre.hpp"
#include "uvector.hpp"
#include "multiloop.hpp"
#include "quadrature.hpp"
#include "math_tools.hpp"
#include "grid.hpp"
#include "smatrix.hpp"
#include "BlockSparseMatrix.hpp"

template<int P, int N>
class PoissonSolver
{
    uniformGrid<N> grid;

    algoim::uvector<algoim::uvector<int, 2>, N> boundary_conditions;

    double c1 = 0.;
    double c2 = 1.-c1;

    double tau_i = 1;
    double tau_d = 1;

    algoim::uvector<smatrix<double, ipow(P,N)>,N> D;

    BlockSparseMatrix<smatrix<double, ipow(P,N)>> L[N], T[N], G[N];

    elem_vec<P,N> rhs, sol;

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
                       const std::function<double(const algoim::uvector<double,N>&)> &func)
    {
        for (algoim::MultiLoop<N> i(0, grid.get_elements_per_dim()); ~i; ++i)
        {
            algoim::uvector<double, N> shift, scale;

            shift = grid.get_nodal_pos(i());
            scale = grid.get_dx();

            auto xfunc = [shift, scale, func](algoim::uvector<double, N> x) {
                return func(shift + x * scale);
            };

            l2_projection_on_reference_element<P,N>(xfunc, projected_func[grid.get_element_id(i())]);
        }
    }

    void project_rhs(const std::function<double(const algoim::uvector<double,N>&)> &func)
    {
        l2_projection(rhs, func);
    }

    void project_sol(const std::function<double(const algoim::uvector<double,N>&)> &func)
    {
        l2_projection(sol, func);
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
        A_ii = 0.;
        A_ij = 0.;
        A_ji = 0.;
        A_jj = 0.;

        constexpr int Q = P;

        algoim::uvector<double, ipow(P,N)> eval_i, eval_j;

        algoim::uvector<double, N> eval_pos_i, eval_pos_j;
        algoim::uvector<double, N-1> pos_Dmo;

        eval_pos_i(dim) = 1.;
        eval_pos_j(dim) = 0.;

        for (algoim::MultiLoop<N-1> i(0, Q); ~i; ++i) {
            double weight = 1.;

            for (int j = 0; j < N - 1; ++j) {
                pos_Dmo(j) = GaussQuad::x(Q, i(j));
                weight *= GaussQuad::w(Q, i(j));
            }

            int t = 0;
            for (int j = 0; j < N; ++j) {
                if (j != dim) {
                    eval_pos_i(j) = pos_Dmo(t);
                    eval_pos_j(j) = pos_Dmo(t);
                    ++t;
                }
            }

//            // look here... MATT
//            for (int j = 0; j < N; ++j)
//                eval_pos_i(j) = (j == dim) ? 1 : pos_Dmo( (j < dim) ? j : j - 1 );

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
            for (algoim::MultiLoop<N> i(0,elements_per_dim); ~i; ++i) {
                algoim::uvector<int, N> element_i, element_j;

                element_i(dim) = (i(dim) - 1 == -1) ? elements_per_dim(dim) - 1 : i(dim) - 1;
                element_j(dim) = i(dim);

                for (int d = 0; d < N; ++d) {
                    if (d != dim) {
                        element_i(d) = i(d);
                        element_j(d) = i(d);
                    }
                }

                compute_lifting_operator_on_ref_face(dim, A_ii, A_ij, A_ji, A_jj);

                L[dim](grid.get_element_id(element_i), grid.get_element_id(element_i)) += (c1-1.) * A_ii / grid.get_dx(dim);
                L[dim](grid.get_element_id(element_i), grid.get_element_id(element_j)) +=      c2 * A_ij / grid.get_dx(dim);
                L[dim](grid.get_element_id(element_j), grid.get_element_id(element_i)) +=     -c1 * A_ji / grid.get_dx(dim);
                L[dim](grid.get_element_id(element_j), grid.get_element_id(element_j)) += (1.-c2) * A_jj / grid.get_dx(dim);

                // construct penalty while we are here
                T[dim](grid.get_element_id(element_i), grid.get_element_id(element_i)) += tau_i * A_ii / grid.get_dx(dim);
                T[dim](grid.get_element_id(element_i), grid.get_element_id(element_j)) += tau_i * A_ij / grid.get_dx(dim);
                T[dim](grid.get_element_id(element_j), grid.get_element_id(element_i)) += tau_i * A_ji / grid.get_dx(dim);
                T[dim](grid.get_element_id(element_j), grid.get_element_id(element_j)) += tau_i * A_jj / grid.get_dx(dim);

            }
        }
    }

    void construct_penalty_operator()
    {
        // maybe we only construct the penalty for specific BCs.
    }

    void construct_gradient_operator()
    {
        construct_broken_gradient_operator();
        construct_lifting_operator_periodic_grid();

        // loop through ever dimension
        for (int dim = 0; dim < N; ++dim) {
            for (int i = 0; i < grid.get_total_elements(); ++i)
            {
                G[dim](i,i) += D(dim) / grid.get_dx()(dim);

                for (auto j : L[dim].row[i])
                {
                    G[dim](i,j) += L[dim](i,j);
                }
            }
        }
    }

    void print_operators_to_file(const std::string& filename) const
    {
        std::ofstream file(filename);

        // print a header
        for (int dim = 0; dim < N; ++dim) {
            file << "D_" << dim << ",";
        }

        for (int dim = 0; dim < N; ++dim) {
            file << "L_" << dim << ",";
        }

        for (int dim = 0; dim < N; ++dim) {
            file << "G_" << dim << ",";
        }

        for (int dim = 0; dim < N; ++dim) {
            file << "T_" << dim << ((dim == N-1) ? "" : ",");
        }
        file << std::endl;

        for (int i = 0; i < grid.get_total_elements()*ipow(P,N); ++i) {
            for (int j = 0; j < grid.get_total_elements()*ipow(P,N); ++j) {

                int e_i = int(i/ipow(P,N));
                int e_j = int(j/ipow(P,N));

                int s_i = i % ipow(P,N);
                int s_j = j % ipow(P,N);

                for (int dim = 0; dim < N; ++dim) {
                    file << std::setprecision(16) << D(dim)(s_i,s_j)/grid.get_dx(dim) << ",";
                }

                for (int dim = 0; dim < N; ++dim) {
                    file << std::setprecision(16) << L[dim](e_i,e_j)(s_i,s_j) << ",";
                }

                for (int dim = 0; dim < N; ++dim) {
                    file << std::setprecision(16) << G[dim](e_i,e_j)(s_i,s_j) << ",";
                }

                for (int dim = 0; dim < N; ++dim) {
                    file << std::setprecision(16) << T[dim](e_i,e_j)(s_i,s_j) << ((dim == N-1) ? "" : ",");
                }

                file << std::endl;
            }
        }
    }

    void print_1d_gradient_to_file(const std::string& filename, int dim) const
    {
        std::ofstream file(filename);
        // print a header
        file << "i," << "j," << "G" << std::endl;
        for (int i = 0; i < grid.get_total_elements(); ++i)
        {
            for (auto j : L[dim].row[i])
            {
                for (int k = 0; k < ipow(P,N); ++k) {
                    for (int l = 0; l < ipow(P,N); ++l) {
                        file << i * ipow(P, N) + k << "," << j * ipow(P, N) + l << "," << G[dim](i, j)(k, l) << std::endl;
                    }
                }
            }
        }
        file.close();
    }

    void print_3d_gradient_to_file(const std::string& filename0, const std::string& filename1, const std::string& filename2) const
    {
        print_1d_gradient_to_file(filename0, 0);
        print_1d_gradient_to_file(filename1, 1);
        print_1d_gradient_to_file(filename2, 2);
    }

    void print_vectors_to_file(const std::string& filename) const
    {
        std::ofstream file(filename);

        // print header
        file << "mass" << "," << "rhs" << "," << "sol" << std::endl;

        for (int elm = 0; elm < grid.get_total_elements(); ++elm) {
            for (int basis_fun = 0; basis_fun < ipow(P, N); ++basis_fun) {
                double dV = 1.;
                for (int dim = 0; dim < N; ++dim) {
                    dV *= grid.get_dx(dim);
                }
                file << std::setprecision(16) << dV << ","
                     << std::setprecision(16) << rhs[elm](basis_fun) << ","
                     << std::setprecision(16) << sol[elm](basis_fun) << std::endl;
            }
        }

        file.close();
    }
};




#endif //LDG_POISSON_POISSON_HPP
