#pragma once

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
#include "multigrid.hpp"

#include "../../eigen/Eigen/Sparse"

// Supported boundary condition types (per dimension, per side).
//   BC_PERIODIC  : periodic pairing of the two domain boundaries (default)
//   BC_DIRICHLET : u = g
//   BC_NEUMANN   : du/dn = g   (natural condition on the auxiliary flux)
//   BC_ROBIN     : alpha*u + beta*du/dn = g
enum BCType { BC_PERIODIC = 0, BC_DIRICHLET = 1, BC_NEUMANN = 2, BC_ROBIN = 3 };

template<int P, int N>
class PoissonSolver
{
    uniformGrid<N> grid;

    // Boundary condition type per dimension (dim) and side (0 = low, 1 = high).
    algoim::uvector<algoim::uvector<int, 2>, N> boundary_conditions;

    // Boundary data g and Robin coefficients per (dim, side).
    std::function<double(const algoim::uvector<double, N>&)> bc_g[N][2];
    double bc_alpha[N][2];
    double bc_beta[N][2];

    double c1 = 0.;
    double c2 = 1.-c1;

    double tau_i = 1;
    double tau_d = 1;

    algoim::uvector<smatrix<double, ipow(P,N)>,N> D;

    BlockSparseMatrix<smatrix<double, ipow(P,N)>> L[N], T[N], G[N];

    elem_vec<P,N> rhs, sol, t_sol;

public:

    PoissonSolver()
    {
        boundary_conditions = 0;
        grid.set_periodic(true);

        for (int dim = 0; dim < N; ++dim) {
            for (int side = 0; side < 2; ++side) {
                bc_alpha[dim][side] = 0.;
                bc_beta[dim][side]  = 1.;
                bc_g[dim][side]     = nullptr;
            }
        }

        construct_broken_gradient_operator();
    }

    void solve_with_Multigrid()
    {
        construct_gradient_operator();
        MultiGrid<P,N> solver(grid, G);

        for (int i = 0; i < prod(grid.get_elements_per_dim()); ++i) {
            rhs[i] *= prod(grid.get_dx());
        }
        sol = solver.solve(rhs);
    }

    // Dispatch to the appropriate solver: geometric multigrid for the fully
    // periodic problem (unchanged), and a direct sparse solve for problems with
    // Dirichlet / Neumann / Robin boundaries.
    void solve()
    {
        if (is_all_periodic())
            solve_with_Multigrid();
        else
            solve_direct();
    }

    // Accumulate the boundary-condition contributions into the assembled system
    // matrix A (Robin surface mass) and right-hand side b:
    //   Dirichlet : b -= dV * G^T (b_G), where b_G is the g-part of the flux
    //               lifting  (n/dx) * \int_face g phi.
    //   Neumann   : b += \int_face g phi ds   (physical surface integral).
    //   Robin     : A(e,e) += (alpha/beta) \int_face phi phi ds
    //               b       += (1/beta)   \int_face g   phi ds.
    void add_boundary_terms(BlockSparseMatrix<smatrix<double, ipow(P,N)>> &A,
                            elem_vec<P,N> &b)
    {
        algoim::uvector<int, N> elements_per_dim = grid.get_elements_per_dim();
        double dV = prod(grid.get_dx());

        // g-part of the Dirichlet flux lifting, one auxiliary field per dim.
        elem_vec<P,N> b_G[N];
        bool have_dirichlet = false;

        for (int dim = 0; dim < N; ++dim) {

            if (grid.is_periodic(dim))
                continue;

            double h_dim         = grid.get_dx(dim);
            double face_measure  = dV / h_dim;   // physical area of the boundary face

            for (int side = 0; side < 2; ++side) {

                int type = boundary_conditions(dim)(side);
                if (type == BC_PERIODIC)
                    continue;

                int fixed  = (side == 0) ? 0 : elements_per_dim(dim) - 1;
                double n   = (side == 1) ? 1. : -1.;

                for (algoim::MultiLoop<N> e(0, elements_per_dim); ~e; ++e) {
                    if (e(dim) != fixed)
                        continue;

                    int eid = grid.get_element_id(e());

                    if (type == BC_DIRICHLET) {
                        // Weakly impose u_hat = g via the flux lifting:
                        // b_G[dim][e] += (n/dx) * \int_face g phi.
                        algoim::uvector<double, ipow(P,N)> load =
                                compute_face_load(dim, side, e(), bc_g[dim][side]);
                        b_G[dim][eid] += (n / h_dim) * load;
                        have_dirichlet = true;

                        // A fixed alternating-flux direction imposes Dirichlet
                        // data on u_hat naturally at only one of the two
                        // boundaries; a consistent penalty  sigma*<u - g, w>
                        // (sigma = tau_d * P^2 / h) restores coercivity at the
                        // other boundary.  It is consistent (u = g at the exact
                        // solution) and preserves the optimal order.
                        double sigma = tau_d * double(P * P) / h_dim;
                        A(eid, eid) += sigma * compute_face_mass(dim, side);
                        b[eid]      += sigma * load;
                    }
                    else if (type == BC_NEUMANN) {
                        algoim::uvector<double, ipow(P,N)> load =
                                compute_face_load(dim, side, e(), bc_g[dim][side]);
                        b[eid] += face_measure * load;
                    }
                    else if (type == BC_ROBIN) {
                        double alpha = bc_alpha[dim][side];
                        double beta  = bc_beta[dim][side];

                        // A(e,e) += (alpha/beta) * face mass
                        A(eid, eid) += (alpha / beta) * face_measure * compute_face_mass(dim, side);

                        // b += (1/beta) * face load
                        algoim::uvector<double, ipow(P,N)> load =
                                compute_face_load(dim, side, e(), bc_g[dim][side]);
                        b[eid] += (1. / beta) * face_measure * load;
                    }
                }
            }
        }

        // Dirichlet right-hand side correction: b -= dV * sum_dim G[dim]^T b_G[dim].
        if (have_dirichlet) {
            int num_elements = grid.get_total_elements();
            for (int dim = 0; dim < N; ++dim) {
                for (int i = 0; i < num_elements; ++i) {
                    for (auto k : G[dim].col[i]) {
                        b[i] -= dV * matvec(G[dim](k, i).transpose(), b_G[dim][k]);
                    }
                }
            }
        }
    }

    // Direct sparse solve of the LDG system for Dirichlet / Neumann / Robin
    // boundary conditions.  Assembles A = sum_dim dV * G[dim]^T G[dim] (plus the
    // Robin surface mass) and solves A sol = dV * rhs + boundary terms with
    // Eigen's SparseLU.  For the singular pure-Neumann problem a mean-zero
    // Lagrange constraint is appended so the (compatible) system is nonsingular
    // and the returned solution has zero mean.
    void solve_direct()
    {
        construct_gradient_operator();

        constexpr int S = ipow(P, N);
        int ne = grid.get_total_elements();
        double dV = prod(grid.get_dx());

        // A = sum_dim dV * G[dim]^T G[dim]
        BlockSparseMatrix<smatrix<double, S>> A;
        for (int dim = 0; dim < N; ++dim) {
            BlockSparseMatrix<smatrix<double, S>> Adim =
                    BlockSparse_matmat<P, N>(G[dim], G[dim], ne, true, dV);
            for (int i = 0; i < ne; ++i)
                for (auto j : Adim.row[i])
                    A(i, j) += Adim(i, j);
        }

        // b = dV * rhs  (mass-scaled right-hand side, matching the MG path)
        elem_vec<P,N> b;
        for (int i = 0; i < ne; ++i)
            b[i] = dV * rhs[i];

        // Fold in the boundary contributions (modifies A and b).
        add_boundary_terms(A, b);

        // Assemble the global sparse system.
        bool singular = is_singular();
        int  ndof     = ne * S;
        int  total    = singular ? ndof + 1 : ndof;

        std::vector<Eigen::Triplet<double>> trips;
        for (int i = 0; i < ne; ++i) {
            for (auto j : A.row[i]) {
                smatrix<double, S> &blk = A(i, j);
                for (int a = 0; a < S; ++a)
                    for (int c = 0; c < S; ++c) {
                        double v = blk(a, c);
                        if (v != 0.)
                            trips.emplace_back(i * S + a, j * S + c, v);
                    }
            }
        }

        Eigen::VectorXd bb = Eigen::VectorXd::Zero(total);
        for (int i = 0; i < ne; ++i)
            for (int a = 0; a < S; ++a)
                bb[i * S + a] = b[i](a);

        if (singular) {
            // Mean-zero constraint: sum_e (constant coeff of element e) = 0.
            for (int i = 0; i < ne; ++i) {
                trips.emplace_back(ndof, i * S + 0, 1.);
                trips.emplace_back(i * S + 0, ndof, 1.);
            }
        }

        Eigen::SparseMatrix<double> Amat(total, total);
        Amat.setFromTriplets(trips.begin(), trips.end());
        Amat.makeCompressed();

        Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
        lu.analyzePattern(Amat);
        lu.factorize(Amat);
        if (lu.info() != Eigen::Success)
            std::cout << "Warning: SparseLU factorization did not succeed." << std::endl;

        Eigen::VectorXd x = lu.solve(bb);

        // Unpack the solution back into the element-wise storage.
        for (int i = 0; i < ne; ++i)
            for (int a = 0; a < S; ++a)
                sol[i](a) = x[i * S + a];
    }

    // L2 error ||sol - t_sol|| (the solver's convention is sol -> +t_sol).  For
    // the pure Neumann problem the solution is only defined up to a constant, so
    // pass subtract_mean = true to compare after removing the domain-average
    // constant mode from both fields.
    double compute_l2_error(bool subtract_mean = false)
    {
        double l2_error = 0;

        int num_elements = prod(grid.get_elements_per_dim());

        double mean = 0.;
        if (subtract_mean) {
            for (int i = 0; i < num_elements; ++i)
                mean += sol[i](0) - t_sol[i](0);
            mean /= num_elements;
        }

        for (int i = 0; i < num_elements; ++i) {
            for (int j = 0; j < ipow(P, N); ++j) {
                double d = sol[i](j) - t_sol[i](j);
                if (subtract_mean && j == 0)
                    d -= mean;
                l2_error += prod(grid.get_dx()) * std::pow(d, 2);
            }
        }

        l2_error = std::sqrt(l2_error);

        std::cout << "L2 Error: " << l2_error << std::endl;

        return l2_error;
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

            // A dimension is only periodic when both of its faces are periodic.
            grid.set_periodic(dim, boundary_conditions(dim)(0) == BC_PERIODIC &&
                                   boundary_conditions(dim)(1) == BC_PERIODIC);
        }
    }

    // Set the boundary condition on a single face (dim, side).  For Dirichlet
    // and Neumann, g is the prescribed data; for Robin, alpha/beta are the
    // coefficients in alpha*u + beta*du/dn = g.
    void set_bc(int dim, int side, int type,
                std::function<double(const algoim::uvector<double, N>&)> g = nullptr,
                double alpha = 0., double beta = 1.)
    {
        boundary_conditions(dim)(side) = type;
        bc_g[dim][side]     = g;
        bc_alpha[dim][side] = alpha;
        bc_beta[dim][side]  = beta;

        grid.set_periodic(dim, boundary_conditions(dim)(0) == BC_PERIODIC &&
                               boundary_conditions(dim)(1) == BC_PERIODIC);
    }

    // Convenience: apply the same boundary condition to every domain face.
    void set_all_bc(int type,
                    std::function<double(const algoim::uvector<double, N>&)> g = nullptr,
                    double alpha = 0., double beta = 1.)
    {
        for (int dim = 0; dim < N; ++dim) {
            set_bc(dim, 0, type, g, alpha, beta);
            set_bc(dim, 1, type, g, alpha, beta);
        }
    }

    // True when the discrete operator is singular (constant nullspace): this
    // happens when no face pins the solution value, i.e. every face is either
    // periodic or Neumann.  Dirichlet and Robin faces both remove the constant
    // nullspace.
    bool is_singular()
    {
        for (int dim = 0; dim < N; ++dim) {
            for (int side = 0; side < 2; ++side) {
                int type = boundary_conditions(dim)(side);
                if (type == BC_DIRICHLET || type == BC_ROBIN)
                    return false;
            }
        }
        return true;
    }

    // True when every dimension is periodic (the original solver path).
    bool is_all_periodic()
    {
        for (int dim = 0; dim < N; ++dim)
            if (boundary_conditions(dim)(0) != BC_PERIODIC ||
                boundary_conditions(dim)(1) != BC_PERIODIC)
                return false;
        return true;
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

    void project_tsol(const std::function<double(const algoim::uvector<double,N>&)> &func)
    {
        l2_projection(t_sol, func);
    }

    void print_sol_on_uniform_grid()
    {
        evaluate_basis_on_uniform_grid(sol, 50, "../../solution.vtk", "solution");
        return;
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

    // Reference-face mass matrix  M_face(a,b) = \int_face phi_a phi_b dxi  for
    // the face of the reference element at (dim, side) (side 0 -> low, 1 -> high).
    smatrix<double, ipow(P,N)> compute_face_mass(int dim, int side)
    {
        smatrix<double, ipow(P,N)> M;

        constexpr int Q = P;

        algoim::uvector<double, N> pos;
        pos(dim) = (side == 0) ? 0. : 1.;

        for (algoim::MultiLoop<N-1> i(0, Q); ~i; ++i) {
            double weight = 1.;

            int t = 0;
            for (int j = 0; j < N; ++j) {
                if (j != dim) {
                    pos(j) = GaussQuad::x(Q, i(t));
                    weight *= GaussQuad::w(Q, i(t));
                    ++t;
                }
            }

            algoim::uvector<double, ipow(P,N)> e = evaluate_basis_as_point<P, N>(pos);
            M += outer_prod(e, e) * weight;
        }

        return M;
    }

    // Reference-face load vector  v(a) = \int_face g(x) phi_a dxi  for the face
    // of element 'elem' at (dim, side), with g evaluated at physical positions.
    algoim::uvector<double, ipow(P,N)> compute_face_load(int dim, int side,
                                                         const algoim::uvector<int, N> &elem,
                                                         const std::function<double(const algoim::uvector<double, N>&)> &g)
    {
        algoim::uvector<double, ipow(P,N)> load = 0.;

        if (!g)
            return load;

        constexpr int Q = P;

        algoim::uvector<double, N> ref_pos;
        ref_pos(dim) = (side == 0) ? 0. : 1.;

        algoim::uvector<double, N> emin  = grid.get_nodal_pos(elem);
        algoim::uvector<double, N> h     = grid.get_dx();

        for (algoim::MultiLoop<N-1> i(0, Q); ~i; ++i) {
            double weight = 1.;

            int t = 0;
            for (int j = 0; j < N; ++j) {
                if (j != dim) {
                    ref_pos(j) = GaussQuad::x(Q, i(t));
                    weight *= GaussQuad::w(Q, i(t));
                    ++t;
                }
            }

            algoim::uvector<double, N> x_phys;
            for (int j = 0; j < N; ++j)
                x_phys(j) = emin(j) + ref_pos(j) * h(j);

            algoim::uvector<double, ipow(P,N)> e = evaluate_basis_as_point<P, N>(ref_pos);
            load += e * (weight * g(x_phys));
        }

        return load;
    }

    // Add the interior-face lifting contributions for a single interior face
    // (alternating flux, u_hat taken from the upper element).
    void add_interior_face_lifting(int dim,
                                   const algoim::uvector<int, N> &element_i,
                                   const algoim::uvector<int, N> &element_j)
    {
        smatrix<double, ipow(P,N)> A_ii, A_ij, A_ji, A_jj;
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

    // Build the lifting operator, honouring the boundary condition on each
    // domain face.  Periodic dimensions reproduce the original wrap-around
    // treatment; for non-periodic dimensions the interior faces use the usual
    // alternating flux and the two domain boundaries are handled per BC type:
    //   Dirichlet -> lift (g - u^int): the -u^int part enters L here, the g
    //                part is deferred to the right-hand side (see below).
    //   Neumann/Robin -> u_hat = u^int, i.e. no lifting at the boundary face.
    void construct_lifting_operator()
    {
        algoim::uvector<int, N> elements_per_dim = grid.get_elements_per_dim();

        for (int dim = 0; dim < N; ++dim) {

            if (grid.is_periodic(dim)) {
                // Original periodic treatment: every face (including the
                // wrapped domain boundary at i(dim) == 0) is an interior face.
                for (algoim::MultiLoop<N> i(0, elements_per_dim); ~i; ++i) {
                    algoim::uvector<int, N> element_i, element_j;

                    element_i(dim) = (i(dim) - 1 == -1) ? elements_per_dim(dim) - 1 : i(dim) - 1;
                    element_j(dim) = i(dim);

                    for (int d = 0; d < N; ++d) {
                        if (d != dim) {
                            element_i(d) = i(d);
                            element_j(d) = i(d);
                        }
                    }

                    add_interior_face_lifting(dim, element_i, element_j);
                }
                continue;
            }

            // Non-periodic dimension.
            for (algoim::MultiLoop<N> i(0, elements_per_dim); ~i; ++i) {

                // Interior face between element (i-1) and element i in 'dim'.
                if (i(dim) >= 1) {
                    algoim::uvector<int, N> element_i = i(), element_j = i();
                    element_i(dim) = i(dim) - 1;
                    element_j(dim) = i(dim);
                    add_interior_face_lifting(dim, element_i, element_j);
                }

                // Lower domain boundary (side 0, outward normal -1).
                if (i(dim) == 0 && boundary_conditions(dim)(0) == BC_DIRICHLET) {
                    // u_hat = g: linear part lifts -u^int.  n = -1, so the
                    // coefficient -n/dx is +1/dx.
                    smatrix<double, ipow(P,N)> Mf = compute_face_mass(dim, 0);
                    L[dim](grid.get_element_id(i()), grid.get_element_id(i())) += Mf / grid.get_dx(dim);
                }

                // Upper domain boundary (side 1, outward normal +1).
                if (i(dim) == elements_per_dim(dim) - 1 && boundary_conditions(dim)(1) == BC_DIRICHLET) {
                    // u_hat = g: linear part lifts -u^int.  n = +1, so the
                    // coefficient -n/dx is -1/dx.
                    smatrix<double, ipow(P,N)> Mf = compute_face_mass(dim, 1);
                    L[dim](grid.get_element_id(i()), grid.get_element_id(i())) += (-1.) * Mf / grid.get_dx(dim);
                }
            }
        }
    }

    void construct_penalty_operator()
    {
        // maybe we only construct the penalty for specific BCs.
    }

    void construct_gradient_operator()
    {
        construct_lifting_operator();

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

    void print_operators_to_file(const std::string& filename)
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

    void print_vectors_to_file(const std::string& filename)
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
