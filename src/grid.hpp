//
// Created by mblomquist on 7/12/23.
//

#ifndef DG_UTILS_GRID_HPP
#define DG_UTILS_GRID_HPP



#include <iostream>
#include <fstream>
#include <unordered_map>
#include "uvector.hpp"
#include "multiloop.hpp"


template<int N>
class uniformGrid
{
    algoim::uvector<int, N> elements_per_dim;

    algoim::uvector<double, N> domain_min;
    algoim::uvector<double, N> domain_max;

    algoim::uvector<double, N> dx;

    bool print_vtk_header = false;

    algoim::uvector<bool, N> is_periodic = false;

    int zOrderMap(algoim::uvector<int,N> elm)
    {
        int id = 0;

        for (int dim = 0; dim < N; ++dim) {
            int lower_dims = 1;
            for (int i = dim; i > 0; --i) {
                lower_dims *= elements_per_dim(i-1);
            }
            id += elm(dim) * lower_dims;
        }

        return id;
    }

    int getNeighbor(algoim::uvector<int,N> elm, int dim, int dir)
    {
        algoim::uvector<int,N> neighbor;

        neighbor = elm;

        neighbor(dim) = elm(dim) + dir;

        int n_elm = -1;

        if (neighbor(dim) >= 0 && neighbor(dim) < elements_per_dim(dim))
            n_elm = zOrderMap(neighbor);

        return n_elm;
    }

    algoim::uvector<int,4*(N-1)> getNodes(algoim::uvector<int,N> elm)
    {
        algoim::uvector<int,4*(N-1)> nodes;

        for (int dim = 2; dim <= N; ++dim)
        {
            nodes(4*(dim-2))   = (N == 3) ? elm(2) + (elements_per_dim(2)+1) * (elm(1) + elm(0) * (elements_per_dim(1)+1)) + (dim-2)*(elements_per_dim(2)+1)*(elements_per_dim(1)+1) : elm(1) + (elements_per_dim(1)+1) * elm(0);
            nodes(4*(dim-2)+1) = nodes(4*(dim-2)) + 1;
            nodes(4*(dim-2)+2) = (N == 3) ? nodes(4*(dim-2)) + elements_per_dim(2)+1 : nodes(4*(dim-2)) + elements_per_dim(1)+1;
            nodes(4*(dim-2)+3) = (N == 3) ? nodes(4*(dim-2)) + elements_per_dim(2)+2 : nodes(4*(dim-2)) + elements_per_dim(1)+2;
        }

        return nodes;
    }

    void update_dx()
    {
        for (int i = 0; i < N; ++i) {
            dx(i) = (domain_max(i) - domain_min(i)) / elements_per_dim(i);
        }
    }

public:

    uniformGrid()
    {
        algoim::uvector<int, N> elements = 1;
        algoim::uvector<int, N> min = 0;
        algoim::uvector<int, N> max = 1;

        uniformGrid<N>(elements, min, max);
    }

    uniformGrid(algoim::uvector<int, N> elements_per_dim_,
                algoim::uvector<double, N> domain_min_,
                algoim::uvector<double, N> domain_max_)
    {
        elements_per_dim = elements_per_dim_;

        domain_min = domain_min_;
        domain_max = domain_max_;

        update_dx();
    }


    void set_domain_min(algoim::uvector<double, N> domain_min_)
    {
        domain_min = domain_min_;
        update_dx();
    }

    void set_domain_max(algoim::uvector<double, N> domain_max_)
    {
        domain_max = domain_max_;
        update_dx();
    }

    void set_elements_per_dim(algoim::uvector<int, N> elements_per_dim_)
    {
        elements_per_dim = elements_per_dim_;
        update_dx();
    }

    int get_element_id(algoim::uvector<int, N> k)
    {
        return zOrderMap(k);
    }

    int get_node_id(algoim::uvector<int, N> node)
    {

        int node_id = node(0);

        for (int dim = 1; dim < N; ++dim) {
            int lower_dims = 1;
            for (int i = dim; i > 0; --i) {
                lower_dims *= (elements_per_dim(dim-1)+1);
            }
            node_id += node(dim) * lower_dims;
        }

        return node_id;
    }

    int get_total_elements()
    {
        int total = 1;

        for (int i = 0; i < N; ++i) {
            total *= elements_per_dim(i);
        }
        return total;
    }

    algoim::uvector<int, N> get_faces_by_dim()
    {
        algoim::uvector<int, N> total_per_dim = 0;

        int total_int = 1;

        for (int dim = 0; dim < N; ++dim) {

            total_int = elements_per_dim(dim)+1;

            for (int idim = 0; idim < N; ++idim) {

                if (idim != dim)
                    total_int *= elements_per_dim(dim);
            }
            total_per_dim(dim) += total_int;
        }

        return total_per_dim;
    }

    int get_total_faces() {
        int total = 0;
        algoim::uvector<int, N> total_per_dim = get_faces_by_dim();

        for (int dim = 0; dim < N; ++dim) {
            total += total_per_dim(dim);
        }
        return total;
    }

    int get_total_nodes()
    {
        int total = 1;

        for (int dim = 0; dim < N; ++dim) {
            total *= (elements_per_dim(dim)+1);
        }
        return total;
    }

    algoim::uvector<double, N> get_xmin()
    {
        return domain_min;
    }

    algoim::uvector<double, N> get_xmax()
    {
        return domain_max;
    }

    algoim::uvector<double, N> get_dx()
    {
        return dx;
    }

    double get_dx(int i)
    {
        if ((i < N) && (i >= 0))
            return dx(i);
        else
            std::cout << "Attempted to access dx out of bounds, returning dx(0)." << std::endl;
            return dx(0);
    }

    double get_volume()
    {
        double vol = 1.;

        for (int i = 0; i < N; ++i) {
            vol *= dx(i);
        }

        return vol;
    }

    algoim::uvector<int, N> get_elements_per_dim()
    {
        return elements_per_dim;
    }

    algoim::uvector<double, N> get_nodal_pos(algoim::uvector<int, N> node)
    {
        algoim::uvector<double, N> nodal_pos;

        for (int i = 0; i < N; ++i) {
            nodal_pos(i) = node(i) * dx(i) + domain_min(i);
        }

        return nodal_pos;
    }

    algoim::uvector<int, N> get_indx_from_pos(algoim::uvector<double, N> pos)
    {
        algoim::uvector<int, N> indx;

        for (int dim = 0; dim < N; ++dim) {
            indx(dim) = floor((pos(dim) - domain_min(dim)) / dx(dim));

            /*
             * note: this function has a bias in the positive direction,
             * so we need to make sure that we don't return an index value
             * greater than the max number of elements.
             */
            if (indx(dim) >= elements_per_dim(dim))
                --indx(dim);
        }

        return indx;
    }
    algoim::uvector<double, N> map_pos_to_ref_element(algoim::uvector<double, N> pos)
    {
        algoim::uvector<double, N> ref_pos = 0.;

        algoim::uvector<int, N> indx = get_indx_from_pos(pos);
        algoim::uvector<double, N> local_min = get_nodal_pos(indx);

        for (int dim = 0; dim < N; ++dim) {
            ref_pos(dim) = (pos(dim) - local_min(dim)) / dx(dim);
        }

        return ref_pos;
    }

    int get_element_from_pos(algoim::uvector<double, N> pos)
    {
        algoim::uvector<int, N> indx = get_indx_from_pos(pos);

        return get_element_id(indx);
    }

    algoim::uvector<double, N> get_element_max(algoim::uvector<int, N> element)
    {
        algoim::uvector<double, N> element_max;

        for (int i = 0; i < N; ++i) {
            element_max(i) = (element(i)+1)*dx(i);
        }

        return element_max;
    }

    void print_grid_to_vtk(const std::string& filename)
    {
        if (N == 2 || N == 3) {

            print_vtk_header = true;

            int numPoints = 1;
            int numCells = 1;

            for (int i = 0; i < N; ++i) {
                numPoints *= (elements_per_dim(i) + 1);
                numCells *= elements_per_dim(i);
            }

            std::ofstream file(filename);

            // Write the VTK header
            file << "# vtk DataFile Version 2.0" << std::endl;
            file << "QuadTree" << std::endl;
            file << "ASCII" << std::endl;
            file << "DATASET UNSTRUCTURED_GRID" << std::endl;

            // Write grid points
            file << "POINTS " << numPoints << " double" << std::endl;
            for (algoim::MultiLoop<N> i(0, elements_per_dim+1); ~i; ++i) {
                for (int j = 0; j < N; ++j) {
                    file << i(j) * dx(j) + domain_min(j) << " ";
                }
                if (N < 3)
                    file << "0.0" << std::endl;
                else
                    file << std::endl;
            }

            // Write the cell connectivity
            int numCellPoints = 4 * (N - 1);
            file << "CELLS " << numCells << " " << numCells * (numCellPoints + 1) << std::endl;
            for (algoim::MultiLoop<N> i(0, elements_per_dim); ~i; ++i)
            {
                file << numCellPoints << " ";

                algoim::uvector<int,4*(N-1)> element_nodes = getNodes(i());

                if (N == 2) {
                    file << element_nodes(0) << " "
                         << element_nodes(1) << " "
                         << element_nodes(3) << " "
                         << element_nodes(2) << " " << std::endl;
                } else {
                    file << element_nodes(0) << " "
                         << element_nodes(1) << " "
                         << element_nodes(3) << " "
                         << element_nodes(2) << " "
                         << element_nodes(4) << " "
                         << element_nodes(5) << " "
                         << element_nodes(7) << " "
                         << element_nodes(6) << " " << std::endl;
                }
            }

            // Write the cell types
            file << "CELL_TYPES " << numCells << std::endl;
            for (int i = 0; i < numCells; i++)
            {
                if (N == 2)
                    file << "9" << std::endl;  // VTK_QUAD cell type
                else
                    file << "12" << std::endl;  // VTK_OCT cell type
            }

            file.close();

        }
        else
        {
            if (N > 3)
                std::cout << "Yeah, I don't really know why you want to print a vtk file for N > 3." << std::endl;
            else
                std::cout << "Yeah, I don't really know why you want to print a vtk file for N < 2." << std::endl;
        }

    }

    void print_grid_point_data_to_vtk(const std::string& filename,
                                const std::string& data_name,
                                const std::vector<double> &data)
    {
        // make sure the header has already been printed
        if (!print_vtk_header)
            print_grid_to_vtk(filename);

        std::ofstream file(filename, std::ios::app);

        int numPoints = 1;

        for (int i = 0; i < N; ++i) {
            numPoints *= (elements_per_dim(i) + 1);
        }

        // Write the point data
        file << "POINT_DATA " << numPoints << std::endl;
        file << "SCALARS " << data_name << " double" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;

        // Write the grid values
        for (int j = 0; j < numPoints; j++) {
            file << data[j] << std::endl;
        }

        file.close();
    }

};

#endif //DG_UTILS_GRID_HPP