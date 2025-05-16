#ifndef SPLISHSPLASH_MESHPROJECTOR_H
#define SPLISHSPLASH_MESHPROJECTOR_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/PLY/PLY_reader.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"

#include <vector>
#include <random>   // std::mt19937, std::uniform_real_distribution

class MeshProjector {
public:
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point = Kernel::Point_3;
    using Mesh = CGAL::Surface_mesh<Point>;

    // typedef concerning AABB tree
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
    using AABB_traits = CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree = CGAL::AABB_tree<AABB_traits>;
    using Inside_tester = CGAL::Side_of_triangle_mesh<Mesh, Kernel>;

private:
    Mesh m_mesh;
    Tree m_tree;
    std::unique_ptr<Inside_tester> m_inside_tester;
    bool m_is_project;
    std::vector<Kernel::Point_3> m_interior_samples;

    static void m_fill_holes_on_mesh(Mesh& mesh);
    bool m_is_inside_mesh(const Point& particle) const;
    std::pair<bool, double> m_project_single_particle(Point& particle) const;
    static Point m_eigen_to_CGALPoint(const Vector3r &vec);
    static Vector3r m_cgalPoint_to_eigen(const Point &p);

public:
    MeshProjector() : m_is_project(true) {}
    bool load_mesh(const std::string& filepath);
    double project_particles(SPH::FluidModel* model) const;
    void move_particles_inside(SPH::FluidModel* model, const unsigned particle_num, const Real radius) const;

    void sample_interior_points(std::size_t num_samples);
    double projection_loss_with_interior(
        const std::vector<Vector3r> &particlePositions,
        const unsigned activate_num,
        bool use_square_error = true) const;

    void reset_is_project() {
        m_is_project = false;
    }
    bool get_is_project() const {
        return m_is_project;
    }
};

#endif
