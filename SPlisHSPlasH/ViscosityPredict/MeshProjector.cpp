#include "MeshProjector.h"

bool MeshProjector::load_mesh(const std::string &filepath) {
    if (!m_mesh.is_empty()) {
        m_mesh.clear();
    }
    if (!CGAL::IO::read_PLY(filepath, m_mesh)) {
        std::cerr << "read mesh PLY error: " << filepath << std::endl;
        return false;
    }
    m_fill_holes_on_mesh(m_mesh);
    m_tree = Tree(faces(m_mesh).first, faces(m_mesh).second, m_mesh);
    m_tree.accelerate_distance_queries();
    m_inside_tester = std::make_unique<Inside_tester>(m_mesh);
    return true;
}

double MeshProjector::project_particles(SPH::FluidModel* model) const {
    double distance_sum{0.0};
    size_t particle_num = model->getAllPosition().size();
    for (size_t i = 0; i <= particle_num; ++i) {
        Point p_point = m_eigen_to_CGALPoint(model->getPosition(i));
        auto [is_projected, distance] = m_project_single_particle(p_point);
        if (is_projected) {
            distance_sum += distance;
            model->setPosition(i, m_cgalPoint_to_eigen(p_point));
            model->setVelocity(i, {0.0,0.0,0.0});
        }
    }
    return distance_sum;
}

void MeshProjector::move_particles_inside(SPH::FluidModel* model, const unsigned particle_num) const {

    // 1 find all the particles in the mesh
    std::vector<unsigned int> inside_indices;
    for (unsigned int i = 0; i < particle_num; ++i) {
        if (m_is_inside_mesh(m_eigen_to_CGALPoint(model->getPosition(i)))) {
            inside_indices.push_back(i);
        }
    }
    if (inside_indices.empty())
        return;

    // 2 change pos & vel if particle is out of mesh
    static std::mt19937 rng(
            static_cast<unsigned int>(
                    std::chrono::system_clock::now().time_since_epoch().count()
            )
    );  // gen random nums
    std::uniform_real_distribution<Real> dist_01(0.0f, 1.0f); // [0,1) distribution
    for (unsigned int i = 0; i < particle_num; ++i)
    {
        if (!m_is_inside_mesh(m_eigen_to_CGALPoint(model->getPosition(i))))
        {
            // 2.1 randomly choose a particle which is inside the mesh
            unsigned int idx_in = inside_indices[rng() % inside_indices.size()];
            Real rho = model->getDensity(idx_in);
            if (rho < 1.0e-8) {
                continue;
            }

            // 2.2 calcu reference radius
            Real m = model->getMass(idx_in);
            Real volume = m / rho;
            Real r_base = std::cbrt(volume);

            // 2.3 gen a random bias
            Vector3f rand_dir;
            while (true) {
                Real rx = dist_01(rng) * 2.0 - 1.0;
                Real ry = dist_01(rng) * 2.0 - 1.0;
                Real rz = dist_01(rng) * 2.0 - 1.0;
                rand_dir = Vector3f(rx, ry, rz);
                if (rand_dir.squaredNorm() <= 1.0)
                    break;
            }

            Real scale = 1.0;
            rand_dir *= (scale * r_base);

            // 2.4 new pos, new vel
            model->setPosition(i, model->getPosition(idx_in) + rand_dir);
            model->setVelocity(i, model->getVelocity(idx_in));
        }
    }
}

void MeshProjector::m_fill_holes_on_mesh(Mesh& mesh) {
    for(auto h : halfedges(mesh)) {
        if(CGAL::is_border(h, mesh)) {
            CGAL::Polygon_mesh_processing::triangulate_hole(mesh, h);
        }
    }
}

bool MeshProjector::m_is_inside_mesh(const Point& particle) const {
    CGAL::Bounded_side side = (*m_inside_tester)(particle);
    return (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY);
}

std::pair<bool, double> MeshProjector::m_project_single_particle(Point& particle) const {
    if (m_is_inside_mesh(particle)) {
        return {true, 0.0};
    }
    Point closest_point = m_tree.closest_point(particle);
    double distance_2 = CGAL::squared_distance(particle, closest_point);
    particle = closest_point;
    return {false, std::sqrt(distance_2)};
}

MeshProjector::Point MeshProjector::m_eigen_to_CGALPoint(const Vector3r &vec)
{
    return {
            static_cast<Real>(vec[0]),
            static_cast<Real>(vec[1]),
            static_cast<Real>(vec[2])
    };
}

Vector3r MeshProjector::m_cgalPoint_to_eigen(const Point &p){
    return {
            static_cast<Real>(p.x()),
            static_cast<Real>(p.y()),
            static_cast<Real>(p.z())
    };
}

void MeshProjector::sample_interior_points(std::size_t num_samples) {
    m_interior_samples.clear();
    if (m_mesh.is_empty() || num_samples == 0) return;

    // Build a point-in-polyhedron tester
    CGAL::Side_of_triangle_mesh<Mesh, Kernel> point_inside(m_mesh);

    // Compute mesh bounding box
    CGAL::Bbox_3 bb = CGAL::Polygon_mesh_processing::bbox(m_mesh);
    std::mt19937_64 gen((unsigned)std::random_device{}());
    std::uniform_real_distribution<double> dx(bb.xmin(), bb.xmax());
    std::uniform_real_distribution<double> dy(bb.ymin(), bb.ymax());
    std::uniform_real_distribution<double> dz(bb.zmin(), bb.zmax());

    // Sample until we have enough interior points
    while (m_interior_samples.size() < num_samples)
    {
        Kernel::Point_3 p(dx(gen), dy(gen), dz(gen));
        if (point_inside(p) == CGAL::ON_BOUNDED_SIDE)
            m_interior_samples.push_back(p);
    }
}


double MeshProjector::projection_loss_with_interior(
    const std::vector<Vector3r> &pts,
    const unsigned activate_num,
    bool use_square_error) const {
    if (pts.empty()) return 0.0;
    // 1) Outside particles loss to mesh surface
    double sum_out = 0.0;
    const int N = static_cast<int>(activate_num);

    // CGAL::Side_of_triangle_mesh<Mesh, Kernel> point_inside(m_mesh);
    // int i = 0;
    // for (auto const &v : pts)
    // {
    //     if (i >= N) break;
    //     ++i;
    //     Kernel::Point_3 p(v.x(), v.y(), v.z());
    //     if (point_inside(p) == CGAL::ON_UNBOUNDED_SIDE)
    //     {
    //         Kernel::Point_3 c = m_tree.closest_point(p);
    //         double sqd = CGAL::squared_distance(p, c);
    //         sum_out += (use_square_error ? sqd : std::sqrt(sqd));
    //     }
    // }

    #pragma omp parallel reduction(+:sum_out)
    {
        CGAL::Side_of_triangle_mesh<Mesh, Kernel> point_inside(m_mesh);
        #pragma omp for schedule(static)
        for (int i = 0; i < N; i++) {
            Kernel::Point_3 p(pts[i].x(), pts[i].y(), pts[i].z());
            if (point_inside(p) == CGAL::ON_UNBOUNDED_SIDE) {
                Kernel::Point_3 c = m_tree.closest_point(p);
                double sqd = CGAL::squared_distance(p, c);
                sum_out += (use_square_error ? sqd : std::sqrt(sqd));
            }
        }
    }

    // 2) Interior samples loss to fluid particles via KD-search
    using Traits = CGAL::Search_traits_3<Kernel>;
    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;
    using Tree3 = Neighbor_search::Tree;

    // Build KD-tree on fluid positions
    std::vector<Kernel::Point_3> fluid_pts;
    fluid_pts.reserve(N);
    for (int i = 0; i < N; ++i) {
        fluid_pts.emplace_back(pts[i].x(), pts[i].y(), pts[i].z());
    }
    Tree3 tree3(fluid_pts.begin(), fluid_pts.end());

    double sum_in = 0.0;
    const int M = m_interior_samples.size();
    // for (auto const &sp : m_interior_samples)
    // {
    //     Neighbor_search search(tree3, sp, 1);
    //     double sqd = search.begin()->second;
    //     sum_in += (use_square_error ? sqd : std::sqrt(sqd));
    // }
    #pragma omp parallel reduction(+:sum_in)
    {
        #pragma omp for schedule(static)
        for (int i =0; i < M; i++) {
            Neighbor_search search(tree3, m_interior_samples[i], 1);
            double sqd = search.begin()->second;
            sum_in += (use_square_error ? sqd : std::sqrt(sqd));
        }
    }
    // Return average loss
    double total = sum_out + sum_in;
    // return total / static_cast<double>(N + M);
    return total;
}
