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

void MeshProjector::move_particles_inside(SPH::FluidModel* model) const {
    unsigned int particle_num = model->getAllPosition().size();

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