//
// Created by huaxi on 2024/12/21.
//

#include "ParticleExporter_PLY.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;
using namespace Utilities;
using namespace std;

ParticleExporter_PLY::ParticleExporter_PLY(SimulatorBase* base) :
        ExporterBase(base)
{
    m_outfile = nullptr;
}

ParticleExporter_PLY::~ParticleExporter_PLY(void)
{
}

void ParticleExporter_PLY::init(const std::string& outputPath)
{
    // define output path for the data
    m_exportPath = FileSystem::normalizePath(outputPath + "/ply");
}

void ParticleExporter_PLY::step(const unsigned int frame)
{
    // check if the exporter is active
    if (!m_active)
        return;

    Simulation* sim = Simulation::getCurrent();
    for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++) {
        FluidModel *model = sim->getFluidModel(i);

        //define the exportFileName
        std::string fileName = "Fluid_" + std::to_string(i) + "_" + std::to_string(frame);
        std::string exportFileName = FileSystem::normalizePath(m_exportPath + "/" + fileName + ".ply");

        // Define the fstream object m_outfile. Open the file
        m_outfile = new std::ofstream(exportFileName, std::ios_base::trunc);
        if (!m_outfile->is_open()) {
            LOG_WARN << "Cannot open a file to save PLY particles.";
            return;
        }

        const unsigned int numParticles = model->numActiveParticles();
        const float particleRadius = sim->getParticleRadius();


        *m_outfile << "ply\n";
        *m_outfile << "format ascii 1.0\n";
        *m_outfile << "comment Exported by SPlisHSPlasH custom exporter\n";

        *m_outfile << "element vertex " << numParticles << "\n";
        *m_outfile << "property float x\n";
        *m_outfile << "property float y\n";
        *m_outfile << "property float z\n";

        *m_outfile << "end_header\n";

        for (unsigned int j = 0; j < numParticles; j++) {
            Vector3f pos = model->getPosition(j);
            *m_outfile << pos.x() << " " << pos.y() << " " << pos.z() << "\n";
        }
    }
    m_outfile->close();
    delete m_outfile;
    m_outfile = nullptr;
}

void ParticleExporter_PLY::reset()
{
}

void ParticleExporter_PLY::setActive(const bool active)
{
    ExporterBase::setActive(active);
    // create output folder
    if (m_active)
        FileSystem::makeDirs(m_exportPath);
}

