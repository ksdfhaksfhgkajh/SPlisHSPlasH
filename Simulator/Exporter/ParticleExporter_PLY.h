//
// Created by huaxi on 2024/12/21.
//

#ifndef SPLISHSPLASH_PARTICLEEXPORTER_PLY_H
#define SPLISHSPLASH_PARTICLEEXPORTER_PLY_H

#include "ExporterBase.h"
#include <fstream>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH {

    class ParticleExporter_PLY : public ExporterBase {
    protected:
        std::string m_exportPath;
        std::ofstream *m_outfile;

    public:
        ParticleExporter_PLY(SimulatorBase *base);

        ParticleExporter_PLY(const ParticleExporter_PLY &) = delete;

        ParticleExporter_PLY &operator=(const ParticleExporter_PLY &) = delete;

        virtual ~ParticleExporter_PLY(void);

        virtual void init(const std::string &outputPath);

        virtual void step(const unsigned int frame);

        virtual void reset();

        virtual void setActive(const bool active);
    };
}
#endif //SPLISHSPLASH_PARTICLEEXPORTER_PLY_H
