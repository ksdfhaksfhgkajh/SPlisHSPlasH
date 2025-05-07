#include "BoundaryModel_Bender2019.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel_Bender2019::BoundaryModel_Bender2019() :
	m_boundaryVolume(),
	m_boundaryXj()
{		
	m_map = nullptr;
	m_maxDist = 0.0;
	m_maxVel = 0.0;
}

BoundaryModel_Bender2019::~BoundaryModel_Bender2019(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int i = 0; i < nModels; i++)
	{
		m_boundaryVolume[i].clear();
		m_boundaryXj[i].clear();
	}
	m_boundaryVolume.clear();
	m_boundaryXj.clear();

	delete m_map;
}


void BoundaryModel_Bender2019::initModel(RigidBodyObject *rbo)
{
	Simulation *sim = Simulation::getCurrent();

	const unsigned int nModels = sim->numberOfFluidModels();
	m_boundaryVolume.resize(nModels);
	m_boundaryXj.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_boundaryVolume[i].resize(fm->numParticles(), 0.0);
		m_boundaryXj[i].resize(fm->numParticles(), Vector3r::Zero());
	}

	if (rbo->isDynamic())
	{
#ifdef _OPENMP
		const int maxThreads = omp_get_max_threads();
#else
		const int maxThreads = 1;
#endif
		m_forcePerThread.resize(maxThreads, Vector3r::Zero());
		m_torquePerThread.resize(maxThreads, Vector3r::Zero());
	}

	m_rigidBody = rbo;
}

void BoundaryModel_Bender2019::reset()
{
	BoundaryModel::reset();

	m_maxDist = 0.0;
	m_maxVel = 0.0;
}

BoundaryModel_Bender2019::BoundaryModelState BoundaryModel_Bender2019::save_state() const
{
	/* ---- 把基类部分存进去 ---- */
	BoundaryModelState s = BoundaryModel::save_state();
	/* ---- 派生类自己的字段 ---- */
	s.max_dist = m_maxDist;
	s.max_vel  = m_maxVel;
	return s;
}

void BoundaryModel_Bender2019::load_state(const BoundaryModelState &s)
{
	/* ---- 先恢复基类 ---- */
	BoundaryModel::load_state(s);    // 隐式地把 s 视作 BoundaryModelState

	/* ---- 再恢复新增字段 ---- */
	m_maxDist = s.max_dist;
	m_maxVel  = s.max_vel;
}