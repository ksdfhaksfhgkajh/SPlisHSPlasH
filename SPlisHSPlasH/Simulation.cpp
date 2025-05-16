#include "Simulation.h"
#include "TimeManager.h"
#include "Utilities/Timing.h"
#include "TimeStep.h"
#include "EmitterSystem.h"
#include "SPlisHSPlasH/WCSPH/TimeStepWCSPH.h"
#include "SPlisHSPlasH/PCISPH/TimeStepPCISPH.h"
#include "SPlisHSPlasH/PBF/TimeStepPBF.h"
#include "SPlisHSPlasH/IISPH/TimeStepIISPH.h"
#include "SPlisHSPlasH/DFSPH/TimeStepDFSPH.h"
#include "SPlisHSPlasH/PF/TimeStepPF.h"
#include "SPlisHSPlasH/ICSPH/TimeStepICSPH.h"
#include "BoundaryModel_Akinci2012.h"
#include "BoundaryModel_Bender2019.h"
#include "BoundaryModel_Koschier2017.h"
#include "Viscosity/ViscosityBase.h"


using namespace SPH;
using namespace std;
using namespace GenParam;

Simulation* Simulation::current = nullptr;
int Simulation::SIM_2D = -1;
int Simulation::PARTICLE_RADIUS = -1;
int Simulation::GRAVITATION = -1;
int Simulation::CFL_METHOD = -1;
int Simulation::CFL_FACTOR = -1;
int Simulation::CFL_MIN_TIMESTEPSIZE = -1;
int Simulation::CFL_MAX_TIMESTEPSIZE = -1;
int Simulation::ENABLE_Z_SORT = -1;
int Simulation::STEPS_PER_Z_SORT = -1;
int Simulation::KERNEL_METHOD = -1;
int Simulation::GRAD_KERNEL_METHOD = -1;
int Simulation::ENUM_KERNEL_CUBIC = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_KERNEL_POLY6 = -1;
int Simulation::ENUM_KERNEL_SPIKY = -1;
int Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_KERNEL_CUBIC_2D = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_GRADKERNEL_POLY6 = -1;
int Simulation::ENUM_GRADKERNEL_SPIKY = -1;
int Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC_2D = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::SIMULATION_METHOD = -1;
int Simulation::ENUM_CFL_NONE = -1;
int Simulation::ENUM_CFL_STANDARD = -1;
int Simulation::ENUM_CFL_ITER = -1;
int Simulation::ENUM_SIMULATION_WCSPH = -1;
int Simulation::ENUM_SIMULATION_PCISPH = -1;
int Simulation::ENUM_SIMULATION_PBF = -1;
int Simulation::ENUM_SIMULATION_IISPH = -1;
int Simulation::ENUM_SIMULATION_DFSPH = -1;
int Simulation::ENUM_SIMULATION_PF = -1;
int Simulation::ENUM_SIMULATION_ICSPH = -1;
int Simulation::BOUNDARY_HANDLING_METHOD = -1;
int Simulation::ENUM_AKINCI2012 = -1;
int Simulation::ENUM_KOSCHIER2017 = -1;
int Simulation::ENUM_BENDER2019 = -1;


Simulation::Simulation () 
{
	m_cflMethod = 1;
	m_cflFactor = 0.5;
	m_cflMinTimeStepSize = 0.0001;
	m_cflMaxTimeStepSize = 0.005;
	m_gravitation = Vector3r(0.0, -9.81, 0.0);

	m_kernelMethod = -1;
	m_gradKernelMethod = -1;

	m_neighborhoodSearch = nullptr;
	m_timeStep = nullptr;
	m_simulationMethod = SimulationMethods::NumSimulationMethods;
	m_simulationMethodChanged = NULL;

	m_simulationIsInitialized = false;
	m_sim2D = false;
	m_enableZSort = true;
	m_stepsPerZSort = 500;
	m_counter = 0;

	m_animationFieldSystem = new AnimationFieldSystem();
	m_boundaryHandlingMethod = static_cast<int>(BoundaryHandlingMethods::Bender2019);

	m_mesh_projector = std::make_unique<MeshProjector>();
	m_prev_state = std::make_unique<SimulationState>();
}

Simulation::~Simulation () 
{
#ifdef USE_DEBUG_TOOLS
	delete m_debugTools;
#endif
	delete m_animationFieldSystem;
	delete m_timeStep;
	delete m_neighborhoodSearch;
	delete TimeManager::getCurrent();

	for (unsigned int i = 0; i < m_fluidModels.size(); i++)
		delete m_fluidModels[i];
	m_fluidModels.clear();

	for (unsigned int i = 0; i < m_boundaryModels.size(); i++)
		delete m_boundaryModels[i];
	m_boundaryModels.clear();

	current = nullptr;
}

Simulation* Simulation::getCurrent()
{
	if (current == nullptr)
	{
		current = new Simulation ();
	}
	return current;
}

void Simulation::setCurrent (Simulation* tm)
{
	current = tm;
}

bool Simulation::hasCurrent()
{
	return (current != nullptr);
}

void Simulation::init(const Real particleRadius, const bool sim2D)
{
	m_sim2D = sim2D;
	initParameters();

	registerNonpressureForces();

	// init kernel
	setParticleRadius(particleRadius);

	setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
	setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);

	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
#ifdef GPU_NEIGHBORHOOD_SEARCH
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius);
#else
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius, false);
#endif
	m_neighborhoodSearch->set_radius(m_supportRadius);
}

void Simulation::deferredInit()
{
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel* fm = getFluidModel(i);
		fm->deferredInit();
	}
}

void Simulation::initParameters()
{
	ParameterObject::initParameters();

	SIM_2D = createBoolParameter("sim2D", "2D Simulation", &m_sim2D);
	setGroup(SIM_2D, "Simulation|Simulation");
	setDescription(SIM_2D, "2D/3D simulation.");
	getParameter(SIM_2D)->setReadOnly(true);

	ENABLE_Z_SORT = createBoolParameter("enableZSort", "Enable z-sort", &m_enableZSort);
	setGroup(ENABLE_Z_SORT, "Simulation|Simulation");
	setDescription(ENABLE_Z_SORT, "Enable z-sort to improve cache hits.");

	STEPS_PER_Z_SORT = createNumericParameter("enableZSort", "Simulation steps per z-sort", &m_stepsPerZSort);
	setGroup(STEPS_PER_Z_SORT, "Simulation|Simulation");
	setDescription(STEPS_PER_Z_SORT, "Number of simulation steps which are performed before a z-sort is applied.");
	static_cast<UnsignedIntParameter*>(getParameter(STEPS_PER_Z_SORT))->setMinValue(1u);

	ParameterBase::GetFunc<Real> getRadiusFct = std::bind(&Simulation::getParticleRadius, this);
	ParameterBase::SetFunc<Real> setRadiusFct = std::bind(&Simulation::setParticleRadius, this, std::placeholders::_1);
	PARTICLE_RADIUS = createNumericParameter("particleRadius", "Particle radius", getRadiusFct, setRadiusFct);
	setGroup(PARTICLE_RADIUS, "Simulation|Simulation");
	setDescription(PARTICLE_RADIUS, "Radius of the fluid particles.");
	getParameter(PARTICLE_RADIUS)->setReadOnly(true);

 	GRAVITATION = createVectorParameter("gravitation", "Gravitation", 3u, m_gravitation.data());
 	setGroup(GRAVITATION, "Simulation|Simulation");
 	setDescription(GRAVITATION, "Vector to define the gravitational acceleration.");

	CFL_METHOD = createEnumParameter("cflMethod", "CFL - method", &m_cflMethod);
	setGroup(CFL_METHOD, "Simulation|CFL");
	setDescription(CFL_METHOD, "CFL method used for adaptive time stepping.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(CFL_METHOD));
	enumParam->addEnumValue("None", ENUM_CFL_NONE);
	enumParam->addEnumValue("CFL", ENUM_CFL_STANDARD);
	enumParam->addEnumValue("CFL - iterations", ENUM_CFL_ITER);

	CFL_FACTOR = createNumericParameter("cflFactor", "CFL - factor", &m_cflFactor);
	setGroup(CFL_FACTOR, "Simulation|CFL");
	setDescription(CFL_FACTOR, "Factor to scale the CFL time step size.");
	static_cast<RealParameter*>(getParameter(CFL_FACTOR))->setMinValue(1e-6);

	CFL_MIN_TIMESTEPSIZE = createNumericParameter("cflMinTimeStepSize", "CFL - min. time step size", &m_cflMinTimeStepSize);
	setGroup(CFL_MIN_TIMESTEPSIZE, "Simulation|CFL");
	setDescription(CFL_MIN_TIMESTEPSIZE, "Min. time step size.");
	static_cast<RealParameter*>(getParameter(CFL_MIN_TIMESTEPSIZE))->setMinValue(1e-9);

	CFL_MAX_TIMESTEPSIZE = createNumericParameter("cflMaxTimeStepSize", "CFL - max. time step size", &m_cflMaxTimeStepSize);
	setGroup(CFL_MAX_TIMESTEPSIZE, "Simulation|CFL");
	setDescription(CFL_MAX_TIMESTEPSIZE, "Max. time step size.");
	static_cast<RealParameter*>(getParameter(CFL_MAX_TIMESTEPSIZE))->setMinValue(1e-6);

	ParameterBase::GetFunc<int> getKernelFct = std::bind(&Simulation::getKernel, this);
	ParameterBase::SetFunc<int> setKernelFct = std::bind(&Simulation::setKernel, this, std::placeholders::_1);
	KERNEL_METHOD = createEnumParameter("kernel", "Kernel", getKernelFct, setKernelFct);
	setGroup(KERNEL_METHOD, "Simulation|Kernel");
	setDescription(KERNEL_METHOD, "Kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_KERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_KERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_KERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_KERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_KERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_KERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_KERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getGradKernelFct = std::bind(&Simulation::getGradKernel, this);
	ParameterBase::SetFunc<int> setGradKernelFct = std::bind(&Simulation::setGradKernel, this, std::placeholders::_1);
	GRAD_KERNEL_METHOD = createEnumParameter("gradKernel", "Gradient of kernel", getGradKernelFct, setGradKernelFct);
	setGroup(GRAD_KERNEL_METHOD, "Simulation|Kernel");
	setDescription(GRAD_KERNEL_METHOD, "Gradient of the kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(GRAD_KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_GRADKERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_GRADKERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_GRADKERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_GRADKERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_GRADKERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getSimulationFct = std::bind(&Simulation::getSimulationMethod, this);
	ParameterBase::SetFunc<int> setSimulationFct = std::bind(&Simulation::setSimulationMethod, this, std::placeholders::_1);
	SIMULATION_METHOD = createEnumParameter("simulationMethod", "Simulation method", getSimulationFct, setSimulationFct);
	setGroup(SIMULATION_METHOD, "Simulation|Simulation");
	setDescription(SIMULATION_METHOD, "Simulation method.");
	enumParam = static_cast<EnumParameter*>(getParameter(SIMULATION_METHOD));
	enumParam->addEnumValue("WCSPH", ENUM_SIMULATION_WCSPH);
	enumParam->addEnumValue("PCISPH", ENUM_SIMULATION_PCISPH);
	enumParam->addEnumValue("PBF", ENUM_SIMULATION_PBF);
	enumParam->addEnumValue("IISPH", ENUM_SIMULATION_IISPH);
	enumParam->addEnumValue("DFSPH", ENUM_SIMULATION_DFSPH);
	enumParam->addEnumValue("Projective Fluids", ENUM_SIMULATION_PF);
	enumParam->addEnumValue("ICSPH", ENUM_SIMULATION_ICSPH);

	BOUNDARY_HANDLING_METHOD = createEnumParameter("boundaryHandlingMethod", "Boundary handling method", &m_boundaryHandlingMethod);
	setGroup(BOUNDARY_HANDLING_METHOD, "Simulation|Simulation");
	setDescription(BOUNDARY_HANDLING_METHOD, "Boundary handling method.");
	enumParam = static_cast<EnumParameter*>(getParameter(BOUNDARY_HANDLING_METHOD));
	enumParam->addEnumValue("Akinci et al. 2012", ENUM_AKINCI2012);
	enumParam->addEnumValue("Koschier and Bender 2017", ENUM_KOSCHIER2017);
	enumParam->addEnumValue("Bender et al. 2019", ENUM_BENDER2019);
	enumParam->setReadOnly(true);
}


void Simulation::setParticleRadius(Real val)
{
	m_particleRadius = val;
	m_supportRadius = static_cast<Real>(4.0) * m_particleRadius;
	initKernels();
}

void Simulation::initKernels()
{
	// init kernel
	Poly6Kernel::setRadius(m_supportRadius);
	SpikyKernel::setRadius(m_supportRadius);
	CubicKernel::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel::setRadius(m_supportRadius);
	PrecomputedCubicKernel::setRadius(m_supportRadius);
	CohesionKernel::setRadius(m_supportRadius);
	AdhesionKernel::setRadius(m_supportRadius);
	CubicKernel2D::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel2D::setRadius(m_supportRadius);
#ifdef USE_AVX
	CubicKernel_AVX::setRadius(m_supportRadius, m_sim2D);
	// 	Poly6Kernel_AVX::setRadius(m_supportRadius);
	// 	SpikyKernel_AVX::setRadius(m_supportRadius);
#endif
}

void Simulation::setGradKernel(int val)
{
	m_gradKernelMethod = val;

	if (!m_sim2D)
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 4))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
			m_gradKernelFct = CubicKernel::gradW;
		else if (m_gradKernelMethod == 1)
			m_gradKernelFct = WendlandQuinticC2Kernel::gradW;
		else if (m_gradKernelMethod == 2)
			m_gradKernelFct = Poly6Kernel::gradW;
		else if (m_gradKernelMethod == 3)
			m_gradKernelFct = SpikyKernel::gradW;
		else if (m_gradKernelMethod == 4)
			m_gradKernelFct = Simulation::PrecomputedCubicKernel::gradW;
	}
	else
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 1))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
			m_gradKernelFct = CubicKernel2D::gradW;
		else if (m_gradKernelMethod == 1)
			m_gradKernelFct = WendlandQuinticC2Kernel2D::gradW;
	}
}

void Simulation::setKernel(int val)
{
	if (val == m_kernelMethod)
		return;

	m_kernelMethod = val;
	if (!m_sim2D)
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 4))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel::W_zero();
			m_kernelFct = CubicKernel::W;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel::W;
		}
		else if (m_kernelMethod == 2)
		{
			m_W_zero = Poly6Kernel::W_zero();
			m_kernelFct = Poly6Kernel::W;
		}
		else if (m_kernelMethod == 3)
		{
			m_W_zero = SpikyKernel::W_zero();
			m_kernelFct = SpikyKernel::W;
		}
		else if (m_kernelMethod == 4)
		{
			m_W_zero = Simulation::PrecomputedCubicKernel::W_zero();
			m_kernelFct = Simulation::PrecomputedCubicKernel::W;
		}
	}
	else
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 1))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel2D::W_zero();
			m_kernelFct = CubicKernel2D::W;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel2D::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel2D::W;
		}
	}
	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		updateBoundaryVolume();
}

void Simulation::updateTimeStepSize()
{
	if (m_cflMethod == 1)
		updateTimeStepSizeCFL();
	else if (m_cflMethod == 2)
	{
		Real h = TimeManager::getCurrent()->getTimeStepSize();
		updateTimeStepSizeCFL();
		const unsigned int iterations = m_timeStep->getValue<unsigned int>(TimeStep::SOLVER_ITERATIONS);
		if (iterations > 10)
			h *= static_cast<Real>(0.9);
		else if (iterations < 5)
			h *= static_cast<Real>(1.1);
		h = min(h, TimeManager::getCurrent()->getTimeStepSize());
		TimeManager::getCurrent()->setTimeStepSize(h);
	}
}

void Simulation::updateTimeStepSizeCFL()
{
	const Real radius = m_particleRadius;
	Real h = TimeManager::getCurrent()->getTimeStepSize();
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.0;
	const Real diameter = static_cast<Real>(2.0)*radius;

	// fluid particles
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < numberOfFluidModels(); fluidModelIndex++)
	{
		FluidModel *fm = getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();
		for (unsigned int i = 0; i < numParticles; i++)
		{
			const Vector3r &vel = fm->getVelocity(i);
			const Vector3r &accel = fm->getAcceleration(i);
			const Real velMag = (vel + accel*h).squaredNorm();
			if (velMag > maxVel)
				maxVel = velMag;
		}
	}

	// boundary particles
	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		{
			BoundaryModel_Akinci2012 *bm = static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(i));
			if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
			{
				for (unsigned int j = 0; j < bm->numberOfParticles(); j++)
				{
					const Vector3r &vel = bm->getVelocity(j);
					const Real velMag = vel.squaredNorm();
					if (velMag > maxVel)
						maxVel = velMag;
				}
			}
		}
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		for (unsigned int boundaryModelIndex = 0; boundaryModelIndex < numberOfBoundaryModels(); boundaryModelIndex++)
		{
			BoundaryModel_Koschier2017 *bm = static_cast<BoundaryModel_Koschier2017*>(getBoundaryModel(boundaryModelIndex));
			if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
			{
				maxVel = std::max(maxVel, bm->getMaxVel());
			}
		}
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		for (unsigned int boundaryModelIndex = 0; boundaryModelIndex < numberOfBoundaryModels(); boundaryModelIndex++)
		{
			BoundaryModel_Bender2019 *bm = static_cast<BoundaryModel_Bender2019*>(getBoundaryModel(boundaryModelIndex));
			if (bm->getRigidBodyObject()->isDynamic() || bm->getRigidBodyObject()->isAnimated())
			{
				maxVel = std::max(maxVel, bm->getMaxVel());
			}
		}
	}

	// avoid division by zero
	if (maxVel < 1.0e-9)
		maxVel = 1.0e-9;

	// Approximate max. time step size 		
	h = m_cflFactor * static_cast<Real>(0.4) * (diameter / (sqrt(maxVel)));

	h = min(h, m_cflMaxTimeStepSize);
	h = max(h, m_cflMinTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

void Simulation::computeNonPressureForces()
{
	START_TIMING("computeNonPressureForces")
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->computeSurfaceTension();
		fm->computeViscosity();
		fm->computeVorticity();
		fm->computeDragForce();
		fm->computeElasticity();
		fm->computeXSPH();
	}
	STOP_TIMING_AVG
}

void Simulation::reset()
{
	// reset fluid models
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->reset();

	// reset boundary models
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->reset();

	if (getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		updateBoundaryVolume();

	if (m_timeStep)
		m_timeStep->reset();

	m_animationFieldSystem->reset();

	m_counter = 0;
	performNeighborhoodSearchSort();

	TimeManager::getCurrent()->setTime(0.0);
}

void Simulation::setSimulationMethod(const int val)
{
	SimulationMethods method = static_cast<SimulationMethods>(val);
	if ((method < SimulationMethods::WCSPH) || (method >= SimulationMethods::NumSimulationMethods))
		method = SimulationMethods::DFSPH;

	if (method == m_simulationMethod)
		return;

	delete m_timeStep;
	m_timeStep = nullptr;

	m_simulationMethod = method;

	if (method == SimulationMethods::WCSPH)
	{
		m_timeStep = new TimeStepWCSPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PCISPH)
	{
		m_timeStep = new TimeStepPCISPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PBF)
	{
		m_timeStep = new TimeStepPBF();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_POLY6);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_SPIKY);
	}
	else if (method == SimulationMethods::IISPH)
	{
		m_timeStep = new TimeStepIISPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::DFSPH)
	{
		m_timeStep = new TimeStepDFSPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::PF)
	{
		m_timeStep = new TimeStepPF();
		m_timeStep->init();

		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::ICSPH)
	{
		m_timeStep = new TimeStepICSPH();
		m_timeStep->init();
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}

	if (m_simulationMethodChanged != nullptr)
		m_simulationMethodChanged();
}



void Simulation::performNeighborhoodSearch()
{
	if (zSortEnabled())
	{
		if (m_counter % m_stepsPerZSort == 0)
		{
			performNeighborhoodSearchSort();
		}
		m_counter++;
	}
	START_TIMING("neighborhood_search");
	m_neighborhoodSearch->find_neighbors();
	STOP_TIMING_AVG;
}

void Simulation::performNeighborhoodSearchSort()
{
	if (!zSortEnabled())
		return;

	m_neighborhoodSearch->z_sort();

	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->performNeighborhoodSearchSort();
	}
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = getBoundaryModel(i);
		bm->performNeighborhoodSearchSort();
	}
	getTimeStep()->performNeighborhoodSearchSort();
#ifdef USE_DEBUG_TOOLS
	m_debugTools->performNeighborhoodSearchSort();
#endif
}

void Simulation::setSimulationMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_simulationMethodChanged = callBackFct;
}

void Simulation::setSimulationInitialized(int val) 
{ 
	m_simulationIsInitialized = val; 
	if (m_simulationIsInitialized)
		deferredInit();
}

void Simulation::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	model->emittedParticles(startIndex);
	m_timeStep->emittedParticles(model, startIndex);
#ifdef USE_DEBUG_TOOLS
	m_debugTools->emittedParticles(model, startIndex);
#endif
}

void Simulation::emitParticles()
{
	START_TIMING("emitParticles");
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->getEmitterSystem()->step();
	}
	STOP_TIMING_AVG
}

void Simulation::animateParticles()
{
	START_TIMING("animateParticles");
	m_animationFieldSystem->step();
	STOP_TIMING_AVG
}

void Simulation::addBoundaryModel(BoundaryModel *bm)
{
	m_boundaryModels.push_back(bm);
}

void Simulation::addFluidModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, unsigned int* fluidObjectIds, const unsigned int nMaxEmitterParticles)
{
	FluidModel *fm = new FluidModel();
	fm->initModel(id, nFluidParticles, fluidParticles, fluidVelocities, fluidObjectIds, nMaxEmitterParticles);
	m_fluidModels.push_back(fm);
}


void Simulation::updateBoundaryVolume()
{
	if (m_neighborhoodSearch == nullptr)
		return;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood

	// Activate only static boundaries
	LOG_INFO << "Initialize boundary volume";
	m_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		if (!getBoundaryModel(i)->getRigidBodyObject()->isDynamic() && !getBoundaryModel(i)->getRigidBodyObject()->isAnimated())
			m_neighborhoodSearch->set_active(i + nFluids, true, true);
	}

	//performNeighborhoodSearchSort();
	m_neighborhoodSearch->find_neighbors();

	// Boundary objects
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		if (!getBoundaryModel(body)->getRigidBodyObject()->isDynamic() && !getBoundaryModel(body)->getRigidBodyObject()->isAnimated())
			static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(body))->computeBoundaryVolume();
	}

	////////////////////////////////////////////////////////////////////////// 
	// Compute boundary psi for all dynamic bodies
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		// Deactivate all
		m_neighborhoodSearch->set_active(false);

		// Only activate next dynamic body
		if (getBoundaryModel(body)->getRigidBodyObject()->isDynamic() || getBoundaryModel(body)->getRigidBodyObject()->isAnimated())
		{
			m_neighborhoodSearch->set_active(body + nFluids, true, true);
			m_neighborhoodSearch->find_neighbors();
			static_cast<BoundaryModel_Akinci2012*>(getBoundaryModel(body))->computeBoundaryVolume();
		}
	}

	// Activate only fluids 
	m_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		for (unsigned int j = 0; j < numberOfFluidModels(); j++)
			m_neighborhoodSearch->set_active(i, j, true);
		for (unsigned int j = numberOfFluidModels(); j < m_neighborhoodSearch->point_sets().size(); j++)
			m_neighborhoodSearch->set_active(i, j, true);
	}
}

void SPH::Simulation::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(m_W_zero);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->saveState(binWriter);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->saveState(binWriter);
	m_timeStep->saveState(binWriter);
}

void SPH::Simulation::loadState(BinaryFileReader &binReader)
{
	binReader.read(m_W_zero);
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->loadState(binReader);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->loadState(binReader);
	m_timeStep->loadState(binReader);
}

void Simulation::particle_project(const unsigned frame_counter) {
	// Only support one phrase fluid
	if (numberOfFluidModels() > 1) return;

	std::string file_path = "../assets/";
	std::string file_name("mesh_");
	file_name += std::to_string(frame_counter) + ".ply";
	file_path += file_name;
	if  (std::filesystem::exists(file_path)) {
		m_mesh_projector->load_mesh(file_path);
		m_mesh_projector->move_particles_inside(getFluidModel(0),
												getFluidModel(0)->numActiveParticles(),
												getSupportRadius());
		getNeighborhoodSearch()->find_neighbors(true);
	} else {
		m_mesh_projector->reset_is_project();
	}
}

bool Simulation::load_mesh(const unsigned frame_counter) {
	std::string file_path = "../assets/";
	std::string file_name("mesh_");
	file_name += std::to_string(frame_counter) + ".ply";
	file_path += file_name;
	if  (std::filesystem::exists(file_path)) {
		m_mesh_projector->load_mesh(file_path);
		return true;
	}
	return false;
}

std::unique_ptr<Simulation::SimulationState> Simulation::save_simulation_state() const {
	auto simulation_state = std::make_unique<SimulationState>();
	simulation_state->fluid_state = getFluidModel(0)->save_state();
	for (unsigned int i = 0; i < numberOfBoundaryModels(); ++i) {
		simulation_state->boundary_states.push_back(getBoundaryModel(i)->save_state());
	}
	simulation_state->time_step_state = m_timeStep->save_state();
	simulation_state->current_time = TimeManager::getCurrent()->getTime();
	simulation_state->current_frame_number = TimeManager::getCurrent()->getFrame();
	simulation_state->counter = m_counter;
	simulation_state->valid_ = true;

	return simulation_state;
}

void Simulation::load_simulation_state(const SimulationState &state) {
	getFluidModel(0)->load_state(*state.fluid_state);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); ++i) {
		getBoundaryModel(i)->load_state(*state.boundary_states[i]);
	}
	m_timeStep->load_state(*state.time_step_state);
	m_counter = state.counter;
	performNeighborhoodSearchSort();
	TimeManager::getCurrent()->setTime(state.current_time);
	TimeManager::getCurrent()->setFrame(state.current_frame_number);
}

void Simulation::load_init_state() {
	TimeManager::getCurrent()->setTime(0.0);
	TimeManager::getCurrent()->setFrame(0);
	getFluidModel(0)->load_init_state();
	for (unsigned int i = 0; i < numberOfBoundaryModels(); ++i)
		getBoundaryModel(i)->reset();
	m_counter = 0;
	if (m_timeStep)
		m_timeStep->reset();
}

double Simulation::calcu_one_step_loss(const double *viscosity,
									   const int viscosity_dimension,
									   const std::unique_ptr<SimulationState> &curr_state,
									   std::unique_ptr<SimulationState> &best_state,
									   double& best_loss,
									   const unsigned frame_interval) {
	const auto fm = getFluidModel(0);
	 if (*m_prev_state) {
	 	load_simulation_state(*m_prev_state);
	 } else {
		load_init_state();
	 }

	fm->getViscosityBase()->setViscosity(static_cast<Real>(viscosity[0])); // change vis
	for (int i = 0; i < frame_interval; ++i) {
		getTimeStep()->step(); // run single step simulation
	}

	m_mesh_projector->sample_interior_points(fm->numActiveParticles());
	const double loss = m_mesh_projector->projection_loss_with_interior(
		fm->getAllPosition(),
		fm->numActiveParticles(),
		true
	);

	if (best_state == nullptr || loss < best_loss) {
		best_state = save_simulation_state();
		best_loss = loss;
	}
	load_simulation_state(*curr_state);

	return loss;
}

void Simulation::viscosity_predict(const unsigned frame_interval) {
	const Real vis_prev = getFluidModel(0)->getViscosityBase()->getViscosity();
	const double sigma = std::max(2.0 * vis_prev, 0.1);
	auto curr_state = save_simulation_state();
	std::unique_ptr<SimulationState> best_state;
	double best_loss = DBL_MAX;

	const std::vector<double> x0{vis_prev};
	libcmaes::CMAParameters<> cmaparams(1, x0.data(), sigma);

	cmaparams.set_ftarget(5.0);	// stop
	cmaparams.set_max_fevals(50);
	// libcmaes::FitFunc fsphere =
	//   std::bind(&Simulation::calcu_one_step_loss,
	// 			this,
	// 			std::placeholders::_1,
	// 			std::placeholders::_2);

	libcmaes::FitFunc fsphere = [this, &curr_state, &best_state, &best_loss, frame_interval]
	(const double *viscosity, const int viscosity_dimension) -> double {
		return calcu_one_step_loss(viscosity,
								   viscosity_dimension,
								   curr_state,
								   best_state,
								   best_loss,
								   frame_interval);
	};

	const libcmaes::CMASolutions cmasols = libcmaes::cmaes<>(fsphere, cmaparams);
	const auto best = cmasols.best_candidate();
	const double vis_opt = best.get_x()[0];
	const double fopt  = best.get_fvalue();

	// const double lamda = exp(-abs(vis_prev - std::max(0.0, vis_opt)));
	// const double vis_accum = (1.0 - lamda) * vis_opt + lamda * vis_prev;
	const double vis_accum = vis_opt;

	if (vis_accum < 0) {
		getFluidModel(0)->getViscosityBase()->setViscosity(0.0);
		std::cout << "best loss = " << fopt << ", viscosity = " << vis_accum << " reset to 0.0" << ", vis_opt = " << vis_opt << std::endl;
		return;
	}
	getFluidModel(0)->getViscosityBase()->setViscosity(static_cast<Real>(vis_accum));
	std::cout << "best loss = " << fopt << ", viscosity = " << vis_accum << ", vis_opt = " << vis_opt << std::endl;

	load_simulation_state(*best_state);
	m_prev_state = std::move(best_state);
}
