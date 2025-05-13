#include "BoundaryModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;


BoundaryModel::BoundaryModel() :
	m_forcePerThread(),
	m_torquePerThread()
{		
}

BoundaryModel::~BoundaryModel(void)
{
	m_forcePerThread.clear();
	m_torquePerThread.clear();

	delete m_rigidBody;
}

void BoundaryModel::reset()
{
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}

void BoundaryModel::getForceAndTorque(Vector3r &force, Vector3r &torque)
{
	force.setZero();
	torque.setZero();
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		force += m_forcePerThread[j];
		torque += m_torquePerThread[j];
	}
}

void BoundaryModel::clearForceAndTorque()
{
	for (int j = 0; j < m_forcePerThread.size(); j++)
	{
		m_forcePerThread[j].setZero();
		m_torquePerThread[j].setZero();
	}
}

std::unique_ptr<BoundaryModel::BoundaryModelState> BoundaryModel::save_state() const
{
	auto s = std::make_unique<BoundaryModelState>();

	// 1) 刚体状态 -------------------------------------------
	s->position        = m_rigidBody->getPosition();
	s->rotation        = m_rigidBody->getRotation();
	s->velocity        = m_rigidBody->getVelocity();
	s->angularVelocity = m_rigidBody->getAngularVelocity();
	s->isDynamic       = m_rigidBody->isDynamic();

	// 2) 力 / 矩 ---------------------------------------------
	s->forcePerThread  = m_forcePerThread;   // 直接拷贝
	s->torquePerThread = m_torquePerThread;

	return s;     // 返回值拷贝或 NRVO
}

void BoundaryModel::load_state(const BoundaryModelState &s)
{
	// 1) 还原刚体状态 ---------------------------------------
	m_rigidBody->setPosition       (s.position);
	m_rigidBody->setRotation       (s.rotation);
	m_rigidBody->setVelocity       (s.velocity);
	m_rigidBody->setAngularVelocity(s.angularVelocity);

	// 2) 线程局部力 / 矩 -------------------------------------
	m_forcePerThread  = s.forcePerThread;
	m_torquePerThread = s.torquePerThread;
}
