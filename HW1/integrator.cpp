#include "integrator.h"
#include <iostream>

#include <vector>
using namespace std;
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review “ODE_basics.pptx” from p.15 - p.16

    int particle_ct = particleSystem.getJellyPointer(0)->getParticleNum();// jelly index = 0
    
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);
        
        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        pos += velocity * delta_t;
        // Integrate velocity , v' = v + at 
        velocity += accel * delta_t;

        // set new position & velocity
        single_p->setAcceleration(accel);
        single_p->setPosition(pos);
        single_p->setVelocity(velocity);
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) { // particleSystem here 

    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review “ODE_basics.pptx” from p.18 - p.19
    int particle_ct = particleSystem.getJellyPointer(0)->getParticleNum();// jelly index = 0
    std::vector<Particle> previous_p;
    
    for (int i = 0; i < particle_ct; i++) {
        previous_p.push_back(particleSystem.getJellyPointer(0)->getParticle(i));
    }

    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);
        
        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        pos += velocity * delta_t;
        // Integrate velocity , v' = v + at 
        velocity += accel * delta_t;
        // set new position & velocity
        single_p->setAcceleration(accel);
        single_p->setPosition(pos);
        single_p->setVelocity(velocity);
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
        
    }

    for (int it = 0; it < 2; it++) { // make orginal data used , as the formula
        particleSystem.computeJellyForce(*particleSystem.getJellyPointer(0));
        for (int i = 0; i < particle_ct; i++) {
            Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

            float delta_t = particleSystem.deltaTime;
            // a' = g + f/m
            Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
            Eigen::Vector3f pos = previous_p[i].getPosition() + single_p->getVelocity() * delta_t;
            Eigen::Vector3f velocity = previous_p[i].getVelocity() + accel * delta_t;
            
            single_p->setAcceleration(accel);
            // set new position & velocity
            single_p->setPosition(pos);
            single_p->setVelocity(velocity);
            // since acceleration is added
            single_p->setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review “ODE_basics.pptx” from p .18 - p .20
    int particle_ct = particleSystem.getJellyPointer(0)->getParticleNum();// jelly index = 0
    std::vector<Particle> previous_p;

    for (int i = 0; i < particle_ct; i++) {
        previous_p.push_back(particleSystem.getJellyPointer(0)->getParticle(i));
    }

    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        pos += (velocity * delta_t)/2;
        // Integrate velocity , v' = v + at 
        velocity += (accel * delta_t)/2;
        // set new position & velocity
        single_p->setAcceleration(accel);
        single_p->setPosition(pos);
        single_p->setVelocity(velocity);
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());

    }

    particleSystem.computeJellyForce(*particleSystem.getJellyPointer(0));
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        float delta_t = particleSystem.deltaTime;
        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = previous_p[i].getPosition() + single_p->getVelocity() * delta_t;
        Eigen::Vector3f velocity = previous_p[i].getVelocity() + accel * delta_t;

        single_p->setAcceleration(accel);
        // set new position & velocity
        single_p->setPosition(pos);
        single_p->setVelocity(velocity);
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review “ODE_basics.pptx” from p.21
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    
    StateStep step;
    int particle_ct = particleSystem.getJellyPointer(0)->getParticleNum();// jelly index = 0
    std::vector<Particle> previous_p;
    std::vector<Particle> new_p;
    
    for (int i = 0; i < particle_ct; i++) {
        previous_p.push_back(particleSystem.getJellyPointer(0)->getParticle(i));
        new_p.push_back(particleSystem.getJellyPointer(0)->getParticle(i));
    }

    // k1
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        step.deltaPos = velocity * delta_t;
        // Integrate velocity , v' = v + at 
        step.deltaVel = accel * delta_t;
        
        // set new position & velocity
        single_p->setAcceleration(accel);
        //single_p->setPosition(pos);
        single_p->setVelocity(previous_p[i].getVelocity() + step.deltaVel);//
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
        new_p[i].addPosition(step.deltaPos / 6.0);
        new_p[i].addVelocity(step.deltaVel / 6.0);
    }
    // k2
    particleSystem.computeJellyForce(*particleSystem.getJellyPointer(0));
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        step.deltaPos = velocity * delta_t;
        // Integrate velocity , v' = v + at 
        step.deltaVel = accel * delta_t;

        // set new position & velocity
        single_p->setAcceleration(accel);
        //single_p->setPosition(pos);
        single_p->setVelocity(previous_p[i].getVelocity() + 0.5 * step.deltaVel);
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
        new_p[i].addPosition(step.deltaPos / 3.0);
        new_p[i].addVelocity(step.deltaVel / 3.0);
    }
    // k3
    particleSystem.computeJellyForce(*particleSystem.getJellyPointer(0));
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        step.deltaPos = velocity * delta_t;
        // Integrate velocity , v' = v + at 
        step.deltaVel = accel * delta_t;

        // set new position & velocity
        single_p->setAcceleration(accel);
        //single_p->setPosition(pos);
        single_p->setVelocity(previous_p[i].getVelocity() + 0.5 * step.deltaVel);//
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
        new_p[i].addPosition(step.deltaPos / 3.0);
        new_p[i].addVelocity(step.deltaVel / 3.0);
    }
    // k4
    particleSystem.computeJellyForce(*particleSystem.getJellyPointer(0));
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);

        // a' = g + f/m
        Eigen::Vector3f accel = particleSystem.gravity + (single_p->getForce() / single_p->getMass());
        Eigen::Vector3f pos = single_p->getPosition();
        Eigen::Vector3f velocity = single_p->getVelocity();
        float delta_t = particleSystem.deltaTime;

        // Integrate position first , x' = x + vt
        step.deltaPos = velocity * delta_t;
        // Integrate velocity , v' = v + at 
        step.deltaVel = accel * delta_t;

        // set new position & velocity
        single_p->setAcceleration(accel);
        //single_p->setPosition(pos);
        single_p->setVelocity(previous_p[i].getVelocity() + step.deltaVel);//
        // since acceleration is added
        single_p->setForce(Eigen::Vector3f::Zero());
        new_p[i].addPosition(step.deltaPos / 6.0);
        new_p[i].addVelocity(step.deltaVel / 6.0);
    }

    // update
    for (int i = 0; i < particle_ct; i++) {
        Particle* single_p = &particleSystem.getJellyPointer(0)->getParticle(i);
        Eigen::Vector3f pos = new_p[i].getPosition();
        Eigen::Vector3f velocity = new_p[i].getVelocity();
        single_p->setPosition(pos);
        single_p->setVelocity(velocity);
        single_p->setForce(Eigen::Vector3f::Zero());
    }
}
}  // namespace simulation
