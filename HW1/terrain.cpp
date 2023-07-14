#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx” from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    Eigen::Vector3f N_force = this->normal / this->normal.norm();
    for (int i = 0; i < jelly.getParticleNum(); i++) {
        // close to the wall & heading in
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f pos = particle.getPosition();
        Eigen::Vector3f velocity = particle.getVelocity();
        //Eigen::Vector3f force = particle.getForce();
        if (abs(N_force.dot(pos - this->position)) < eEPSILON && N_force.dot(velocity) < 0 && (pos - this->hole_position).norm() > this->hole_radius){
           
            Eigen::Vector3f vn = velocity.dot(N_force) * N_force;
            Eigen::Vector3f vt = velocity - vn;
            Eigen::Vector3f v_new = -coefResist * vn + vt;
            particle.setVelocity(v_new);

            if (N_force.dot(particle.getForce()) < 0) {
                Eigen::Vector3f contact_f = -(N_force.dot(particle.getForce()) * N_force);
                particle.addForce(contact_f);

                Eigen::Vector3f friction_f = -coefFriction * (-N_force.dot(particle.getForce()) * vt);
                particle.addForce(friction_f);
            }   
        }
    }
}
// BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are jelly.getParticleNum() particles.
    //   2. See TODOs in Jelly::computeInternalForce and functions in particle.h if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f velocity = particle.getVelocity();
        Eigen::Vector3f g_force = Eigen::Vector3f(0, -9.8f, 0);
        Eigen::Vector3f normal = (particle.getPosition() - position).normalized();
        Eigen::Vector3f N_force = normal / normal.norm();

        float hole_radius = this->radius / sqrt(2);
        float hole_dist = pow(particle.getPosition()(0) - HOLE_X, 2) + pow(particle.getPosition()(2) - HOLE_Z, 2);// since position = (HOLE_X, radius / sqrt(2), HOLE_Z)

        if (hole_dist <= pow(hole_radius, 2) && particle.getPosition()(1) <= 0) { // if in range of "half sphere",bowl && position.y < plain.y

            if (this->radius - (particle.getPosition() - position).norm() < eEPSILON) {  // distance < epsilon => collision
                //Eigen::Vector3f N_force = (pos - this->position) / (pos - this->position).norm();
                if (N_force.dot(velocity) > 0) {
                    Eigen::Vector3f vn = velocity.dot(N_force) * N_force;
                    Eigen::Vector3f vt = velocity - vn;
                    // bowl collision
                    Eigen::Vector3f v_new = -coefResist * (vn * (particle.getMass() - this->mass) / (particle.getMass() + this->mass)) + vt;
                    particle.setVelocity(v_new);

                    Eigen::Vector3f contact_f = -(N_force.dot(particle.getForce()) * N_force);
                    particle.addForce(contact_f);

                    Eigen::Vector3f friction_f = -coefFriction * (-N_force.dot(particle.getForce()) * vt);//
                    particle.addForce(friction_f);

                    float N_bowl = (particle.getMass() * -g_force).dot(normal);
                    particle.addVelocity(N_bowl * normal);
                }
            }
        }
    }
}
}  // namespace simulation
