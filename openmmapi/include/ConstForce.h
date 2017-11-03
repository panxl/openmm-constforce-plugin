#ifndef OPENMM_CONSTFORCE_H_
#define OPENMM_CONSTFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportConstForce.h"

using namespace OpenMM;

namespace ConstForcePlugin {

class OPENMM_EXPORT_CONSTFORCE ConstForce : public OpenMM::Force {
public:
    /**
     * Create an ConstForce.
     */
    ConstForce();
    /**
     * Get the number of particles.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Add a particle term to the force field.
     *
     * @param particle     the index of the particle this term is applied to
     * @param parameters   the list of parameters for the new force term
     * @return the index of the particle term that was added
     */
    int addParticle(int particle, const Vec3& force = Vec3(0.0, 0.0, 0.0));
    /**
     * Get the force and particle index of a constant force.
     * 
     * @param index     the index of the constant force
     * @param particle  the index of the particle
     * @param pforce    the value of the constant force
     */
    void getParticleForce(int index, int& particle, Vec3& pforce) const;
    /**
     * Set the force and particle index of a constant force.
     * 
     * @param index     the index of the constant force
     * @param particle  the index of the particle
     * @param pforce    the value of the constant force
     */
    void setParticleForce(int index, int particle, Vec3 pforce);
    /**
     * Get the energy
     * @return the energy
     */
    double getEnergy() const;
    /**
     * Set the energy
     */
    void setEnergy(double input_energy);
    /**
     * Update the constant forces in a Context to match those stored in this Force object.
     */
    void updateForceInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    std::vector<ParticleInfo> particles;
    double energy;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ConstForce::ParticleInfo {
public:
    int particle;
    Vec3 pforce;
    ParticleInfo() : particle(-1) {
    }
    ParticleInfo(int particle, const Vec3& pforce) : particle(particle), pforce(pforce) {
    }
};

} // namespace ConstForcePlugin

#endif /*OPENMM_CONSTFORCE_H_*/
