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

#include "ReferenceConstForceKernels.h"
#include "ConstForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/ReferencePlatform.h"

using namespace ConstForcePlugin;
using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

void ReferenceCalcConstForceKernel::initialize(const System& system, const ConstForce& force) {
    // Initialize constant force.
    
    int numParticles = force.getNumParticles();
    particle.resize(numParticles);
    pforce.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
        force.getParticleForce(i, particle[i], pforce[i]);
}

double ReferenceCalcConstForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& force = extractForces(context);
    int numParticles = particle.size();
    double energy = 0;
    
    // Assign the constant force.
    
    for (int i = 0; i < numParticles; i++) {
        int p = particle[i];
        force[p] += pforce[i];
    }
    return energy;
}

void ReferenceCalcConstForceKernel::copyForceToContext(ContextImpl& context, const ConstForce& force) {
    if (force.getNumParticles() != particle.size())
        throw OpenMMException("updateForceInContext: The number of particles has changed");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int p;
        force.getParticleForce(i, p, pforce[i]);
        if (p != particle[i])
            throw OpenMMException("updateForceInContext: A particle index has changed");
    }
}
