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

/**
 * This tests the Reference implementation of ConstForce.
 */

#include "ConstForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace ConstForcePlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerConstForceReferenceKernelFactories();

void testForce() {
    // Create a system with one particle.
    
    System system;
    vector<Vec3> positions(1);
    system.addParticle(1.0);
    positions[0] = Vec3(0.0, 0.0, 0.0);
    ConstForce* force = new ConstForce();
    Vec3 testForce = Vec3(-0.1, 0.0, 0.1);
    force->addParticle(0, testForce);
    system.addForce(force);
    
    // Compute the forces and energy.

    VerletIntegrator integ(1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    
    // See if the energy and forces are correct.
    
    double expectedEnergy = 0;
    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);

    for (int i = 0; i < 3; i++) {
        ASSERT_EQUAL_TOL(state.getForces()[0][i], testForce[i], 1e-5);
    }

    // Modify the parameters.
    
    const Vec3 testForce2 = Vec3(0.2, -0.1, 0.3);
    expectedEnergy = 1.0;
    force->setParticleForce(0, 0, testForce2);
    force->setEnergy(expectedEnergy);
    force->updateForceInContext(context);
    state = context.getState(State::Energy | State::Forces);

    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);

    for (int i = 0; i < 3; i++) {
        ASSERT_EQUAL_TOL(state.getForces()[0][i], testForce2[i], 1e-5);
    }
}

int main() {
    try {
        registerConstForceReferenceKernelFactories();
        testForce();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
