/* -------------------------------------------------------------------------- *
 *                                OpenMMExample                                 *
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

#include "ConstForceProxy.h"
#include "ConstForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace ConstForcePlugin;
using namespace OpenMM;
using namespace std;

ConstForceProxy::ConstForceProxy() : SerializationProxy("ConstForce") {
}

void ConstForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ConstForce& force = *reinterpret_cast<const ConstForce*>(object);
    SerializationNode& constforces = node.createChildNode("ConstForces");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int particle;
        Vec3 pforce;
        force.getParticleForce(i, particle, pforce);
        constforces.createChildNode("Force").setIntProperty("p", particle).setDoubleProperty("fx", pforce[0]).setDoubleProperty("fy", pforce[1]).setDoubleProperty("fz", pforce[2]);
    }
}

void* ConstForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    ConstForce* force = new ConstForce();
    try {
        const SerializationNode& constforces = node.getChildNode("ConstForces");
        for (int i = 0; i < (int) constforces.getChildren().size(); i++) {
            const SerializationNode& constforce = constforces.getChildren()[i];
            force->addParticle(constforce.getIntProperty("p"), Vec3(constforce.getIntProperty("fx"), constforce.getDoubleProperty("fy"), constforce.getDoubleProperty("fz")));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
