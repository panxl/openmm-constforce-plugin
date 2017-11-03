%module constforceplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/* Convert C++ (int&, Vec3&) object to python list */
%typemap(argout, fragment="Vec3_to_PyVec3") (int& particle, Vec3& pforce) {
    PyObject* pyInt = PyInt_FromLong(*$1);
    PyObject* pyVec = Vec3_to_PyVec3(*$2);
    $result = Py_BuildValue("[N, N]", pyInt, pyVec);
}

%typemap(in, numinputs=0) (int& particle, Vec3& pforce) (int tempP, Vec3 tempF) {
    $1 = &tempP;
    $2 = &tempF;
}

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "ConstForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/
%pythonappend ConstForcePlugin::ConstForce::getParticleForce(int index, int& particle, Vec3& pforce) const %{
    val[1] = unit.Quantity(val[1], unit.kilojoule_per_mole/unit.nanometer)
%}

%pythonappend ConstForcePlugin::ConstForce::getEnergy() const %{
    val = unit.Quantity(val, unit.kilojoule_per_mole)
%}

namespace ConstForcePlugin {

class ConstForce : public OpenMM::Force {
public:
    ConstForce();

    int getNumParticles() const;

    int addParticle(int particle, Vec3 force = Vec3(0.0, 0.0, 0.0));

    void setParticleForce(int index, int particle, Vec3 pforce);

    double getEnergy() const;

    void setEnergy(double input_energy);

    void updateForceInContext(OpenMM::Context& context);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply int& OUTPUT {int& particle};
    %apply Vec3& OUTPUT {Vec3& pforce};
    void getParticleForce(int index, int& particle, Vec3& pforce) const;
    %clear int& particle;
    %clear Vec3& pforce;
};

}
