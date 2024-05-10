#ifndef ELECTRODES_HH
#define ELECTRODES_HH

#include <mutex>

// ROOT
#include "Math/Point3D.h"

// local
#include "geomodel.hh"
#include "fields.hh"

using namespace ROOT::Math;

//***********************************
// Charge signal class
// to be used as an interface
// to the algorithm.
//***********************************
class Electrode {
 private:
  // pointer to geometry for constructing fields
  GeometryModel* gm;
  // pointer to comsol fields for constructing fields
  ComsolFields* femfields;

  Fields* field; // specific for each electrode, constructed at creation
  std::mutex mtx;

 protected:

 public:
  // Constructor
  // using geometry data
  Electrode(ComsolFields* fem, GeometryModel* g);
  
  // Default destructor
  ~Electrode();

  // access
  void initfields(); // out of constructor - takes time.

  XYZPoint getFieldValue(bool& analytic, XYZPoint& p);
};
#endif
