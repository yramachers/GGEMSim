#ifndef GEOMODEL_HH
#define GEOMODEL_HH

#include <vector>
#include <string>

// ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"

// local

//***********************************
// Charge signal class
// to be used as an interface
// to the algorithm.
//***********************************
class GeometryModel {
 private:
  TGeoManager* geom;

 protected:  


 public:

  // Constructor
  // from file: geometry
  GeometryModel() = default; 
  
  // Default destructor
  ~GeometryModel();

  // Methods
  void init(std::string fn);
  int whereami(double xv, double yv, double zv); // int coding of regions
};
#endif
