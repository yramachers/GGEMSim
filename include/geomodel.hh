#ifndef GEOMODEL_HH
#define GEOMODEL_HH

#include <string>

// ROOT includes
#include "TGeoManager.h"

//***********************************
// Geometry class
//***********************************
class GeometryModel {
 private:
  TGeoManager* geom;

 protected:  


 public:

  // Constructor
  GeometryModel() = default; 
  
  // Default destructor
  ~GeometryModel();

  // Methods
  void init(std::string fn); // read GDML file
  int whereami(double xv, double yv, double zv); // int coding of regions
};
#endif
