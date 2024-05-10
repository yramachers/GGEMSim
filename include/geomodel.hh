#ifndef GEOMODEL_HH
#define GEOMODEL_HH

#include <vector>

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

  // all plates
  std::vector<TGeoNode*> plates; // stores electrodes as TGeoNodes
  

 protected:  
  void fill_plates();


 public:

  // Constructor
  // from file: geometry
  GeometryModel(const char* filename); 
  
  // Default destructor
  ~GeometryModel();

  // Methods
  int whereami(double xv, double yv, double zv); // int coding of regions

  // geometry get/set

  // access geometry data
  std::vector<TGeoNode*> electrodes() {return plates;}

};
#endif
