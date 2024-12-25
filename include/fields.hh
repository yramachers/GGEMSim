#ifndef FIELDS_HH
#define FIELDS_HH

#include <vector>
#include <string>
#include <mutex>

// ROOT includes
#include "TKDTree.h"
#include "Math/Point3D.h"

// us
#include "geomodel.hh"

using namespace ROOT::Math;

//***********************************
// Field map classes
//***********************************
class ComsolFields {
 private:
  // scaling the weighting field for two electrodes
  double bias;

  std::mutex mtx;

  // which ROOT file to read the weighting field
  std::vector<XYZPoint> coords;
  std::vector<XYZPoint> dmap;

 protected:

 public:
  // Constructor
  ComsolFields() = default; // from file
  
  // Default destructor
  ~ComsolFields() = default;

  // Methods
  void read_fields(std::string fname);
  void setBias(double b) {bias = b;};
  void clear();
  std::vector<XYZPoint> positions() {return coords;}
  std::vector<XYZPoint> driftmap() {return dmap;}
};


class Fields {
 private:
  // access to geometry
  GeometryModel* gm;
  
  // container for field coordinates here
  TKDTreeID* coordinates;

  // storage container
  double* allx;
  double* allz;
  double* alldx;
  double* alldz;

 protected:

 public:
  // Constructor
  Fields(GeometryModel* g); // from file
  
  // Default destructor
  ~Fields();

  // Methods
  void prepare_fields(ComsolFields& fem);

  // return field values in [V/m]
  XYZPoint getFieldValue(XYZPoint& p, bool& analytic);  
};
#endif
