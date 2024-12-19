#ifndef FIELDS_HH
#define FIELDS_HH

#include <vector>
#include <string>

// ROOT includes
#include "TKDTree.h"
#include "Math/Point3D.h"


//local
#include "geomodel.hh"

using namespace ROOT::Math;

//***********************************
// Field map classes
//***********************************
class ComsolFields {
 private:
  // scaling the weighting field for two electrodes
  double bias;
  // which ROOT file to read the weighting field
  std::string fname;
  std::vector<XYZPoint> coords;
  std::vector<XYZPoint> dmap;

 protected:

 public:
  // Constructor
  ComsolFields(std::string fname); // from file
  
  // Default destructor
  ~ComsolFields() = default;

  // Methods
  void read_fields();
  void setBias(double b) {bias = b;};
  std::vector<XYZPoint> positions() {return coords;}
  std::vector<XYZPoint> driftmap() {return dmap;}
};


class Fields {
 private:
  // pointer to geometry for asking
  GeometryModel* gm;
  // container for field coordinates here
  TKDTreeID* coordinates;
  // storage container
  double* allx;
  double* allz;
  double* alldx;
  double* alldz;

 protected:
  void prepare_fields(ComsolFields* fem);
  XYZPoint getFieldValue(XYZPoint& p, bool& analytic);  

 public:
  // Constructor
  Fields(ComsolFields* fem, GeometryModel* gm); // from file
  
  // Default destructor
  ~Fields();

  // Methods
  // return field values in [V/m]
  XYZPoint getDriftField(XYZPoint& p, bool& analytic);
};
#endif
