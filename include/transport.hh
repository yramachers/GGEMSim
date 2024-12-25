#ifndef TRANSPORT_HH
#define TRANSPORT_HH

#include <list>
#include <vector>
#include <random>
#include <mutex>

// ROOT
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

//local
#include "fields.hh"
#include "pmodel.hh"

using namespace ROOT::Math;

//***********************************
// Charge signal class
// to be used as an interface
// to the algorithm.
//***********************************
class Transport {
 private:
  int photon_number;
  int ion_number;
  double density;
  std::list<XYZPoint> charges;
  std::vector<XYZPoint> chargeStore;
  std::vector<XYZPoint> photonStore;
  std::mutex mtx;

  std::mt19937 generator;
  std::uniform_real_distribution<double> rnd(0.0, 1.0);
  std::uniform_real_distribution<double> rndHalf(0.01, 0.49);  

  Physics_Model* pm;
 
  // used by task function
  void book_charge(XYZPoint q);
  void book_photon(XYZPoint q);
  void addToGammas(int g);
  void addToIons(int i);
  void kin_factor(XYZVector& v0, bool momentum_flag);


 protected:
  bool taskfunction(Fields& fd, XYZPoint q, double en);

 public:
  // Constructor
  Transport(int seed);
  
  // Default destructor
  ~Transport();

  // Methods
  // preparation, required input from main()
  // otherwise no transport possible
  // work on this electrode id with charges
  int transport(Fields& fd, std::list<XYZPoint>& q, double energy); 

  int getPhotons() {return photon_number;}
  int getIons() {return ion_number;}
  double getDensity() {return density;}
  void setDensity(double d) {density = d;};
  std::vector<XYZPoint> allphotons() {return photonStore;}
  std::vector<XYZPoint> allcharges() {return chargeStore;}
};
#endif
