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
// Transport class
// to be used as an interface
// to the algorithm.
//***********************************
class Transport {
 private:
  bool zero_latch;

  int photon_number;
  int half_counter;
  int ion_number;
  int anode_number;
  double density;
  std::list<XYZPoint> charges;
  std::vector<XYZPoint> chargeStore;
  std::vector<XYZPoint> photonStore;
  std::mutex mtx;

  std::mt19937 generator;
  std::uniform_real_distribution<double> rnd;

  Physics_Model* pm;
  Fields* fd;
  
  // used by task functions
  void addAnode();
  void addHalf();
  void book_charge(XYZPoint q, bool flag);
  void book_photon(XYZPoint q);
  void addToGammas(int g);
  void addToIons(int i);
  void kin_factor(XYZVector& v0, bool momentum_flag);


 protected:
  bool singletask(XYZPoint q, double en);
  bool multitask(XYZPoint q, double en);

 public:
  // Constructor
  Transport(Fields* f, int seed);
  
  // Default destructor
  ~Transport();

  // Methods
  int single_transport(std::list<XYZPoint>& q, double energy); 
  int multi_transport(std::list<XYZPoint>& q, double energy); 

  int getPhotons() {return photon_number;}
  int getHalfCounter() {return half_counter;}
  int getIons() {return ion_number;}
  double getDensity() {return density;}
  void setDensity(double d) {density = d;};
  std::vector<XYZPoint> allphotons() {return photonStore;}
  std::vector<XYZPoint> allcharges() {return chargeStore;}
  inline void clear_counters() {
    half_counter=photon_number=ion_number=anode_number=0;
    zero_latch = false;
  }
};
#endif
