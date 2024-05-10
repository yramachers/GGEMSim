#ifndef TRANSPORT_HH
#define TRANSPORT_HH

#include <list>
#include <vector>
#include <mutex>

// ROOT
#include "TRandom3.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"

//local
#include "electrode.hh"

using namespace ROOT::Math;

// helper struct
struct charge_t {
  XYZPoint location;
  double charge;
  int chargeID; // distinguish e- (1) and gamma (0)
};


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
  std::list<charge_t> charges;
  std::vector<XYZPoint> chargeStore;
  std::vector<XYZPoint> photonStore;
  std::mutex mtx;
  TRandom3* rnd;

  // used by task function
  void book_charge(charge_t q);
  void book_photon(charge_t q);
  void addToGammas(int g);
  void addToIons(int i);
  void find_phaseshift(double e, double* phaseshift);
  double time_update(double tau);
  double angle_function2(double energy);
  double cross_section(double energy, bool &momentum_flag, int &inel_flag);
  double stauffer_elastic_cs(double energy);
  double stauffer_momentum_cs(double energy);
  double func_a(int l);
  double func_b(int l);
  double inelastic_cs(double energy);
  double inelastic_cs_mason_newell(double energy);
  double emp_fit(double x, double* par);
  double* calc_phaseshift1to6(double energy);
  XYZVector speed_update(int charge, XYZPoint dfield, double time);
  XYZVector d_update(XYZVector v0, double time);
  XYZVector kin_factor2(XYZVector v0, bool momentum_flag);


 protected:
  bool run(Electrode* electrode, double en);
  bool taskfunction(Electrode* electrode, charge_t q, double en);

 public:
  // Constructor
  Transport(int seed);
  
  // Default destructor
  ~Transport();

  // Methods
  // preparation, required input from main()
  // otherwise no transport possible
  // work on this electrode id with charges
  int transport(double energy, Electrode* electrode, std::list<charge_t> q); 
  int getPhotons() {return photon_number;}
  int getIons() {return ion_number;}
  double getDensity() {return density;}
  void setDensity(double d) {density = d;};
  std::vector<XYZPoint> allphotons() {return photonStore;}
  std::vector<XYZPoint> allcharges() {return chargeStore;}
};
#endif
