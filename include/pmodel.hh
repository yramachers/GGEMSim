#ifndef PMODEL_HH
#define PMODEL_HH

#include <array>
// #include <mutex>

// ROOT
#include "TRandom3.h"
#include "TF1.h"
#include "TGraph.h"

using namespace ROOT::Math;

// Physics Model
// cross section calculation
class Physics_Model {

private:
  // operation functions
  void find_phaseshift(double e, double* phaseshift);
  double stauffer_elastic_cs(double energy);
  double stauffer_momentum_cs(double energy);
  double func_a(int l);
  double func_b(int l);
  double inelastic_cs(double energy);
  double inelastic_cs_mason_newell(double energy);
  double emp_fit(double x, double* par);
  std::array<double,6> calc_phaseshift1to6(double energy);

  // data members
  TRandom3* rnd;
  TF1* afunc;
  TF1* inel;
  TGraph* data0;
  TGraph* data1;
  TGraph* data2;
  TGraph* data3;
  
public:
  // Constructor
  Physics_Model();
  
  // Destructor
  ~Physics_Model();

  // interface
  double cross_section(double energy, bool &momentum_flag, int &inel_flag);
  double angle_function2(double energy);

};
#endif
