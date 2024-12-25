#ifndef PMODEL_HH
#define PMODEL_HH

#include <array>
#include <random>
#include <mutex>

// ROOT
#include "TMath.h"
#include "TGraph.h" // for interpolation

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
  double inelsum(double x, double* par);
  double emp_fit(double x, double* par);
  std::array<double,6> calc_phaseshift1to6(double energy);

  // data members
  TGraph* data0;
  TGraph* data1;
  TGraph* data2;
  TGraph* data3;

  std::mutex mtx;
  std::ranlux24      generator;
  std::uniform_real_distribution<double> rnd;

  
public:
  // Constructor
  Physics_Model(int seed);
  
  // Destructor
  ~Physics_Model();

  // interface
  double cross_section(double energy, bool &momentum_flag, int &inel_flag);
  double angle_function(double energy);

};

// angle functor
class angleFunction {
private:
  double par[5];

public:
  angleFunction(double* params)
  {
    par[0] = params[0];
    par[1] = params[1];
    par[2] = params[2];
    par[3] = params[3];
    par[4] = params[4];
  } // paramter constructor

  double operator () (double x) // Unary for piecewise_linear distribution
  {
    double xx = TMath::Cos(x);
    double k = TMath::Sqrt(par[0]/13.605); // units of 1/a0
    double eta0 = par[1];
    double eta1 = par[2];
    double eta2 = par[3];
    double eta3 = par[4];
    
    double sum1 = 0.0;
    double sum2 = 0.0;
    double p0 = 1.0;
    double p1 = xx;
    double p2 = 0.5*(3.0*xx*xx - 1.0);
    double p3 = 0.5*(5.0*xx*xx*xx - 3.0*xx);
    double alpha = 11.08; // units of a0^3
    double dummy1, dummy2, result;
    
    sum1 += TMath::Cos(eta0)*TMath::Sin(eta0)*p0;
    sum1 += 3.0*TMath::Cos(eta1)*TMath::Sin(eta1)*p1;
    sum1 += 5.0*TMath::Cos(eta2)*TMath::Sin(eta2)*p2;
    sum1 += 7.0*TMath::Cos(eta3)*TMath::Sin(eta3)*p3;
    dummy1 = 1/k * sum1;
    
    sum2 += p1/5.0;
    sum2 += p2/21.0;
    sum2 += p3/45.0;
    dummy2 = TMath::Pi()*alpha*k*(1.0/3.0 -  0.5*TMath::Sin(x/2.0) - sum2);
    dummy1 += dummy2;
    
    sum1 = 0.0;
    sum1 += TMath::Sin(eta0)*TMath::Sin(eta0)*p0;
    sum1 += 3.0*TMath::Sin(eta1)*TMath::Sin(eta1)*p1;
    sum1 += 5.0*TMath::Sin(eta2)*TMath::Sin(eta2)*p2;
    sum1 += 7.0*TMath::Sin(eta3)*TMath::Sin(eta3)*p3;
    dummy2 = 1/k * sum1;
    
    result = dummy1*dummy1 + dummy2*dummy2;
    return result*2.80028561e-21; // convert from Bohr radius^2 to SI
  }

};
#endif
