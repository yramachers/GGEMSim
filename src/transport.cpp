// us
#include "transport.hh"
#include "thread_pool.hpp"

// standard includes
#include <vector>
#include <iostream>
#include <future>
#include <functional>

// ROOT includes
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"


//*******
// Transport
//*******
Transport::Transport(int seed) {
  photon_number = 0;
  ion_number = 0;
  density = 1.6903; // [kg/m^3] fix NTP (295K) argon gas density
  charges.clear();
  chargeStore.clear();
  photonStore.clear();
  rnd = new TRandom3(seed);
}


Transport::~Transport() {
  delete rnd;
}


// calculate a signal on electrode for any charges in region of interest
int Transport::transport(double energy, Electrode* electrode, std::list<charge_t> q) {

  if (q.size()<1) {
    std::cout << "Error: container of charges is empty" << std::endl;
    return false;
  }
  // got all charges as initial input
  charges.clear(); // copy to data member
  for (charge_t cc : q) {
    charges.push_front(cc); // insert from front
  }
  //  if (!charges.empty()) { // don't run on empty charges
  std::cout << "Transport: number of charges in box: " << charges.size() << std::endl;
  int counter = 0;
  while (!run(electrode, energy) && counter<1000) { // runs first charge
  // next attempt with initial charge
    charges.clear();
    for (charge_t cc : q) {
      charges.push_front(cc); // insert from front
    }
    counter++;
  }
  return counter+1;
}


bool Transport::taskfunction(Electrode* electrode, charge_t q, double en) {
  // have a charge and info about all fields for each thread

  //Init
  XYZVector speed;
  double energy;
  double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
  double time_sum = 0.0;
  bool momentum_flag = kTRUE;
  int inel_flag = 0;
  int photon_counter = 0;
  int ion_counter = 0;
  int check = 0;
  int real_collisions = 0;
  
  XYZVector distance_sum, distance_step;

  distance_step.SetXYZ(0.0,0.0,0.0);
  speed.SetXYZ(0.0,0.0,0.0);
  
  double e_mass = 0.511e-3; // [GeV/c^2]
  double argon_mass = 39.948*1.0735; // [GeV/c^2]
  double mumass_eV = 1.0e9 * e_mass * argon_mass / (e_mass + argon_mass);
  double prob;
  double localdensity = density * 6.023e26 / 39.948;// convert to number density [m^-3]
  // for E=2.12e8, gives E/N = 10Td = 1.e-16 Vcm^2
  
  double time_step, running_time;
  double kmax = 2.e-12; // constant for null coll. method
  double kv;
  double tau = 1/(localdensity * kmax);
  
  double init_energy;
  double speed_start, tangle;
  bool stopflag = kFALSE; 

  // speed vector init/ positive z-direction
  init_energy = 1.e-9 * en; // [GeV] 
  speed_start = TMath::Sqrt(2.0*init_energy / e_mass * c2);
  tangle = TMath::Pi()*rnd->Uniform(0.01,0.49);// isotropic, not full pi
  speed.SetX(speed_start*TMath::Cos(tangle));
  speed.SetY(0.0);
  speed.SetZ(speed_start*TMath::Sin(tangle));
  
  time_sum = running_time = 0.0;
  real_collisions = 0;

  double x, y, z;
  double xe, ye, ze;
  int elcharge;
  XYZPoint exyz; // Drift field
  XYZPoint point;
  bool analytic = false; // default False

  point = q.location; // start location, XYZPoint object; [cm] from ROOT

  // starting XYZVector from point
  distance_sum.SetXYZ(point.x()*0.01,point.y()*0.01,point.z()*0.01); // [cm]->[m]

  elcharge = q.charge; // -1: e-
  charge_t cc;

  exyz = electrode->getFieldValue(analytic,point); // [V/m]

  // transport loop
  while (!analytic) { 
	
    // prepare and update
    time_step = time_update(tau);
    running_time += time_step;
    // keep track of total time
    time_sum += time_step;

    // vector addition stepwise turns velocity vector
    speed += speed_update(elcharge,exyz,time_step);
    
    // CMS system energy
    energy = 0.5*mumass_eV*speed.Mag2()/c2; // non-rel. energy in [eV]
    // artificially raise the cross section 
    kv = speed.R() * cross_section(energy,momentum_flag,inel_flag);

    if (inel_flag>0) { // was ionization
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      if (rnd->Rndm()<0.1) { // 10% recombination prob
      	analytic = true; // Stop
      }
      else { // produce an electron
      	cc.location = point; // last known collision location
	      cc.chargeID = 1; // was an electron
	      cc.charge = -1; //
	      book_charge(cc); // store in object container
	      ion_counter++;
// 	std::cout << "inelastic collision electron booked at energy " << energy << std::endl;
// 	std::cout << "last field values " << exyz.xc() << " "<< exyz.yc() << " " << exyz.zc() << std::endl;
      }
    }
    else if (inel_flag<0) { // was excitation
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      cc.location = point; // last known collision location
      cc.chargeID = 1; // was an electron
      cc.charge = -1; //
      book_photon(cc); // store in object container
      photon_counter++; // count photons
//       std::cout << "inelastic collision photon booked at energy " << energy << std::endl;
//       std::cout << "last field values " << exyz.xc() << " "<< exyz.yc() << " " << exyz.zc() << std::endl;
    }
    
    if (kv>=kmax) {
      std::cout << "kmax too small" << std::endl;
      break;
    }
    
    // random number collision decision
    prob = rnd->Rndm();
	    
    // collision decision
    if (prob <= (kv/kmax)) {

      // book position of collision
      distance_step = d_update(speed,running_time);
      distance_sum += distance_step; // in [m]
      point.SetXYZ(distance_sum.X()*100.0,distance_sum.Y()*100.0,distance_sum.Z()*100.0); // [cm]
//       if (point.zc() < -0.085)  // particle out of hole 1.6mm in z
// 	analytic = 1; // Stop

      // new speed from elastic collision kinematics
      speed = kin_factor2(speed,momentum_flag);
      
      // check geometry and fields
      exyz = electrode->getFieldValue(analytic,point);
      //     std::cout << "analytic bool " << analytic << std::endl;
      //     std::cout << "in transport: x,z field values " << exyz.xc() << " " << exyz.zc() << std::endl;
      //     std::cout << "in transport: x,z coordinates " << point.xc() << " " << point.zc() << std::endl;
//      std::cout << "collision at energy " << energy << std::endl;
//      std::cout << "speed X: " << speed.X() << " Z: " << speed.Z() << std::endl;
//      std::cout << "time between coll " << running_time << std::endl;

      // reset system
      running_time = 0.0;
    }
    if (time_sum>=1.0e-7) { // particle got stuck
      analytic = true; // Stop
      std::cout << "STUCK: time = " << time_sum << std::endl;
    }

  }
  // one charge done, book counters to global
  addToGammas(photon_counter);
  addToIons(ion_counter);
  return true;
}

bool Transport::run(Electrode* electrode, double en) {
  // First, prepare electrode object for transport
  electrode->initfields(); // ready to transport

  unsigned int nthreads = std::thread::hardware_concurrency();
  if (nthreads>4) nthreads = 4; // limit max CPU number
  std::vector<std::future<bool> > results; 
  thread_pool* pool = new thread_pool(nthreads); // task pool

  charge_t q;
  int counter = 0;

  while (!charges.empty()) { // stop when refilling stopped
//    std::cout << "from threads, charge basket size = " << charges.size() << std::endl;

    // empty charges and store tasks in blocks of nthreads
    for (int n=0;n<nthreads && !charges.empty();n++) { // drain charges basket
      q = charges.front(); // get front element of std::list
      results.push_back(pool->async(std::function<bool(Electrode*, charge_t, double)>(std::bind(&Transport::taskfunction, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3)), electrode, q, en)); // tasks
      charges.pop_front(); // remove first charge from list
      counter++; // counts tasks/electrons launched
    }

      //    std::future<bool> status = pool.async(std::function<bool(Electrode*, charge_t)>(std::bind(&Transport::taskfunction, this, std::placeholders::_1, std::placeholders::_2)), electrode, q); // tasks in pool
    // drain task pool
    for (std::future<bool>& status : results) 
      status.get(); // wait for completion before the next round
    
      //	counter++; // counts tasks launched
    //    }
    // all tasks from pool finished - clear it. Next batch of charges in pool.
    results.clear();
  }

  //  std::cout << "from threads, total task counter = " << counter << std::endl;
  // charge loop finished
  delete pool;
  if (getPhotons()>0 || counter>1)
    return true;
  else
    return false;
}


void Transport::book_charge(charge_t q) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  charges.push_back(q); // total charge list to be filled/drained in threads
  chargeStore.push_back(q.location); // permanent storage
  return;
}


void Transport::book_photon(charge_t q) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  photonStore.push_back(q.location); // permanent storage
  return;
}


void Transport::addToGammas(int g) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  photon_number += g;
  return;
}


void Transport::addToIons(int i) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  ion_number += i;
  return;
}

double Transport::time_update(double tau)
{
    return -tau*TMath::Log(rnd->Rndm());
}

XYZVector Transport::speed_update(int charge, XYZPoint dfield, double time)
{
  XYZVector Efield(charge*dfield.x(), charge*dfield.y(), charge*dfield.z()); // in [V/m]
  double eoverm = 1.759e11; // Coulomb / kg
  XYZVector v = eoverm * Efield * time;
  return v;
}

XYZVector Transport::d_update(XYZVector v0, double time)
{
    XYZVector dstep = v0*time; // acceleration done already in speed_update
    return dstep;
}

XYZVector Transport::kin_factor2(XYZVector v0, bool momentum_flag)
{
  Polar3D<double> vel(v0);
    double theta0 = v0.Theta();
    double phi0 = v0.Phi();
    double energy, transfer;
    double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
    double argon_mass = 39.95*1.0735; // [GeV/c^2]
    double e_mass = 0.511e-3; // [GeV/c^2]
    double reduced_mass = (4.0*argon_mass*e_mass)/((argon_mass + e_mass)*(argon_mass + e_mass));
    double mumass_eV = 1.0e9 * e_mass * argon_mass / (e_mass + argon_mass);

    energy = 0.5 * mumass_eV * v0.Mag2() / c2; // non-rel. energy in [eV]
    double azimuth = TMath::TwoPi()*rnd->Rndm();
    double theta = angle_function2(energy);
    double phi = 0.5*TMath::ASin((argon_mass + e_mass)/argon_mass*TMath::Sin(theta));
    transfer = TMath::Sqrt((1.0 - reduced_mass * TMath::Cos(phi)*TMath::Cos(phi)));

    if (momentum_flag) {
      vel.SetTheta(theta+theta0); // relative to old theta
      vel.SetPhi(azimuth+phi0);
      XYZVector res(vel);
      return res; // return XYZVector type
    }
    else {
      vel.SetR(vel.R() * transfer); // magnitude change
      XYZVector res(vel);
      return res;
    }
}

//****************************
// theory from JPhysB16,4023
//****************************
double Transport::angle_function2(double energy)
{
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  double dcs(double* x, double* par);
  if (energy <= 0.02) // trial correction
    return TMath::Pi()*rnd->Rndm();//isotropic for <0.02 eV
  TF1* ff = new TF1("name",dcs,0.0,TMath::Pi(),5);
  double p[5];
  double phaseshift[4];
  find_phaseshift(energy,phaseshift);
  
  p[0] = energy; 
  for (Int_t i=0;i<4;i++) {
    p[i+1] = phaseshift[i];
  }
  ff->SetParameters(p);
  double theta = ff->GetRandom();
  //  cout << "theta=" << theta << " at [eV] " << energy << endl;
  delete ff;
  return theta;
}

double dcs(double* x, double* par)
{
  double xx = TMath::Cos(x[0]);
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
  dummy2 = TMath::Pi()*alpha*k*(1.0/3.0 -  0.5*TMath::Sin(x[0]/2.0) - sum2);
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

void Transport::find_phaseshift(double e, double* phaseshift)
{
    if (e<0.02) e = 0.02; // minimum energy value fixed

    TGraph* data0 = new TGraph();
    TGraph* data1 = new TGraph();
    TGraph* data2 = new TGraph();
    TGraph* data3 = new TGraph();
    double energy[11] = {0.02,0.06,0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0,3.0};
    double values[44] = {4.426e-2,3.092e-3,5.836e-4,1.537e-4,5.363e-2,8.155e-3,1.577e-3,5.21e-4,5.133e-2,1.219e-2,2.686e-3,8.937e-4,3.255e-2,1.93e-2,5.516e-3,1.831e-3,8.595e-3,2.307e-2,8.379e-3,2.756e-3,-1.639e-2,2.46e-2,1.128e-2,3.667e-3,-4.122e-2,2.444e-2,1.419e-2,4.574e-3,-8.931e-2,2.04e-2,2.016e-2,6.374e-3,-1.566e-1,8.199e-3,2.927e-2,9.059e-3,-3.48e-1,-5.612e-3,6.446e-2,1.788e-2,-5.057e-1,-1.329e-1,1.111e-1,2.669e-2};
    
    Int_t m = 11, n = 4;
    Int_t i;
    for (i=0;i<m;i++) {
    	data0->SetPoint(i,energy[i],values[i*n]);
	    data1->SetPoint(i,energy[i],values[1+i*n]);
	    data2->SetPoint(i,energy[i],values[2+i*n]);
	    data3->SetPoint(i,energy[i],values[3+i*n]);
    }
    phaseshift[0] = data0->Eval(e);
    phaseshift[1] = data1->Eval(e);
    phaseshift[2] = data2->Eval(e);
    phaseshift[3] = data3->Eval(e);

    delete data3;
    delete data2;
    delete data1;
    delete data0;

}


// check new cross sections better than Eachran, Stauffer (1983)
double Transport::cross_section(double energy, bool &momentum_flag, int &inel_flag)
{
  double cs_el, cs_inel,ratio;
  double cs_etrans;
  double cs_ptrans;
  
  inel_flag = 0; // default case
  momentum_flag = kTRUE; // default case
  
  // purely elastic scattering possible
  if (energy>0.0 && energy<=11.55) {
    cs_ptrans = stauffer_momentum_cs(energy);
    cs_etrans = stauffer_elastic_cs(energy);
    ratio = cs_ptrans / (cs_etrans+cs_ptrans);

    if (energy>4.5 && energy<=11.55) { // overestimated sigma
      cs_ptrans *= 0.73;
      cs_etrans *= 0.73;
    }
    if (rnd->Rndm() <  ratio) momentum_flag = kTRUE;
    else momentum_flag = kFALSE; // energy transfer cs
    return (momentum_flag ? cs_ptrans : cs_etrans);
  }
  // otherwise decide
  else {	
    cs_ptrans = stauffer_momentum_cs(energy);
    cs_el = stauffer_elastic_cs(energy);
    cs_inel = inelastic_cs(energy);
    ratio = cs_inel / (cs_el+cs_inel);
    if (rnd->Rndm() < ratio) {
      inel_flag = -1; // photon
      if (energy > 15.76) inel_flag = 1; // electron
      //      std::cout << " inelastic sigma = " << cs_inel << " ratio = " << ratio << std::endl;
      return 0.0; // stops further transport after inel
    }
    // elastic as normal
    ratio = cs_ptrans / cs_etrans;
    if (rnd->Rndm() <  ratio) momentum_flag = kTRUE;
    else momentum_flag = kFALSE; // energy transfer cs
    return (momentum_flag ? cs_ptrans+cs_inel : cs_etrans+cs_inel);
  }
}


//****************************
//****************************
// Stauffer formulae
//****************************
double Transport::stauffer_elastic_cs(double energy)
{
    int l;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    double sigma, factor, sum = 0.0;

    double* delta = calc_phaseshift1to6(energy);

    factor = 4.0*TMath::Pi()/(k*k);
    for (l=0;l<6;l++) 
    	sum += (2.0*l+1)*TMath::Sin(delta[l])*TMath::Sin(delta[l]);

    sigma = factor * sum;
    
    delete [] delta;
    // unit conversion Bohr radius^2 to SI
    sigma *= 2.80028561e-21;

    return sigma;
}

double Transport::stauffer_momentum_cs(double energy)
{
    int l;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    double sigma, factor, sum = 0.0;

    double* delta = calc_phaseshift1to6(energy);

    factor = 4.0*TMath::Pi()/(k*k);
    for (l=0;l<5;l++) 
    	sum += (l+1)*TMath::Sin(delta[l]-delta[l+1])*TMath::Sin(delta[l]-delta[l+1]);
    
    sigma = factor * sum;
    
    delete [] delta;
    // unit conversion Bohr radius^2 to SI
    sigma *= 2.80028561e-21;

    return sigma;
}

double* Transport::calc_phaseshift1to6(double energy)
{
    int l;
    double beta_argon = 10.755;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    double* delta = new double [6];
    // Energy input in [eV]
    // output cross section in m^2
    double h[5][4] = {{0.150575e1,0.0,0.0,0.0},{0.671099e1,0.225686e1,0.319699,0.107419},{-0.103246e3,-0.805571e1,-0.471566,0.0},{0.11987e2,0.213514e1,0.662807,-0.183263},{0.0,0.0,0.0,0.315322}};


    double d[4][4] = {{0.122028e2,-0.80281e-2,-0.145269e1,0.0},{0.429028e2,0.425572e1,0.201048,-0.177524e1},{0.547674e1,-0.677844,0.584853,0.293544e1},{0.0,0.0,0.0,-0.497051}};

    double p4[2] = {0.324336e-2,-0.867549e-3};
    double p5[2] = {0.382582e-3,0.0};
    
    double N[4], D[4];
    N[0] = h[0][0]*k + h[1][0]*k*k + h[2][0]*k*k*k + h[3][0]*k*k*k*k;
    for (l=1;l<4;l++) 
    	N[l] = h[1][l]*k*k + h[2][l]*k*k*k + h[3][l]*k*k*k*k + h[4][l]*k*k*k*k*k;
    
    D[0] = 1.0 + d[0][0]*k + d[1][0]*k*k + d[2][0]*k*k*k + d[3][0]*k*k*k*k;;
    D[1] = 1.0 + d[0][1]*k + d[1][1]*k*k + d[2][1]*k*k*k*TMath::Log(k);
    for (l=2;l<4;l++) 
    	D[l] = 1.0 + d[0][l]*k + d[1][l]*k*k + d[2][l]*k*k*k + d[3][l]*k*k*k*k;

    for (l=0;l<4;l++) 
	    delta[l] = N[l] / D[l];
    
    delta[4] = k*k*beta_argon*func_a(4) + k*k*k*k*beta_argon*beta_argon*func_b(4) + p4[0]*k*k*k*k*k*k + p4[1]*k*k*k*k*k*k*k;

    delta[5] = k*k*beta_argon*func_a(5) + k*k*k*k*beta_argon*beta_argon*func_b(5) + p5[0]*k*k*k*k*k*k;
    return delta;
}

double Transport::func_a(int l)
{
    double dummy = (2.0*l-1)*(2.0*l+1)*(2.0*l+3);
    return TMath::Pi()/dummy;
}

double Transport::func_b(int l)
{
    double dummy = (2.0*l-1)*(2.0*l+1)*(2.0*l+3);
    double dummy1 = (2.0*l-3)*dummy*dummy*dummy*(2.0*l+5);
    double dummy2 = (2.0*l+1)*(2.0*l+1);
    return TMath::Pi()*(15.0*dummy2*dummy2 - 140.0*dummy + 128.0)/dummy1;
}

//**************************
//** Inelastic from fits in NIM A242, 327
//**************************

double Transport::inelastic_cs(double energy)
{
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  double inelsum(double* x, double* par);
  int flag;

  double bins[6] = {11.55,11.62,11.72,11.83,14.0,15.76};

  if (energy <= bins[0]) return 0.0;
  if (energy <= bins[5]) return inelastic_cs_mason_newell(energy);

  double array[7][4] = {{-0.4087e-2,0.5321e-3,-0.1543e-4,0.0},{-0.1775e-2,0.2037e-3,-0.4394e-5,0.0},{-0.6713e-3,0.8463e-4,-0.2335e-5,0.0},{-0.5562e-2,0.7999e-3,-0.3438e-4,0.551e-6},{-0.2844e-1,0.3025e-2,-0.7096e-4,0.0},{0.681e1,-0.9458,0.3259e-1,0.0},{-9.884,0.1894,0.0,0.0}};

  if (energy>bins[5] && energy <=17.5) flag = 6;
  if (energy>17.5) flag = 7;

  TF1* func = new TF1("inelsum",inelsum,11.55,50.0,4);
  double par[4], result = 0.0;

  for (int j=1;j<flag;j++) { // all previous
  	par[0] = array[j-1][0];
  	par[1] = array[j-1][1];
	  par[2] = array[j-1][2];
	  par[3] = array[j-1][3];
	  func->SetParameters(par);
	  result += func->Eval(bins[j]-1.0e-3);
  }
  par[0] = array[flag-1][0]; // + actual
  par[1] = array[flag-1][1];
  par[2] = array[flag-1][2];
  par[3] = array[flag-1][3];
  func->SetParameters(par);
  result += func->Eval(energy);
  // unit conversion pi*Bohr radius^2 to SI
  result *= 8.797356696e-21;
  result += inelastic_cs_mason_newell(energy); // excitation + ionization
  //
  delete func;
  return result;        
}

double inelsum(double* x, double* par)
{
    double result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
    return result;
}

double Transport::inelastic_cs_mason_newell(double energy)
{
    double par[7] = {34.7531,-10.2377,1.10632,-0.0562047,0.00147617,-1.94974e-5,1.02707e-7};
    double result = emp_fit(energy,par);
    return result*1.0e-21; // SI units
}

double Transport::emp_fit(double x, double* par)
{
    // pol6 fits by 190.0/20 chi2/ndf up to 50 eV, not further
    return par[0]+par[1]*x+par[2]*x*x+par[3]*x*x*x+par[4]*x*x*x*x+par[5]*x*x*x*x*x+par[6]*x*x*x*x*x*x;
}


