// us
#include "transport.hh"

#include "thread_pool.hpp"

// standard includes
#include <iostream>
#include <cmath>
#include <future>
#include <functional>

// ROOT includes
#include "TMath.h"


//*******
// Transport
//*******
Transport::Transport(Fields* f,int seed) :
  fd(f) // pointer copy
{
  photon_number = 0;
  ion_number = 0;
  anode_number = 0;
  density = 1.6903; // [kg/m^3] fix NTP (295K) argon gas density
  charges.clear();
  chargeStore.clear();
  photonStore.clear();

  pm = new Physics_Model(seed+1); // seed+1 in pm
  generator.seed(seed); // use given seed

}

Transport::~Transport() {
  delete pm;
}


// transport a single charge from a definite location
int Transport::multi_transport(std::list<XYZPoint>& q, double energy) {

  if (q.size()<1) {
    std::cout << "Error: container of charges is empty" << std::endl;
    return 0;
  }

  unsigned int nthreads = std::thread::hardware_concurrency();
  if (nthreads>4) nthreads = 4; // limit max CPU number

  // got all charges as initial input
  charges.clear(); // copy to data member
  for (XYZPoint cc : q) {
    charges.push_front(cc); // insert from front
  }
  
  std::vector<std::future<bool> > results; 
  thread_pool* pool = new thread_pool(nthreads); // task pool

  XYZPoint eloc;
  while (!charges.empty()) { // stop when refilling stopped

    // empty charges and store tasks in blocks of nthreads
    for (int n=0;n<nthreads && !charges.empty();n++) { // drain charges basket
      eloc = charges.front(); // get front element of std::list
      results.push_back(pool->async(std::function<bool(XYZPoint, double)>(std::bind(&Transport::multitask, this, std::placeholders::_1, std::placeholders::_2)), eloc, energy)); // tasks
      charges.pop_front(); // remove first charge from list
    }

    // drain task pool
    for (std::future<bool>& status : results) 
      status.get(); // wait for completion before the next round

    // all tasks from pool finished - clear it. Next batch of charges in pool.
    results.clear();
  }

  // charge loop finished
  delete pool;
  return anode_number; // e- on anode from starter population
}


// transport a single charge from a definite location
int Transport::single_transport(std::list<XYZPoint>& q, double energy) {

  if (q.size()<1) {
    std::cout << "Error: container of charges is empty" << std::endl;
    return 0;
  }

  unsigned int nthreads = std::thread::hardware_concurrency();
  if (nthreads>4) nthreads = 4; // limit max CPU number

  // got all charges as initial input
  charges.clear(); // copy to data member
  for (XYZPoint cc : q) {
    charges.push_front(cc); // insert from front
  }
  
  std::vector<std::future<bool> > results; 
  thread_pool* pool = new thread_pool(nthreads); // task pool

  XYZPoint eloc;
  while (!charges.empty()) { // stop when refilling stopped

    // empty charges and store tasks in blocks of nthreads
    for (int n=0;n<nthreads && !charges.empty();n++) { // drain charges basket
      eloc = charges.front(); // get front element of std::list
      results.push_back(pool->async(std::function<bool(XYZPoint, double)>(std::bind(&Transport::singletask, this, std::placeholders::_1, std::placeholders::_2)), eloc, energy)); // tasks
      charges.pop_front(); // remove first charge from list
    }

    // drain task pool
    for (std::future<bool>& status : results) 
      status.get(); // wait for completion before the next round

    // all tasks from pool finished - clear it. Next batch of charges in pool.
    results.clear();
  }

  // charge loop finished
  delete pool;
  return anode_number; // e- on anode from starter population
}


bool Transport::singletask(XYZPoint point, double en) {
  // have a charge and info about all fields for each thread

  //Init, local storage
  XYZVector speed;
  XYZVector distance_sum, distance_step;
  XYZPoint exyz; // Drift field

  double energy;
  double prob;
  double time_step, running_time;
  double init_energy;
  double speed_start, tangle;
  double kv;
  double x, y, z;
  double xe, ye, ze;
  double csection;

  // init constants
  double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
  double eoverm = 1.759e11; // Coulomb / kg
  double time_sum = 0.0;
  bool momentum_flag = kTRUE;
  int inel_flag = 0;
  int photon_counter = 0;
  int ion_counter = 0;
  
  distance_step.SetXYZ(0.0,0.0,0.0);
  speed.SetXYZ(0.0,0.0,0.0);
  std::uniform_real_distribution<double> rndHalf(0.01, 0.49);  
  
  double e_mass = 0.511e-3; // [GeV/c^2]
  double argon_mass = 39.948*1.0735; // [GeV/c^2]
  double mumass_eV = 1.0e9 * e_mass * argon_mass / (e_mass + argon_mass);
  double localdensity = density * 6.023e26 / 39.948;// convert to number density [m^-3]
  // for E=2.12e8, gives E/N = 10Td = 1.e-16 Vcm^2
  double kmax = 2.e-12; // constant for null coll. method
  double tau = 1/(localdensity * kmax);
  bool stopflag = kFALSE; 

  // speed vector init/ positive z-direction
  init_energy = 1.e-9 * en; // [GeV] 
  speed_start = TMath::Sqrt(2.0*init_energy / e_mass * c2);
  tangle = TMath::Pi()*rndHalf(generator);// isotropic, not full pi
  speed.SetX(speed_start*TMath::Cos(tangle));
  speed.SetY(0.0);
  speed.SetZ(speed_start*TMath::Sin(tangle));
  
  time_sum = running_time = 0.0;
  int geovol = 0; // geometry volume encoding
  bool analytic = false; // default False

  // starting XYZVector from point
  distance_sum.SetXYZ(point.x()*0.01,point.y()*0.01,point.z()*0.01); // [cm]->[m]
  exyz = fd->getFieldValue(point,geovol,analytic); // [V/m]

  // transport loop
  while (!analytic) { 
    // prepare and update
    time_step = -tau * TMath::Log(rnd(generator));
    // free-flight time
    running_time += time_step;
    // keep track of total time
    time_sum += time_step;

    // vector addition stepwise turns velocity vector
    XYZVector Efield(-1*exyz.x(), -1*exyz.y(), -1*exyz.z()); // in [V/m]
    speed += eoverm * Efield * time_step; // in-place update

    // CMS system energy
    energy = 0.5*mumass_eV*speed.Mag2()/c2; // non-rel. energy in [eV]

    csection = pm->cross_section(energy, momentum_flag, inel_flag); // local storage
    kv = speed.R() * csection;
    
    if (inel_flag>0) { // was ionization
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      if (rnd(generator)<0.1) { // 10% recombination prob
      	analytic = true; // Stop
      }
      else { // produce an electron
	book_charge(point); // store in object container
	ion_counter++;
	// std::cout << "inelastic collision electron booked at energy " << energy << std::endl;
	// std::cout << "last field values " << exyz.x() << " "<< exyz.y() << " " << exyz.z() << std::endl;
      }
    }
    else if (inel_flag<0) { // was excitation
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      book_photon(point); // store in object container
      photon_counter++; // count photons
      // std::cout << "inelastic collision photon booked at energy " << energy << std::endl;
      // std::cout << "last field values " << exyz.x() << " "<< exyz.y() << " " << exyz.z() << std::endl;
    }
    
    if (kv>=kmax) {
      std::cout << "kmax too small" << std::endl;
      break;
    }
    
    // random number collision decision
    prob = rnd(generator);
	    
    // collision decision
    if (prob <= (kv/kmax)) {

      // set position of collision
      distance_step = speed * running_time;
      distance_sum += distance_step; // in [m]
      point.SetXYZ(distance_sum.x()*100.0,distance_sum.y()*100.0,distance_sum.z()*100.0); // [cm]

      // new speed from elastic collision kinematics
      kin_factor(speed, momentum_flag);
      
      // check geometry and fields
      exyz = fd->getFieldValue(point,geovol,analytic);
      // std::cout << "analytic bool " << analytic << std::endl;
      // std::cout << "in transport: x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
      // std::cout << "in transport: x,z coordinates " << point.x() << " " << point.z() << std::endl;
      //      std::cout << "collision at energy " << energy << std::endl;
      //      std::cout << "speed X: " << speed.X() << " Z: " << speed.Z() << std::endl;
      // std::cout << "time between coll " << running_time << std::endl;
      if (geovol>1) // e- stopped, now count
	addAnode(); // as hit on anode
      
      // reset system
      running_time = 0.0;
    }
    if (time_sum>=1.0e-7) { // particle got stuck
      analytic = true; // Stop
      std::cout << "STUCK: time = " << time_sum << std::endl;
      std::cout << "STUCK: last x,z coordinates " << point.x() << " " << point.z() << std::endl;
      std::cout << "STUCK: last x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
    }
    if (point.z()>=0.4) { // particle escapes holes, z > 3xhole diameter
      analytic = true; // Stop
      std::cout << "ESCAPE: last x,z coordinates [cm] " << point.x() << " " << point.z() << std::endl;
      std::cout << "ESCAPE: last x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
    }
  }
  // one charge done, book counters to global
  addToGammas(photon_counter);
  addToIons(ion_counter);
  return true;
}


bool Transport::multitask(XYZPoint point, double en) {
  // single difference: do not store interaction locations, just count.
  // have a charge and info about all fields for each thread

  //Init, local storage
  XYZVector speed;
  XYZVector distance_sum, distance_step;
  XYZPoint exyz; // Drift field

  double energy;
  double prob;
  double time_step, running_time;
  double init_energy;
  double speed_start, tangle;
  double kv;
  double x, y, z;
  double xe, ye, ze;
  double csection;

  // init constants
  double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
  double eoverm = 1.759e11; // Coulomb / kg
  double time_sum = 0.0;
  bool momentum_flag = kTRUE;
  int inel_flag = 0;
  int photon_counter = 0;
  int ion_counter = 0;
  
  distance_step.SetXYZ(0.0,0.0,0.0);
  speed.SetXYZ(0.0,0.0,0.0);
  std::uniform_real_distribution<double> rndHalf(0.01, 0.49);  
  
  double e_mass = 0.511e-3; // [GeV/c^2]
  double argon_mass = 39.948*1.0735; // [GeV/c^2]
  double mumass_eV = 1.0e9 * e_mass * argon_mass / (e_mass + argon_mass);
  double localdensity = density * 6.023e26 / 39.948;// convert to number density [m^-3]
  // for E=2.12e8, gives E/N = 10Td = 1.e-16 Vcm^2
  double kmax = 2.e-12; // constant for null coll. method
  double tau = 1/(localdensity * kmax);
  bool stopflag = kFALSE; 

  // speed vector init/ positive z-direction
  init_energy = 1.e-9 * en; // [GeV] 
  speed_start = TMath::Sqrt(2.0*init_energy / e_mass * c2);
  tangle = TMath::Pi()*rndHalf(generator);// isotropic, not full pi
  speed.SetX(speed_start*TMath::Cos(tangle));
  speed.SetY(0.0);
  speed.SetZ(speed_start*TMath::Sin(tangle));
  
  time_sum = running_time = 0.0;
  int geovol = 0; // geometry volume encoding
  bool analytic = false; // default False

  // starting XYZVector from point
  distance_sum.SetXYZ(point.x()*0.01,point.y()*0.01,point.z()*0.01); // [cm]->[m]
  exyz = fd->getFieldValue(point,geovol,analytic); // [V/m]

  // transport loop
  while (!analytic) { 
    // prepare and update
    time_step = -tau * TMath::Log(rnd(generator));
    // free-flight time
    running_time += time_step;
    // keep track of total time
    time_sum += time_step;

    // vector addition stepwise turns velocity vector
    XYZVector Efield(-1*exyz.x(), -1*exyz.y(), -1*exyz.z()); // in [V/m]
    speed += eoverm * Efield * time_step; // in-place update

    // CMS system energy
    energy = 0.5*mumass_eV*speed.Mag2()/c2; // non-rel. energy in [eV]

    csection = pm->cross_section(energy, momentum_flag, inel_flag); // local storage
    kv = speed.R() * csection;
    
    if (inel_flag>0) { // was ionization
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      if (rnd(generator)<0.1) { // 10% recombination prob
      	analytic = true; // Stop
      }
      else { // produce an electron
	ion_counter++;
	// std::cout << "inelastic collision electron booked at energy " << energy << std::endl;
	// std::cout << "last field values " << exyz.x() << " "<< exyz.y() << " " << exyz.z() << std::endl;
      }
    }
    else if (inel_flag<0) { // was excitation
      speed.SetXYZ(0.0,0.0,-1.0); // inelastic takes energy off e-
      kv = 0.0;
      photon_counter++; // count photons
      // std::cout << "inelastic collision photon booked at energy " << energy << std::endl;
      // std::cout << "last field values " << exyz.x() << " "<< exyz.y() << " " << exyz.z() << std::endl;
    }
    
    if (kv>=kmax) {
      std::cout << "kmax too small" << std::endl;
      break;
    }
    
    // random number collision decision
    prob = rnd(generator);
	    
    // collision decision
    if (prob <= (kv/kmax)) {

      // set position of collision
      distance_step = speed * running_time;
      distance_sum += distance_step; // in [m]
      point.SetXYZ(distance_sum.x()*100.0,distance_sum.y()*100.0,distance_sum.z()*100.0); // [cm]

      // new speed from elastic collision kinematics
      kin_factor(speed, momentum_flag);
      
      // check geometry and fields
      exyz = fd->getFieldValue(point,geovol,analytic);
      // std::cout << "analytic bool " << analytic << std::endl;
      // std::cout << "in transport: x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
      // std::cout << "in transport: x,z coordinates " << point.x() << " " << point.z() << std::endl;
      //      std::cout << "collision at energy " << energy << std::endl;
      //      std::cout << "speed X: " << speed.X() << " Z: " << speed.Z() << std::endl;
      // std::cout << "time between coll " << running_time << std::endl;
      if (geovol>1) // e- stopped, now count
	addAnode(); // as hit on anode
      
      // reset system
      running_time = 0.0;
    }
    if (time_sum>=1.0e-7) { // particle got stuck
      analytic = true; // Stop
      std::cout << "STUCK: time = " << time_sum << std::endl;
      std::cout << "STUCK: last x,z coordinates " << point.x() << " " << point.z() << std::endl;
      std::cout << "STUCK: last x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
    }
    if (point.z()>=0.4) { // particle escapes holes, z > 3xhole diameter
      analytic = true; // Stop
      std::cout << "ESCAPE: last x,z coordinates [cm] " << point.x() << " " << point.z() << std::endl;
      std::cout << "ESCAPE: last x,z field values " << exyz.x() << " " << exyz.z() << std::endl;
    }
  }
  // one charge done, book counters to global
  addToGammas(photon_counter);
  addToIons(ion_counter);
  return true;
}


void Transport::book_charge(XYZPoint q) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  charges.push_back(q); // total charge list to be filled/drained in threads
  chargeStore.push_back(q); // permanent storage
  return;
}


void Transport::book_photon(XYZPoint q) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  photonStore.push_back(q); // permanent storage
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

void Transport::addAnode() {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  anode_number++;
  return;
}

void Transport::kin_factor(XYZVector& v0, bool momentum_flag)
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
  double azimuth = TMath::TwoPi()*rnd(generator);
  double theta = pm->angle_function(energy);
  double phi = 0.5*TMath::ASin((argon_mass)/(argon_mass+ e_mass) * TMath::Sin(theta));
  transfer = TMath::Sqrt((1.0 - reduced_mass * TMath::Cos(phi)*TMath::Cos(phi)));
  
  if (momentum_flag) {
    vel.SetTheta(theta+theta0); // relative to old theta
    vel.SetPhi(azimuth+phi0);
  }
  else {
    vel.SetR(vel.R() * transfer); // magnitude change
  }
  XYZVector res(vel);
  v0 = res; // in-place update
}
