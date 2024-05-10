// *********************************
// THGEM: transport through hole
//**********************************

#include <list>
#include <iostream>
#include <string>
#include <cmath>

// us
#include "transport.hh"
#include "electrode.hh"
#include "fields.hh"
#include "geomodel.hh"
#include "CLI11.hpp"


int main(int argc, char** argv) {
  // function declare
  void signal_calculation(int nstart, int seed, double en, double bias, double xstart, double ystart, double zstart);

    // command line interface
  CLI::App app{"thgem signal calculator"};
  int    seed = 1234;
  int    nstart = 1;
  double en = 1.0;
  double bias = 400.0;
  double xs = 0.055;
  double ys = 0.0;
  double zs = 0.087;
  double dens = 1.6903;
  
  app.add_option("-x,--xstart", xs, "<x start position [cm]> Default: 0.055");
  app.add_option("-y,--ystart", zs, "<y start position [cm]> Default: 0.0");
  app.add_option("-z,--zstart", xs, "<z start position [cm]> Default: 0.087");
  app.add_option("-b,--bias", bias, "<bias value [V]]> Default: 400.0");
  app.add_option("-e,--energy", en, "<initial energy [eV]> Default: 1.0");
  app.add_option("-s,--seed", seed, "<random seed integer> Default: 1234");
  app.add_option("-n,--nstart", nstart, "<number of sims> Default: 1");

  CLI11_PARSE(app, argc, argv);

  //run the code
  signal_calculation(nstart,seed,en,bias,xs,ys,zs);
  
  return 0;
}



void signal_calculation(int nstart, int seed, double en, double bias, double xstart, double ystart, double zstart) {
  //----------------------------------------------------------
  // Hit data
  // THGEM in x,y plane
  // z drift direction to THGEM
  // |  x, 
  // | /
  // |/___ y
  // 
  // MGEM [cm]: +-0.08 shift of plates, plate 2x2 full sides, thickness +- 0.005
  // single charge start on top of THGEM, off right in x from central hole
  // 0.7 mm to x, 1 mum above plate
  charge_t hit;
  XYZPoint loc(xstart, ystart, zstart); // [cm] unit from root geometry
  double qq = -1.0; // [e]
  hit.location = loc;
  hit.charge = qq;
  hit.chargeID = 0;
  std::list<charge_t> hits;
  hits.push_back(hit); // let's have the one

  // distribute rest on line
  double radius = sqrt(xstart*xstart + ystart*ystart);
  double length = 2.0*radius;
  double ival = length/nstart;
  for (int j=1;j < nstart;j++) {
    loc.SetXYZ(-radius + j*ival, 0.0, zstart);
    hit.location = loc;
    hits.push_back(hit); 
  }
    
  //----------------------------------------------------------
  // Geometry
  const char* gfname = "MGEMgeometry.gdml";
  GeometryModel* gmodel = new GeometryModel(gfname);

  //----------------------------------------------------------
  // FEM fields from file
  // hard-coded field map files
  ComsolFields* fem = new ComsolFields("Mgemweight.root");
  fem->setBias(bias);

  //----------------------------------------------------------
  // Transport
  Transport* transportation = new Transport(seed);
  // setting up

  //----------------------------------------------------------
  // transport start
  //----------------------------------------------------------
  fem->read_fields();
  Electrode* anode = new Electrode(fem, gmodel);

  int attempts = transportation->transport(en, anode, hits);

  std::cout << "Total excitation photons counted: " << transportation->getPhotons() << std::endl;
  std::cout << "Total ionization electrons counted: " << transportation->getIons() << std::endl;

  // output for awk filter
  std::cout << "** " << transportation->getPhotons() << std::endl;
  std::cout << "*** " << transportation->getIons() << std::endl;

  delete transportation;
  delete fem;
  delete gmodel;

  return;
}


