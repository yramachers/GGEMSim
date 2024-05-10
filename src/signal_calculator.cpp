// *********************************
// THGEM: transport through hole
//**********************************

#include <list>
#include <iostream>
#include <string>

// us
#include "transport.hh"
#include "electrode.hh"
#include "fields.hh"
#include "geomodel.hh"
#include "CLI11.hpp"

// ROOT
#include "TFile.h"
#include "TNtuple.h"
#include "TParameter.h"


int main(int argc, char** argv) {
  // function declare
  void signal_calculation(int seed, double en, double bias, double xstart, double ystart, double zstart, std::string fname);

  // command line interface
  CLI::App app{"thgem signal calculator"};
  int    seed = 1234;
  double en = 1.0;
  double bias = 400.0;
  double xs = 0.055;
  double ys = 0.0;
  double zs = 0.087;
  std::string outputFileName = "avalanche.root";
  
  app.add_option("-x,--xstart", xs, "<x start position [cm]> Default: 0.055");
  app.add_option("-y,--ystart", zs, "<y start position [cm]> Default: 0.0");
  app.add_option("-z,--zstart", xs, "<z start position [cm]> Default: 0.087");
  app.add_option("-b,--bias", bias, "<bias value [V]]> Default: 400.0");
  app.add_option("-e,--energy", en, "<initial energy [eV]> Default: 1.0");
  app.add_option("-s,--seed", seed, "<random seed integer> Default: 1234");
  app.add_option("-o,--outfile", outputFileName, "<output file name> Default: avalanche.root");

  CLI11_PARSE(app, argc, argv);

  //run the code
  signal_calculation(seed,en,bias,xs,ys,zs,outputFileName);
  
  return 0;
}



void signal_calculation(int seed, double en, double bias, double xstart, double ystart, double zstart, std::string fname) {
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

  std::cout << "attempt: " << attempts << " from 1000" << std::endl;
  std::cout << "Total excitation photons counted: " << transportation->getPhotons() << std::endl;
  std::cout << "Total ionization electrons counted: " << transportation->getIons() << std::endl;

  //----------------------------------------------------------
  // to storage
  //----------------------------------------------------------
  // metainfo
  TParameter<double> xpar("xstart",xstart);
  TParameter<double> ypar("ystart",ystart);
  TParameter<double> zpar("zstart",zstart);
  TParameter<double> bpar("bias",bias);
  TParameter<double> epar("initenergy",en);
  // file
  TFile ff(fname.c_str(),"RECREATE");
  TNtuple* ntcharge = new TNtuple("charge","Ionization charge locations","cx:cy:cz");
  TNtuple* ntgamma = new TNtuple("gamma","photon locations","px:py:pz");
  std::vector<XYZPoint> ac = transportation->allcharges();
  std::vector<XYZPoint> ap = transportation->allphotons();
  for (unsigned int i=0;i<ac.size();i++)
    ntcharge->Fill(ac.at(i).x(),ac.at(i).y(),ac.at(i).z());
  for (unsigned int i=0;i<ap.size();i++)
    ntgamma->Fill(ap.at(i).x(),ap.at(i).y(),ap.at(i).z());

  // store metainfo in both ntuples
  ntcharge->GetUserInfo()->Add(&xpar);
  ntcharge->GetUserInfo()->Add(&ypar);
  ntcharge->GetUserInfo()->Add(&zpar);
  ntcharge->GetUserInfo()->Add(&bpar);
  ntcharge->GetUserInfo()->Add(&epar);
  ntgamma->GetUserInfo()->Add(&xpar);
  ntgamma->GetUserInfo()->Add(&ypar);
  ntgamma->GetUserInfo()->Add(&zpar);
  ntgamma->GetUserInfo()->Add(&bpar);
  ntgamma->GetUserInfo()->Add(&epar);

  ntcharge->Write();
  ntgamma->Write();
  ff.Close();

  delete transportation;
  delete fem;
  delete gmodel;

  return;
}


