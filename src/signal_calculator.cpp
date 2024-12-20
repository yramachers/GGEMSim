// *********************************
// THGEM: transport through hole
//**********************************

#include <vector>
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
  void signal_calculation(int, const XYZPoint&, const std::vector<double>&, const std::vector<std::string>&);

  // command line interface
  CLI::App app{"thgem signal calculator"};
  int    seed = 1234;
  double en = 0.2; // [eV]
  double bias = 600.0; // [V]
  double xs = 0.06; // [cm]
  double ys = 0.0;  // [cm]
  double zs = 0.151; // [cm]
  std::vector<std::string> names;
  std::vector<double> config;
  std::string outputFileName = "avalanche.root";
  std::string gdmlName = "DGGEM_5_2_5.gdml";
  std::string fieldName = "2DGGEM_5_2_5-1V.root";
  
  app.add_option("-x,--xstart", xs, "<x start position [cm]> Default: 0.06");
  app.add_option("-y,--ystart", zs, "<y start position [cm]> Default: 0.0");
  app.add_option("-z,--zstart", xs, "<z start position [cm]> Default: 0.151");
  app.add_option("-b,--bias", bias, "<bias value [V]]> Default: 600.0");
  app.add_option("-e,--energy", en, "<initial energy [eV]> Default: 0.2");
  app.add_option("-s,--seed", seed, "<random seed integer> Default: 1234");
  app.add_option("-o,--outfile", outputFileName, "<output file name> Default: avalanche.root");
  app.add_option("-g,--gdmlfile", gdmlName, "<gdml file name> Default: DGGEM_5_2_5.gdml");
  app.add_option("-f,--fieldfile", fieldName, "<field file name> Default: 2DGGEM_5_2_5-1V.root");

  CLI11_PARSE(app, argc, argv);

  //run the code
  config.push_back(en);
  config.push_back(bias);
  names.push_back(outputFileName);
  names.push_back(gdmlName);
  names.push_back(fieldName);
  XYZPoint loc(xs, ys, zs); // [cm] unit from root geometry
  
  signal_calculation(seed, loc, config, names);
  
  return 0;
}



void signal_calculation(int seed, const XYZPoint& loc, const std::vector<double>& conf, const std::vector<std::string>& names) {
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
  double qq = -1.0; // [e]
  hit.location = loc;
  hit.charge = qq;
  hit.chargeID = 0;
  std::list<charge_t> hits;
  hits.push_back(hit); // let's have the one

  //----------------------------------------------------------
  // Geometry
  GeometryModel* gmodel = new GeometryModel(names.at(1));

  //----------------------------------------------------------
  // FEM fields from file
  // hard-coded field map files
  ComsolFields* fem = new ComsolFields(names.at(2));
  fem->setBias(conf.at(1));

  //----------------------------------------------------------
  // Transport
  Transport* transportation = new Transport(seed);
  // setting up

  //----------------------------------------------------------
  // transport start
  //----------------------------------------------------------
  fem->read_fields();
  Electrode* anode = new Electrode(fem, gmodel);

  int attempts = transportation->transport(conf.at(0), anode, hits);

  std::cout << "attempt: " << attempts << " from 1000" << std::endl;
  std::cout << "Total excitation photons counted: " << transportation->getPhotons() << std::endl;
  std::cout << "Total ionization electrons counted: " << transportation->getIons() << std::endl;

  //----------------------------------------------------------
  // to storage
  //----------------------------------------------------------
  // metainfo
  TParameter<double> xpar("xstart",loc.x());
  TParameter<double> ypar("ystart",loc.y());
  TParameter<double> zpar("zstart",loc.z());
  TParameter<double> bpar("bias",conf.at(1));
  TParameter<double> epar("initenergy",conf.at(0));
  // file
  TFile ff(names.at(0).c_str(),"RECREATE");
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


