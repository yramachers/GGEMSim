// *********************************
// GGEM: transport through holes
//**********************************

#include <vector>
#include <list>
#include <iostream>
#include <string>

// us
#include "transport.hh"
#include "fields.hh"
#include "geomodel.hh"
#include "CLI11.hpp"

// ROOT
#include "TFile.h"
#include "TNtupleD.h"
#include "TParameter.h"


int main(int argc, char** argv) {

  // command line interface
  CLI::App app{"ggem multi charge scan"};
  int    seed = 12345;
  int    nsims = 1;
  double en = 0.2; // [eV]
  double bias = 600.0; // [V]
  double xs = 0.06; // [cm]
  double ys = 0.0;  // [cm]
  double zs = 0.151; // [cm]
  std::string outputFileName = "avalanche.root";
  std::string gdmlName = "DGGEM_5_2_5.gdml";
  std::string fieldName = "2DGGEM_5_2_5-1V.root";
  
  app.add_option("-x,--xstart", xs, "<x start position [cm]> Default: 0.06");
  app.add_option("-y,--ystart", ys, "<y start position [cm]> Default: 0.0");
  app.add_option("-z,--zstart", zs, "<z start position [cm]> Default: 0.151");
  app.add_option("-b,--bias", bias, "<bias value [V]]> Default: 600.0");
  app.add_option("-e,--energy", en, "<initial energy [eV]> Default: 0.2");
  app.add_option("-s,--seed", seed, "<random seed integer> Default: 12345");
  app.add_option("-n,--nsims", nsims, "<number of simulations> Default: 1");
  app.add_option("-o,--outfile", outputFileName, "<output file name> Default: avalanche.root");
  app.add_option("-g,--gdmlfile", gdmlName, "<gdml file name> Default: DGGEM_5_2_5.gdml");
  app.add_option("-f,--fieldfile", fieldName, "<field file name> Default: 2DGGEM_5_2_5-1V.root");

  CLI11_PARSE(app, argc, argv);

  // set up charge population
  std::list<XYZPoint> hits;

  //----------------------------------------------------------
  // Geometry
  GeometryModel gmodel;
  gmodel.init(gdmlName);
  
  //----------------------------------------------------------
  // FEM fields from file
  ComsolFields fem;
  fem.setBias(bias); // set first
  fem.read_fields(fieldName); /// then read and prepare

  Fields field(&gmodel);
  field.prepare_fields(fem); // clears fem vectors at end
  
  //----------------------------------------------------------
  // Transport
  Transport transportation(&field, seed);

  //----------------------------------------------------------
  // transport loop
  //----------------------------------------------------------
  XYZPoint loc(xs, ys, zs); // [cm] unit from root geometry

  // metainfo
  TParameter<double> bpar("bias",bias);
  // file
  TFile ff(outputFileName.data(),"RECREATE");
  TNtupleD* ntcharge = new TNtupleD("charge","charge counters","x,y,z,energy,nions,nphot,nanode");

  // store metainfo in ntuple
  ntcharge->GetUserInfo()->Add(&bpar);
  // Add loops over locations and energies
  for (int count;count<nsims;++count) { // statistics loop
    hits.push_front(loc); // let's have the one

    int anode_count = transportation.multi_transport(hits, en);
    ntcharge->Fill(xs,ys,zs,en,transportation.getIons(),transportation.getPhotons(),anode_count);
    hits.clear();
    transportation.clear_counters();
  }

  ntcharge->Write();
  ff.Close();

  return 0;
}

