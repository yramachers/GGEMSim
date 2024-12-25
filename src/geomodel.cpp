// us
#include "geomodel.hh"

// standard includes
#include <iostream>

// ROOT includes
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TString.h"


//****************
// Geometry Model
//****************
// Constructors
void GeometryModel::init(std::string fname) {
  // import geometry from file
  // closes geometry but ID's can be different to initial building
  TGeoManager* g = new TGeoManager("dummy","");
  geom = g->Import(fname.data());

  if (!geom) {
    std::cout << "Could not import " << std::endl;
    return;
  }

}


// default Destructor
GeometryModel::~GeometryModel() {
  if (geom) delete geom;
}


int GeometryModel::whereami(double xv, double yv, double zv) {
  // Talk to geometry
  geom->SetCurrentPoint(xv,yv,zv);
  TGeoNode* nd = geom->FindNode();
  TGeoVolume* vol = nd->GetVolume();
  
  TString region(vol->GetName()); // changed from GetNumber()
  //  std::cout << "Geometry model: region = " << region << std::endl;
  //  std::cout << "Geometry model: coords: " << xv << " " << yv << " " << zv << std::endl;

  // the four considered cases as answers to Where Am I
  if (region.Contains("Comsol")) 
    return 1; // drifting region, field map

  else if (region.Contains("Anode"))
    return 2; // stop, add to anode counter
  
  else if (region.Contains("Plate")) 
    return -1; // stop
  
  else if (region.Contains("World")) 
    return -1; // out, stop transport
  
  else {
    std::cout << "Geometry model: undefined place" << std::endl;
    std::cout << "Geometry model: coords: " << xv << " " << yv << " " << zv << std::endl;
    std::cout << "Geometry model: name: " << region << std::endl;
    return -1; // should never get here, stop
  }
}

