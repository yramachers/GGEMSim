// us
#include "geomodel.hh"

// standard includes
#include <iostream>

// ROOT includes
#include "TGeoVolume.h"
#include "TString.h"
#include "TObjArray.h"


//****************
// Geometry Model
//****************
// Constructors
GeometryModel::GeometryModel(const char* filename) {

  // import geometry from file
  // closes geometry but ID's can be different to initial building
  TGeoManager* g = new TGeoManager("dummy","");
  geom = g->Import(filename);

  if (!geom) {
    std::cout << "Could not import " << std::endl;
    return;
  }
  else 
    fill_plates(); // only constructed when geometry is known from file
}


// default Destructor
GeometryModel::~GeometryModel() {
  plates.clear();
  if (geom) delete geom;
}



void GeometryModel::fill_plates() {
  TGeoVolume* vol = geom->FindVolumeFast("Comsol"); // hard-wired region name
  TObjArray* lon = vol->GetNodes(); // should all be in comsol region
  
  TString name;
  
  for (int n=0;n<lon->GetEntries();n++) {
    name = lon->At(n)->GetName();
    if (name.Contains("Plate"))  // hard-wired plate name
      plates.push_back((TGeoNode*)lon->At(n));
  }
  std::cout << "Geometry model; got "<< plates.size() << " electrodes" << std::endl;

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

  else if (region.Contains("Plate")) {
    return -1; // stop
  }
  else if (region.Contains("World")) {
    return -1; // out, stop transport
  }
  else {
    std::cout << "Geometry model: undefined place" << std::endl;
    std::cout << "Geometry model: coords: " << xv << " " << yv << " " << zv << std::endl;
    std::cout << "Geometry model: name: " << region << std::endl;
    return -1; // should never get here, stop
  }
}

