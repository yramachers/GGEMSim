// Build a double GGEM in ROOT geometry and store as GDML file.
// single hole to single hole assumption.
// could be 0.8 mm diameter twice or 1.2 mm on second plate
//
#include <iostream>
// ROOT items
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoTorus.h"
#include "TGeoVolume.h"
#include "TGeoCompositeShape.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"

int main() {
  void DGGEM08mm();

  // simply call geometry function
  DGGEM08mm();
  return 0;
}


void DGGEM08mm()
{
  TGeoManager* geom = new TGeoManager("ggemtest","GGEM geometry");

  // Materials irrelevant but required for GDML file
  TGeoMaterial* mat = new TGeoMaterial("Vacuum",0,0,1.e-20);
  TGeoMaterial* air = new TGeoMaterial("Ar",40,18,1.e-20);
  TGeoMaterial *matAl = new TGeoMaterial("Al",26,13,2.7);

  TGeoMedium* med = new TGeoMedium("Vacuum",1,mat);
  TGeoMedium* al = new TGeoMedium("Al",2,matAl);
  TGeoMedium* airmed = new TGeoMedium("Ar",3,air);

  // units [cm]
  // note hard-coded names for volumes required
  double worldhalf = 2.0; // [cm]
  double halfside = 0.5; // [cm]
  double halfheight = 0.01; // height=200mum plate thickness
  double radius = 0.04; // r=0.4 mm
  double circlerad = 0.01; // 100 mum
  double shift = 0.025; // 0.5 mm plate distance
  double halftransfer = 0.06; // 1.2 mm transfer gap
  
  // make shape components
  TGeoBBox *sbox  = new TGeoBBox("P",halfside,halfside,halfheight);
  TGeoTube *stub  = new TGeoTube("H",0,radius,halfheight);


  TGeoVolume* world = geom->MakeBox("World",med,worldhalf+0.1,worldhalf+0.1,worldhalf+0.1); // larger than rest
  TGeoVolume* comsol = geom->MakeBox("Comsol",airmed,worldhalf,worldhalf,worldhalf); // field volume

  TGeoRotation* ident = new TGeoRotation("id",0,0,0); // nothing
  ident->RegisterYourself();

  TGeoCompositeShape *cs = new TGeoCompositeShape("cs",
						  "(P-H:id)");
  TGeoVolume *plate = new TGeoVolume("Plate",cs,al);


  geom->SetTopVolume(world);

  // plate volume shifted
  world->AddNode(comsol,1);

  comsol->AddNode(plate,1,new TGeoTranslation(0.0,0.0,2*shift+2*halfheight+halftransfer));
  comsol->AddNode(plate,2,new TGeoTranslation(0.0,0.0,halfheight+halftransfer));
  comsol->AddNode(plate,3,new TGeoTranslation(0.0,0.0,-halfheight-halftransfer));
  comsol->AddNode(plate,4,new TGeoTranslation(0.0,0.0,-2*shift-2*halfheight-halftransfer));

  world->SetLineColor(1);
  comsol->SetLineColor(1);
  plate->SetLineColor(2);
  geom->SetTopVisible();

  geom->CloseGeometry();
  // *** GEOMETRY closed ***

  //  geom->Export("DGGEM_08_12.gdml");
  geom->SetVisLevel(4);
  geom->SetVisOption(0);
  world->Draw();

}


