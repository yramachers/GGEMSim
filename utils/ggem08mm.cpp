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
  void GGEM08mm();

  // simply call geometry function
  GGEM08mm();
  return 1;
}


void GGEM08mm()
{
  TGeoManager* geom = new TGeoManager("ggemtest","GGEM geometry");

  TGeoMaterial* mat = new TGeoMaterial("Vacuum",0,0,1.e-20);
  TGeoMaterial* air = new TGeoMaterial("Ar",40,18,1.e-20);
  TGeoMaterial *matAl = new TGeoMaterial("Al",26,13,2.7);

  TGeoMedium* med = new TGeoMedium("Vacuum",1,mat);
  TGeoMedium* al = new TGeoMedium("Al",2,matAl);
  TGeoMedium* airmed = new TGeoMedium("Ar",3,air);

  // units [cm]
  // note hard-coded names for volumes required
  double worldhalf = 5.0; // 1cm
  double halfside = 1.0; // 1cm
  double halfheight = 0.0175; // height=350mum
  double radius = 0.0575; // 0.575mm
  double circlerad = 0.0175; // 175mum
  double shift = 0.0575; // 0.575mm
  //  double diag = 0.141*0.8; // sqrt(2)*0.8mm diam for pitch

  // make shape components
  TGeoBBox *sbox  = new TGeoBBox("P",halfside,halfside,halfheight);
  TGeoTube *stub  = new TGeoTube("H",0,radius,halfheight);
  TGeoTorus* edge = new TGeoTorus("Edge",radius, 0.0, circlerad, 0.0, 360.0); // circular hole edge


  TGeoVolume* world = geom->MakeBox("World",med,worldhalf+0.1,worldhalf+0.1,worldhalf+0.1); // larger than rest
  TGeoVolume* comsol = geom->MakeBox("Comsol",airmed,worldhalf,worldhalf,worldhalf); // field volume

  TGeoRotation* ident = new TGeoRotation("id",0,0,0); // nothing
  ident->RegisterYourself();

//   TGeoTranslation* mv1 = new TGeoTranslation("mv1",diag,diag,0.0);
//   mv1->RegisterYourself();
//   TGeoTranslation* mv2 = new TGeoTranslation("mv2",diag,-diag,0.0);
//   mv2->RegisterYourself();
//   TGeoTranslation* mv3 = new TGeoTranslation("mv3",-diag,diag,0.0);
//   mv3->RegisterYourself();
//   TGeoTranslation* mv4 = new TGeoTranslation("mv4",-diag,-diag,0.0);
//   mv4->RegisterYourself();

  // create a composite plate with 5 holes
//   TGeoCompositeShape *cs = new TGeoCompositeShape("cs",
// 						  "(P-H:id-H:mv1-H:mv2-H:mv3-H:mv4)");
  TGeoCompositeShape *cs = new TGeoCompositeShape("cs",
						  "(P-H:id)+Edge:id");
  TGeoVolume *plate = new TGeoVolume("Plate",cs,al);


  geom->SetTopVolume(world);

  // plate volume shifted
  world->AddNode(comsol,1);

  comsol->AddNode(plate,1,new TGeoTranslation(0.0,0.0,shift));
  comsol->AddNode(plate,2,new TGeoTranslation(0.0,0.0,-shift));

  world->SetLineColor(1);
  comsol->SetLineColor(1);
  plate->SetLineColor(2);
  geom->SetTopVisible();

  geom->CloseGeometry();
  // *** GEOMETRY closed ***

  geom->Export("GGEMgeometry.gdml");

  //  world->Draw();

}


