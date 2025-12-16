#include <iostream>

// us
#include "fields.hh"

// ROOT includes
#include "TFile.h"
#include "TNtupleD.h"
#include "TMath.h"
#include "Math/Vector3D.h"


// as 2D data in x,y from comsol
void ComsolFields::read_fields(std::string fn) {
  coords.clear();
  dmap.clear();
  
  XYZPoint p;
  XYZPoint pe;

  TFile* ffd = new TFile(fn.data(),"read");
  TNtupleD* ntd = (TNtupleD*)ffd->Get("weighting");
  int entries = ntd->GetEntries();

  double x,y;
  double wx,wy;

  ntd->SetBranchAddress("x",&x);
  ntd->SetBranchAddress("y",&y);

  ntd->SetBranchAddress("wx",&wx);
  ntd->SetBranchAddress("wy",&wy);


  // comsol y-coord becomes z-coord in geometry
  for (int i=0;i<entries;i++){
    ntd->GetEntry(i);
    p.SetXYZ(x*1.0e2, 0.0, y*1.0e2); // [m]->[cm]
    coords.push_back(p);
    pe.SetXYZ(bias*wx, 0.0, bias*wy); // 2 electrode weighting to drift
    dmap.push_back(pe);
  }
  std::cout << "in Comsol Fields: read field entries " << entries << std::endl;
  // all done and in memory
  ffd->Close();
}


void ComsolFields::clear()
{
  coords.clear();
  dmap.clear();  
}


// read 3D data from comsol
void ComsolFields3D::read_fields(std::string fn) {
  coords.clear();
  dmap.clear();
  
  XYZPoint p;
  XYZPoint pe;

  TFile* ffd = new TFile(fn.data(),"read");
  TNtupleD* ntd = (TNtupleD*)ffd->Get("weighting");
  int entries = ntd->GetEntries();

  double x,y,z;
  double wx,wy,wz;

  ntd->SetBranchAddress("x",&x);
  ntd->SetBranchAddress("y",&y);
  ntd->SetBranchAddress("z",&z);

  ntd->SetBranchAddress("wx",&wx);
  ntd->SetBranchAddress("wy",&wy);
  ntd->SetBranchAddress("wz",&wz);


  // comsol y-coord becomes z-coord in geometry
  for (int i=0;i<entries;i++){
    ntd->GetEntry(i);
    p.SetXYZ(x*1.0e2, y*1.0e2, z*1.0e2); // [m]->[cm]
    coords.push_back(p);
    pe.SetXYZ(wx, wy, wz); // drift electric field strength V/m
    dmap.push_back(pe);
  }
  std::cout << "in Comsol Fields: read field entries " << entries << std::endl;
  // all done and in memory
  ffd->Close();
}


void ComsolFields3D::clear()
{
  coords.clear();
  dmap.clear();  
}


// Fields class
Fields:: Fields(GeometryModel* g) :
  gm(g) // copy pointer
{}


Fields::~Fields() {
  if (allx) { // all set together
    delete [] allx ;
    delete [] allz ;
    delete [] alldx;
    delete [] alldz;
  }
  if (coordinates) delete coordinates;
}

void Fields::prepare_fields(ComsolFields& fem) {
  std::vector<XYZPoint> cdata = fem.positions();
  int nentries = cdata.size();

  coordinates = new TKDTreeID(nentries,2,1);
  
  allx = new double [nentries];
  allz = new double [nentries];

  for (int i=0;i<nentries;i++){
    allx[i] = cdata[i].x();
    allz[i] = cdata[i].z();
  }
  coordinates->SetData(0,allx);
  coordinates->SetData(1,allz);
  coordinates->Build();
  // KDTree built, all in memory

  alldx = new double [nentries];
  alldz = new double [nentries];

  std::vector<XYZPoint> ddata = fem.driftmap();
  for (int i=0;i<nentries;i++){
    alldx[i] = ddata[i].x();
    alldz[i] = ddata[i].z();
  }  
  //  std::cout << "in Fields::prepare fields finished." << std::endl;
  // all done and in memory
  fem.clear();
}


XYZPoint Fields::getFieldValue(XYZPoint& p, int& gv, bool& analytic)
{
  // thread access protection
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  
  // ask the geometry
  double xv = p.x();
  double yv = p.y();
  double zv = p.z();
  double rad = TMath::Sqrt(xv*xv+yv*yv); // x-y-plane
  double angle = TMath::ATan2(yv,xv); // x-y-plane

  int value = gm->whereami(xv,yv,zv); // single access point to geomodel
  XYZPoint triplet;

  if (value==2) {
    gv = 2; // anode encoding
    triplet.SetXYZ(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  else if (value==1) { // comsol region
    gv = 1; // drift region code
    double point[2];
    double dist[5]; // check on nearest 5 neighbours in grid
    int indx[5];
    
    XYZVector fieldvec;
    XYZVector sumvec;
    std::vector<XYZVector> nnvec;
    
    double dsum = 0.0;
    point[0] = rad;  // relative to origin x
    point[1] = zv;  // relative to z-axis origin

    coordinates->FindNearestNeighbors(point,5,indx,dist);
    for (int j=0;j<5;j++) {
      fieldvec.SetXYZ(alldx[indx[j]], 0.0, alldz[indx[j]]);
      nnvec.push_back(fieldvec);
      dsum += dist[j];
    }

    double denom = 0.0;
    for (int j=0;j<5;j++) denom += (1.0-dist[j]/dsum);
    
    sumvec.SetXYZ(0.,0.,0.);
    for (int j=0;j<5;j++) {
      fieldvec = nnvec.at(j)*((1.0-dist[j]/dsum)/denom);
      sumvec += fieldvec;
    }
    // getting the proportions right between x and y field components 
    triplet.SetXYZ(sumvec.X()*fabs(TMath::Cos(angle)),sumvec.X()*fabs(TMath::Sin(angle)),sumvec.Z());
//     std::cout << "in Fields: point coordinates " << xv << " " << yv << " " << zv << std::endl;
//     std::cout << "in Fields: average field value " << sumvec.X() << " " << sumvec.Y() << " " << sumvec.Z() << std::endl;
    nnvec.clear();
  }
  else {
    // outside anything relevant, stop transport.
    gv = 0;
    triplet.SetXYZ(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  return triplet;
  
}


// Fields3D class
Fields3D::Fields3D(GeometryModel* g) : gm(g) // copy pointer
{}


Fields3D::~Fields3D() {
  if (allx) { // all set together
    delete [] allx;
    delete [] ally;
    delete [] allz;
    delete [] alldx;
    delete [] alldy;
    delete [] alldz;
  }
  if (coordinates) delete coordinates;
}

void Fields3D::prepare_fields(ComsolFields3D& fem) {
  std::vector<XYZPoint> cdata = fem.positions();
  int nentries = cdata.size();

  coordinates = new TKDTreeID(nentries,3,1); // 3D points
  
  allx = new double [nentries];
  ally = new double [nentries];
  allz = new double [nentries];

  for (int i=0;i<nentries;i++){
    allx[i] = cdata[i].x();
    ally[i] = cdata[i].y();
    allz[i] = cdata[i].z();
  }
  coordinates->SetData(0,allx);
  coordinates->SetData(1,ally);
  coordinates->SetData(2,allz);
  coordinates->Build();
  // KDTree built, all in memory

  alldx = new double [nentries];
  alldy = new double [nentries];
  alldz = new double [nentries];

  std::vector<XYZPoint> ddata = fem.driftmap();
  for (int i=0;i<nentries;i++){
    alldx[i] = ddata[i].x();
    alldy[i] = ddata[i].y();
    alldz[i] = ddata[i].z();
  }  
  //  std::cout << "in Fields::prepare fields finished." << std::endl;
  // all done and in memory
  fem.clear();
}


XYZPoint Fields3D::getFieldValue(XYZPoint& p, int& gv, bool& analytic)
{
  // thread access protection
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  
  // ask the geometry
  double xv = p.x();
  double yv = p.y();
  double zv = p.z();
  double angle = TMath::ATan2(yv,xv); // x-y-plane

  int value = gm->whereami(xv,yv,zv); // single access point to geomodel
  XYZPoint triplet;

  if (value==2) {
    gv = 2; // anode encoding
    triplet.SetXYZ(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  else if (value==1) { // comsol region
    gv = 1; // drift region code
    double point[3];
    double dist[7]; // check on nearest 5 neighbours in grid
    int indx[7];
    
    XYZVector fieldvec;
    XYZVector sumvec;
    std::vector<XYZVector> nnvec;
    
    double dsum = 0.0;
    point[0] = xv;
    point[1] = yv;
    point[2] = zv;  // relative to z-axis origin
    // in 3D look for minimum 6+1 neighbours incl itself
    coordinates->FindNearestNeighbors(point,7,indx,dist);
    for (int j=0;j<7;j++) {
      fieldvec.SetXYZ(alldx[indx[j]], alldy[indx[j]], alldz[indx[j]]);
      nnvec.push_back(fieldvec);
      dsum += dist[j];
    }

    double denom = 0.0;
    for (int j=0;j<7;j++) denom += (1.0-dist[j]/dsum);
    
    sumvec.SetXYZ(0.,0.,0.);
    for (int j=0;j<7;j++) {
      fieldvec = nnvec.at(j)*((1.0-dist[j]/dsum)/denom);
      sumvec += fieldvec;
    }
    // getting the proportions right between x and y field components 
    triplet.SetXYZ(sumvec.X()*fabs(TMath::Cos(angle)),sumvec.X()*fabs(TMath::Sin(angle)),sumvec.Z());
//     std::cout << "in Fields: point coordinates " << xv << " " << yv << " " << zv << std::endl;
//     std::cout << "in Fields: average field value " << sumvec.X() << " " << sumvec.Y() << " " << sumvec.Z() << std::endl;
    nnvec.clear();
  }
  else {
    // outside anything relevant, stop transport.
    gv = 0;
    triplet.SetXYZ(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  return triplet;
  
}
