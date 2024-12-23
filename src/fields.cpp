#include <iostream>

// us
#include "fields.hh"

// ROOT includes
#include "TFile.h"
#include "TNtupleD.h"
#include "TMath.h"
#include "Math/Vector3D.h"


// come now as 2D data in x,y from comsol
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


// Fields class
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
  //  std::cout << "in Fields::prepare fields." << std::endl;

  coordinates = new TKDTreeID(nentries,2,1);
  
  allx = new double [nentries];
  allz = new double [nentries];

  for (int i=0;i<nentries;i++){
    // no transf needed, requests come as vectors
    allx[i] = cdata[i].x();
    allz[i] = cdata[i].z();
  }
  coordinates->SetData(0,allx);
  coordinates->SetData(1,allz);
  coordinates->Build();
  //  std::cout << "in Fields::KDTree built." << std::endl;
  // KDTree built, all in memory

  alldx = new double [nentries];
  alldz = new double [nentries];

  // no scaling if Comsol fields are at 500V/cm
  std::vector<XYZPoint> ddata = fem.driftmap();
  for (int i=0;i<nentries;i++){
    alldx[i] = ddata[i].x();
    alldz[i] = ddata[i].z();
  }  
  //  std::cout << "in Fields::prepare fields finished." << std::endl;
  // all done and in memory
  fem.clear();
}


XYZPoint Fields::getFieldValue(GeometryModel& gm, XYZPoint& p, bool& analytic) {

  // common routine to ask for field value
  // bool drift decides between Drift field: drift=True
  // or weighting field values: drift=False

  // ask the geometry
  double xv = p.x();
  double yv = p.y();
  double zv = p.z();
  double rad = TMath::Sqrt(xv*xv+yv*yv); // x-y-plane
  double angle = TMath::ATan2(yv,xv); // x-y-plane
  int value = gm.whereami(xv,yv,zv);
  //  std::cout << "in Fields::answer to whereami: " << value << std::endl;
  XYZPoint triplet;
  
  if (value==1) { // comsol region
    double point[2];
    double dist[4]; // check on nearest 4 neighbours in grid
    int indx[4];
    
    XYZVector fieldvec;
    XYZVector sumvec;
    std::vector<XYZVector> nnvec;
    
    double dsum = 0.0;
    point[0] = rad;  // relative to origin x
    point[1] = zv;  // relative to z-axis origin

    //    std::cout << "in Fields: point coordinates " << xv << " " << yv << " " << zv << std::endl;
        
    coordinates->FindNearestNeighbors(point,4,indx,dist);
    for (int j=0;j<4;j++) {
      fieldvec.SetXYZ(alldx[indx[j]], 0.0, alldz[indx[j]]);
      // 	std::cout << "in Fields: nearest coords: " << allx[indx[j]] << " " << ally[indx[j]] << " " << allz[indx[j]] << std::endl;
      // 	std::cout << "in Fields: Drift field value: " << alldx[indx[j]] << " " << alldy[indx[j]] << " " << alldz[indx[j]] << std::endl;
      nnvec.push_back(fieldvec);
      dsum += dist[j];
    }

    double denom = 0.0;
    for (int j=0;j<4;j++) denom += (1.0-dist[j]/dsum);
    
    sumvec.SetXYZ(0.,0.,0.);
    for (int j=0;j<4;j++) {
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
    triplet.SetXYZ(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  return triplet;
  
}
