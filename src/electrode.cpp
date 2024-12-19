// us
#include "electrode.hh"

//****************
// Electrode Model
//****************
// Constructors
Electrode::Electrode(ComsolFields* fem, GeometryModel* g)
  : gm(g), femfields(fem)
{
  field = NULL;
}


// default Destructor
Electrode::~Electrode() {
  if (field) delete field;
}


void Electrode::initfields() {
  field = new Fields(femfields, gm); // create from file + geometry info
}


XYZPoint Electrode::getFieldValue(bool& analytic, XYZPoint& p) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  XYZPoint triplet;
  
  triplet = field->getDriftField(p, analytic);

  return triplet;
}


