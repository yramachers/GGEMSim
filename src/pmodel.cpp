#include "pmodel.hh"

// ROOT includes
#include "TMath.h"

Physics_Model::Physics_Model() {
  generator.seed(rd()); // use random seed
  
  // local declarations
  double inelsum(double* x, double* par);
  inel  = new TF1("inelsum",inelsum,11.55,50.0,4);

  // phaseshift calculation, set interpolator data
  data0 = new TGraph(11);
  data1 = new TGraph(11);
  data2 = new TGraph(11);
  data3 = new TGraph(11);
  double energy[11] = {0.02,0.06,0.1,0.2,0.3,0.4,0.5,0.7,1.0,2.0,3.0};
  double values[44] = {4.426e-2,3.092e-3,5.836e-4,1.537e-4,5.363e-2,8.155e-3,1.577e-3,5.21e-4,5.133e-2,1.219e-2,2.686e-3,8.937e-4,3.255e-2,1.93e-2,5.516e-3,1.831e-3,8.595e-3,2.307e-2,8.379e-3,2.756e-3,-1.639e-2,2.46e-2,1.128e-2,3.667e-3,-4.122e-2,2.444e-2,1.419e-2,4.574e-3,-8.931e-2,2.04e-2,2.016e-2,6.374e-3,-1.566e-1,8.199e-3,2.927e-2,9.059e-3,-3.48e-1,-5.612e-3,6.446e-2,1.788e-2,-5.057e-1,-1.329e-1,1.111e-1,2.669e-2};
    
  int m = 11, n = 4;
  int i;
  for (i=0;i<m;i++) {
    data0->SetPoint(i,energy[i],values[i*n]);
    data1->SetPoint(i,energy[i],values[1+i*n]);
    data2->SetPoint(i,energy[i],values[2+i*n]);
    data3->SetPoint(i,energy[i],values[3+i*n]);
  }
  
}


Physics_Model::~Physics_Model() {
  delete inel;
  delete data0;
  delete data1;
  delete data2;
  delete data3;
}


//****************************
// theory from JPhysB16,4023
//****************************
double Physics_Model::angle_function(TRandom3& rnd, double energy)
{
  using pld_type = std::piecewise_linear_distribution<double>;

  //  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  if (energy <= 0.02) // trial correction
    return TMath::Pi()*rnd.Rndm();//isotropic for <0.02 eV
  double p[5];
  double phaseshift[4];
  find_phaseshift(energy,phaseshift);
  
  p[0] = energy; 
  for (Int_t i=0;i<4;i++) {
    p[i+1] = phaseshift[i];
  }
  pld_type adist(360,0.0,TMath::Pi(),angleFunction(p)); // unary functor in header
  double theta = adist(generator); // random angle theta
  //  cout << "theta=" << theta << " at [eV] " << energy << endl;
  return theta;
}

void Physics_Model::find_phaseshift(double e, double* phaseshift)
{
    if (e<0.02) e = 0.02; // minimum energy value fixed

    phaseshift[0] = data0->Eval(e);
    phaseshift[1] = data1->Eval(e);
    phaseshift[2] = data2->Eval(e);
    phaseshift[3] = data3->Eval(e);
}


// check new cross sections better than Eachran, Stauffer (1983)
double Physics_Model::cross_section(TRandom3& rnd, double energy, bool &momentum_flag, int &inel_flag)
{
  double cs_el, cs_inel,ratio;
  double cs_etrans;
  double cs_ptrans;
  
  inel_flag = 0; // default case
  momentum_flag = kTRUE; // default case
  
  // purely elastic scattering possible
  if (energy>0.0 && energy<=11.55) {
    cs_ptrans = stauffer_momentum_cs(energy);
    cs_etrans = stauffer_elastic_cs(energy);
    ratio = cs_ptrans / (cs_etrans+cs_ptrans);

    if (energy>4.5 && energy<=11.55) { // overestimated sigma
      cs_ptrans *= 0.73;
      cs_etrans *= 0.73;
    }
    if (rnd.Rndm() <  ratio) momentum_flag = kTRUE;
    else momentum_flag = kFALSE; // energy transfer cs
    return (momentum_flag ? cs_ptrans : cs_etrans);
  }
  // otherwise decide
  else {	
    cs_ptrans = stauffer_momentum_cs(energy);
    cs_el = stauffer_elastic_cs(energy);
    cs_inel = inelastic_cs(energy);
    ratio = cs_inel / (cs_el+cs_inel);
    if (rnd.Rndm() < ratio) {
      inel_flag = -1; // photon
      if (energy > 15.76) inel_flag = 1; // electron
      //      std::cout << " inelastic sigma = " << cs_inel << " ratio = " << ratio << std::endl;
      return 0.0; // stops further transport after inel
    }
    // elastic as normal
    ratio = cs_ptrans / cs_el;
    if (rnd.Rndm() <  ratio) momentum_flag = kTRUE;
    else momentum_flag = kFALSE; // energy transfer cs
    return (momentum_flag ? cs_ptrans+cs_inel : cs_el+cs_inel);
  }
}


//****************************
//****************************
// Stauffer formulae
//****************************
double Physics_Model::stauffer_elastic_cs(double energy)
{
    int l;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    double sigma, factor, sum = 0.0;

    std::array<double,6> delta = calc_phaseshift1to6(energy);

    factor = 4.0*TMath::Pi()/(k*k);
    for (l=0;l<6;l++) 
    	sum += (2.0*l+1)*TMath::Sin(delta[l])*TMath::Sin(delta[l]);

    sigma = factor * sum;
    
    // unit conversion Bohr radius^2 to SI
    sigma *= 2.80028561e-21;

    return sigma;
}

double Physics_Model::stauffer_momentum_cs(double energy)
{
    int l;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    double sigma, factor, sum = 0.0;

    std::array<double,6> delta = calc_phaseshift1to6(energy);

    factor = 4.0*TMath::Pi()/(k*k);
    for (l=0;l<5;l++) 
    	sum += (l+1)*TMath::Sin(delta[l]-delta[l+1])*TMath::Sin(delta[l]-delta[l+1]);
    
    sigma = factor * sum;
    
    // unit conversion Bohr radius^2 to SI
    sigma *= 2.80028561e-21;

    return sigma;
}

std::array<double,6> Physics_Model::calc_phaseshift1to6(double energy)
{
    int l;
    double beta_argon = 10.755;
    double k = TMath::Sqrt(energy/13.6025); // k in [amu]
    std::array<double,6> delta;
    // Energy input in [eV]
    // output cross section in m^2
    double h[5][4] = {{0.150575e1,0.0,0.0,0.0},{0.671099e1,0.225686e1,0.319699,0.107419},{-0.103246e3,-0.805571e1,-0.471566,0.0},{0.11987e2,0.213514e1,0.662807,-0.183263},{0.0,0.0,0.0,0.315322}};


    double d[4][4] = {{0.122028e2,-0.80281e-2,-0.145269e1,0.0},{0.429028e2,0.425572e1,0.201048,-0.177524e1},{0.547674e1,-0.677844,0.584853,0.293544e1},{0.0,0.0,0.0,-0.497051}};

    double p4[2] = {0.324336e-2,-0.867549e-3};
    double p5[2] = {0.382582e-3,0.0};
    
    double N[4], D[4];
    N[0] = h[0][0]*k + h[1][0]*k*k + h[2][0]*k*k*k + h[3][0]*k*k*k*k;
    for (l=1;l<4;l++) 
    	N[l] = h[1][l]*k*k + h[2][l]*k*k*k + h[3][l]*k*k*k*k + h[4][l]*k*k*k*k*k;
    
    D[0] = 1.0 + d[0][0]*k + d[1][0]*k*k + d[2][0]*k*k*k + d[3][0]*k*k*k*k;;
    D[1] = 1.0 + d[0][1]*k + d[1][1]*k*k + d[2][1]*k*k*k*TMath::Log(k);
    for (l=2;l<4;l++) 
    	D[l] = 1.0 + d[0][l]*k + d[1][l]*k*k + d[2][l]*k*k*k + d[3][l]*k*k*k*k;

    for (l=0;l<4;l++) 
	    delta[l] = N[l] / D[l];
    
    delta[4] = k*k*beta_argon*func_a(4) + k*k*k*k*beta_argon*beta_argon*func_b(4) + p4[0]*k*k*k*k*k*k + p4[1]*k*k*k*k*k*k*k;

    delta[5] = k*k*beta_argon*func_a(5) + k*k*k*k*beta_argon*beta_argon*func_b(5) + p5[0]*k*k*k*k*k*k;
    return delta;
}

double Physics_Model::func_a(int l)
{
    double dummy = (2.0*l-1)*(2.0*l+1)*(2.0*l+3);
    return TMath::Pi()/dummy;
}

double Physics_Model::func_b(int l)
{
    double dummy = (2.0*l-1)*(2.0*l+1)*(2.0*l+3);
    double dummy1 = (2.0*l-3)*dummy*dummy*dummy*(2.0*l+5);
    double dummy2 = (2.0*l+1)*(2.0*l+1);
    return TMath::Pi()*(15.0*dummy2*dummy2 - 140.0*dummy + 128.0)/dummy1;
}

//**************************
//** Inelastic from fits in NIM A242, 327
//**************************

double Physics_Model::inelastic_cs(double energy)
{
  //  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  int flag;

  double bins[6] = {11.55,11.62,11.72,11.83,14.0,15.76};

  if (energy <= bins[0]) return 0.0;
  if (energy <= bins[5]) return inelastic_cs_mason_newell(energy);

  double array[7][4] = {{-0.4087e-2,0.5321e-3,-0.1543e-4,0.0},{-0.1775e-2,0.2037e-3,-0.4394e-5,0.0},{-0.6713e-3,0.8463e-4,-0.2335e-5,0.0},{-0.5562e-2,0.7999e-3,-0.3438e-4,0.551e-6},{-0.2844e-1,0.3025e-2,-0.7096e-4,0.0},{0.681e1,-0.9458,0.3259e-1,0.0},{-9.884,0.1894,0.0,0.0}};

  if (energy>bins[5] && energy <=17.5) flag = 6;
  if (energy>17.5) flag = 7;

  double par[4], result = 0.0;

  for (int j=1;j<flag;j++) { // all previous
  	par[0] = array[j-1][0];
  	par[1] = array[j-1][1];
	par[2] = array[j-1][2];
	par[3] = array[j-1][3];
	inel->SetParameters(par);
	result += inel->Eval(bins[j]-1.0e-3);
  }
  par[0] = array[flag-1][0]; // + actual
  par[1] = array[flag-1][1];
  par[2] = array[flag-1][2];
  par[3] = array[flag-1][3];
  inel->SetParameters(par);
  result += inel->Eval(energy);
  // unit conversion pi*Bohr radius^2 to SI
  result *= 8.797356696e-21;
  result += inelastic_cs_mason_newell(energy); // excitation + ionization
  //
  return result;        
}

double inelsum(double* x, double* par)
{
    double result = par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
    return result;
}

double Physics_Model::inelastic_cs_mason_newell(double energy)
{
    double par[7] = {34.7531,-10.2377,1.10632,-0.0562047,0.00147617,-1.94974e-5,1.02707e-7};
    double result = emp_fit(energy,par);
    return result*1.0e-21; // SI units
}

double Physics_Model::emp_fit(double x, double* par)
{
    // pol6 fits by 190.0/20 chi2/ndf up to 50 eV, not further
    return par[0]+par[1]*x+par[2]*x*x+par[3]*x*x*x+par[4]*x*x*x*x+par[5]*x*x*x*x*x+par[6]*x*x*x*x*x*x;
}


