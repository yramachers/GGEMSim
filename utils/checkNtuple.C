void check(const char* fn) {
  TFile ff(fn, "READ");
  TNtupleD* nt = (TNtupleD*)ff.Get("charge");
  double nanode = 0.;
  double nions = 0.;
  double nphot = 0.;
  nt->SetBranchAddress("nanode",&nanode);
  nt->SetBranchAddress("nions",&nions);
  nt->SetBranchAddress("nphot",&nphot);
  for (int i=0;i<nt->GetEntries();++i) {
    nt->GetEntry(i);
    std::cout << "charges: " << nions << " photons: " << nphot << " on anode: " << nanode << std::endl;
  }
  ff.Close();
}
