#include "TFile.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH2.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "JetscapeAnalysis.h"

//-----------------------------------------------------------------
// Default constructor
JetscapeAnalysis::JetscapeAnalysis():
  fHistos(),
  fEventID(0),
  fCrossSection(0),
  fJetR()
{
  fHistos = new THashList();
  fHistos->SetName("histos");
}

//-----------------------------------------------------------------
JetscapeAnalysis::~JetscapeAnalysis()
{
}

//-----------------------------------------------------------------
void JetscapeAnalysis::Init()
{
  
  // Create histograms
  
  // Event histograms
  CreateTH1("hCrossSection", "hCrossSection", 10000, 0, 100);
  
  // Hadron histograms
  CreateTH1("hHadronN", "hHadronN", 1000, 0, 1000);
  CreateTH1("hHadronPt", "hHadronPt", 2000, 0, 200);
  CreateTH2("hHadronEtaPhi", "hHadronEtaPhi", 100, -5, 5, 100, -6.28, 6.28);
  
  // Jet histograms
  for (auto jetR : fJetR) {
    
    std::string histname = FormJetHistoName("N", jetR);
    CreateTH1(histname.c_str(), histname.c_str(), 1000, 0, 1000);

    histname = FormJetHistoName("Pt", jetR);
    CreateTH1(histname.c_str(), histname.c_str(), 2000, 0, 200);

    histname = FormJetHistoName("EtaPhi", jetR);
    CreateTH2(histname.c_str(), histname.c_str(), 100, -5, 5, 100, -6.28, 6.28);
    
  }
  
}

//-----------------------------------------------------------------
void JetscapeAnalysis::AnalyzeEvent(const HepMC::GenEvent &event) {
  
  // Print and store basic event info
  GetEventInfo(event);
  
  // Get list of hadrons from the event, and fill some histograms
  const std::vector<HepMC::GenParticlePtr> hadrons = GetHadrons(event);
  FillHadronHistograms(hadrons);
  
  // Perform jet finding
  double jetR = 0.2;
  double minJetPt = 1.;
  double absEtaMax = 1.;
  
  std::vector<fastjet::PseudoJet> fjHadrons = FillFastjetConstituents(hadrons);
  for (auto jetR : fJetR) {
    
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    fastjet::ClusterSequence cs(fjHadrons, jetDef);
    std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(minJetPt));
                        
    // Apply selection on jet: Pt, eta
    std::vector<fastjet::PseudoJet> jets_accepted = GetAcceptedJets(jets, minJetPt, absEtaMax);
    
    // Fill some jet histograms
    FillJetHistograms(jets_accepted, jetR);
    
    jets.clear();
    jets_accepted.clear();
  }

  // Clean up
  fjHadrons.clear();
  
  // How to identify the final state partons?
  
  // From Manual: “Note, however, that HepMC always requires incoming and outgoing edges for every vertex,
  // whereas JETSCAPE graphs normally start and end with a vertex; the HepMC output therefore duplicates some edges.”
  
  // Particles and vertices in HepMC3 are stored in topological order. This means that when creating vertices,
  // incoming particles must have id lower than any of the outgoing particles.
  
  /*
   for (auto particle : event.particles()) {
   int pdg = particle->pdg_id();
   HepMC::FourVector momentum = particle->momentum();
   }
   */
  
}

//-----------------------------------------------------------------
// Create output file and write histograms
void JetscapeAnalysis::WriteOutput()
{
  
  // Fill cross-section with last event's value, which is most accurate
  FillTH1("hCrossSection", fCrossSection/(1e9));

  // Create output file
  TFile* f = new TFile("AnalysisResults.root", "RECREATE");
  
  // Write all histograms in fHistos list
  TIter next(fHistos);
  TObject* obj = 0;
  while ((obj = next())) {
    obj->Write();
  }
  
  f->Close();
}

//-----------------------------------------------------------------
void JetscapeAnalysis::GetEventInfo(const HepMC::GenEvent &event) {
  
  // Get cross-section
  std::shared_ptr<HepMC::GenCrossSection> crossSection = event.attribute<HepMC::GenCrossSection>("GenCrossSection");
  fCrossSection = crossSection->cross_section;
  double xsec_error = crossSection->cross_section_error;
  
  // Print event number
  if (fEventID % 100 == 0) {
    std::cout << "Event: " << fEventID << std::endl;
  }
  
  // Print some basic info for first event
  if (fEventID == 0) {
    
    std::cout << "xsec: " << fCrossSection << " +/- " << xsec_error << "pb" << std::endl;
    
    // Get heavy ion attributes
    std::shared_ptr<HepMC::GenHeavyIon> heavyIon = event.attribute<HepMC::GenHeavyIon>("GenHeavyIon");
    int nColl = heavyIon->Ncoll;
    int nPart = heavyIon->Npart_proj;
    double eventPlaneAngle = heavyIon->event_plane_angle;
    std::cout << "NColl = " << nColl << ", NPart = " << nPart << ", EP-angle = " << eventPlaneAngle << std::endl;
    
    //HepMC::Print::listing(event);
    //HepMC::Print::content(event);
  }
  
  fEventID++;
  
}

//-----------------------------------------------------------------
// Get list of hadrons.
// Final state hadrons (from jet + bulk) are stored as outgoing particles in a disjoint vertex with t = 100
std::vector<HepMC::GenParticlePtr> JetscapeAnalysis::GetHadrons(const HepMC::GenEvent &event) {
  
  std::vector<HepMC::GenParticlePtr> hadrons;
  for (auto vertex : event.vertices()) {
    
    double vertexTime = vertex->position().t();
    if ( abs(vertexTime - 100) < 1e-3 ) {
      hadrons = vertex->particles_out();
    }
    
  }
  
  return hadrons;
  
}

//-----------------------------------------------------------------
// Fill hadron histograms
void JetscapeAnalysis::FillHadronHistograms(const std::vector<HepMC::GenParticlePtr> hadrons) {
  
  // Loop through hadrons
  for (auto hadron : hadrons) {
    
    // Fill some basic hadron info
    int pid = hadron->pid();
    
    HepMC::FourVector momentum = hadron->momentum();
    double pt = momentum.pt();
    double eta = momentum.eta();
    double phi = momentum.phi(); // [-pi, pi]
    
    FillTH1("hHadronPt", pt);
    FillTH2("hHadronEtaPhi", eta, phi);
    
  }
  FillTH1("hHadronN", hadrons.size());
  
}

//-----------------------------------------------------------------
// Fill hadrons into vector of fastjet pseudojets
std::vector<fastjet::PseudoJet> JetscapeAnalysis::FillFastjetConstituents(const std::vector<HepMC::GenParticlePtr> hadrons) {
  
  std::vector<fastjet::PseudoJet> fjConstituents;
  for (auto hadron : hadrons) {
    
    HepMC::FourVector momentum = hadron->momentum();
    double px = momentum.px();
    double py = momentum.py();
    double pz = momentum.pz();
    double e = momentum.e();
    fjConstituents.push_back(fastjet::PseudoJet(px, py, pz, e));
    
  }
  return fjConstituents;
}


//-----------------------------------------------------------------
// Apply jet selection
std::vector<fastjet::PseudoJet> JetscapeAnalysis::GetAcceptedJets(const std::vector<fastjet::PseudoJet> jets, double minPt, double absEtaMax) {
  
  fastjet::Selector select_pt = fastjet::SelectorPtMin(minPt);
  fastjet::Selector select_eta = fastjet::SelectorAbsEtaMax(absEtaMax);
  fastjet::Selector selection = select_pt && select_eta;
  return selection(jets);
  
}

//-----------------------------------------------------------------
// Fill jet histograms
void JetscapeAnalysis::FillJetHistograms(const std::vector<fastjet::PseudoJet> jets, double jetR) {
  
  std::string histname = FormJetHistoName("N", jetR);
  FillTH1(histname.c_str(), jets.size());
  
  for (auto jet : jets) {
    
    double jetPt = jet.pt();
    double jetEta = jet.eta();
    double jetPhi = jet.phi(); // [0, 2pi]
    
    histname = FormJetHistoName("Pt", jetR);
    FillTH1(histname.c_str(), jetPt);
    
    histname = FormJetHistoName("EtaPhi", jetR);
    FillTH2(histname.c_str(), jetEta, jetPhi);
    
    // Loop through jet constituents
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    for (auto constituent : constituents ) {
      double constituentPt = constituent.pt();
    }
  
  }

}

//-----------------------------------------------------------------
std::string JetscapeAnalysis::FormJetHistoName(const char* title, double jetR) {
  
  std::string hname = "hJet0";
  std::string jetR_str = GetRLabel(jetR);
  hname.append(jetR_str);
  hname.append("_");
  hname.append(title);
  return hname;
  
}

//-----------------------------------------------------------------
std::string JetscapeAnalysis::GetRLabel(double jetR) {
  
  int jetR_int = jetR*100;
  std::string jetR_str = std::to_string(jetR_int);
  return jetR_str;
  
}

//-----------------------------------------------------------------
TH1* JetscapeAnalysis::CreateTH1(const char* name, const char* title, int nbins, double xmin, double xmax)
{
  if(fHistos->FindObject(name)){
    Printf("Histogram already exists: %s", name);
    return 0;
  }
  
  TH1* h = new TH1D(name, title, nbins, xmin, xmax);
  fHistos->Add(h);
  return h;
}

//-----------------------------------------------------------------
TH2* JetscapeAnalysis::CreateTH2(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax)
{
  if(fHistos->FindObject(name)){
    Printf("Histogram already exists: %s", name);
    return 0;
  }
  
  TH2* h = new TH2D(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  fHistos->Add(h);
  return h;
}

//-----------------------------------------------------------------
void JetscapeAnalysis::FillTH1(const char *name, double x) {
  
  TH1* hist = dynamic_cast<TH1*>(fHistos->FindObject(name));
  if(!hist){
    Printf("Histogram Fill not found: %s", name);
    return;
  }
  hist->Fill(x);
}

//-----------------------------------------------------------------
void JetscapeAnalysis::FillTH2(const char *name, double x, double y) {
  
  TH2* hist = dynamic_cast<TH2*>(fHistos->FindObject(name));
  if(!hist){
    Printf("Histogram Fill not found: %s", name);
    return;
  }
  hist->Fill(x, y);
}
