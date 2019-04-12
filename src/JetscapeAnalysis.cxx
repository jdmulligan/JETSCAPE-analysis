#include "TFile.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "JetscapeAnalysis.h"

//-----------------------------------------------------------------
// Default constructor
JetscapeAnalysis::JetscapeAnalysis():
  fEventID(0),
  fCrossSection(0),
  hCrossSection(nullptr),
  hHadronN(nullptr),
  hHadronPt(nullptr),
  hHadronEtaPhi(nullptr),
  hJetN(nullptr),
  hJetPt(nullptr),
  hJetEtaPhi(nullptr)
{
}

//-----------------------------------------------------------------
JetscapeAnalysis::~JetscapeAnalysis()
{
}

//-----------------------------------------------------------------
void JetscapeAnalysis::Init()
{
  // Create histograms
  
  hCrossSection = new TH1F("hCrossSection", "hCrossSection", 10000, 0, 100);
  
  hHadronN = new TH1F("hHadronN", "hHadronN", 1000, 0, 1000);
  hHadronPt = new TH1F("hHadronPt", "hHadronPt", 2000, 0, 200);
  hHadronEtaPhi = new TH2F("hHadronEtaPhi", "hHadronEtaPhi", 100, -5, 5, 100, -6.28, 6.28);
  
  hJetN = new TH1F("hJetN", "hJetN", 1000, 0, 1000);
  hJetPt = new TH1F("hJetPt", "hJetPt", 2000, 0, 200);
  hJetEtaPhi = new TH2F("hJetEtaPhi", "hJetEtaPhi", 100, -5, 5, 100, -6.28, 6.28);
  
}

//-----------------------------------------------------------------
void JetscapeAnalysis::AnalyzeEvent(HepMC::GenEvent &event) {
  
  // Get some event info
  GetEventInfo(event);
  
  // Get list of hadrons
  std::vector<HepMC::GenParticlePtr> hadrons = GetHadrons(event);
  
  // Loop through hadrons, do jet finding
  double jetR = 0.2;
  double minJetPt = 1.;
  fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
  std::vector<fastjet::PseudoJet> fjHadrons;
  for (auto hadron : hadrons) {
    
    // Fill some basic hadron info
    int pid = hadron->pid();
    
    HepMC::FourVector momentum = hadron->momentum();
    double pt = momentum.pt();
    double eta = momentum.eta();
    double phi = momentum.phi(); // [-pi, pi]
    
    hHadronPt->Fill(pt);
    hHadronEtaPhi->Fill(eta, phi);
    
    // Fill fastjet constituents
    double px = momentum.px();
    double py = momentum.py();
    double pz = momentum.pz();
    double e = momentum.e();
    fjHadrons.push_back(fastjet::PseudoJet(px, py, pz, e));
    
  }
  hHadronN->Fill(hadrons.size());
  
  // Jet finding
  fastjet::ClusterSequence cs(fjHadrons, jetDef);
  std::vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(minJetPt));
                        
  // Jet selection
  fastjet::Selector select_pt = fastjet::SelectorPtMin(1.);
  fastjet::Selector select_eta = fastjet::SelectorAbsEtaMax(1.);
  fastjet::Selector selection = select_pt && select_eta;
  std::vector<fastjet::PseudoJet> jets_accepted = selection(jets);

  // Loop through jets
  hJetN->Fill(jets_accepted.size());
  for (auto jet : jets_accepted) {
    double jetPt = jet.pt();
    double jetEta = jet.eta();
    double jetPhi = jet.phi(); // [0, 2pi]
    
    hJetPt->Fill(jetPt);
    hJetEtaPhi->Fill(jetEta, jetPhi);
    
    // Loop through jet constituents
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    for (auto constituent : constituents ) {
      double constituentPt = constituent.pt();
    }
    
  }

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

  TFile* f = new TFile("AnalysisResults.root", "RECREATE");

  // Fill cross-section with last event's value, which is most accurate
  hCrossSection->Fill(fCrossSection/(1e9));
  hCrossSection->Write();
  
  hHadronN->Write();
  hHadronPt->Write();
  hHadronEtaPhi->Write();
  
  hJetN->Write();
  hJetPt->Write();
  hJetEtaPhi->Write();

  f->Close();
}

//-----------------------------------------------------------------
void JetscapeAnalysis::GetEventInfo(HepMC::GenEvent &event) {
  
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
std::vector<HepMC::GenParticlePtr> JetscapeAnalysis::GetHadrons(HepMC::GenEvent &event) {
  
  std::vector<HepMC::GenParticlePtr> hadrons;
  for (auto vertex : event.vertices()) {
    
    double vertexTime = vertex->position().t();
    if ( abs(vertexTime - 100) < 1e-3 ) {
      hadrons = vertex->particles_out();
    }
    
  }
  
  return hadrons;
  
}
