#include "TFile.h"

#include "JetscapeAnalysis.h"

//-----------------------------------------------------------------
// Default constructor
JetscapeAnalysis::JetscapeAnalysis():
  fEventID(0),
  fCrossSection(0),
  hCrossSection(nullptr),
  hHadronN(nullptr),
  hHadronPt(nullptr),
  hHadronEtaPhi(nullptr)
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
  
}

//-----------------------------------------------------------------
void JetscapeAnalysis::AnalyzeEvent(HepMC::GenEvent &event) {
  
  // Get some event info
  GetEventInfo(event);
  
  // Get list of hadrons
  std::vector<HepMC::GenParticlePtr> hadrons = GetHadrons(event);
  
  // Loop through hadrons, do stuff
  for (auto hadron : hadrons) {
    
    int pid = hadron->pid();
    
    HepMC::FourVector momentum = hadron->momentum();
    double pt = momentum.pt();
    double eta = momentum.eta();
    double phi = momentum.phi();
    
    hHadronPt->Fill(pt);
    hHadronEtaPhi->Fill(eta, phi);
    
  }
  hHadronN->Fill(hadrons.size());
  
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
    std::cout << "Event: " << fEventID << std::endl;
    std::cout << "xsec: " << fCrossSection << " +/- " << xsec_error << "pb" << std::endl;
    // Note that PYTHIA improves its accuracy by Monte Carlo integration in the course of the run, so the values associated with the last generated event should be the most accurate ones.
    
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
