// Author: James Mulligan (james.mulligan@berkeley.edu)
/*
 *  Macro for analysis of jetscape events in HepMC3 format.
 */

#include "HepMC/GenEvent.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/Print.h"

// Declare helper functions
void PrintEventInfo(HepMC::GenEvent &event);
std::vector<HepMC::GenParticlePtr> GetHadrons(HepMC::GenEvent &event);

//-----------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
    
  HepMC::ReaderAscii reader(argv[1]);

  int eventID = 0;
  while (!reader.failed()) {

    // Read event
    printf("Event %d \n", eventID);
    HepMC::GenEvent event(HepMC::Units::GEV,HepMC::Units::MM);
    reader.read_event(event);
    if (reader.failed()) {
      break;
    }
    
    // Print some event info
    if (eventID == 0) {
      PrintEventInfo(event);
    }
    
    // Get list of hadrons
    std::vector<HepMC::GenParticlePtr> hadrons = GetHadrons(event);
    
    // Loop through hadrons, do stuff
    for (auto hadron : hadrons) {
      
      int pid = hadron->pid();
      
      HepMC::FourVector momentum = hadron->momentum();
      double pt = momentum.pt();
      double eta = momentum.eta();
      double phi = momentum.phi();
      
      if (abs(eta) < 1) {
        printf("pt = %f, eta = %f, phi = %f \n", pt, eta, phi);
      }
      
    }
    
    
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
    
    eventID++;
    
  }
  
  reader.close();
  
  return 0;
  
}

//-----------------------------------------------------------------
void PrintEventInfo(HepMC::GenEvent &event) {
  
  // Get cross-section
  std::shared_ptr<HepMC::GenCrossSection> crossSection = event.attribute<HepMC::GenCrossSection>("GenCrossSection");
  double xsec = crossSection->cross_section;
  double xsec_error = crossSection->cross_section_error;
  std::cout << "xsec: " << xsec << " +/- " << xsec_error << std::endl;
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

//-----------------------------------------------------------------
// Get list of hadrons.
// Final state hadrons (from jet + bulk) are stored as outgoing particles in a disjoint vertex with t = 100
std::vector<HepMC::GenParticlePtr> GetHadrons(HepMC::GenEvent &event) {
  
  std::vector<HepMC::GenParticlePtr> hadrons;
  for (auto vertex : event.vertices()) {
    
    double vertexTime = vertex->position().t();
    if ( abs(vertexTime - 100) < 1e-3 ) {
      hadrons = vertex->particles_out();
    }
    
  }
  
  if (hadrons.empty()) {
    printf("No hadrons \n");
  }
  else {
    printf("Found %d hadrons!\n", hadrons.size());
  }
  
  return hadrons;
  
}
