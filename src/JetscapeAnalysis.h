#ifndef JETSCAPEANALYSIS_H
#define JETSCAPEANALYSIS_H

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"

#include "TH1.h"
#include "TH2.h"
  
class JetscapeAnalysis {

  public:
  
    // Default constructor
    JetscapeAnalysis();
    
    // Destructor
    virtual ~JetscapeAnalysis();
  
    // Main analysis functions, called by steering macro
    void Init();
    void AnalyzeEvent(HepMC::GenEvent &event);
    void WriteOutput();
  
  protected:
  
    // Helper functions
    void GetEventInfo(HepMC::GenEvent &event);
    std::vector<HepMC::GenParticlePtr> GetHadrons(HepMC::GenEvent &event);
  
    // Data members
    unsigned int fEventID;
    float fCrossSection;
  
    // Histograms
    TH1F* hCrossSection;
  
    TH1F* hHadronN;
    TH1F* hHadronPt;
    TH2F* hHadronEtaPhi;
  
    TH1F* hJetN;
    TH1F* hJetPt;
    TH2F* hJetEtaPhi;  
    
};

#endif
