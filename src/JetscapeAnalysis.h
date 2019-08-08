#ifndef JETSCAPEANALYSIS_H
#define JETSCAPEANALYSIS_H

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "fastjet/PseudoJet.hh"

#include "THashList.h"

class TH1;
class TH2;

class JetscapeAnalysis {

  public:
  
    // Default constructor
    JetscapeAnalysis(int bin, std::string outputDirBin);
    
    // Destructor
    virtual ~JetscapeAnalysis();
  
    // Main analysis functions, called by steering macro
    void Init();
    void AnalyzeEvent(const HepMC3::GenEvent &event);
    void WriteOutput();
  
    // Setters
    void SetJetR(std::vector<double> r) { fJetR = r; }
    void SetMinJetPt(double d) { fMinJetPt = d; }
    void SetAbsJetEtaMax(double d) {fAbsJetEtaMax = d; }

  protected:
  
    // Helper functions
    void                                   GetEventInfo(const HepMC3::GenEvent &event);
    std::vector<HepMC3::ConstGenParticlePtr>     GetHadrons(const HepMC3::GenEvent &event);
    void                                   FillHadronHistograms(const std::vector<HepMC3::ConstGenParticlePtr> hadrons);
  
    std::vector<fastjet::PseudoJet>        FillFastjetConstituents(const std::vector<HepMC3::ConstGenParticlePtr> hadrons);
    std::vector<fastjet::PseudoJet>        GetAcceptedJets(const std::vector<fastjet::PseudoJet> jets);
    void                                   FillJetHistograms(const std::vector<fastjet::PseudoJet> jets, double jetR);
  
    // Data members
    std::string                            fOutputDirBin;
    unsigned int                           fEventID;
    unsigned int                           fPtHatBin;
    float                                  fCrossSection;
    std::vector<double>                    fJetR;
    double                                 fMinJetPt;
    double                                 fAbsJetEtaMax;
  
    // Histogram management (based loosely on THistManager http://alidoc.cern.ch/AliPhysics/master/class_t_hist_manager.html)
    THashList *fHistos;                   ///< List of histograms
    std::string GetRLabel(double jetR);
    std::string FormJetHistoName(const char* title, double jetR);
    TH1* CreateTH1(const char *name, const char *title, int nbins, double xmin, double xmax);
    TH2* CreateTH2(const char* name, const char* title, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax);
    void FillTH1(const char *name, double x, double weight=1.);
    void FillTH2(const char *name, double x, double y);
  
};

#endif
