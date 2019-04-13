// Author: James Mulligan (james.mulligan@berkeley.edu)
/*
 *  Macro for analysis of jetscape events in HepMC3 format.
 */

#include "HepMC/GenEvent.h"
#include "HepMC/ReaderAscii.h"

#include "JetscapeAnalysis.h"

//-----------------------------------------------------------------
// Main function
int main(int argc, char** argv)
{
  
  // Read command line arguments
  int bin = atoi(argv[1]);
  std::string outputDirBin = argv[2];
  
  // Configure analysis task
  JetscapeAnalysis* analyzer = new JetscapeAnalysis(bin, outputDirBin);
  analyzer->SetJetR({0.2, 0.4});
  
  // Initialize analysis task
  analyzer->Init();
  
  // Read HepMC file
  std::string hepmcFile = outputDirBin.append("test_out.hepmc");
  HepMC::ReaderAscii reader(hepmcFile.c_str());
  
  // Loop over HepMC events, and call analysis task to process them
  while (!reader.failed()) {

    // Read event
    HepMC::GenEvent event(HepMC::Units::GEV,HepMC::Units::MM);
    reader.read_event(event);
    if (reader.failed()) {
      break;
    }
    
    analyzer->AnalyzeEvent(event);
    
  }
  reader.close();
  
  // Write analysis task output to ROOT file
  analyzer->WriteOutput();
  
  return 0;
  
}
