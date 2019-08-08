// Author: James Mulligan (james.mulligan@berkeley.edu)
/*
 *  Macro for analysis of jetscape events in HepMC3 format.
 */

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"

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
  HepMC3::ReaderAscii reader(hepmcFile.c_str());
  
  // Loop over HepMC events, and call analysis task to process them
  while (!reader.failed()) {

    // Read event
    HepMC3::GenEvent event(HepMC3::Units::GEV,HepMC3::Units::MM);
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
