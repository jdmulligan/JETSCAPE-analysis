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
    
  HepMC::ReaderAscii reader(argv[1]);
  
  // Initialize analysis task
  JetscapeAnalysis* analyzer = new JetscapeAnalysis();
  analyzer->Init();

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
