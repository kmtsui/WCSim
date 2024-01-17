#include "WCSimPrimaryGeneratorMessenger.hh"
#include "WCSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ios.hh"

WCSimPrimaryGeneratorMessenger::WCSimPrimaryGeneratorMessenger(WCSimPrimaryGeneratorAction* pointerToAction)
:myAction(pointerToAction)
{
  mydetDirectory = new G4UIdirectory("/mygen/");
  mydetDirectory->SetGuidance("WCSim detector control commands.");

  genCmd = new G4UIcmdWithAString("/mygen/generator",this);
  genCmd->SetGuidance("Select primary generator.");

  genCmd->SetGuidance(" Available generators : muline, gun, laser, gps, rootracker, gamma-conversion");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("muline");
  genCmd->SetCandidates("muline gun laser gps rootracker, gamma-conversion");

  fileNameCmd = new G4UIcmdWithAString("/mygen/vecfile",this);
  fileNameCmd->SetGuidance("Select the file of vectors.");
  fileNameCmd->SetGuidance(" Enter the file name of the vector file");
  fileNameCmd->SetParameterName("fileName",true);
  fileNameCmd->SetDefaultValue("inputvectorfile");

  genSet = false;

  //C. Vilela: Adding PMTPoisson for generating photoelectrons directly on PMTs according to a Poisson distribution.
  poisCmd = new G4UIcmdWithABool("/mygen/pmtPoisson",this);
  poisCmd->SetGuidance("Flag for generating photoelectrons directly on PMTs according to a Poisson distribution. These PE's will be generated in addition to light produce by any particles generated. Set dark rate to 0 and do not generate any particles for events with only Poisson PE's.");
  poisCmd->SetGuidance("Set poisson mean with /mygen/poissonMean");
  poisCmd->SetParameterName("pmtPoisson", true);

  poisMeanCmd = new G4UIcmdWithADouble("/mygen/poissonMean",this);
  poisMeanCmd->SetGuidance("Set Poisson mean to be used with /mygen/pmtPoisson. Defaults to 1.");
  poisMeanCmd->SetParameterName("poissonMean", true);
  poisMeanCmd->SetDefaultValue(1);

  multiVtxFileNameCmd = new G4UIcmdWithAString("/mygen/gps/multiVtxFile",this);
  multiVtxFileNameCmd->SetGuidance("Select the file of for multi-vertex file.");
  multiVtxFileNameCmd->SetGuidance(" Enter the file name of the txt file");
  multiVtxFileNameCmd->SetParameterName("fileName",true);
  multiVtxFileNameCmd->SetDefaultValue("multiVtxFile");
}

WCSimPrimaryGeneratorMessenger::~WCSimPrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
  delete multiVtxFileNameCmd;
}

void WCSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  {
    genSet = true;
    if (newValue == "muline")
    {
      myAction->SetMulineEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetRootrackerEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
    }
    else if ( newValue == "gun")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(true);
      myAction->SetRootrackerEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
    }
    else if ( newValue == "laser")   //T. Akiri: Addition of laser
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetRootrackerEvtGenerator(false);
      myAction->SetLaserEvtGenerator(true);
      myAction->SetGPSEvtGenerator(false);
    }
    else if ( newValue == "gps")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetRootrackerEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetNeedConversion(false);
    }
    else if ( newValue == "rootracker")   //M. Scott: Addition of Rootracker events
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetRootrackerEvtGenerator(true);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(false);
    }
    else if ( newValue == "gamma-conversion")
    {
      myAction->SetMulineEvtGenerator(false);
      myAction->SetRootrackerEvtGenerator(false);
      myAction->SetGunEvtGenerator(false);
      myAction->SetLaserEvtGenerator(false);
      myAction->SetGPSEvtGenerator(true);
      myAction->SetNeedConversion(true);
    }
  }

  if( command == fileNameCmd )
  {
    if(genSet){
        if(myAction->IsUsingRootrackerEvtGenerator()){
            myAction->OpenRootrackerFile(newValue);
        }
        else{
            myAction->OpenVectorFile(newValue);
        }
        G4cout << "Input vector file set to " << newValue << G4endl;
    }
    else{
        G4cout << "Generator has not been set, guessing input vector file is NOT in the Rootracker format - this will crash if you are using a Rootracker input file" << G4endl;
        G4cout << "Please put the '/mygen/generator' command above the '/mygen/vecfile' command in the mac file." << G4endl;
    }
  }

  if( command == poisCmd )
    {
      if ( poisCmd->GetNewBoolValue(newValue) ){
	myAction->SetPoissonPMT(true);
	G4cout << "Running with PoissonPMT flag. Photoelectrons will be generated directly on the PMTs according to a Poisson distribuition. Any hits resulting from physics generated elsewhere will be discarded !!!" << G4endl;
      }
      else myAction->SetPoissonPMT(false);
    }
  
  if( command == poisMeanCmd )
    {
      myAction->SetPoissonPMTMean(poisMeanCmd->GetNewDoubleValue(newValue));
      G4cout << "PoissonPMT mean set to: " << poisMeanCmd->GetNewDoubleValue(newValue) << G4endl;
    }

  if( command == multiVtxFileNameCmd )
  {
    G4cout << "Multi-vertex input set to: " << newValue << G4endl;
    myAction->SetMultiVertexInput(newValue);
  }
}

G4String WCSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->IsUsingMulineEvtGenerator())
      { cv = "muline"; }
    else if(myAction->IsUsingGunEvtGenerator())
      { cv = "gun"; }
    else if(myAction->IsUsingLaserEvtGenerator())
      { cv = "laser"; }   //T. Akiri: Addition of laser
    else if(myAction->IsUsingGPSEvtGenerator())
      { cv = myAction->NeedsConversion() ? "gamma-conversion" : "gps"; }
    else if(myAction->IsUsingRootrackerEvtGenerator())
      { cv = "rootracker"; }   //M. Scott: Addition of Rootracker events
  }
  
  return cv;
}

