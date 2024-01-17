#include "WCSimPrimaryGeneratorAction.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimPrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4EventManager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <fstream>
#include <vector>
#include <string>

//#include "WCSimEnumerations.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

#include "TRandom3.h"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using std::vector;
using std::string;
using std::fstream;

vector<string> tokenize( string separators, string input );

inline vector<string> readInLine(fstream& inFile, int lineSize, char* inBuf)
{
    // Read in line break it up into tokens
    inFile.getline(inBuf,lineSize);
    return tokenize(" $", inBuf);
}

inline float atof( const string& s ) {return std::atof( s.c_str() );}
inline int   atoi( const string& s ) {return std::atoi( s.c_str() );}

WCSimPrimaryGeneratorAction::WCSimPrimaryGeneratorAction(
					  WCSimDetectorConstruction* myDC)
  :myDetector(myDC), vectorFileName("")
{
  //T. Akiri: Initialize GPS to allow for the laser use 
  MyGPS = new G4GeneralParticleSource();

  // Initialize to zero
  mode = UNKNOWN;    //0;
  vtxvol = 0;
  vtx = G4ThreeVector(0.,0.,0.);
  nuEnergy = 0.;
  _counterRock=0; // counter for generated in Rock
  _counterCublic=0; // counter generated
  
  //---Set defaults. Do once at beginning of session.
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.0));

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  particleGun->
    SetParticleDefinition(particleTable->FindParticle(particleName="mu+"));

  particleGun->
    SetParticlePosition(G4ThreeVector(0.*m,0.*m,0.*m));
    
  messenger = new WCSimPrimaryGeneratorMessenger(this);
  useMulineEvt = true;
  useGunEvt    = false;
  useLaserEvt  = false;
  useGPSEvt    = false;
  useRootrackerEvt = false;
  
  fEvNum = 0;
  fInputRootrackerFile = NULL;
  fNEntries = 1;
  
  needConversion = false;
  foundConversion = true;

  useMultiVertex = false;
}

WCSimPrimaryGeneratorAction::~WCSimPrimaryGeneratorAction()
{
    if (IsGeneratingVertexInRock()){
        G4cout << "Fraction of Rock volume is : " << G4endl;
        G4cout << " Random number generated in Rock / in Cublic = " 
            << _counterRock << "/" << _counterCublic 
            << " = " << _counterRock/(G4double)_counterCublic << G4endl;
    }
    inputFile.close();

    if(useRootrackerEvt) delete fRooTrackerTree;

    delete particleGun;
    delete MyGPS;   //T. Akiri: Delete the GPS variable
    delete messenger;
}

void WCSimPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

    // We will need a particle table
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    // Temporary kludge to turn on/off vector text format 
    G4bool useNuanceTextFormat = true;

  // Do for every event

  if (useMulineEvt)
  { 

    if ( !inputFile.is_open() )
    {
      G4cout << "Set a vector file using the command /mygen/vecfile name"
	     << G4endl;
      exit(-1);
    }

    //
    // Documentation describing the nuance text format can be found here: 
    // http://neutrino.phy.duke.edu/nuance-format/
    //
    // The format must be strictly adhered to for it to be processed correctly.
    // The lines and their meanings from begin through info are fixed, and then
    // a variable number of tracks may follow.
    //
    if (useNuanceTextFormat)
      {
	const int lineSize=100;
	char      inBuf[lineSize];
	vector<string> token(1);
	
	token = readInLine(inputFile, lineSize, inBuf);
	  
        if (token.size() == 0) 
	  {
	    G4cout << "end of nuance vector file!" << G4endl;
	  }
	else if (token[0] != "begin")
	  {
	    G4cout << "unexpected line begins with " << token[0] << G4endl;
	  }
	else   // normal parsing begins here
	  {
	    // Read the nuance line (ignore value now)

	    token = readInLine(inputFile, lineSize, inBuf);
	    mode = BEAM;       //atoi(token[1]);    //break backwards compatibility, should deprecate Nuance though

	    // Read the Vertex line
	    token = readInLine(inputFile, lineSize, inBuf);
	    vtx = G4ThreeVector(atof(token[1])*cm,
				atof(token[2])*cm,
				atof(token[3])*cm);
	    
            // true : Generate vertex in Rock , false : Generate vertex in WC tank
            SetGenerateVertexInRock(false);

	    // Next we read the incoming neutrino and target
	    
	    // First, the neutrino line

	    token=readInLine(inputFile, lineSize, inBuf);
	    beampdg = atoi(token[1]);
	    beamenergy = atof(token[2])*MeV;
	    beamdir = G4ThreeVector(atof(token[3]),
				    atof(token[4]),
				    atof(token[5]));

	    // Now read the target line

	    token=readInLine(inputFile, lineSize, inBuf);
	    targetpdg = atoi(token[1]);
	    targetenergy = atof(token[2])*MeV;
	    targetdir = G4ThreeVector(atof(token[3]),
				      atof(token[4]),
				      atof(token[5]));

	    // Read the info line, basically a dummy
	    token=readInLine(inputFile, lineSize, inBuf);
	    G4cout << "Vector File Record Number " << token[2] << G4endl;
            vecRecNumber = atoi(token[2]);
	    
	    // Now read the outgoing particles
	    // These we will simulate.


	    while ( token=readInLine(inputFile, lineSize, inBuf),
		    token[0] == "track" )
	      {
		// We are only interested in the particles
		// that leave the nucleus, tagged by "0"


		if ( token[6] == "0")
		  {
		    G4int pdgid = atoi(token[1]);
		    G4double energy = atof(token[2])*MeV;
		    G4ThreeVector dir = G4ThreeVector(atof(token[3]),
						      atof(token[4]),
						      atof(token[5]));
		    particleGun->
		      SetParticleDefinition(particleTable->
					    FindParticle(pdgid));
		    G4double mass = 
		      particleGun->GetParticleDefinition()->GetPDGMass();

		    G4double ekin = energy - mass;

		    particleGun->SetParticleEnergy(ekin);
		    //G4cout << "Particle: " << pdgid << " KE: " << ekin << G4endl;
		    particleGun->SetParticlePosition(vtx);
		    particleGun->SetParticleMomentumDirection(dir);
		    particleGun->GeneratePrimaryVertex(anEvent);
		  }
	      }
	  }
      }
    else 
      {    // old muline format  
	inputFile >> nuEnergy >> energy >> xPos >> yPos >> zPos 
		  >> xDir >> yDir >> zDir;
	
	G4double random_z = ((myDetector->GetWaterTubePosition())
			     - .5*(myDetector->GetWaterTubeLength()) 
			     + 1.*m + 15.0*m*G4UniformRand())/m;
	zPos = random_z;
	G4ThreeVector vtx = G4ThreeVector(xPos, yPos, random_z);
	G4ThreeVector dir = G4ThreeVector(xDir,yDir,zDir);

	particleGun->SetParticleEnergy(energy*MeV);
	particleGun->SetParticlePosition(vtx);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->GeneratePrimaryVertex(anEvent);
      }
  }

    else if (useRootrackerEvt)
    { 
        if ( !fInputRootrackerFile->IsOpen() )
        {
            G4cout << "Set a Rootracker vector file using the command /mygen/vecfile name"
                << G4endl;
            return;
        }

        //Generate 1 event

        //Load event from file
        if(fEvNum == 0){
            fSettingsTree->SetBranchAddress("NuBeamAng",&fNuBeamAng);
            fSettingsTree->SetBranchAddress("DetRadius",&fNuPrismRadius);
            fSettingsTree->SetBranchAddress("NuIdfdPos",fNuPlanePos);
        }
        if (fEvNum<fNEntries){
            fRooTrackerTree->GetEntry(fEvNum);
            fSettingsTree->GetEntry(fEvNum);
            fEvNum++;
        }
        else{
            G4cout << "End of File" << G4endl; 
            return; 
        }

        // Get the neutrino direction
        xDir=fTmpRootrackerVtx->StdHepP4[0][0];
        yDir=fTmpRootrackerVtx->StdHepP4[0][1];
        zDir=fTmpRootrackerVtx->StdHepP4[0][2];

        // Calculate offset from neutrino generation plane to centre of nuPRISM detector (in metres)
        double z_offset = fNuPlanePos[2]/100.0;
        double y_offset = 0;//(fNuPrismRadius/zDir)*yDir;
      	double x_offset = fNuPlanePos[0]/100.0;

        //Subtract offset to get interaction position in WCSim coordinates
        xPos = fTmpRootrackerVtx->EvtVtx[0] - x_offset;
        yPos = fTmpRootrackerVtx->EvtVtx[1] - y_offset;
        zPos = fTmpRootrackerVtx->EvtVtx[2] - z_offset;

        //Check if event is outside detector; skip to next event if so; keep
        //loading events until one is found within the detector or there are
        //no more interaction to simulate for this event.
        //The current neut vector files do not correspond directly to the detector dimensions, so only keep those events within the detector
        while (sqrt(pow(xPos,2)+pow(zPos,2))*m > (myDetector->GetWCIDDiameter()/2.) || (abs(yPos*m - myDetector->GetWCIDVerticalPosition()) > (myDetector->GetWCIDHeight()/2.))){
            //Load another event
            if (fEvNum<fNEntries){
                fRooTrackerTree->GetEntry(fEvNum);
                G4cout << "Skipped event# " << fEvNum - 1 << " (event vertex outside detector)" << G4endl;
                fEvNum++;
            }
            else{
                G4cout << "End of File" << G4endl; 
                return; 
            }
            //Convert coordinates
            xPos = fTmpRootrackerVtx->EvtVtx[0] - x_offset;
            yPos = fTmpRootrackerVtx->EvtVtx[1] - y_offset; 
            zPos = fTmpRootrackerVtx->EvtVtx[2] - z_offset;
        } 

        //Generate particles
        //i = 0 is the neutrino
        //i = 1 is the target nucleus
        //i = 2 is the target nucleon
        //i > 2 are the outgoing particles

        // First simulate the incoming neutrino
        xDir=fTmpRootrackerVtx->StdHepP4[0][0];
        yDir=fTmpRootrackerVtx->StdHepP4[0][1];
        zDir=fTmpRootrackerVtx->StdHepP4[0][2];

        double momentumGeV=sqrt((xDir*xDir)+(yDir*yDir)+(zDir*zDir))*GeV;
        double momentum=sqrt((xDir*xDir)+(yDir*yDir)+(zDir*zDir));

        G4ThreeVector vtx = G4ThreeVector(xPos*m, yPos*m, zPos*m);
        G4ThreeVector dir = G4ThreeVector(-xDir, -yDir, -zDir);

        dir = dir*(momentumGeV/momentum);

        particleGun->SetParticleDefinition(particleTable->FindParticle(fTmpRootrackerVtx->StdHepPdgTemp[0]));
        double kin_energy = momentumGeV;//fabs(fTmpRootrackerVtx->StdHepP4[i][3])*GeV - particleGun->GetParticleDefinition()->GetPDGMass();
        particleGun->SetParticleEnergy(kin_energy);
        particleGun->SetParticlePosition(vtx);
        particleGun->SetParticleMomentumDirection(dir);
        // Will want to include some beam time structure at some point, but not needed at the moment since we only simulate 1 interaction per events
        // particleGun->SetParticleTime(time);
        particleGun->GeneratePrimaryVertex(anEvent);  //Place vertex in stack

        // Now simulate the outgoing particles
        for (int i = 0; i < fTmpRootrackerVtx->StdHepN; i++){

            // Skip the initial neutrino, target nucleus and target nucleon
            if( i < 3){
                int pdg = abs(fTmpRootrackerVtx->StdHepPdgTemp[i]);
                if(pdg > 100000 || pdg == 12 || pdg == 14 || pdg == 2112 || pdg == 2212){
                    continue;
                }
            }

            xDir=fTmpRootrackerVtx->StdHepP4[i][0];
            yDir=fTmpRootrackerVtx->StdHepP4[i][1];
            zDir=fTmpRootrackerVtx->StdHepP4[i][2];

            momentumGeV=sqrt((xDir*xDir)+(yDir*yDir)+(zDir*zDir))*GeV;
            momentum=sqrt((xDir*xDir)+(yDir*yDir)+(zDir*zDir));

            vtx.setX(xPos*m);
            vtx.setY(yPos*m);
            vtx.setZ(zPos*m);

            dir.setX(xDir);
            dir.setY(yDir);
            dir.setZ(zDir);

            dir = dir*(momentumGeV/momentum);

            particleGun->SetParticleDefinition(particleTable->FindParticle(fTmpRootrackerVtx->StdHepPdgTemp[i]));

            double kin_energy = fabs(fTmpRootrackerVtx->StdHepP4[i][3])*GeV - particleGun->GetParticleDefinition()->GetPDGMass();

            particleGun->SetParticleEnergy(kin_energy);
            particleGun->SetParticlePosition(vtx);
            particleGun->SetParticleMomentumDirection(dir);
            // Will want to include some beam time structure at some point, but not needed at the moment since we only simulate 1 interaction per events
            // particleGun->SetParticleTime(time);
            particleGun->GeneratePrimaryVertex(anEvent);  //Place vertex in stack
        }
    }

  else if (useGunEvt)
  {      // manual gun operation
    particleGun->GeneratePrimaryVertex(anEvent);

    G4ThreeVector P  =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
    G4ThreeVector vtx=anEvent->GetPrimaryVertex()->GetPosition();
    G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
    G4int pdg        =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();

    G4ThreeVector dir  = P.unit();
    G4double E         = std::sqrt((P.dot(P))+(m*m));

    mode            = PARTICLEGUN;

//     particleGun->SetParticleEnergy(E);
//     particleGun->SetParticlePosition(vtx);
//     particleGun->SetParticleMomentumDirection(dir);
    SetVtx(vtx);
    SetBeamEnergy(E);
    SetBeamDir(dir);
    SetBeamPDG(pdg);
  }
  else if (useLaserEvt)
    {
      //T. Akiri: Create the GPS LASER event
      MyGPS->GeneratePrimaryVertex(anEvent);
      
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vtx =anEvent->GetPrimaryVertex()->GetPosition();
      G4double m       =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass(); // will be 0 for photon anyway, but for other gps particles not
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      //G4double E         = std::sqrt((P.dot(P)));
      G4double E         = std::sqrt((P.dot(P))+(m*m));
      //std::cout << "Energy " << E << " eV " << std::endl;

      mode            = LASER; //actually could also be particle gun here. Gps and laser will be separate soon!!

      SetVtx(vtx);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);
    }
  else if (useGPSEvt)
    {
      if (useMultiVertex) // minimal workaround to use multi-vertex in the same run
      {
        G4int evtID = anEvent->GetEventID();
        G4int vtxID = evtID%vertexEntries.size(); // loop over vertex entry if necessary
        //G4cout<<"evtID = "<<evtID<<" vertexEntries.size() = "<<vertexEntries.size()<<" vtxID = "<<vtxID<<" vertexEntries[vtxID].pdgId = "<<vertexEntries[vtxID].pdgId<<G4endl;
        if (vertexEntries[vtxID].pdgId==0) // special treatment for optical photon
          MyGPS->GetCurrentSource()->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton"));
        else MyGPS->GetCurrentSource()->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(vertexEntries[vtxID].pdgId));

        G4ThreeVector position(vertexEntries[vtxID].xpos*m,vertexEntries[vtxID].ypos*m,vertexEntries[vtxID].zpos*m);
        MyGPS->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
        MyGPS->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);

        G4ThreeVector direction(vertexEntries[vtxID].xdir,vertexEntries[vtxID].ydir,vertexEntries[vtxID].zdir);
        MyGPS->GetCurrentSource()->GetAngDist()->SetAngDistType("planar");
        MyGPS->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(direction);

        MyGPS->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
        MyGPS->GetCurrentSource()->GetEneDist()->SetMonoEnergy(vertexEntries[vtxID].ke*MeV);

        MyGPS->GetCurrentSource()->SetNumberOfParticles(vertexEntries[vtxID].nPart);
      }

      MyGPS->GeneratePrimaryVertex(anEvent);
      G4ThreeVector P   =anEvent->GetPrimaryVertex()->GetPrimary()->GetMomentum();
      G4ThreeVector vertex =anEvent->GetPrimaryVertex()->GetPosition();
      G4double mass        =anEvent->GetPrimaryVertex()->GetPrimary()->GetMass();
      G4int pdg         =anEvent->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
      
      G4ThreeVector dir  = P.unit();
      G4double E         = std::sqrt((P.dot(P))+(mass*mass));
      
      SetVtx(vertex);
      SetBeamEnergy(E);
      SetBeamDir(dir);
      SetBeamPDG(pdg);

      //G4cout<<"vtx = "<<vertex.x()<<" "<<vertex.y()<<" "<<vertex.z()<<" pdg = "<<pdg<<" KE = "<<E-mass<<" dir = "<<dir.x()<<" "<<dir.y()<<" "<<dir.z()<<" "<<G4endl;

      if(needConversion) {
          G4PrimaryParticle *primaryParticle = anEvent->GetPrimaryVertex(0)->GetPrimary(0);
          G4ThreeVector tmpDir(1, 0, 0);
          G4ThreeVector tmpPos(0, 0, 0);
          particleGun->SetParticleDefinition(primaryParticle->GetG4code());
          particleGun->SetParticleEnergy(primaryParticle->GetKineticEnergy());
          particleGun->SetParticlePosition(tmpPos);
          particleGun->SetParticleMomentumDirection(tmpDir);
          foundConversion = false;
          while (!foundConversion) {
              G4Event *tmpEvent = new G4Event(-1);
              particleGun->GeneratePrimaryVertex(tmpEvent);
              G4EventManager::GetEventManager()->ProcessOneEvent(tmpEvent);
              delete tmpEvent;
          }
          G4ThreeVector newDir = primaryParticle->GetMomentumDirection();
          if (!newDir.isParallel(tmpDir, 1e-5)) {
              G4ThreeVector rotationAxis = tmpDir.cross(newDir);
              rotationAxis = rotationAxis / rotationAxis.mag();
              double rotationAngle = acos(newDir.dot(tmpDir));
              for (int i = 0; i < 2; i++) {
                  conversionProductMomentum[i].rotate(rotationAngle, rotationAxis);
              }
          }
          primaryParticle->SetParticleDefinition(conversionProductParticle[0]);
          primaryParticle->SetMomentum(conversionProductMomentum[0].getX(),
                                       conversionProductMomentum[0].getY(),
                                       conversionProductMomentum[0].getZ());
          G4PrimaryParticle *secondProduct = new G4PrimaryParticle(conversionProductParticle[1],
                                                                   conversionProductMomentum[1].getX(),
                                                                   conversionProductMomentum[1].getY(),
                                                                   conversionProductMomentum[1].getZ());
          anEvent->GetPrimaryVertex(0)->SetPrimary(secondProduct);
      }
    }
}

void WCSimPrimaryGeneratorAction::SaveOptionsToOutput(WCSimRootOptions * wcopt)
{
  if(useMulineEvt)
    wcopt->SetVectorFileName(vectorFileName);
  else
    wcopt->SetVectorFileName("");
  wcopt->SetGeneratorType(GetGeneratorTypeString());
}

G4String WCSimPrimaryGeneratorAction::GetGeneratorTypeString()
{
  if(useMulineEvt)
    return "muline";
  else if(useGunEvt)
    return "gun";
  else if(useGPSEvt)
    return "gps";
  else if(useLaserEvt)
    return "laser";
  else if(useRootrackerEvt)
    return "rooTrackerEvt";
  return "";
}

// Returns a vector with the tokens
vector<string> tokenize( string separators, string input ) 
{
  std::size_t startToken = 0, endToken; // Pointers to the token pos
  vector<string> tokens;  // Vector to keep the tokens
  
  if( separators.size() > 0 && input.size() > 0 ) 
    {

        while( startToken < input.size() )
        {
            // Find the start of token
            startToken = input.find_first_not_of( separators, startToken );

            // If found...
            if( startToken != input.npos ) 
            {
                // Find end of token
                endToken = input.find_first_of( separators, startToken );
                if( endToken == input.npos )
                    // If there was no end of token, assign it to the end of string
                    endToken = input.size();

                // Extract token
                tokens.push_back( input.substr( startToken, endToken - startToken ) );

                // Update startToken
                startToken = endToken;
            }
        }
    }

    return tokens;
}

void WCSimPrimaryGeneratorAction::OpenRootrackerFile(G4String fileName)
{
    if (fInputRootrackerFile) fInputRootrackerFile->Delete();

    fInputRootrackerFile = TFile::Open(fileName.data());
    if (!fInputRootrackerFile){ 
        G4cout << "Cannot open: " << fileName << G4endl; 
        exit(1);
    }

    fRooTrackerTree = (TTree*) fInputRootrackerFile->Get("nRooTracker");
    fSettingsTree = (TTree*) fInputRootrackerFile->Get("Settings");
    if (!fRooTrackerTree){
        G4cout << "File: " << fileName << " does not contain a Rootracker nRooTracker tree - please check you intend to process Rootracker events" << G4endl;
        exit(1);
    }
    fNEntries=fRooTrackerTree->GetEntries();

    fTmpRootrackerVtx = new NRooTrackerVtx();
    SetupBranchAddresses(fTmpRootrackerVtx); //link fTmpRootrackerVtx and current input file

    fSettingsTree->SetBranchAddress("NuBeamAng",&fNuBeamAng);
    fSettingsTree->SetBranchAddress("DetRadius",&fNuPrismRadius);
    fSettingsTree->SetBranchAddress("NuIdfdPos",fNuPlanePos);

}

void WCSimPrimaryGeneratorAction::SetupBranchAddresses(NRooTrackerVtx* nrootrackervtx){

    // Set up branch address for rooTrackerVertex tree in nuPRISM files
    fRooTrackerTree->SetBranchAddress("EvtCode",        &(nrootrackervtx->EvtCode) );
    fRooTrackerTree->SetBranchAddress("EvtNum",         &(nrootrackervtx->EvtNum) );
    fRooTrackerTree->SetBranchAddress("EvtXSec",        &(nrootrackervtx->EvtXSec) );
    fRooTrackerTree->SetBranchAddress("EvtDXSec",       &(nrootrackervtx->EvtDXSec) );
    fRooTrackerTree->SetBranchAddress("EvtWght",        &(nrootrackervtx->EvtWght) );
    fRooTrackerTree->SetBranchAddress("EvtProb",        &(nrootrackervtx->EvtProb) );
    fRooTrackerTree->SetBranchAddress("EvtVtx",          (nrootrackervtx->EvtVtx) );
    fRooTrackerTree->SetBranchAddress("StdHepN",        &(nrootrackervtx->StdHepN) );
    fRooTrackerTree->SetBranchAddress("StdHepPdg",       (nrootrackervtx->StdHepPdgTemp) );
    fRooTrackerTree->SetBranchAddress("StdHepStatus",    (nrootrackervtx->StdHepStatusTemp) );
    fRooTrackerTree->SetBranchAddress("StdHepX4",        (nrootrackervtx->StdHepX4) );
    fRooTrackerTree->SetBranchAddress("StdHepP4",        (nrootrackervtx->StdHepP4) );
    fRooTrackerTree->SetBranchAddress("StdHepPolz",      (nrootrackervtx->StdHepPolz) );
    fRooTrackerTree->SetBranchAddress("StdHepFd",        (nrootrackervtx->StdHepFdTemp) );
    fRooTrackerTree->SetBranchAddress("StdHepLd",        (nrootrackervtx->StdHepLdTemp) );
    fRooTrackerTree->SetBranchAddress("StdHepFm",        (nrootrackervtx->StdHepFmTemp) );
    fRooTrackerTree->SetBranchAddress("StdHepLm",        (nrootrackervtx->StdHepLmTemp) );

    // NEUT > v5.0.7 && MCP > 1 (>10a)
    fRooTrackerTree->SetBranchAddress("NEnvc",          &(nrootrackervtx->NEnvc)    );
    fRooTrackerTree->SetBranchAddress("NEipvc",          (nrootrackervtx->NEipvcTemp)   );
    fRooTrackerTree->SetBranchAddress("NEpvc",           (nrootrackervtx->NEpvc)    );
    fRooTrackerTree->SetBranchAddress("NEiorgvc",        (nrootrackervtx->NEiorgvcTemp) );
    fRooTrackerTree->SetBranchAddress("NEiflgvc",        (nrootrackervtx->NEiflgvcTemp) );
    fRooTrackerTree->SetBranchAddress("NEicrnvc",        (nrootrackervtx->NEicrnvcTemp) );

    fRooTrackerTree->SetBranchAddress("NEnvert",        &(nrootrackervtx->NEnvert)  );
    fRooTrackerTree->SetBranchAddress("NEposvert",       (nrootrackervtx->NEposvert) );
    fRooTrackerTree->SetBranchAddress("NEiflgvert",      (nrootrackervtx->NEiflgvertTemp) );
    fRooTrackerTree->SetBranchAddress("NEnvcvert",      &(nrootrackervtx->NEnvcvert) );
    fRooTrackerTree->SetBranchAddress("NEdirvert",       (nrootrackervtx->NEdirvert) );
    fRooTrackerTree->SetBranchAddress("NEabspvert",      (nrootrackervtx->NEabspvertTemp) );
    fRooTrackerTree->SetBranchAddress("NEabstpvert",     (nrootrackervtx->NEabstpvertTemp) );
    fRooTrackerTree->SetBranchAddress("NEipvert",        (nrootrackervtx->NEipvertTemp)  );
    fRooTrackerTree->SetBranchAddress("NEiverti",        (nrootrackervtx->NEivertiTemp)  );
    fRooTrackerTree->SetBranchAddress("NEivertf",        (nrootrackervtx->NEivertfTemp)  );
    // end Rootracker > v5.0.7 && MCP > 1 (>10a)

    // NEUT > v5.1.2 && MCP >= 5
    fRooTrackerTree->SetBranchAddress("NEcrsx",          &(nrootrackervtx->NEcrsx)    );
    fRooTrackerTree->SetBranchAddress("NEcrsy",          &(nrootrackervtx->NEcrsy)    );
    fRooTrackerTree->SetBranchAddress("NEcrsz",          &(nrootrackervtx->NEcrsz)    );
    fRooTrackerTree->SetBranchAddress("NEcrsphi",        &(nrootrackervtx->NEcrsphi)    );

    fRooTrackerTree->SetBranchAddress("NuParentDecMode", &(nrootrackervtx->NuParentDecMode));
    fRooTrackerTree->SetBranchAddress("NuParentPdg",     &(nrootrackervtx->NuParentPdg));
    fRooTrackerTree->SetBranchAddress("NuNg",            &(nrootrackervtx->NuNg));
    fRooTrackerTree->SetBranchAddress("NuGp",             (nrootrackervtx->NuGp));
    fRooTrackerTree->SetBranchAddress("NuGv",             (nrootrackervtx->NuGv));
    fRooTrackerTree->SetBranchAddress("NuGpid",           (nrootrackervtx->NuGpid));
    fRooTrackerTree->SetBranchAddress("NuGmat",           (nrootrackervtx->NuGmat));
    fRooTrackerTree->SetBranchAddress("NuGmec",           (nrootrackervtx->NuGmec));
    fRooTrackerTree->SetBranchAddress("NuGdistc",         (nrootrackervtx->NuGdistc));
    fRooTrackerTree->SetBranchAddress("NuGdistal",        (nrootrackervtx->NuGdistal));

}
void WCSimPrimaryGeneratorAction::CopyRootrackerVertex(NRooTrackerVtx* nrootrackervtx){
    nrootrackervtx->Copy(fTmpRootrackerVtx);
    nrootrackervtx->TruthVertexID = -999;
}

void WCSimPrimaryGeneratorAction::SetMultiVertexInput(G4String fileName)
{
  useMultiVertex = true;
  vertexEntries.clear();

  std::ifstream Data(fileName.data(),std::ios_base::in);
  if (!Data)
  {
    G4cout<<"Vertex data file "<<fileName<<" could not be opened --> Exiting..."<<G4endl;
    exit(-1);
  }
  else
    G4cout<<"Vertex data file "<<fileName<<" is opened to read ..."<<G4endl;

  std::string str, tmp;
	G4int Column=0;
	while (std::getline(Data, str)) {
		if (str=="#DATASTART") break;
	}
	std::ifstream::pos_type SavePoint = Data.tellg();
	std::getline(Data, str);
	std::istringstream stream(str);
	while (std::getline(stream,tmp,' ')) Column++;
	if (Column!=9)
  {
    G4cerr<<"Number of column = "<<Column<<" which is not equal to 9. "<<G4endl;
    G4cerr<<"Inappropriate input --> Exiting..."<<G4endl;
    exit(-1);
  }
  Data.seekg(SavePoint);

  G4int nVerticesRead = 0;
  G4int pdgId, nPart;
  G4double xpos, ypos, zpos, xdir, ydir, zdir, ke;
  while (!Data.eof())
  {
    Data>>pdgId>>xpos>>ypos>>zpos>>xdir>>ydir>>zdir>>ke>>nPart;
    VertexEntry ve;
    ve.pdgId = pdgId;
    ve.xpos = xpos;
    ve.ypos = ypos;
    ve.zpos = zpos;
    ve.xdir = xdir;
    ve.ydir = ydir;
    ve.zdir = zdir;
    ve.ke = ke;
    ve.nPart = nPart;
    nVerticesRead++;

    vertexEntries.push_back(ve);

    //G4cout<<"pdgId = "<<pdgId<<" "<<xpos<<" "<<ypos<<" "<<zpos<<G4endl;
  }
  Data.close();

  if (nVerticesRead==0)
  {
    G4cerr<<"Number of vertices = 0. "<<G4endl;
    G4cerr<<"Inappropriate input --> Exiting..."<<G4endl;
    exit(-1);
  }
  else
    G4cout<<nVerticesRead<<" vertices are read into memory"<<G4endl;
}