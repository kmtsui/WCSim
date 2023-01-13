//  -*- mode:c++; tab-width:4;  -*-
#include "WCSimDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalSurface.hh"
#include "G4UserLimits.hh"
#include "G4ReflectionFactory.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "WCSimTuningParameters.hh" //jl145

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//#define DEBUG

/***********************************************************
 *
 * This file contains the functions which construct a 
 * cylindrical WC detector.  It used by both the SK and 
 * LBNE WC detector modes.  It is called in the Construct()
 * method in WCSimDetectorConstruction.cc.
 *
 * Sourcefile for the WCSimDetectorConstruction class
 *
 ***********************************************************/


G4LogicalVolume* WCSimDetectorConstruction::ConstructSinglemPMTWorld()
{
    G4cout << "**** Building SinglemPMTWorld Detector ****" << G4endl;

  //-----------------------------------------------------
  // Positions
  //-----------------------------------------------------

  debugMode = false;

  //SetHyperK_HybridmPMT10PCGeometry();
  
  WCPosition=0.;//Set the WC tube offset to zero

  WCIDRadius = WCIDDiameter/2.;
  // the number of regular cell in phi direction:
  WCBarrelRingNPhi     = (G4int)(WCBarrelNumPMTHorizontal/WCPMTperCellHorizontal); 
  // the part of one ring, that is covered by the regular cells: 
  totalAngle  = 2.0*pi*rad*(WCBarrelRingNPhi*WCPMTperCellHorizontal/WCBarrelNumPMTHorizontal) ;
  // angle per regular cell:
  dPhi        =  totalAngle/ WCBarrelRingNPhi;
  // it's height:
  barrelCellHeight  = (WCIDHeight-2.*WCBarrelPMTOffset)/WCBarrelNRings;
  // the height of all regular cells together:
  mainAnnulusHeight = WCIDHeight -2.*WCBarrelPMTOffset -2.*barrelCellHeight;
  
  //TF: has to change for mPMT vessel:

  G4cout << WCIDRadius << ", " << WCPMTExposeHeight << ", " << mPMT_vessel_cyl_height << ", " << mPMT_vessel_radius << G4endl;
#ifdef DEBUG
	G4cout << "HYBRID = " << hybrid << G4endl;
#endif
  
  if(hybrid){
	if(mPMT_vessel_cyl_height + mPMT_vessel_radius < 1.*mm)
	  innerAnnulusRadius = WCIDRadius - std::max(WCPMTExposeHeight,WCPMTExposeHeight2)-1.*mm;
	else
	  innerAnnulusRadius = WCIDRadius - std::max(WCPMTExposeHeight,(mPMT_vessel_cyl_height + mPMT_vessel_radius)) -1.*mm;
  }
  else{
	if(mPMT_vessel_cyl_height + mPMT_vessel_radius < 1.*mm)
	  innerAnnulusRadius = WCIDRadius - WCPMTExposeHeight -1.*mm;
	else
	  innerAnnulusRadius = WCIDRadius - (mPMT_vessel_cyl_height + mPMT_vessel_radius) -1.*mm;
  }
  //TF: need to add a Polyhedra on the other side of the outerAnnulusRadius for the OD
  outerAnnulusRadius = WCIDRadius + WCBlackSheetThickness + 1.*mm;//+ Stealstructure etc.

  // the radii are measured to the center of the surfaces
  // (tangent distance). Thus distances between the corner and the center are bigger.
  //BQ: Updated with new HK OD size (2020/12/06). Simply assume no tyvek thickness or dead space.
  WCODTyvekSheetThickness = 0.*m;
  WCODDeadSpace = 0.*m;
  WCODHeightWaterDepth = 2*m + 0.6*m;//OD + support structure
  WCODLateralWaterDepth = 1*m + 0.6*m;//OD + support structure
  WCLength    = WCIDHeight + 2*(WCODHeightWaterDepth + WCBlackSheetThickness + WCODDeadSpace + WCODTyvekSheetThickness);
  WCRadius    = (outerAnnulusRadius + WCODLateralWaterDepth)/cos(dPhi/2.) ;
  //WCLength    = WCIDHeight + (WCODHeight+WCOuterStructure)*2;	//jl145 - reflects top veto blueprint, cf. Farshid Feyzi
  //WCRadius    = (WCIDDiameter/2. + WCBlackSheetThickness + (WCODRadius+WCODHeight)))/cos(dPhi/2.) ; // BQ: Updated with new HK OD size (2020/12/06) 
 
  // now we know the extend of the detector and are able to tune the tolerance
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(WCLength > WCRadius ? WCLength : WCRadius);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  //Decide if adding Gd
  water = "Water";
  if (WCAddGd)
  {
	  water = "Doped Water";
	  if(!G4Material::GetMaterial("Doped Water", false)) AddDopedWater();
  }

  // Optionally switch on the checkOverlaps. Once the geometry is debugged, switch off for 
  // faster running. Checking ALL overlaps is crucial. No overlaps are allowed, otherwise
  // G4Navigator will be confused and result in wrong photon tracking and resulting yields.
  // ToDo: get these options from .mac
  checkOverlaps = true; 
  checkOverlapsPMT = false; 
  // Optionally place parts of the detector. Very useful for visualization and debugging 
  // geometry overlaps in detail.
  placeBarrelPMTs = true;
  placeCapPMTs = true;
  placeBorderPMTs = true; 

  //-----------------------------------------------------
  // Volumes
  //-----------------------------------------------------
  // The water barrel is placed in an tubs of air
  
  G4Tubs* solidWC = new G4Tubs("WC",
			       0.0*m,
			       WCRadius+2.*m, 
			       .5*WCLength+4.2*m,	//jl145 - per blueprint
			       0.*deg,
			       360.*deg);
  
  G4LogicalVolume* logicWC = 
    new G4LogicalVolume(solidWC,
			G4Material::GetMaterial("Air"),
			"WC",
			0,0,0);
 
 
   G4VisAttributes* showColor = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
   logicWC->SetVisAttributes(showColor);

   logicWC->SetVisAttributes(G4VisAttributes::Invisible); //amb79
  
  //-----------------------------------------------------
  // everything else is contained in this water tubs
  //-----------------------------------------------------
  G4Tubs* solidWCBarrel = new G4Tubs("WCBarrel",
				     0.0*m,
				     WCRadius+1.*m, // add a bit of extra space
				     .5*WCLength,  //jl145 - per blueprint
				     0.*deg,
				     360.*deg);
  
  G4LogicalVolume* logicWCBarrel = 
    new G4LogicalVolume(solidWCBarrel,
			G4Material::GetMaterial(water),
			"WCBarrel",
			0,0,0);

    G4VPhysicalVolume* physiWCBarrel = 
    new G4PVPlacement(0,
		      G4ThreeVector(0.,0.,0.),
		      logicWCBarrel,
		      "WCBarrel",
		      logicWC,
		      false,
			  0,
			  checkOverlaps); 

  G4cout << "Just a single water volume" << G4endl;

// This volume needs to made invisible to view the blacksheet and PMTs with RayTracer
  if (Vis_Choice == "RayTracer")
   {logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 

  else
   {//{if(!debugMode)
    //logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 
   }
 
  //-----------------------------------------------------
  // The PMT
  //-----------------------------------------------------

  ////////////J Felde: The PMT logical volume is now constructed separately 
  // from any specific detector geometry so that any geometry can use the same definition. 
  // K.Zbiri: The PMT volume and the PMT glass are now put in parallel. 
  // The PMT glass is the sensitive volume in this new configuration.

  // TF: Args are set to properties of the class which is somehow global (see the ConstructDetector.hh)
  //     They are set in the WCSimDetectorConfigs and are property of the PMT.
  G4LogicalVolume* logicWCPMT;
  if(nID_PMTs<=1) logicWCPMT = ConstructPMT(WCPMTName, WCIDCollectionName,"tank",nID_PMTs);
  else logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName,"tank",nID_PMTs);
  if(!logicWCPMT){
    G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
    return NULL; 
  }
	G4LogicalVolume* logicWCPMT2;
#ifdef DEBUG
	G4cout << "HYBRID2 = " << hybrid << G4endl;
#endif
	if(hybrid){
	  G4cout<<"First type of PMT LV is created. Now creating the LV for the second type of PMT"<<G4endl;
	  if(nID_PMTs2<=1) logicWCPMT2 = ConstructPMT(WCPMTName2, WCIDCollectionName2,"tankPMT2",nID_PMTs2);
	  else logicWCPMT2 = ConstructMultiPMT(WCPMTName2, WCIDCollectionName2,"tankPMT2",nID_PMTs2);
	  if(!logicWCPMT2){
		G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
		return NULL; 
	  }
	}

  //G4LogicalVolume* logicWCPMT = ConstructPMT(WCPMTName, WCIDCollectionName);
  G4String pmtname = "WCMultiPMT";


  /*These lines of code will give color and volume to the PMTs if it hasn't been set in WCSimConstructPMT.cc.
I recommend setting them in WCSimConstructPMT.cc. 
If used here, uncomment the SetVisAttributes(WClogic) line, and comment out the SetVisAttributes(G4VisAttributes::Invisible) line.*/
  
  G4VisAttributes* WClogic 
      = new G4VisAttributes(G4Colour(0.4,0.0,0.8));
     WClogic->SetForceSolid(true);
	 WClogic->SetForceAuxEdgeVisible(true);

    //logicWCPMT->SetVisAttributes(WClogic);
	 //logicWCPMT->SetVisAttributes(G4VisAttributes::Invisible);


    ///////////////   Single PMT placement
  G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
  G4VPhysicalVolume* physiWCBarrelPMT =
	new G4PVPlacement(WCPMTRotation,              // its rotation
					  G4ThreeVector(0,0,0), 
					  logicWCPMT2,                // its logical volume
					  pmtname,//"WCPMT",             // its name
					  logicWCBarrel,         // its mother volume
					  false,                     // no boolean operations
					  0,
					  checkOverlaps);   

  return logicWC;
}
