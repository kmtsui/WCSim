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

#include "G4VSolid.hh" // following 4 .hh required for CADMesh.
#include "G4Trd.hh"
#include "G4NistManager.hh"
#include "CADMesh.hh"


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


G4LogicalVolume* WCSimDetectorConstruction::ConstructCylinder()
{
    G4cout << "**** Building Cylindrical Detector ****" << G4endl;

  //-----------------------------------------------------
  // Positions
  //-----------------------------------------------------

  debugMode = false;
  
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
  if(vessel_cyl_height + vessel_radius < 1.*mm)
	innerAnnulusRadius = WCIDRadius - WCPMTExposeHeight-1.*mm;
  else
	innerAnnulusRadius = WCIDRadius - vessel_cyl_height - vessel_radius -1.*mm;
  
  //TF: need to add a Polyhedra on the other side of the outerAnnulusRadius for the OD
  outerAnnulusRadius = WCIDRadius + WCBlackSheetThickness + 1.*mm;//+ Stealstructure etc.

  // the radii are measured to the center of the surfaces
  // (tangent distance). Thus distances between the corner and the center are bigger.
  WCLength    = WCIDHeight + 2*2.3*m;	//jl145 - reflects top veto blueprint, cf. Farshid Feyzi
  WCRadius    = (WCIDDiameter/2. + WCBlackSheetThickness + 1.5*m)/cos(dPhi/2.) ; // ToDo: OD 
 
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

// This volume needs to made invisible to view the blacksheet and PMTs with RayTracer
  if (Vis_Choice == "RayTracer")
   {logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 

  else
   {//{if(!debugMode)
    //logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 
   }
	
  //----------------- example CADMesh -----------------------
  //                 L.Anthony 23/06/2022 
  //---------------------------------------------------------
  /*
  auto shape_bunny = CADMesh::TessellatedMesh::FromSTL("/path/to/WCSim//bunny.stl");
  
  // set scale
  shape_bunny->SetScale(1);
  G4ThreeVector posBunny = G4ThreeVector(0*m, 0*m, 0*m);
  
  // make new shape a solid
  G4VSolid* solid_bunny = shape_bunny->GetSolid();
  
  G4LogicalVolume* bunny_logical =                                   //logic name
	new G4LogicalVolume(solid_bunny,                                 //solid name
						G4Material::GetMaterial("StainlessSteel"),          //material
						"bunny");                                     //objects name
  // rotate if necessary
  G4RotationMatrix* bunny_rot = new G4RotationMatrix; // Rotates X and Z axes only
  bunny_rot->rotateX(270*deg);                 
  bunny_rot->rotateY(45*deg);
  
  G4VPhysicalVolume* physiBunny =
	new G4PVPlacement(bunny_rot,                       //no rotation
					  posBunny,                    //at position
					  bunny_logical,           //its logical volume
					  "bunny",                //its name
					  logicWCBarrel,                //its mother  volume
					  false,                   //no boolean operation
					  0);                       //copy number
  //                      checkOverlaps);          //overlaps checking
  
  new G4LogicalSkinSurface("BunnySurface",bunny_logical,ReflectorSkinSurface);
  */
	
  //-----------------------------------------------------
  // Form annular section of barrel to hold PMTs 
  //----------------------------------------------------

 
  G4double mainAnnulusZ[2] = {-mainAnnulusHeight/2., mainAnnulusHeight/2};
  G4double mainAnnulusRmin[2] = {innerAnnulusRadius, innerAnnulusRadius};
  G4double mainAnnulusRmax[2] = {outerAnnulusRadius, outerAnnulusRadius};

  G4Polyhedra* solidWCBarrelAnnulus = new G4Polyhedra("WCBarrelAnnulus",
                                                   0.*deg, // phi start
                                                   totalAngle, 
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   2,
                                                   mainAnnulusZ,
                                                   mainAnnulusRmin,
                                                   mainAnnulusRmax);
  
  G4LogicalVolume* logicWCBarrelAnnulus = 
    new G4LogicalVolume(solidWCBarrelAnnulus,
			G4Material::GetMaterial(water),
			"WCBarrelAnnulus",
			0,0,0);
  // G4cout << *solidWCBarrelAnnulus << G4endl; 
  G4VPhysicalVolume* physiWCBarrelAnnulus = 
    new G4PVPlacement(0,
		      G4ThreeVector(0.,0.,0.),
		      logicWCBarrelAnnulus,
		      "WCBarrelAnnulus",
		      logicWCBarrel,
		      false,
		      0,
			  checkOverlaps);
if(!debugMode)
   logicWCBarrelAnnulus->SetVisAttributes(G4VisAttributes::Invisible); //amb79

  //-----------------------------------------------------
  // Subdivide the BarrelAnnulus into rings
  //-----------------------------------------------------
  G4double RingZ[2] = {-barrelCellHeight/2.,
                        barrelCellHeight/2.};

  G4Polyhedra* solidWCBarrelRing = new G4Polyhedra("WCBarrelRing",
                                                   0.*deg,//+dPhi/2., // phi start
                                                   totalAngle, //phi end
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   2,
                                                   RingZ,
                                                   mainAnnulusRmin,
                                                   mainAnnulusRmax);

  G4LogicalVolume* logicWCBarrelRing = 
    new G4LogicalVolume(solidWCBarrelRing,
			G4Material::GetMaterial(water),
			"WCBarrelRing",
			0,0,0);

  G4VPhysicalVolume* physiWCBarrelRing = 
    new G4PVReplica("WCBarrelRing",
		    logicWCBarrelRing,
		    logicWCBarrelAnnulus,
		    kZAxis,
		    (G4int)WCBarrelNRings-2,
		    barrelCellHeight);

if(!debugMode)
  {G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
	tmpVisAtt->SetForceWireframe(true);// This line is used to give definition to the rings in OGLSX Visualizer
	logicWCBarrelRing->SetVisAttributes(tmpVisAtt);
	//If you want the rings on the Annulus to show in OGLSX, then comment out the line below.
	logicWCBarrelRing->SetVisAttributes(G4VisAttributes::Invisible);
  }
else {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
        tmpVisAtt->SetForceWireframe(true);
        logicWCBarrelRing->SetVisAttributes(tmpVisAtt);
  }

  //-----------------------------------------------------
  // Subdivisions of the BarrelRings are cells
  //------------------------------------------------------


  G4Polyhedra* solidWCBarrelCell = new G4Polyhedra("WCBarrelCell",
                                                   -dPhi/2.+0.*deg, // phi start
                                                   dPhi, //total Phi
                                                   1, //NPhi-gon
                                                   2,
                                                   RingZ,
                                                   mainAnnulusRmin,
                                                   mainAnnulusRmax); 
  //G4cout << *solidWCBarrelCell << G4endl; 
  G4LogicalVolume* logicWCBarrelCell = 
    new G4LogicalVolume(solidWCBarrelCell,
			G4Material::GetMaterial(water),
			"WCBarrelCell",
			0,0,0);

  G4VPhysicalVolume* physiWCBarrelCell = 
    new G4PVReplica("WCBarrelCell",
		    logicWCBarrelCell,
		    logicWCBarrelRing,
		    kPhi,
		    (G4int)WCBarrelRingNPhi,
		    dPhi,
                    0.); 

  if(!debugMode)
  	{G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	logicWCBarrelCell->SetVisAttributes(tmpVisAtt); 
	//If you want the columns on the Annulus to show in OGLSX, then comment out the line below.
  	logicWCBarrelCell->SetVisAttributes(G4VisAttributes::Invisible);
	}
  else {
  	G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceWireframe(true);
  	logicWCBarrelCell->SetVisAttributes(tmpVisAtt);
  }

  //-----------------------------------------------------------
  // The Blacksheet, a daughter of the cells containing PMTs,
  // and also some other volumes to make the edges light tight
  //-----------------------------------------------------------

  //-------------------------------------------------------------
  // add barrel blacksheet to the normal barrel cells 
  // ------------------------------------------------------------
  G4double annulusBlackSheetRmax[2] = {(WCIDRadius+WCBlackSheetThickness),
                                        WCIDRadius+WCBlackSheetThickness};
  G4double annulusBlackSheetRmin[2] = {(WCIDRadius),
                                        WCIDRadius};

  G4Polyhedra* solidWCBarrelCellBlackSheet = new G4Polyhedra("WCBarrelCellBlackSheet",
                                                   -dPhi/2., // phi start
                                                   dPhi, //total phi
                                                   1, //NPhi-gon
                                                   2,
                                                   RingZ,
                                                   annulusBlackSheetRmin,
                                                   annulusBlackSheetRmax);

  logicWCBarrelCellBlackSheet =
    new G4LogicalVolume(solidWCBarrelCellBlackSheet,
                        G4Material::GetMaterial("Blacksheet"),
                        "WCBarrelCellBlackSheet",
                          0,0,0);

   G4VPhysicalVolume* physiWCBarrelCellBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCBarrelCellBlackSheet,
                      "WCBarrelCellBlackSheet",
                      logicWCBarrelCell,
                      false,
                      0,
                      checkOverlaps);

   
   G4LogicalBorderSurface * WaterBSBarrelCellSurface 
	 = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
								  physiWCBarrelCell,
								  physiWCBarrelCellBlackSheet, 
								  OpWaterBSSurface);
   
   new G4LogicalSkinSurface("BSBarrelCellSkinSurface",logicWCBarrelCellBlackSheet,
							BSSkinSurface);

 // Change made here to have the if statement contain the !debugmode to be consistent
 // This code gives the Blacksheet its color. 

if (Vis_Choice == "RayTracer"){

   G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
     WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
	 WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown
      if(!debugMode)
        logicWCBarrelCellBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
      else
		logicWCBarrelCellBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);}

else {

   G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2));
      if(!debugMode)
        logicWCBarrelCellBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
      else
        logicWCBarrelCellBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}

 //-----------------------------------------------------------
 // add extra tower if necessary
 // ---------------------------------------------------------
 
  // we have to declare the logical Volumes 
  // outside of the if block to access it later on 
  G4LogicalVolume* logicWCExtraTowerCell;
  G4LogicalVolume* logicWCExtraBorderCell;
  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

    // as the angles between the corners of the main annulus 
    // and the corners extra tower are different, we need to adjust the 
    // tangent distance the surfaces of the extra tower. Otherwise
    // the corners of the main annulus and the extra tower would 
    // not match. 
    G4double extraTowerRmin[2];
    G4double extraTowerRmax[2];
    for(int i = 0; i < 2 ; i++){
      extraTowerRmin[i] = mainAnnulusRmin[i] != 0 ? mainAnnulusRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
      extraTowerRmax[i] = mainAnnulusRmax[i] != 0 ? mainAnnulusRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
    }
    G4Polyhedra* solidWCExtraTower = new G4Polyhedra("WCextraTower",
  			 totalAngle-2.*pi,//+dPhi/2., // phi start
			 2.*pi -  totalAngle // total angle.
			 -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m),
			   // we need this little Gap between the extra tower and the main annulus
			   // to avoid a shared surface. Without the gap the photons stuck
			   // at this place for mare than 25 steps and the howl simulation
			   // crashes.
	         	 1, //NPhi-gon
		         2,
			 mainAnnulusZ,
			 extraTowerRmin,
		         extraTowerRmax);

    G4LogicalVolume* logicWCExtraTower = 
      new G4LogicalVolume(solidWCExtraTower,
			  G4Material::GetMaterial(water),
			  "WCExtraTower",
			  0,0,0);
    G4VPhysicalVolume* physiWCExtraTower = 
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCExtraTower,
			"WCExtraTower",
			logicWCBarrel,
			false,
			0,
			checkOverlaps);

 

    logicWCExtraTower->SetVisAttributes(G4VisAttributes::Invisible);
  //-------------------------------------------
  // subdivide the extra tower into cells  
  //------------------------------------------

    G4Polyhedra* solidWCExtraTowerCell = new G4Polyhedra("WCExtraTowerCell",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
			   1, //NPhi-gon
			   2,
			   RingZ,
			   extraTowerRmin,
			   extraTowerRmax); 
    //G4cout << * solidWCExtraTowerCell << G4endl;
    logicWCExtraTowerCell = 
      new G4LogicalVolume(solidWCExtraTowerCell,
			  G4Material::GetMaterial(water),
			  "WCExtraTowerCell",
			  0,0,0);
    G4VPhysicalVolume* physiWCTowerCell = 
      new G4PVReplica("extraTowerCell",
		      logicWCExtraTowerCell,
		      logicWCExtraTower,
		      kZAxis,
		      (G4int)WCBarrelNRings-2,
		      barrelCellHeight);
    
    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	//logicWCExtraTowerCell->SetVisAttributes(tmpVisAtt); 
	logicWCExtraTowerCell->SetVisAttributes(G4VisAttributes::Invisible);
	//TF vis.

    //---------------------------------------------
    // add blacksheet to this cells
    //--------------------------------------------

    G4double towerBSRmin[2];
    G4double towerBSRmax[2];
    for(int i = 0; i < 2; i++){
      towerBSRmin[i] = annulusBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
      towerBSRmax[i] = annulusBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
    }
    G4Polyhedra* solidWCTowerBlackSheet = new G4Polyhedra("WCExtraTowerBlackSheet",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
		           1, //NPhi-gon
			   2,
			   RingZ,
			   towerBSRmin,
			   towerBSRmax);
    //G4cout << * solidWCTowerBlackSheet << G4endl;
    logicWCTowerBlackSheet =
      new G4LogicalVolume(solidWCTowerBlackSheet,
			  G4Material::GetMaterial("Blacksheet"),
			  "WCExtraTowerBlackSheet",
			    0,0,0);

     G4VPhysicalVolume* physiWCTowerBlackSheet =
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCTowerBlackSheet,
			"WCExtraTowerBlackSheet",
			logicWCExtraTowerCell,
			false,
			0,
			checkOverlaps);


	 G4LogicalBorderSurface * WaterBSTowerCellSurface 
	   = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
									physiWCTowerCell,
									physiWCTowerBlackSheet, 
									OpWaterBSSurface);

	 new G4LogicalSkinSurface("BSTowerCellSkinSurface",logicWCTowerBlackSheet,
							  BSSkinSurface);

// These lines add color to the blacksheet in the extratower. If using RayTracer, comment the first chunk and use the second. The Blacksheet should be green.

  if (Vis_Choice == "OGLSX"){

   G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color

	if(!debugMode)
	  {logicWCTowerBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);}
	else
	  {logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}}
  
  if (Vis_Choice == "RayTracer"){
   
    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
     WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
	 WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown

	if(!debugMode)
	  logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
	else
	  logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}

  else {

   G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color

	if(!debugMode)
	  {logicWCTowerBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);}
	else
	  {logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}}
  

}

  //jl145------------------------------------------------
  // Add top veto volume
  //-----------------------------------------------------

  G4bool WCTopVeto = (WCSimTuningParams->GetTopVeto());

  G4LogicalVolume* logicWCTopVeto;

  if(WCTopVeto){

	  G4double WCTyvekThickness = 1.0*mm; //completely made up

	  G4VSolid* solidWCTopVeto;
	  solidWCTopVeto =
			new G4Tubs(			"WCTopVeto",
								0.0*m,
								WCIDRadius + WCTyvekThickness,
								0.5*m + WCTyvekThickness,
								0.*deg,
								360.*deg);

	  logicWCTopVeto = 
			new G4LogicalVolume(solidWCTopVeto,
								G4Material::GetMaterial(water),
								"WCTopVeto",
								0,0,0);

	  G4VPhysicalVolume* physiWCTopVeto =
			new G4PVPlacement(	0,
								G4ThreeVector(0.,0.,WCIDHeight/2
													+1.0*m),
								logicWCTopVeto,
								"WCTopVeto",
								logicWCBarrel,
								false,0,
								checkOverlaps);

	  //Add the top veto Tyvek
	  //-----------------------------------------------------

	  G4VSolid* solidWCTVTyvek;
	  solidWCTVTyvek =
			new G4Tubs(			"WCTVTyvek",
								0.0*m,
								WCIDRadius,
								WCTyvekThickness/2,
								0.*deg,
								360.*deg);


	  G4LogicalVolume* logicWCTVTyvek =
			new G4LogicalVolume(solidWCTVTyvek,
								G4Material::GetMaterial("Tyvek"),
								"WCTVTyvek",
								0,0,0);

	  //Bottom
	  G4VPhysicalVolume* physiWCTVTyvekBot =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,-0.5*m
													-WCTyvekThickness/2),
								logicWCTVTyvek,
		               			"WCTVTyvekBot",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceBot =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceBot",
										physiWCTopVeto,
										physiWCTVTyvekBot,
										OpWaterTySurface);

	  //Top
	  G4VPhysicalVolume* physiWCTVTyvekTop =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,0.5*m
													+WCTyvekThickness/2),
								logicWCTVTyvek,
		               			"WCTVTyvekTop",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceTop =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceTop",
										physiWCTopVeto,
										physiWCTVTyvekTop,
										OpWaterTySurface);

	  //Side
	  G4VSolid* solidWCTVTyvekSide;
	  solidWCTVTyvekSide =
			new G4Tubs(			"WCTVTyvekSide",
								WCIDRadius,
								WCIDRadius + WCTyvekThickness,
								0.5*m + WCTyvekThickness,
								0.*deg,
								360.*deg);


	  G4LogicalVolume* logicWCTVTyvekSide =
			new G4LogicalVolume(solidWCTVTyvekSide,
								G4Material::GetMaterial("Tyvek"),
								"WCTVTyvekSide",
								0,0,0);

	  G4VPhysicalVolume* physiWCTVTyvekSide =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,0.),
								logicWCTVTyvekSide,
		               			"WCTVTyvekSide",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceSide =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceSide",
										physiWCTopVeto,
										physiWCTVTyvekSide,
										OpWaterTySurface);

  }

  //
  //
  //jl145------------------------------------------------


  //////////// M Fechner : I need to  declare the PMT  here in order to
  // place the PMTs in the truncated cells
  //-----------------------------------------------------
  // The PMT
  //-----------------------------------------------------

  ////////////J Felde: The PMT logical volume is now constructed separately 
  // from any specific detector geometry so that any geometry can use the same definition. 
  // K.Zbiri: The PMT volume and the PMT glass are now put in parallel. 
  // The PMT glass is the sensitive volume in this new configuration.

  // TF: Args are set to properties of the class which is somehow global (see the ConstructDetector.hh)
  //     They are set in the WCSimDetectorConfigs and are property of the PMT.
  G4LogicalVolume* logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName);
  if(!logicWCPMT){
    G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
    return NULL; 
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

  //jl145------------------------------------------------
  // Add top veto PMTs
  //-----------------------------------------------------

  if(WCTopVeto){

	  G4double WCTVPMTSpacing = (WCSimTuningParams->GetTVSpacing())*cm;
	  G4double WCTVEdgeLimit = WCCapEdgeLimit;
	  G4int TVNCell = WCTVEdgeLimit/WCTVPMTSpacing + 2;

	  int icopy = 0;

	  for ( int i = -TVNCell ; i <  TVNCell; i++) {
		for (int j = -TVNCell ; j <  TVNCell; j++)   {

		  G4double xoffset = i*WCTVPMTSpacing + WCTVPMTSpacing*0.5;
		  G4double yoffset = j*WCTVPMTSpacing + WCTVPMTSpacing*0.5;

		  G4ThreeVector cellpos =
		  		G4ThreeVector(	xoffset, yoffset, -0.5*m);

		  if ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCTVEdgeLimit) {

		    G4VPhysicalVolume* physiCapPMT =
		    		new G4PVPlacement(	0,						// no rotation
		    							cellpos,				// its position
		    							logicWCPMT,				// its logical volume
		    							pmtname,//"WCPMT",				// its name 
		    							logicWCTopVeto,			// its mother volume
		    							false,					// no boolean os
		    							icopy,					// every PMT need a unique id.
										checkOverlapsPMT);

		    icopy++;
		  }
		}
	  }

	  G4double WCTVEfficiency = icopy*WCPMTRadius*WCPMTRadius/((WCIDRadius)*(WCIDRadius));
	  G4cout << "Total on top veto: " << icopy << "\n";
	  G4cout << "Coverage was calculated to be: " << WCTVEfficiency << "\n";

  }

  //
  //
  //jl145------------------------------------------------


    ///////////////   Barrel PMT placement
  G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;

  //for multiPMT, for single PMT need to test defaults
  if(orientation == PERPENDICULAR)
	WCPMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
  else if(orientation == VERTICAL)
	WCPMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
  else if(orientation == HORIZONTAL)
	WCPMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
  
  


  
  



  G4double barrelCellWidth = 2.*WCIDRadius*tan(dPhi/2.);
  G4double horizontalSpacing   = barrelCellWidth/WCPMTperCellHorizontal;
  G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

  if(placeBarrelPMTs){

  for(G4double i = 0; i < WCPMTperCellHorizontal; i++){
    for(G4double j = 0; j < WCPMTperCellVertical; j++){
      G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius,
						 -barrelCellWidth/2.+(i+0.5)*horizontalSpacing,
						 -barrelCellHeight/2.+(j+0.5)*verticalSpacing);

      G4VPhysicalVolume* physiWCBarrelPMT =
	new G4PVPlacement(WCPMTRotation,              // its rotation
			  PMTPosition, 
			  logicWCPMT,                // its logical volume
			  pmtname,//"WCPMT",             // its name
			  logicWCBarrelCell,         // its mother volume
			  false,                     // no boolean operations
			  (int)(i*WCPMTperCellVertical+j),
			  checkOverlapsPMT);                       
      
   // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
     // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
		// this is still the case.
    }
  }
  //-------------------------------------------------------------
  // Add PMTs in extra Tower if necessary
  //------------------------------------------------------------


  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

    G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
    WCPMTRotation->rotateY(90.*deg);
    WCPMTRotation->rotateX((2*pi-totalAngle)/2.);//align the PMT with the Cell
                                                 
    G4double towerWidth = WCIDRadius*tan(2*pi-totalAngle);

    G4double horizontalSpacing   = towerWidth/(WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal);
    G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

    for(G4double i = 0; i < (WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal); i++){
      for(G4double j = 0; j < WCPMTperCellVertical; j++){
	G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.),
				towerWidth/2.-(i+0.5)*horizontalSpacing,
			       -barrelCellHeight/2.+(j+0.5)*verticalSpacing);
	PMTPosition.rotateZ(-(2*pi-totalAngle)/2.); // align with the symmetry 
	                                            //axes of the cell 

	G4VPhysicalVolume* physiWCBarrelPMT =
	  new G4PVPlacement(WCPMTRotation,             // its rotation
			    PMTPosition, 
			    logicWCPMT,                // its logical volume
			    "WCPMT",             // its name
			    logicWCExtraTowerCell,         // its mother volume
			    false,                     // no boolean operations
			    (int)(i*WCPMTperCellVertical+j),
			    checkOverlapsPMT);                       
	
		// logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
		// daughter volumes to the PMTs (e.g. a acryl cover) you have to check, if
		// this is still the case.
      }
    }

  }

  }//end if placeBarrelPMTs

  G4LogicalVolume* logicTopCapAssembly = ConstructCaps(-1);
  G4LogicalVolume* logicBottomCapAssembly = ConstructCaps(1);

 // These lines make the large cap volume invisible to view the caps blacksheets. Need to make invisible for 
 // RayTracer
  if (Vis_Choice == "RayTracer"){
	logicBottomCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
	logicTopCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
  } else {
	logicBottomCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
	logicTopCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
  }

  G4VPhysicalVolume* physiTopCapAssembly =
  new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(mainAnnulusHeight/2.+ capAssemblyHeight/2.)),
                  logicTopCapAssembly,
                  "TopCapAssembly",
                  logicWCBarrel,
                  false, 0,
				  checkOverlaps);

  G4VPhysicalVolume* physiBottomCapAssembly =
    new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(-mainAnnulusHeight/2.- capAssemblyHeight/2.)),
                  logicBottomCapAssembly,
                  "BottomCapAssembly",
                  logicWCBarrel,
                  false, 0,
				  checkOverlaps);

  return logicWC;
}


G4LogicalVolume* WCSimDetectorConstruction::ConstructCaps(G4int zflip)
{

  capAssemblyHeight = (WCIDHeight-mainAnnulusHeight)/2+1*mm+WCBlackSheetThickness;

  G4Tubs* solidCapAssembly = new G4Tubs("CapAssembly",
							0.0*m,
							outerAnnulusRadius/cos(dPhi/2.), 
							capAssemblyHeight/2,
							0.*deg,
							360.*deg);

  G4LogicalVolume* logicCapAssembly =
    new G4LogicalVolume(solidCapAssembly,
                        G4Material::GetMaterial(water),
                        "CapAssembly",
                        0,0,0);

  
  //----------------------------------------------------
  // extra rings for the top and bottom of the annulus
  //---------------------------------------------------
  G4double borderAnnulusZ[3] = {(-barrelCellHeight/2.-(WCIDRadius-innerAnnulusRadius))*zflip, 
								-barrelCellHeight/2.*zflip,
								barrelCellHeight/2.*zflip};
  G4double borderAnnulusRmin[3] = { WCIDRadius, innerAnnulusRadius, innerAnnulusRadius};
  G4double borderAnnulusRmax[3] = {outerAnnulusRadius, outerAnnulusRadius,outerAnnulusRadius};
  G4Polyhedra* solidWCBarrelBorderRing = new G4Polyhedra("WCBarrelBorderRing",
                                                   0.*deg, // phi start
                                                   totalAngle,
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   3,
                                                   borderAnnulusZ,
                                                   borderAnnulusRmin,
                                                   borderAnnulusRmax);
  G4LogicalVolume* logicWCBarrelBorderRing =
    new G4LogicalVolume(solidWCBarrelBorderRing,
                        G4Material::GetMaterial(water),
                        "WCBarrelRing",
                        0,0,0);
  //G4cout << *solidWCBarrelBorderRing << G4endl;

  G4VPhysicalVolume* physiWCBarrelBorderRing =
    new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(capAssemblyHeight/2.- barrelCellHeight/2.)*zflip),
                  logicWCBarrelBorderRing,
                  "WCBarrelBorderRing",
                  logicCapAssembly,
                  false, 0,
				  checkOverlaps);


                  
  if(!debugMode){ 

    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	//logicWCBarrelBorderRing->SetVisAttributes(tmpVisAtt); 
    logicWCBarrelBorderRing->SetVisAttributes(G4VisAttributes::Invisible); 
	//TF vis.
  }

  //----------------------------------------------------
  // Subdivide border rings into cells
  // --------------------------------------------------
  G4Polyhedra* solidWCBarrelBorderCell = new G4Polyhedra("WCBarrelBorderCell",
                                                   -dPhi/2., // phi start
                                                   dPhi, //total angle 
                                                   1, //NPhi-gon
                                                   3,
                                                   borderAnnulusZ,
                                                   borderAnnulusRmin,
                                                   borderAnnulusRmax);

  G4LogicalVolume* logicWCBarrelBorderCell =
    new G4LogicalVolume(solidWCBarrelBorderCell, 
                        G4Material::GetMaterial(water),
                        "WCBarrelBorderCell", 
                        0,0,0);
  //G4cout << *solidWCBarrelBorderCell << G4endl;
  G4VPhysicalVolume* physiWCBarrelBorderCell =
    new G4PVReplica("WCBarrelBorderCell",
                    logicWCBarrelBorderCell,
                    logicWCBarrelBorderRing,
                    kPhi,
                    (G4int)WCBarrelRingNPhi,
                    dPhi,
                    0.);

// These lines of code below will turn the border rings invisible. 

// used for RayTracer
 if (Vis_Choice == "RayTracer"){

  if(!debugMode){
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
        tmpVisAtt->SetForceSolid(true);
        logicWCBarrelBorderCell->SetVisAttributes(tmpVisAtt);
		logicWCBarrelBorderCell->SetVisAttributes(G4VisAttributes::Invisible);}
  else {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
        tmpVisAtt->SetForceWireframe(true);
        logicWCBarrelBorderCell->SetVisAttributes(tmpVisAtt);
		logicWCBarrelBorderCell->SetVisAttributes(G4VisAttributes::Invisible);}}

// used for OGLSX
 else {

  if(!debugMode)
        {
		  G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
		  tmpVisAtt->SetForceSolid(true);
		  logicWCBarrelBorderCell->SetVisAttributes(tmpVisAtt);
		  logicWCBarrelBorderCell->SetVisAttributes(G4VisAttributes::Invisible);
		}
  else {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
        //tmpVisAtt->SetForceWireframe(true);
		tmpVisAtt->SetForceSolid(true);
        logicWCBarrelBorderCell->SetVisAttributes(tmpVisAtt);}}


  //------------------------------------------------------------
  // add blacksheet to the border cells.
  // We can use the same logical volume as for the normal
  // barrel cells.
  // ---------------------------------------------------------


   G4VPhysicalVolume* physiWCBarrelBorderCellBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCBarrelCellBlackSheet,
                      "WCBarrelCellBlackSheet",
                      logicWCBarrelBorderCell,
                      false,
                      0,
					  checkOverlaps);

  G4LogicalBorderSurface * WaterBSBarrelBorderCellSurface
    = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
                                 physiWCBarrelBorderCell,
                                 physiWCBarrelBorderCellBlackSheet,
                                 OpWaterBSSurface);

  // we have to declare the logical Volumes 
  // outside of the if block to access it later on 
  G4LogicalVolume* logicWCExtraTowerCell;
  G4LogicalVolume* logicWCExtraBorderCell;
  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){
    //----------------------------------------------
    // also the extra tower need special cells at the 
    // top and th bottom.
    // (the top cell is created later on by reflecting the
    // bottom cell) 
    //---------------------------------------------
    G4double extraBorderRmin[3];
    G4double extraBorderRmax[3];
    for(int i = 0; i < 3; i++){
      extraBorderRmin[i] = borderAnnulusRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
      extraBorderRmax[i] = borderAnnulusRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
    } 
    G4Polyhedra* solidWCExtraBorderCell = new G4Polyhedra("WCspecialBarrelBorderCell",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //total phi
			   1, //NPhi-gon
			   3,
			   borderAnnulusZ,
			   extraBorderRmin,
			   extraBorderRmax);

    logicWCExtraBorderCell =
      new G4LogicalVolume(solidWCExtraBorderCell, 
			  G4Material::GetMaterial(water),
			  "WCspecialBarrelBorderCell", 
			  0,0,0);
    //G4cout << *solidWCExtraBorderCell << G4endl;

    G4VPhysicalVolume* physiWCExtraBorderCell =
      new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(capAssemblyHeight/2.- barrelCellHeight/2.)*zflip),
                  logicWCExtraBorderCell,
                  "WCExtraTowerBorderCell",
                  logicCapAssembly,
                  false, 0,
				  checkOverlaps);

    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	logicWCExtraBorderCell->SetVisAttributes(tmpVisAtt); 
	logicWCExtraBorderCell->SetVisAttributes(G4VisAttributes::Invisible);
	//TF vis.


    G4VPhysicalVolume* physiWCExtraBorderBlackSheet =
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCTowerBlackSheet,
			"WCExtraTowerBlackSheet",
			logicWCExtraBorderCell,
			false,
			0,
			checkOverlaps);

    G4LogicalBorderSurface * WaterBSExtraBorderCellSurface 
      = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
				   physiWCExtraBorderCell,
				   physiWCExtraBorderBlackSheet, 
				   OpWaterBSSurface);

  }
 //------------------------------------------------------------
 // add caps
 // -----------------------------------------------------------
  //crucial to match with borderAnnulusZ
  G4double capZ[4] = { (-WCBlackSheetThickness-1.*mm)*zflip,
					   (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))*zflip,
					   (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))*zflip,
					   WCBarrelPMTOffset*zflip} ;

  G4double capRmin[4] = {  0. , 0., 0., 0.} ;
  G4double capRmax[4] = {outerAnnulusRadius, outerAnnulusRadius,  WCIDRadius, innerAnnulusRadius};
  G4VSolid* solidWCCap;
  if(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal){
    solidWCCap
      = new G4Polyhedra("WCCap",
			0.*deg, // phi start
			totalAngle, //phi end
			(int)WCBarrelRingNPhi, //NPhi-gon
			4, // 2 z-planes
			capZ, //position of the Z planes
			capRmin, // min radius at the z planes
			capRmax// max radius at the Z planes
			);
  } else {
  // if there is an extra tower, the cap volume is a union of
  // to polyhedra. We have to unite both parts, because there are 
  // PMTs that are on the border between both parts.
   G4Polyhedra* mainPart 
      = new G4Polyhedra("WCCapMainPart",
			0.*deg, // phi start
			totalAngle, //phi end
			(int)WCBarrelRingNPhi, //NPhi-gon
			4, // 2 z-planes
			capZ, //position of the Z planes
			capRmin, // min radius at the z planes
			capRmax// max radius at the Z planes
			);
    G4double extraCapRmin[4]; 
    G4double extraCapRmax[4]; 
    for(int i = 0; i < 4 ; i++){
      extraCapRmin[i] = capRmin[i] != 0. ?  capRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
      extraCapRmax[i] = capRmax[i] != 0. ? capRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
    }
    G4Polyhedra* extraSlice 
      = new G4Polyhedra("WCCapExtraSlice",
			totalAngle-2.*pi, // phi start
			2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //total phi 
			// fortunately there are no PMTs an the gap!
			1, //NPhi-gon
			4, //  z-planes
			capZ, //position of the Z planes
			extraCapRmin, // min radius at the z planes
			extraCapRmax// max radius at the Z planes
			);
     solidWCCap =
        new G4UnionSolid("WCCap", mainPart, extraSlice);

     //G4cout << *solidWCCap << G4endl;
   
  }
  // G4cout << *solidWCCap << G4endl;
  G4LogicalVolume* logicWCCap = 
    new G4LogicalVolume(solidWCCap,
			G4Material::GetMaterial(water),
			"WCCapPolygon",
			0,0,0);

  G4VPhysicalVolume* physiWCCap =
    new G4PVPlacement(0,                           // no rotation
		      G4ThreeVector(0.,0.,(-capAssemblyHeight/2.+1*mm+WCBlackSheetThickness)*zflip),     // its position
		      logicWCCap,          // its logical volume
		      "WCCap",             // its name
		      logicCapAssembly,                  // its mother volume
		      false,                       // no boolean operations
		      0,                          // Copy #
			  checkOverlaps);


// used for RayTracer
 if (Vis_Choice == "RayTracer"){
  if(!debugMode){  
    G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(1,0.5,0.5));
	tmpVisAtt2->SetForceSolid(true);
	logicWCCap->SetVisAttributes(tmpVisAtt2);
    logicWCCap->SetVisAttributes(G4VisAttributes::Invisible);

  } else{
	
	G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(0.6,0.5,0.5));
	tmpVisAtt2->SetForceSolid(true);
    logicWCCap->SetVisAttributes(tmpVisAtt2);}}

// used for OGLSX
 else{
  if(!debugMode){  

    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	//logicWCCap->SetVisAttributes(tmpVisAtt); 
    logicWCCap->SetVisAttributes(G4VisAttributes::Invisible);
	//TF vis.

  } else
	{G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(.6,0.5,0.5));
    tmpVisAtt2->SetForceWireframe(true);
    logicWCCap->SetVisAttributes(tmpVisAtt2);}}

  //---------------------------------------------------------------------
  // add cap blacksheet
  // -------------------------------------------------------------------
  
 G4double capBlackSheetZ[4] = {-WCBlackSheetThickness*zflip, 0., 0., (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius)) *zflip};
  G4double capBlackSheetRmin[4] = {0., 0., WCIDRadius, WCIDRadius};
  G4double capBlackSheetRmax[4] = {WCIDRadius+WCBlackSheetThickness, 
                                   WCIDRadius+WCBlackSheetThickness,
								   WCIDRadius+WCBlackSheetThickness,
								   WCIDRadius+WCBlackSheetThickness};
  G4VSolid* solidWCCapBlackSheet;
  if(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal){
    solidWCCapBlackSheet
      = new G4Polyhedra("WCCapBlackSheet",
			0.*deg, // phi start
			totalAngle, //total phi
			WCBarrelRingNPhi, //NPhi-gon
			4, //  z-planes
			capBlackSheetZ, //position of the Z planes
			capBlackSheetRmin, // min radius at the z planes
			capBlackSheetRmax// max radius at the Z planes
			);
    // G4cout << *solidWCCapBlackSheet << G4endl;
  } else { 
    // same as for the cap volume
     G4Polyhedra* mainPart
      = new G4Polyhedra("WCCapBlackSheetMainPart",
			0.*deg, // phi start
			totalAngle, //phi end
			WCBarrelRingNPhi, //NPhi-gon
			4, //  z-planes
			capBlackSheetZ, //position of the Z planes
			capBlackSheetRmin, // min radius at the z planes
			capBlackSheetRmax// max radius at the Z planes
			);
     G4double extraBSRmin[4];
     G4double extraBSRmax[4];
     for(int i = 0; i < 4 ; i++){
       extraBSRmin[i] = capBlackSheetRmin[i] != 0. ? capBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
       extraBSRmax[i] = capBlackSheetRmax[i] != 0. ? capBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
     }
     G4Polyhedra* extraSlice
      = new G4Polyhedra("WCCapBlackSheetextraSlice",
			totalAngle-2.*pi, // phi start
			2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //
			WCBarrelRingNPhi, //NPhi-gon
			4, //  z-planes
			capBlackSheetZ, //position of the Z planes
			extraBSRmin, // min radius at the z planes
			extraBSRmax// max radius at the Z planes
			);
    
     solidWCCapBlackSheet =
        new G4UnionSolid("WCCapBlackSheet", mainPart, extraSlice);
  }
  G4LogicalVolume* logicWCCapBlackSheet =
    new G4LogicalVolume(solidWCCapBlackSheet,
			G4Material::GetMaterial("Blacksheet"),
			"WCCapBlackSheet",
			0,0,0);
  G4VPhysicalVolume* physiWCCapBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCCapBlackSheet,
                      "WCCapBlackSheet",
                      logicWCCap,
                      false,
                      0,
					  checkOverlaps);

  G4LogicalBorderSurface * WaterBSBottomCapSurface 
	= new G4LogicalBorderSurface("WaterBSCapPolySurface",
								 physiWCCap,physiWCCapBlackSheet,
								 OpWaterBSSurface);
  
  new G4LogicalSkinSurface("BSBottomCapSkinSurface",logicWCCapBlackSheet,
						   BSSkinSurface);
  
  G4VisAttributes* WCCapBlackSheetVisAtt 
      = new G4VisAttributes(G4Colour(0.9,0.2,0.2));
    
// used for OGLSX
 if (Vis_Choice == "OGLSX"){

   G4VisAttributes* WCCapBlackSheetVisAtt 
   = new G4VisAttributes(G4Colour(0.9,0.2,0.2));
 
	if(!debugMode)
        logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    else
        logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);}

// used for RayTracer (makes the caps blacksheet yellow)
 if (Vis_Choice == "RayTracer"){

   G4VisAttributes* WCCapBlackSheetVisAtt 
   = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

	if(!debugMode)
       //logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible); //Use this line if you want to make the blacksheet on the caps invisible to view through
	   logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);
    else
        logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);}

  //---------------------------------------------------------
  // Add top and bottom PMTs
  // -----------------------------------------------------
  
 G4LogicalVolume* logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName);
 //G4LogicalVolume* logicWCPMT = ConstructPMT(WCPMTName, WCIDCollectionName);
 G4String pmtname = "WCMultiPMT";
 
 // If using RayTracer and want to view the detector without caps, comment out the top and bottom PMT's
  G4double xoffset;
  G4double yoffset;
  G4int    icopy = 0;

  G4RotationMatrix* WCCapPMTRotation = new G4RotationMatrix;
  //if mPMT: perp to wall
  if(orientation == PERPENDICULAR){
	if(zflip==-1){
	  WCCapPMTRotation->rotateY(180.*deg); 
	}
  }
  else if (orientation == VERTICAL)
	WCCapPMTRotation->rotateY(90.*deg);
  else if (orientation == HORIZONTAL)
	WCCapPMTRotation->rotateX(90.*deg);

  // loop over the cap
  if(placeCapPMTs){
  G4int CapNCell = WCCapEdgeLimit/WCCapPMTSpacing + 2;
  for ( int i = -CapNCell ; i <  CapNCell; i++) {
    for (int j = -CapNCell ; j <  CapNCell; j++)   {

      // Jun. 04, 2020 by M.Shinoki
      // For IWCD (NuPRISM_mPMT Geometry)
      xoffset = i*WCCapPMTSpacing + WCCapPMTSpacing*0.5;
      yoffset = j*WCCapPMTSpacing + WCCapPMTSpacing*0.5;
      // For WCTE (NuPRISMBeamTest_mPMT Geometry)
      if (isNuPrismBeamTest || isNuPrismBeamTest_16cShort){
      	xoffset = i*WCCapPMTSpacing;
      	yoffset = j*WCCapPMTSpacing;
      }
      G4ThreeVector cellpos = G4ThreeVector(xoffset, yoffset, 0);     
      //      G4double WCBarrelEffRadius = WCIDDiameter/2. - WCCapPMTSpacing;
      //      double comp = xoffset*xoffset + yoffset*yoffset 
      //	- 2.0 * WCBarrelEffRadius * sqrt(xoffset*xoffset+yoffset*yoffset)
      //	+ WCBarrelEffRadius*WCBarrelEffRadius;
      //      if ( (comp > WCPMTRadius*WCPMTRadius) && ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCCapEdgeLimit) ) {
	  if ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCCapEdgeLimit) 

		// for debugging boundary cases: 
		// &&  ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) > (WCCapEdgeLimit-100)) ) 
		{
		

	G4VPhysicalVolume* physiCapPMT =
	  new G4PVPlacement(WCCapPMTRotation,
			    cellpos,                   // its position
			    logicWCPMT,                // its logical volume
			    pmtname, // its name 
			    logicWCCap,         // its mother volume
			    false,                 // no boolean os
				icopy,               // every PMT need a unique id.
				checkOverlapsPMT);
	

 // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
    // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
	// this is still the case.

	icopy++;
      }
    }
  }


  G4cout << "total on cap: " << icopy << "\n";
  G4cout << "Coverage was calculated to be: " << (icopy*WCPMTRadius*WCPMTRadius/(WCIDRadius*WCIDRadius)) << "\n";
  }//end if placeCapPMTs

    ///////////////   Barrel PMT placement
  G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
  if(orientation == PERPENDICULAR)
	WCPMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
  else if(orientation == VERTICAL)
	WCPMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
  else if(orientation == HORIZONTAL)
	WCPMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
  

  G4double barrelCellWidth = 2.*WCIDRadius*tan(dPhi/2.);
  G4double horizontalSpacing   = barrelCellWidth/WCPMTperCellHorizontal;
  G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

  if(placeBorderPMTs){

  for(G4double i = 0; i < WCPMTperCellHorizontal; i++){
    for(G4double j = 0; j < WCPMTperCellVertical; j++){
      G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius,
						 -barrelCellWidth/2.+(i+0.5)*horizontalSpacing,
						 (-barrelCellHeight/2.+(j+0.5)*verticalSpacing)*zflip);

     G4VPhysicalVolume* physiWCBarrelBorderPMT =
	new G4PVPlacement(WCPMTRotation,                      // its rotation
			  PMTPosition,
			  logicWCPMT,                // its logical volume
			  pmtname,             // its name
			  logicWCBarrelBorderCell,         // its mother volume
			  false,                     // no boolean operations
			  (int)(i*WCPMTperCellVertical+j),
			  checkOverlapsPMT); 

   // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
     // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
		// this is still the case.
	}
  }
  //-------------------------------------------------------------
  // Add PMTs in extra Tower if necessary
  //------------------------------------------------------------


  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

    G4RotationMatrix* WCPMTRotation = new G4RotationMatrix;
    WCPMTRotation->rotateY(90.*deg);
    WCPMTRotation->rotateX((2*pi-totalAngle)/2.);//align the PMT with the Cell
                                                 
    G4double towerWidth = WCIDRadius*tan(2*pi-totalAngle);

    G4double horizontalSpacing   = towerWidth/(WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal);
    G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

    for(G4double i = 0; i < (WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal); i++){
      for(G4double j = 0; j < WCPMTperCellVertical; j++){
	G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.),
				towerWidth/2.-(i+0.5)*horizontalSpacing,
			       (-barrelCellHeight/2.+(j+0.5)*verticalSpacing)*zflip);
	PMTPosition.rotateZ(-(2*pi-totalAngle)/2.); // align with the symmetry 
	                                            //axes of the cell 
	
	G4VPhysicalVolume* physiWCBarrelBorderPMT =
	  new G4PVPlacement(WCPMTRotation,                          // its rotation
			    PMTPosition,
			    logicWCPMT,                // its logical volume
			    "WCPMT",             // its name
			    logicWCExtraBorderCell,         // its mother volume
			    false,                     // no boolean operations
				(int)(i*WCPMTperCellVertical+j),
			    checkOverlapsPMT);

		// logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
		// daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
		// this is still the case.
      }
    }
  
  }
  }//end if placeBorderPMTs
return logicCapAssembly;

}


G4LogicalVolume* WCSimDetectorConstruction::ConstructCylinderNoReplica()
{
    G4cout << "**** Building Cylindrical Detector without using replica ****" << G4endl;

  //-----------------------------------------------------
  // Positions
  //-----------------------------------------------------

  debugMode = false;
  
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
  if(vessel_cyl_height + vessel_radius < 1.*mm)
	innerAnnulusRadius = WCIDRadius - WCPMTExposeHeight-1.*mm;
  else
	innerAnnulusRadius = WCIDRadius - vessel_cyl_height - vessel_radius -1.*mm;
  
  //TF: need to add a Polyhedra on the other side of the outerAnnulusRadius for the OD
  outerAnnulusRadius = WCIDRadius + WCBlackSheetThickness + 1.*mm;//+ Stealstructure etc.

  // the radii are measured to the center of the surfaces
  // (tangent distance). Thus distances between the corner and the center are bigger.
  WCLength    = WCIDHeight + 2*2.3*m;	//jl145 - reflects top veto blueprint, cf. Farshid Feyzi
  WCRadius    = (WCIDDiameter/2. + WCBlackSheetThickness + 1.5*m)/cos(dPhi/2.) ; // ToDo: OD 
 
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

// This volume needs to made invisible to view the blacksheet and PMTs with RayTracer
  if (Vis_Choice == "RayTracer")
   {logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 

  else
   {//{if(!debugMode)
    //logicWCBarrel->SetVisAttributes(G4VisAttributes::Invisible);} 
   }
	
  //----------------- example CADMesh -----------------------
  //                 L.Anthony 23/06/2022 
  //---------------------------------------------------------
  /*
  auto shape_bunny = CADMesh::TessellatedMesh::FromSTL("/path/to/WCSim//bunny.stl");
  
  // set scale
  shape_bunny->SetScale(1);
  G4ThreeVector posBunny = G4ThreeVector(0*m, 0*m, 0*m);
  
  // make new shape a solid
  G4VSolid* solid_bunny = shape_bunny->GetSolid();
  
  G4LogicalVolume* bunny_logical =                                   //logic name
	new G4LogicalVolume(solid_bunny,                                 //solid name
						G4Material::GetMaterial("StainlessSteel"),          //material
						"bunny");                                     //objects name
  // rotate if necessary
  G4RotationMatrix* bunny_rot = new G4RotationMatrix; // Rotates X and Z axes only
  bunny_rot->rotateX(270*deg);                 
  bunny_rot->rotateY(45*deg);
  
  G4VPhysicalVolume* physiBunny =
	new G4PVPlacement(bunny_rot,                       //no rotation
					  posBunny,                    //at position
					  bunny_logical,           //its logical volume
					  "bunny",                //its name
					  logicWCBarrel,                //its mother  volume
					  false,                   //no boolean operation
					  0);                       //copy number
  //                      checkOverlaps);          //overlaps checking
  
  new G4LogicalSkinSurface("BunnySurface",bunny_logical,ReflectorSkinSurface);
  */
	
  //-----------------------------------------------------
  // Form annular section of barrel to hold PMTs 
  //----------------------------------------------------

  const int nLayers = 100; // just a large number to avoid compile warning
  // Divide into WCBarrelNRings-2 layers
  G4double mainAnnulusZ[nLayers];
  G4double mainAnnulusRmin[nLayers];
  G4double mainAnnulusRmax[nLayers];
  for (int i=0;i<(G4int)WCBarrelNRings-1;i++)
  {
    mainAnnulusZ[i] =  -mainAnnulusHeight/2 + mainAnnulusHeight*i/((G4int)WCBarrelNRings-2.);
    G4double dr = GetRadiusChange(mainAnnulusZ[i]); // radius change at end-point
    mainAnnulusRmin[i] = innerAnnulusRadius+dr;
    mainAnnulusRmax[i] = outerAnnulusRadius+dr;
  }

  G4Polyhedra* solidWCBarrelAnnulus = new G4Polyhedra("WCBarrelAnnulus",
                                                   0.*deg, // phi start
                                                   totalAngle, 
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   (G4int)WCBarrelNRings-1,
                                                   mainAnnulusZ,
                                                   mainAnnulusRmin,
                                                   mainAnnulusRmax);
  
  G4LogicalVolume* logicWCBarrelAnnulus = 
    new G4LogicalVolume(solidWCBarrelAnnulus,
			G4Material::GetMaterial(water),
			"WCBarrelAnnulus",
			0,0,0);
  // G4cout << *solidWCBarrelAnnulus << G4endl; 
  G4VPhysicalVolume* physiWCBarrelAnnulus = 
    new G4PVPlacement(0,
		      G4ThreeVector(0.,0.,0.),
		      logicWCBarrelAnnulus,
		      "WCBarrelAnnulus",
		      logicWCBarrel,
		      false,
		      0,
			  checkOverlaps);
  if(!debugMode)
    logicWCBarrelAnnulus->SetVisAttributes(G4VisAttributes::Invisible); //amb79

  if(!debugMode)
  {
    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
    tmpVisAtt->SetForceWireframe(true);// This line is used to give definition to the rings in OGLSX Visualizer
    logicWCBarrelAnnulus->SetVisAttributes(tmpVisAtt);
    //If you want the rings on the Annulus to show in OGLSX, then comment out the line below.
    logicWCBarrelAnnulus->SetVisAttributes(G4VisAttributes::Invisible);
  }
  else {
    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(0,0.5,1.));
    tmpVisAtt->SetForceWireframe(true);
    logicWCBarrelAnnulus->SetVisAttributes(tmpVisAtt);
  }

  //-----------------------------------------------------------
  // The Blacksheet, a daughter of the cells containing PMTs,
  // and also some other volumes to make the edges light tight
  //-----------------------------------------------------------

  //-------------------------------------------------------------
  // add barrel blacksheet to the normal barrel cells 
  // ------------------------------------------------------------
  G4double annulusBlackSheetRmax[nLayers];
  G4double annulusBlackSheetRmin[nLayers];
  for (int i=0;i<(G4int)WCBarrelNRings-1;i++)
  {
    G4double dr = GetRadiusChange(mainAnnulusZ[i]);
    annulusBlackSheetRmin[i] = WCIDRadius+dr;
    annulusBlackSheetRmax[i] = WCIDRadius+WCBlackSheetThickness+dr;
  }

  G4Polyhedra* solidWCBarrelBlackSheet = new G4Polyhedra("WCBarrelAnnulusBlackSheet",
                                                   0, // phi start
                                                   totalAngle, //total phi
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   (G4int)WCBarrelNRings-1,
                                                   mainAnnulusZ,
                                                   annulusBlackSheetRmin,
                                                   annulusBlackSheetRmax);

  logicWCBarrelAnnulusBlackSheet =
    new G4LogicalVolume(solidWCBarrelBlackSheet,
                        G4Material::GetMaterial("Blacksheet"),
                        "WCBarrelAnnulusBlackSheet",
                          0,0,0);

   G4VPhysicalVolume* physiWCBarrelAnnulusBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCBarrelAnnulusBlackSheet,
                      "WCBarrelAnnulusBlackSheet",
                      logicWCBarrelAnnulus,
                      false,
                      0,
                      checkOverlaps);

   
   G4LogicalBorderSurface * WaterBSBarrelSurface 
	 = new G4LogicalBorderSurface("WaterBSBarrelAnnulusSurface",
								  physiWCBarrelAnnulus,
								  physiWCBarrelAnnulusBlackSheet, 
								  OpWaterBSSurface);
   
   new G4LogicalSkinSurface("BSBarrelAnnulusSkinSurface",logicWCBarrelAnnulusBlackSheet,
							BSSkinSurface);

  // Change made here to have the if statement contain the !debugmode to be consistent
  // This code gives the Blacksheet its color. 

  if (Vis_Choice == "RayTracer"){

    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
    WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
	  WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown
    if(!debugMode)
    {
      logicWCBarrelAnnulusBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }
    else
    {
      logicWCBarrelAnnulusBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    }
  }

  else {

    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2));
    if(!debugMode)
    {
      logicWCBarrelAnnulusBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    }
    else
    {
      logicWCBarrelAnnulusBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }
  }

  //-----------------------------------------------------------
  // add extra tower if necessary
  // ---------------------------------------------------------
  G4LogicalVolume* logicWCExtraTower;
  G4double towerBSRmin[nLayers];
  G4double towerBSRmax[nLayers];
  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

    // as the angles between the corners of the main annulus 
    // and the corners extra tower are different, we need to adjust the 
    // tangent distance the surfaces of the extra tower. Otherwise
    // the corners of the main annulus and the extra tower would 
    // not match. 
    G4double extraTowerRmin[nLayers];
    G4double extraTowerRmax[nLayers];
    for(int i = 0; i < (G4int)WCBarrelNRings-1 ; i++){
      extraTowerRmin[i] = mainAnnulusRmin[i] != 0 ? mainAnnulusRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
      extraTowerRmax[i] = mainAnnulusRmax[i] != 0 ? mainAnnulusRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
    }
    G4Polyhedra* solidWCExtraTower = new G4Polyhedra("WCextraTower",
  			totalAngle-2.*pi,//+dPhi/2., // phi start
			  2.*pi -  totalAngle // total angle.
			  -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m),
        // we need this little Gap between the extra tower and the main annulus
        // to avoid a shared surface. Without the gap the photons stuck
        // at this place for mare than 25 steps and the howl simulation
        // crashes.
	      1, //NPhi-gon
		    (G4int)WCBarrelNRings-1,
			  mainAnnulusZ,
			  extraTowerRmin,
		    extraTowerRmax);

    logicWCExtraTower = 
      new G4LogicalVolume(solidWCExtraTower,
			  G4Material::GetMaterial(water),
			  "WCExtraTower",
			  0,0,0);
    G4VPhysicalVolume* physiWCExtraTower = 
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCExtraTower,
			"WCExtraTower",
			logicWCBarrel,
			false,
			0,
			checkOverlaps);

    logicWCExtraTower->SetVisAttributes(G4VisAttributes::Invisible);
    
    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	//logicWCExtraTower->SetVisAttributes(tmpVisAtt); 
	  //TF vis.

    //---------------------------------------------
    // add blacksheet to this cells
    //--------------------------------------------

    for(int i = 0; i < (G4int)WCBarrelNRings-1; i++){
      towerBSRmin[i] = annulusBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
      towerBSRmax[i] = annulusBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
    }
    G4Polyhedra* solidWCTowerBlackSheet = new G4Polyhedra("WCExtraTowerBlackSheet",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
		     1, //NPhi-gon
			   (G4int)WCBarrelNRings-1,
			   mainAnnulusZ,
			   towerBSRmin,
			   towerBSRmax);
    //G4cout << * solidWCTowerBlackSheet << G4endl;
    logicWCTowerBlackSheet =
      new G4LogicalVolume(solidWCTowerBlackSheet,
			  G4Material::GetMaterial("Blacksheet"),
			  "WCExtraTowerBlackSheet",
			    0,0,0);

    G4VPhysicalVolume* physiWCTowerBlackSheet =
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCTowerBlackSheet,
			"WCExtraTowerBlackSheet",
			logicWCExtraTower,
			false,
			0,
			checkOverlaps);


    G4LogicalBorderSurface * WaterBSTowerCellSurface 
      = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
                    physiWCExtraTower,
                    physiWCTowerBlackSheet, 
                    OpWaterBSSurface);

    new G4LogicalSkinSurface("BSTowerSkinSurface",logicWCTowerBlackSheet,
                  BSSkinSurface);

    // These lines add color to the blacksheet in the extratower. If using RayTracer, comment the first chunk and use the second. The Blacksheet should be green.
  
    if (Vis_Choice == "RayTracer"){
    
      G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
        = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
      WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
      WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown

      if(!debugMode)
        logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
      else
        logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }

    else {

      G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
        = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color

      if(!debugMode)
        {logicWCTowerBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);}
      else
        {logicWCTowerBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}
    }
  

  }

  //jl145------------------------------------------------
  // Add top veto volume
  //-----------------------------------------------------

  G4bool WCTopVeto = (WCSimTuningParams->GetTopVeto());

  G4LogicalVolume* logicWCTopVeto;

  if(WCTopVeto){

	  G4double WCTyvekThickness = 1.0*mm; //completely made up

	  G4VSolid* solidWCTopVeto;
	  solidWCTopVeto =
			new G4Tubs(			"WCTopVeto",
								0.0*m,
								WCIDRadius + WCTyvekThickness + GetRadiusChange(WCIDHeight/2),
								0.5*m + WCTyvekThickness,
								0.*deg,
								360.*deg);

	  logicWCTopVeto = 
			new G4LogicalVolume(solidWCTopVeto,
								G4Material::GetMaterial(water),
								"WCTopVeto",
								0,0,0);

	  G4VPhysicalVolume* physiWCTopVeto =
			new G4PVPlacement(	0,
								G4ThreeVector(0.,0.,WCIDHeight/2
													+1.0*m),
								logicWCTopVeto,
								"WCTopVeto",
								logicWCBarrel,
								false,0,
								checkOverlaps);

	  //Add the top veto Tyvek
	  //-----------------------------------------------------

	  G4VSolid* solidWCTVTyvek;
	  solidWCTVTyvek =
			new G4Tubs(			"WCTVTyvek",
								0.0*m,
								WCIDRadius + GetRadiusChange(WCIDHeight/2),
								WCTyvekThickness/2,
								0.*deg,
								360.*deg);


	  G4LogicalVolume* logicWCTVTyvek =
			new G4LogicalVolume(solidWCTVTyvek,
								G4Material::GetMaterial("Tyvek"),
								"WCTVTyvek",
								0,0,0);

	  //Bottom
	  G4VPhysicalVolume* physiWCTVTyvekBot =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,-0.5*m
													-WCTyvekThickness/2),
								logicWCTVTyvek,
		               			"WCTVTyvekBot",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceBot =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceBot",
										physiWCTopVeto,
										physiWCTVTyvekBot,
										OpWaterTySurface);

	  //Top
	  G4VPhysicalVolume* physiWCTVTyvekTop =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,0.5*m
													+WCTyvekThickness/2),
								logicWCTVTyvek,
		               			"WCTVTyvekTop",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceTop =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceTop",
										physiWCTopVeto,
										physiWCTVTyvekTop,
										OpWaterTySurface);

	  //Side
	  G4VSolid* solidWCTVTyvekSide;
	  solidWCTVTyvekSide =
			new G4Tubs(			"WCTVTyvekSide",
								WCIDRadius+ GetRadiusChange(WCIDHeight/2),
								WCIDRadius + WCTyvekThickness+ GetRadiusChange(WCIDHeight/2),
								0.5*m + WCTyvekThickness,
								0.*deg,
								360.*deg);


	  G4LogicalVolume* logicWCTVTyvekSide =
			new G4LogicalVolume(solidWCTVTyvekSide,
								G4Material::GetMaterial("Tyvek"),
								"WCTVTyvekSide",
								0,0,0);

	  G4VPhysicalVolume* physiWCTVTyvekSide =
			new G4PVPlacement(	0,
		                  		G4ThreeVector(0.,0.,0.),
								logicWCTVTyvekSide,
		               			"WCTVTyvekSide",
		          				logicWCTopVeto,
				 				false,0,
								checkOverlaps);

	  G4LogicalBorderSurface * WaterTyTVSurfaceSide =
			new G4LogicalBorderSurface(	"WaterTyTVSurfaceSide",
										physiWCTopVeto,
										physiWCTVTyvekSide,
										OpWaterTySurface);

  }

  //
  //
  //jl145------------------------------------------------


  //////////// M Fechner : I need to  declare the PMT  here in order to
  // place the PMTs in the truncated cells
  //-----------------------------------------------------
  // The PMT
  //-----------------------------------------------------

  ////////////J Felde: The PMT logical volume is now constructed separately 
  // from any specific detector geometry so that any geometry can use the same definition. 
  // K.Zbiri: The PMT volume and the PMT glass are now put in parallel. 
  // The PMT glass is the sensitive volume in this new configuration.

  // TF: Args are set to properties of the class which is somehow global (see the ConstructDetector.hh)
  //     They are set in the WCSimDetectorConfigs and are property of the PMT.
  G4LogicalVolume* logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName);
  if(!logicWCPMT){
    G4cerr << "Overlapping PMTs in multiPMT" << G4endl;
    return NULL; 
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

  //jl145------------------------------------------------
  // Add top veto PMTs
  //-----------------------------------------------------

  if(WCTopVeto){

	  G4double WCTVPMTSpacing = (WCSimTuningParams->GetTVSpacing())*cm;
	  G4double WCTVEdgeLimit = WCCapEdgeLimit;
	  G4int TVNCell = WCTVEdgeLimit/WCTVPMTSpacing + 2;

	  int icopy = 0;

	  for ( int i = -TVNCell ; i <  TVNCell; i++) {
		for (int j = -TVNCell ; j <  TVNCell; j++)   {

		  G4double xoffset = i*WCTVPMTSpacing + WCTVPMTSpacing*0.5;
		  G4double yoffset = j*WCTVPMTSpacing + WCTVPMTSpacing*0.5;

		  G4ThreeVector cellpos =
		  		G4ThreeVector(	xoffset, yoffset, -0.5*m);

		  if ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCTVEdgeLimit) {

		    G4VPhysicalVolume* physiCapPMT =
		    		new G4PVPlacement(	0,						// no rotation
		    							cellpos,				// its position
		    							logicWCPMT,				// its logical volume
		    							pmtname,//"WCPMT",				// its name 
		    							logicWCTopVeto,			// its mother volume
		    							false,					// no boolean os
		    							icopy,					// every PMT need a unique id.
										checkOverlapsPMT);

		    icopy++;
		  }
		}
	  }

	  G4double WCTVEfficiency = icopy*WCPMTRadius*WCPMTRadius/((WCIDRadius)*(WCIDRadius));
	  G4cout << "Total on top veto: " << icopy << "\n";
	  G4cout << "Coverage was calculated to be: " << WCTVEfficiency << "\n";

  }

  //
  //
  //jl145------------------------------------------------


  ///////////////   Barrel PMT placement

  if(placeBarrelPMTs){

    G4int copyNo = 0;
    for (int iz=0;iz<(G4int)WCBarrelNRings-2;iz++)
    {
      G4double z_offset = -mainAnnulusHeight/2.+barrelCellHeight*(iz+0.5);

      // slight different cell width at different z
      G4double barrelCellWidth = (annulusBlackSheetRmin[iz+1]+annulusBlackSheetRmin[iz])*tan(dPhi/2.);
      G4double horizontalSpacing   = barrelCellWidth/WCPMTperCellHorizontal;
      G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;
      
      for (int iphi=0;iphi<(G4int)WCBarrelRingNPhi;iphi++)
      {
        G4double phi_offset = (iphi+0.5)*dPhi;

        G4RotationMatrix* PMTRotation = new G4RotationMatrix;
        if(orientation == PERPENDICULAR)
          PMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
        else if(orientation == VERTICAL)
          PMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
        else if(orientation == HORIZONTAL)
          PMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
        PMTRotation->rotateX(-phi_offset);//align the PMT with the Cell
        // inclined barrel wall
        G4double dth = atan((annulusBlackSheetRmin[iz+1]-annulusBlackSheetRmin[iz])/(mainAnnulusZ[iz+1]-mainAnnulusZ[iz]));
        if(orientation == PERPENDICULAR)
          PMTRotation->rotateY(-dth); 
        else if(orientation == VERTICAL)
          PMTRotation->rotateY(-dth); 

        for(G4double i = 0; i < WCPMTperCellHorizontal; i++){
          for(G4double j = 0; j < WCPMTperCellVertical; j++){
            G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius,
                  -barrelCellWidth/2.+(i+0.5)*horizontalSpacing + G4RandGauss::shoot(0,pmtPosVar),
                  -barrelCellHeight/2.+(j+0.5)*verticalSpacing + z_offset + G4RandGauss::shoot(0,pmtPosVar));

            // ID radius is changed
            G4double newR = annulusBlackSheetRmin[iz]+(annulusBlackSheetRmin[iz+1]-annulusBlackSheetRmin[iz])*(PMTPosition.z()-mainAnnulusZ[iz])/(mainAnnulusZ[iz+1]-mainAnnulusZ[iz]);
            PMTPosition.setX(newR);

            PMTPosition.rotateZ(phi_offset);  // align with the symmetry 
	                                            //axes of the cell 

            G4VPhysicalVolume* physiWCBarrelPMT =
              new G4PVPlacement(PMTRotation,              // its rotation
                    PMTPosition, 
                    logicWCPMT,                // its logical volume
                    pmtname,//"WCPMT",             // its name
                    logicWCBarrelAnnulus,         // its mother volume
                    false,                     // no boolean operations
                    copyNo++,
                    checkOverlapsPMT);                       
          }
        }
      }
    }

    //-------------------------------------------------------------
    // Add PMTs in extra Tower if necessary
    //------------------------------------------------------------


    if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

      copyNo = 0;
      for (int iz=0;iz<(G4int)WCBarrelNRings-2;iz++)
      {
        G4double z_offset = -mainAnnulusHeight/2.+barrelCellHeight*(iz+0.5);

        G4double towerWidth = (annulusBlackSheetRmin[iz+1]+annulusBlackSheetRmin[iz])/2.*tan(2*pi-totalAngle);
        G4double horizontalSpacing   = towerWidth/(WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal);
        G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

        G4RotationMatrix* PMTRotation = new G4RotationMatrix;
        if(orientation == PERPENDICULAR)
          PMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
        else if(orientation == VERTICAL)
          PMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
        else if(orientation == HORIZONTAL)
          PMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
        PMTRotation->rotateX((2*pi-totalAngle)/2.);//align the PMT with the Cell
        G4double dth = atan((towerBSRmin[iz+1]-towerBSRmin[iz])/(mainAnnulusZ[iz+1]-mainAnnulusZ[iz]));
        if(orientation == PERPENDICULAR)
          PMTRotation->rotateY(-dth); 
        else if(orientation == VERTICAL)
          PMTRotation->rotateY(-dth); 

        for(G4double i = 0; i < (WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal); i++){
          for(G4double j = 0; j < WCPMTperCellVertical; j++){
            G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.),
                  towerWidth/2.-(i+0.5)*horizontalSpacing + G4RandGauss::shoot(0,pmtPosVar),
                      -barrelCellHeight/2.+(j+0.5)*verticalSpacing+z_offset + G4RandGauss::shoot(0,pmtPosVar));
            G4double newR = towerBSRmin[iz]+(towerBSRmin[iz+1]-towerBSRmin[iz])*(PMTPosition.z()-mainAnnulusZ[iz])/(mainAnnulusZ[iz+1]-mainAnnulusZ[iz]);
            PMTPosition.setX(newR);
            PMTPosition.rotateZ(-(2*pi-totalAngle)/2.); // align with the symmetry 
                                                        //axes of the cell 

            G4VPhysicalVolume* physiWCBarrelPMT =
              new G4PVPlacement(PMTRotation,             // its rotation
                    PMTPosition, 
                    logicWCPMT,                // its logical volume
                    "WCPMT",             // its name
                    logicWCExtraTower,         // its mother volume
                    false,                     // no boolean operations
                    copyNo++,
                    checkOverlapsPMT);                       
            
              // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
              // daughter volumes to the PMTs (e.g. a acryl cover) you have to check, if
              // this is still the case.
          }
        }

      }
    }
  
  }//end if placeBarrelPMTs

  G4LogicalVolume* logicTopCapAssembly = ConstructCapsNoReplica(-1);
  G4LogicalVolume* logicBottomCapAssembly = ConstructCapsNoReplica(1);

 // These lines make the large cap volume invisible to view the caps blacksheets. Need to make invisible for 
 // RayTracer
  if (Vis_Choice == "RayTracer"){
    logicBottomCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
    logicTopCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
  } else {
    logicBottomCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
    logicTopCapAssembly->SetVisAttributes(G4VisAttributes::Invisible);
  }

  G4VPhysicalVolume* physiTopCapAssembly =
  new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(mainAnnulusHeight/2.+capAssemblyHeight/2.)),
                  logicTopCapAssembly,
                  "TopCapAssembly",
                  logicWCBarrel,
                  false, 0,
				  checkOverlaps);

  G4VPhysicalVolume* physiBottomCapAssembly =
    new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(-mainAnnulusHeight/2.-capAssemblyHeight/2.)),
                  logicBottomCapAssembly,
                  "BottomCapAssembly",
                  logicWCBarrel,
                  false, 0,
				  checkOverlaps);

  return logicWC;
}


G4LogicalVolume* WCSimDetectorConstruction::ConstructCapsNoReplica(G4int zflip)
{

  capAssemblyHeight = (WCIDHeight-mainAnnulusHeight)/2+1*mm+WCBlackSheetThickness;

  G4Tubs* solidCapAssembly = new G4Tubs("CapAssembly",
							0.0*m,
              // use the largest radius in cap region
							(outerAnnulusRadius+std::max(GetRadiusChange(-zflip*WCIDHeight/2),GetRadiusChange(-zflip*mainAnnulusHeight/2)))/cos(dPhi/2.), 
							capAssemblyHeight/2,
							0.*deg,
							360.*deg);

  G4LogicalVolume* logicCapAssembly =
    new G4LogicalVolume(solidCapAssembly,
                        G4Material::GetMaterial(water),
                        "CapAssembly",
                        0,0,0);

  
  //----------------------------------------------------
  // extra rings for the top and bottom of the annulus
  //---------------------------------------------------
  G4double borderAnnulusZ[3] = {(-barrelCellHeight/2.-(WCIDRadius-innerAnnulusRadius))*zflip, 
								-barrelCellHeight/2.*zflip,
								barrelCellHeight/2.*zflip};
  // Avoid overflowing into cap
  if (WCIDRadius-innerAnnulusRadius>WCBarrelPMTOffset) borderAnnulusZ[0] = (-barrelCellHeight/2.-WCBarrelPMTOffset)*zflip;
  // Get the radius at global z 
  G4double borderAnnulusRmin[3] = { WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight+(WCIDRadius-innerAnnulusRadius))),
                                    innerAnnulusRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight)), 
                                    innerAnnulusRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2)) };
  G4double borderAnnulusRmax[3] = { outerAnnulusRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight+(WCIDRadius-innerAnnulusRadius))), 
                                    outerAnnulusRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight)),
                                    outerAnnulusRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2)) };
  G4Polyhedra* solidWCBarrelBorderRing = new G4Polyhedra("WCBarrelBorderRing",
                                                   0.*deg, // phi start
                                                   totalAngle,
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   3,
                                                   borderAnnulusZ,
                                                   borderAnnulusRmin,
                                                   borderAnnulusRmax);
  G4LogicalVolume* logicWCBarrelBorderRing =
    new G4LogicalVolume(solidWCBarrelBorderRing,
                        G4Material::GetMaterial(water),
                        "WCBarrelRing",
                        0,0,0);
  //G4cout << *solidWCBarrelBorderRing << G4endl;

  G4VPhysicalVolume* physiWCBarrelBorderRing =
    new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(capAssemblyHeight/2.- barrelCellHeight/2.)*zflip),
                  logicWCBarrelBorderRing,
                  "WCBarrelBorderRing",
                  logicCapAssembly,
                  false, 0,
				  checkOverlaps);
           
  // These lines of code below will turn the border rings invisible. 
  // used for RayTracer
  if (Vis_Choice == "RayTracer"){
    if(!debugMode){
      G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
      tmpVisAtt->SetForceSolid(true);
      logicWCBarrelBorderRing->SetVisAttributes(tmpVisAtt);
      logicWCBarrelBorderRing->SetVisAttributes(G4VisAttributes::Invisible);
    }
    else {
      G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
      tmpVisAtt->SetForceWireframe(true);
      logicWCBarrelBorderRing->SetVisAttributes(tmpVisAtt);
      logicWCBarrelBorderRing->SetVisAttributes(G4VisAttributes::Invisible);
    }
  }
  // used for OGLSX
  else {
    if(!debugMode)
    {
        G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
        tmpVisAtt->SetForceSolid(true);
        logicWCBarrelBorderRing->SetVisAttributes(tmpVisAtt);
        logicWCBarrelBorderRing->SetVisAttributes(G4VisAttributes::Invisible);
    }
    else {
      G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
      //tmpVisAtt->SetForceWireframe(true);
      tmpVisAtt->SetForceSolid(true);
      logicWCBarrelBorderRing->SetVisAttributes(tmpVisAtt);
    }
  }

  //------------------------------------------------------------
  // add blacksheet to the border cells.
  // ---------------------------------------------------------
  // the part between the barrel ring and cap was missing before
  G4double annulusBlackSheetRmax[3] = { WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight+(WCIDRadius-innerAnnulusRadius))) + WCBlackSheetThickness,
                                        WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight)) + WCBlackSheetThickness,
                                        WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2)) + WCBlackSheetThickness};
  G4double annulusBlackSheetRmin[3] = { WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight+(WCIDRadius-innerAnnulusRadius))),
                                        WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2+barrelCellHeight)),
                                        WCIDRadius + GetRadiusChange(-zflip*(mainAnnulusHeight/2))};
  G4Polyhedra* solidWCBarrelBlackSheet = new G4Polyhedra("WCBarrelBorderBlackSheet",
                                                   0, // phi start
                                                   totalAngle, //total phi
                                                   (G4int)WCBarrelRingNPhi, //NPhi-gon
                                                   3,
                                                   borderAnnulusZ,
                                                   annulusBlackSheetRmin,
                                                   annulusBlackSheetRmax);

  logicWCBarrelBorderBlackSheet =
    new G4LogicalVolume(solidWCBarrelBlackSheet,
                        G4Material::GetMaterial("Blacksheet"),
                        "WCBarrelBorderBlackSheet",
                          0,0,0);

  G4VPhysicalVolume* physiWCBarrelBorderBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCBarrelBorderBlackSheet,
                      "WCBarrelBorderBlackSheet",
                      logicWCBarrelBorderRing,
                      false,
                      0,
                      checkOverlaps);

   
  G4LogicalBorderSurface * WaterBSBarrelSurface 
	 = new G4LogicalBorderSurface("WaterBSBarrelBorderSurface",
								  physiWCBarrelBorderRing,
								  physiWCBarrelBorderBlackSheet, 
								  OpWaterBSSurface);
   
  new G4LogicalSkinSurface("BSBarrelBorderSkinSurface",logicWCBarrelBorderBlackSheet,
							BSSkinSurface);

  // Change made here to have the if statement contain the !debugmode to be consistent
  // This code gives the Blacksheet its color. 

  if (Vis_Choice == "RayTracer"){

    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
    WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
	  WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown
    if(!debugMode)
    {
      logicWCBarrelBorderBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }
    else
    {
      logicWCBarrelBorderBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    }
  }

  else {

    G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
      = new G4VisAttributes(G4Colour(0.2,0.9,0.2));
    if(!debugMode)
    {
      logicWCBarrelBorderBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
    }
    else
    {
      logicWCBarrelBorderBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }
  }

  // we have to declare the logical Volumes 
  // outside of the if block to access it later on 
  G4LogicalVolume* logicWCExtraBorderCell;
  G4double towerBSRmin[3];
  G4double towerBSRmax[3];
  if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){
    //----------------------------------------------
    // also the extra tower need special cells at the 
    // top and bottom.
    // (the top cell is created later on by reflecting the
    // bottom cell) 
    //---------------------------------------------
    G4double extraBorderRmin[3];
    G4double extraBorderRmax[3];
    for(int i = 0; i < 3; i++){
      extraBorderRmin[i] = borderAnnulusRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
      extraBorderRmax[i] = borderAnnulusRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
    } 
    G4Polyhedra* solidWCExtraBorderCell = new G4Polyhedra("WCspecialBarrelBorderCell",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //total phi
			   1, //NPhi-gon
			   3,
			   borderAnnulusZ,
			   extraBorderRmin,
			   extraBorderRmax);

    logicWCExtraBorderCell =
      new G4LogicalVolume(solidWCExtraBorderCell, 
			  G4Material::GetMaterial(water),
			  "WCspecialBarrelBorderCell", 
			  0,0,0);
    //G4cout << *solidWCExtraBorderCell << G4endl;

    G4VPhysicalVolume* physiWCExtraBorderCell =
      new G4PVPlacement(0,
                  G4ThreeVector(0.,0.,(capAssemblyHeight/2.- barrelCellHeight/2.)*zflip),
                  logicWCExtraBorderCell,
                  "WCExtraTowerBorderCell",
                  logicCapAssembly,
                  false, 0,
				          checkOverlaps);

    G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
  	tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
  	logicWCExtraBorderCell->SetVisAttributes(tmpVisAtt); 
	  logicWCExtraBorderCell->SetVisAttributes(G4VisAttributes::Invisible);
	  //TF vis.

    // add black sheet to the ExtraBorder cell
    for(int i = 0; i < 3; i++){
      towerBSRmin[i] = annulusBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
      towerBSRmax[i] = annulusBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.);
    }
    G4Polyhedra* solidWCExtraBorderBlackSheet = new G4Polyhedra("WCExtraBorderBlackSheet",
			   totalAngle-2.*pi,//+dPhi/2., // phi start
			   2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //phi end
		     1, //NPhi-gon
			   3,
			   borderAnnulusZ,
			   towerBSRmin,
			   towerBSRmax);
    logicWCExtraBorderBlackSheet =
      new G4LogicalVolume(solidWCExtraBorderBlackSheet,
			  G4Material::GetMaterial("Blacksheet"),
			  "WCExtraBorderBlackSheet",
			    0,0,0);

    G4VPhysicalVolume* physiWCExtraBorderBlackSheet =
      new G4PVPlacement(0,
			G4ThreeVector(0.,0.,0.),
			logicWCExtraBorderBlackSheet,
			"WCExtraBorderBlackSheet",
			logicWCExtraBorderCell,
			false,
			0,
			checkOverlaps);

    G4LogicalBorderSurface * WaterBSExtraBorderCellSurface 
      = new G4LogicalBorderSurface("WaterBSBarrelCellSurface",
				   physiWCExtraBorderCell,
				   physiWCExtraBorderBlackSheet, 
				   OpWaterBSSurface);

    new G4LogicalSkinSurface("BSExtraBorderSkinSurface",logicWCExtraBorderBlackSheet,
                  BSSkinSurface);

    // These lines add color to the blacksheet in the extratower. If using RayTracer, comment the first chunk and use the second. The Blacksheet should be green.
  
    if (Vis_Choice == "RayTracer"){
    
      G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
        = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color
      WCBarrelBlackSheetCellVisAtt->SetForceSolid(true); // force the object to be visualized with a surface
      WCBarrelBlackSheetCellVisAtt->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown

      if(!debugMode)
        logicWCExtraBorderBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
      else
        logicWCExtraBorderBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);
    }

    else {

      G4VisAttributes* WCBarrelBlackSheetCellVisAtt 
        = new G4VisAttributes(G4Colour(0.2,0.9,0.2)); // green color

      if(!debugMode)
        {logicWCExtraBorderBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);}
      else
        {logicWCExtraBorderBlackSheet->SetVisAttributes(WCBarrelBlackSheetCellVisAtt);}
    }

  }
 //------------------------------------------------------------
 // add caps
 // -----------------------------------------------------------
  //crucial to match with borderAnnulusZ
  //need to get correct radius within barrel
  G4double capZ[5] = { (-WCBlackSheetThickness-1.*mm)*zflip,
                        0,
                        (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))*zflip,
                        (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))*zflip,
                        WCBarrelPMTOffset*zflip} ;
  // ensure squential z 
  if(WCBarrelPMTOffset < (WCIDRadius-innerAnnulusRadius)) { capZ[2] = 0; capZ[3] = 0;}

  G4double capRmin[5] = {  0. , 0., 0., 0., 0.} ;
  G4double capRmax[5] = { outerAnnulusRadius + GetRadiusChange(-zflip*WCIDHeight/2), 
                          outerAnnulusRadius + GetRadiusChange(-zflip*WCIDHeight/2),  
                          outerAnnulusRadius + GetRadiusChange(-zflip*(WCIDHeight/2-(WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius)))),  
                          WCIDRadius + GetRadiusChange(-zflip*(WCIDHeight/2-(WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius)))), 
                          innerAnnulusRadius+ GetRadiusChange(-zflip*(WCIDHeight/2-WCBarrelPMTOffset))};
  G4VSolid* solidWCCap;
  if(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal){
    solidWCCap
      = new G4Polyhedra("WCCap",
			0.*deg, // phi start
			totalAngle, //phi end
			(int)WCBarrelRingNPhi, //NPhi-gon
			5, // 2 z-planes
			capZ, //position of the Z planes
			capRmin, // min radius at the z planes
			capRmax// max radius at the Z planes
			);
  } else {
    // if there is an extra tower, the cap volume is a union of
    // to polyhedra. We have to unite both parts, because there are 
    // PMTs that are on the border between both parts.
    G4Polyhedra* mainPart 
      = new G4Polyhedra("WCCapMainPart",
      0.*deg, // phi start
      totalAngle, //phi end
      (int)WCBarrelRingNPhi, //NPhi-gon
      5, // 2 z-planes
      capZ, //position of the Z planes
      capRmin, // min radius at the z planes
      capRmax// max radius at the Z planes
      );
    G4double extraCapRmin[5]; 
    G4double extraCapRmax[5]; 
    for(int i = 0; i < 5 ; i++){
      extraCapRmin[i] = capRmin[i] != 0. ?  capRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
      extraCapRmax[i] = capRmax[i] != 0. ? capRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
    }
    G4Polyhedra* extraSlice 
      = new G4Polyhedra("WCCapExtraSlice",
			totalAngle-2.*pi, // phi start
			2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //total phi 
			// fortunately there are no PMTs an the gap!
			1, //NPhi-gon
			5, //  z-planes
			capZ, //position of the Z planes
			extraCapRmin, // min radius at the z planes
			extraCapRmax// max radius at the Z planes
			);
     solidWCCap =
        new G4UnionSolid("WCCap", mainPart, extraSlice);

     //G4cout << *solidWCCap << G4endl;
   
  }
  // G4cout << *solidWCCap << G4endl;
  G4LogicalVolume* logicWCCap = 
    new G4LogicalVolume(solidWCCap,
			G4Material::GetMaterial(water),
			"WCCapPolygon",
			0,0,0);

  G4VPhysicalVolume* physiWCCap =
    new G4PVPlacement(0,                           // no rotation
		      G4ThreeVector(0.,0.,(-capAssemblyHeight/2.+1*mm+WCBlackSheetThickness)*zflip),     // its position
		      logicWCCap,          // its logical volume
		      "WCCap",             // its name
		      logicCapAssembly,                  // its mother volume
		      false,                       // no boolean operations
		      0,                          // Copy #
			  checkOverlaps);


  // used for RayTracer
  if (Vis_Choice == "RayTracer"){
    if(!debugMode){  
      G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(1,0.5,0.5));
      tmpVisAtt2->SetForceSolid(true);
      logicWCCap->SetVisAttributes(tmpVisAtt2);
      logicWCCap->SetVisAttributes(G4VisAttributes::Invisible);

    } else{
    
      G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(0.6,0.5,0.5));
      tmpVisAtt2->SetForceSolid(true);
      logicWCCap->SetVisAttributes(tmpVisAtt2);
    }
  }

  // used for OGLSX
  else{
    if(!debugMode){  

      G4VisAttributes* tmpVisAtt = new G4VisAttributes(G4Colour(1.,0.5,0.5));
      tmpVisAtt->SetForceSolid(true);// This line is used to give definition to the cells in OGLSX Visualizer
      //logicWCCap->SetVisAttributes(tmpVisAtt); 
      logicWCCap->SetVisAttributes(G4VisAttributes::Invisible);
    //TF vis.

    } else {
      G4VisAttributes* tmpVisAtt2 = new G4VisAttributes(G4Colour(.6,0.5,0.5));
      tmpVisAtt2->SetForceWireframe(true);
      logicWCCap->SetVisAttributes(tmpVisAtt2);
    }
  }

  //---------------------------------------------------------------------
  // add cap blacksheet
  // -------------------------------------------------------------------
  
  G4double capBlackSheetZ[4] = {-WCBlackSheetThickness*zflip, 0., 0., (WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius)) *zflip};
  // ensure squential z 
  if(WCBarrelPMTOffset < (WCIDRadius-innerAnnulusRadius)) capBlackSheetZ[3] = 0;
  G4double capBlackSheetRmin[4] = { 0., 
                                    0., 
                                    WCIDRadius + GetRadiusChange(-zflip*WCIDHeight/2), 
                                    WCIDRadius + GetRadiusChange(-zflip*(WCIDHeight/2-(WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))))};
  G4double capBlackSheetRmax[4] = { WCIDRadius+WCBlackSheetThickness + GetRadiusChange(-zflip*WCIDHeight/2), 
                                    WCIDRadius+WCBlackSheetThickness + GetRadiusChange(-zflip*WCIDHeight/2),
								                    WCIDRadius+WCBlackSheetThickness + GetRadiusChange(-zflip*WCIDHeight/2),
								                    WCIDRadius+WCBlackSheetThickness + GetRadiusChange(-zflip*(WCIDHeight/2-(WCBarrelPMTOffset - (WCIDRadius-innerAnnulusRadius))))};
  G4VSolid* solidWCCapBlackSheet;
  if(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal){
    solidWCCapBlackSheet
      = new G4Polyhedra("WCCapBlackSheet",
			0.*deg, // phi start
			totalAngle, //total phi
			WCBarrelRingNPhi, //NPhi-gon
			4, //  z-planes
			capBlackSheetZ, //position of the Z planes
			capBlackSheetRmin, // min radius at the z planes
			capBlackSheetRmax// max radius at the Z planes
			);
    // G4cout << *solidWCCapBlackSheet << G4endl;
  } else { 
    // same as for the cap volume
    G4Polyhedra* mainPart
      = new G4Polyhedra("WCCapBlackSheetMainPart",
			0.*deg, // phi start
			totalAngle, //phi end
			WCBarrelRingNPhi, //NPhi-gon
			4, //  z-planes
			capBlackSheetZ, //position of the Z planes
			capBlackSheetRmin, // min radius at the z planes
			capBlackSheetRmax// max radius at the Z planes
			);
    G4double extraBSRmin[4];
    G4double extraBSRmax[4];
    for(int i = 0; i < 4 ; i++){
      extraBSRmin[i] = capBlackSheetRmin[i] != 0. ? capBlackSheetRmin[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
      extraBSRmax[i] = capBlackSheetRmax[i] != 0. ? capBlackSheetRmax[i]/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.) : 0.;
    }
    G4Polyhedra* extraSlice
    = new G4Polyhedra("WCCapBlackSheetextraSlice",
      totalAngle-2.*pi, // phi start
      2.*pi -  totalAngle -G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/(10.*m), //
      WCBarrelRingNPhi, //NPhi-gon
      4, //  z-planes
      capBlackSheetZ, //position of the Z planes
      extraBSRmin, // min radius at the z planes
      extraBSRmax// max radius at the Z planes
      );
  
    solidWCCapBlackSheet =
      new G4UnionSolid("WCCapBlackSheet", mainPart, extraSlice);
  }
  G4LogicalVolume* logicWCCapBlackSheet =
    new G4LogicalVolume(solidWCCapBlackSheet,
			G4Material::GetMaterial("Blacksheet"),
			"WCCapBlackSheet",
			0,0,0);
  G4VPhysicalVolume* physiWCCapBlackSheet =
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,0.),
                      logicWCCapBlackSheet,
                      "WCCapBlackSheet",
                      logicWCCap,
                      false,
                      0,
					  checkOverlaps);

  G4LogicalBorderSurface * WaterBSBottomCapSurface 
	= new G4LogicalBorderSurface("WaterBSCapPolySurface",
								 physiWCCap,physiWCCapBlackSheet,
								 OpWaterBSSurface);
  
  new G4LogicalSkinSurface("BSBottomCapSkinSurface",logicWCCapBlackSheet,
						   BSSkinSurface);
  
  G4VisAttributes* WCCapBlackSheetVisAtt 
      = new G4VisAttributes(G4Colour(0.9,0.2,0.2));
    
  // used for OGLSX
  if (Vis_Choice == "OGLSX"){

    G4VisAttributes* WCCapBlackSheetVisAtt 
    = new G4VisAttributes(G4Colour(0.9,0.2,0.2));
  
    if(!debugMode)
          logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible);
      else
          logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);}

  // used for RayTracer (makes the caps blacksheet yellow)
  if (Vis_Choice == "RayTracer"){

    G4VisAttributes* WCCapBlackSheetVisAtt 
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

    if(!debugMode)
        //logicWCCapBlackSheet->SetVisAttributes(G4VisAttributes::Invisible); //Use this line if you want to make the blacksheet on the caps invisible to view through
      logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);
      else
          logicWCCapBlackSheet->SetVisAttributes(WCCapBlackSheetVisAtt);}

  //---------------------------------------------------------
  // Add top and bottom PMTs
  // -----------------------------------------------------
    
  G4LogicalVolume* logicWCPMT = ConstructMultiPMT(WCPMTName, WCIDCollectionName);
  //G4LogicalVolume* logicWCPMT = ConstructPMT(WCPMTName, WCIDCollectionName);
  G4String pmtname = "WCMultiPMT";
 
  // If using RayTracer and want to view the detector without caps, comment out the top and bottom PMT's
  G4double xoffset;
  G4double yoffset;
  G4int    icopy = 0;

  G4RotationMatrix* WCCapPMTRotation = new G4RotationMatrix;
  //if mPMT: perp to wall
  if(orientation == PERPENDICULAR){
    if(zflip==-1){
      WCCapPMTRotation->rotateY(180.*deg); 
    }
  }
  else if (orientation == VERTICAL)
	WCCapPMTRotation->rotateY(90.*deg);
  else if (orientation == HORIZONTAL)
	WCCapPMTRotation->rotateX(90.*deg);

  // loop over the cap
  if(placeCapPMTs){
    G4int CapNCell = WCCapEdgeLimit/WCCapPMTSpacing + 2;
    for ( int i = -CapNCell ; i <  CapNCell; i++) {
      for (int j = -CapNCell ; j <  CapNCell; j++)   {

        // Jun. 04, 2020 by M.Shinoki
        // For IWCD (NuPRISM_mPMT Geometry)
        xoffset = i*WCCapPMTSpacing + WCCapPMTSpacing*0.5 + G4RandGauss::shoot(0,pmtPosVar);
        yoffset = j*WCCapPMTSpacing + WCCapPMTSpacing*0.5 + G4RandGauss::shoot(0,pmtPosVar);
        // For WCTE (NuPRISMBeamTest_mPMT Geometry)
        if (isNuPrismBeamTest || isNuPrismBeamTest_16cShort){
          xoffset = i*WCCapPMTSpacing + G4RandGauss::shoot(0,pmtPosVar);
          yoffset = j*WCCapPMTSpacing + G4RandGauss::shoot(0,pmtPosVar);
        }
        G4ThreeVector cellpos = G4ThreeVector(xoffset, yoffset, 0);     
        //      G4double WCBarrelEffRadius = WCIDDiameter/2. - WCCapPMTSpacing;
        //      double comp = xoffset*xoffset + yoffset*yoffset 
        //	- 2.0 * WCBarrelEffRadius * sqrt(xoffset*xoffset+yoffset*yoffset)
        //	+ WCBarrelEffRadius*WCBarrelEffRadius;
        //      if ( (comp > WCPMTRadius*WCPMTRadius) && ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCCapEdgeLimit) ) {
        if ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) < WCCapEdgeLimit) 

        // for debugging boundary cases: 
        // &&  ((sqrt(xoffset*xoffset + yoffset*yoffset) + WCPMTRadius) > (WCCapEdgeLimit-100)) ) 
        {
        

          G4VPhysicalVolume* physiCapPMT =
            new G4PVPlacement(WCCapPMTRotation,
                  cellpos,                   // its position
                  logicWCPMT,                // its logical volume
                  pmtname, // its name 
                  logicWCCap,         // its mother volume
                  false,                 // no boolean os
                icopy,               // every PMT need a unique id.
                checkOverlapsPMT);
          

        // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
            // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
          // this is still the case.

          icopy++;
        }
      }
    }


    G4cout << "total on cap: " << icopy << "\n";
    G4cout << "Coverage was calculated to be: " << (icopy*WCPMTRadius*WCPMTRadius/(WCIDRadius*WCIDRadius)) << "\n";
  }//end if placeCapPMTs

  ///////////////   Barrel PMT placement

  if(placeBorderPMTs){

    G4double barrelCellWidth = (annulusBlackSheetRmin[1]+annulusBlackSheetRmin[2])*tan(dPhi/2.);
    G4double horizontalSpacing   = barrelCellWidth/WCPMTperCellHorizontal;
    G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

    G4int copyNo = 0;
    for (int iphi=0;iphi<(G4int)WCBarrelRingNPhi;iphi++)
    {
      G4double phi_offset = (iphi+0.5)*dPhi;
      G4double dth = atan((borderAnnulusRmin[2]-borderAnnulusRmin[1])/(borderAnnulusZ[2]-borderAnnulusZ[1]));
      G4RotationMatrix* PMTRotation = new G4RotationMatrix;
      if(orientation == PERPENDICULAR)
        PMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
      else if(orientation == VERTICAL)
        PMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
      else if(orientation == HORIZONTAL)
        PMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
      PMTRotation->rotateX(-phi_offset);//align the PMT with the Cell
      if(orientation == PERPENDICULAR)
        PMTRotation->rotateY(-dth); 
      else if(orientation == VERTICAL)
        PMTRotation->rotateY(-dth); 

      for(G4double i = 0; i < WCPMTperCellHorizontal; i++){
        for(G4double j = 0; j < WCPMTperCellVertical; j++){
          G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius,
                -barrelCellWidth/2.+(i+0.5)*horizontalSpacing + G4RandGauss::shoot(0,pmtPosVar),
                (-barrelCellHeight/2.+(j+0.5)*verticalSpacing)*zflip + G4RandGauss::shoot(0,pmtPosVar));
          G4double newR = annulusBlackSheetRmin[1]+(annulusBlackSheetRmin[2]-annulusBlackSheetRmin[1])*(PMTPosition.z()-borderAnnulusZ[1])/(borderAnnulusZ[2]-borderAnnulusZ[1]);
          PMTPosition.setX(newR);
          PMTPosition.rotateZ(phi_offset);  // align with the symmetry axes of the cell 

          G4VPhysicalVolume* physiWCBarrelBorderPMT =
            new G4PVPlacement(PMTRotation,                      // its rotation
                  PMTPosition,
                  logicWCPMT,                // its logical volume
                  pmtname,             // its name
                  logicWCBarrelBorderRing,         // its mother volume
                  false,                     // no boolean operations
                  copyNo++,
                  checkOverlapsPMT); 

          // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
          // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
          // this is still the case.
        }
      }
    }
  
    //-------------------------------------------------------------
    // Add PMTs in extra Tower if necessary
    //------------------------------------------------------------
    if(!(WCBarrelRingNPhi*WCPMTperCellHorizontal == WCBarrelNumPMTHorizontal)){

      G4double dth = atan((towerBSRmin[2]-towerBSRmin[1])/(borderAnnulusZ[2]-borderAnnulusZ[1]));
      G4RotationMatrix* PMTRotation = new G4RotationMatrix;
      if(orientation == PERPENDICULAR)
        PMTRotation->rotateY(90.*deg); //if mPMT: perp to wall
      else if(orientation == VERTICAL)
        PMTRotation->rotateY(0.*deg); //if mPMT: vertical/aligned to wall
      else if(orientation == HORIZONTAL)
        PMTRotation->rotateX(90.*deg); //if mPMT: horizontal to wall
      PMTRotation->rotateX((2*pi-totalAngle)/2.);//align the PMT with the Cell
      if(orientation == PERPENDICULAR)
        PMTRotation->rotateY(-dth); 
      else if(orientation == VERTICAL)
        PMTRotation->rotateY(-dth); 
                                                  
      G4double towerWidth = (annulusBlackSheetRmin[1]+annulusBlackSheetRmin[2])/2.*tan(2*pi-totalAngle);

      G4double horizontalSpacing   = towerWidth/(WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal);
      G4double verticalSpacing     = barrelCellHeight/WCPMTperCellVertical;

      for(G4double i = 0; i < (WCBarrelNumPMTHorizontal-WCBarrelRingNPhi*WCPMTperCellHorizontal); i++){
        for(G4double j = 0; j < WCPMTperCellVertical; j++){
          G4ThreeVector PMTPosition =  G4ThreeVector(WCIDRadius/cos(dPhi/2.)*cos((2.*pi-totalAngle)/2.),
                towerWidth/2.-(i+0.5)*horizontalSpacing + G4RandGauss::shoot(0,pmtPosVar),
                    (-barrelCellHeight/2.+(j+0.5)*verticalSpacing)*zflip + G4RandGauss::shoot(0,pmtPosVar));
          G4double newR = towerBSRmin[1]+(towerBSRmin[2]-towerBSRmin[1])*(PMTPosition.z()-borderAnnulusZ[1])/(borderAnnulusZ[2]-borderAnnulusZ[1]);
          PMTPosition.setX(newR);
          PMTPosition.rotateZ(-(2*pi-totalAngle)/2.); // align with the symmetry axes of the cell 
          
          G4VPhysicalVolume* physiWCBarrelBorderPMT =
            new G4PVPlacement(PMTRotation,                          // its rotation
                  PMTPosition,
                  logicWCPMT,                // its logical volume
                  "WCPMT",             // its name
                  logicWCExtraBorderCell,         // its mother volume
                  false,                     // no boolean operations
                (int)(i*WCPMTperCellVertical+j),
                  checkOverlapsPMT);

            // logicWCPMT->GetDaughter(0),physiCapPMT is the glass face. If you add more 
            // daugter volumes to the PMTs (e.g. a acryl cover) you have to check, if
            // this is still the case.
        }
      }
    
    }
  }//end if placeBorderPMTs

  return logicCapAssembly;

}
