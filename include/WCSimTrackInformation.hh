#ifndef WCSimTrackInformation_h
#define WCSimTrackInformation_h 1


#include "globals.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

// Maximilien Fechner, december 2004
// Information class for flagging the secondaries
// I'm interested in (namely gammas from pi0s and secondaries
// from muon decay 
// TF: Also gamma's from neutron capture and oxygen de-excitation and 
// electrons from pion decay are very relevant!!
class WCSimTrackInformation : public G4VUserTrackInformation {
private:
  G4bool saveit;
  G4int  primaryParentID;
  G4float  photonStartTime;
  G4ThreeVector  photonStartPos;
  G4ThreeVector  photonStartDir;

  // Tracking photon reflection and scattering history
  G4int nReflection;
  G4int nRayScattering;
  G4int nMieScattering;

public:
  WCSimTrackInformation() : saveit(false), primaryParentID(-99), nReflection(0), nRayScattering(0), nMieScattering(0) {}  //TF: initialize to value with NO meaning instead of DN
  WCSimTrackInformation(const WCSimTrackInformation* aninfo) {
      saveit = aninfo->saveit;
      primaryParentID = aninfo->primaryParentID;
      photonStartTime = aninfo->photonStartTime;
      photonStartPos = aninfo->photonStartPos;
      photonStartDir = aninfo->photonStartDir;

      nReflection = 0;
      nRayScattering = 0;
      nMieScattering = 0;
  }
  virtual ~WCSimTrackInformation() {}
  WCSimTrackInformation(const G4Track* );
  
  G4bool isSaved() { return saveit;}
  void WillBeSaved(G4bool choice) { saveit = choice;}

  void SetPrimaryParentID(G4int i) { primaryParentID = i;}
  void SetPhotonStartTime(G4float time) { photonStartTime = time;}
  void SetPhotonStartPos(const G4ThreeVector &pos) { photonStartPos = pos;}
  void SetPhotonStartDir(const G4ThreeVector &dir) { photonStartDir = dir;}
  G4int GetPrimaryParentID() {return primaryParentID;}
  G4float GetPhotonStartTime() {return photonStartTime;}
  G4ThreeVector GetPhotonStartPos() {return photonStartPos;}
  G4ThreeVector GetPhotonStartDir() {return photonStartDir;}

  void ResetPhotonHistory() { nReflection=0; nRayScattering=0; nMieScattering=0; }
  void AddRaySct() { nRayScattering++; }
  void AddMieSct() { nMieScattering++; }
  void AddReflec() { nReflection++; }

  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator ==(const WCSimTrackInformation& right) const
  {return (this==&right);}

  void Print() const;

};

extern G4Allocator<WCSimTrackInformation> aWCSimTrackInfoAllocator;

inline void* WCSimTrackInformation::operator new(size_t)
{ void* aTrackInfo;
 aTrackInfo = (void*)aWCSimTrackInfoAllocator.MallocSingle();
 return aTrackInfo;
}

inline void WCSimTrackInformation::operator delete(void *aTrackInfo)
{ aWCSimTrackInfoAllocator.FreeSingle((WCSimTrackInformation*)aTrackInfo);}


#endif
