#ifndef MCTRACKCOLLECTIONANAALG_H
#define MCTRACKCOLLECTIONANAALG_H

/*!
 * Title:   MCTrackCollectionAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg to put properties of collection of MCTracks in a tree.
 * 
 */

class TTree;

namespace sim{

  class MCTrack;

  class MCTrackCollectionAnaAlg{

  public:

    MCTrackCollectionAnaAlg(){};

    void SetOutputTree(TTree*, bool fill=true);

    void FillTree(unsigned int, unsigned int,
		  const std::vector<sim::MCTrack>&);
    
  private:

    TTree* fTree;
    bool   fFillTree;

    unsigned int fRun;
    unsigned int fEvent;

    unsigned int fNMCTracks;
    unsigned int fDParticle;
    float        fDParticleFraction;

    int          fDParticleOrigin;
    int          fDParticlePdgCode;
    unsigned int fDParticleTrackId;
    float        fDParticleStartY;
    float        fDParticleStartZ;
    float        fDParticleStartX;
    float        fDParticleStartE;
    float        fDParticleEndY;
    float        fDParticleEndZ;
    float        fDParticleEndX;
    float        fDParticleEndE;
    int          fDParticleMotherPdgCode;
    unsigned int fDParticleMotherTrackId;
    int          fDParticleAncestorPdgCode;
    unsigned int fDParticleAncestorTrackId;

    float        fCollectionY;
    float        fCollectionZ;
    float        fCollectionX;
    float        fCollectionRMSY;
    float        fCollectionRMSZ;
    float        fCollectionRMSX;
    float        fCollectionEnergy;
    float        fMinX;
    float        fMaxX;
    float        fMinY;
    float        fMaxY;
    float        fMinZ;
    float        fMaxZ;

    void         FillDominantParticleInfo(const sim::MCTrack&);
    
  };
  
}

#endif
