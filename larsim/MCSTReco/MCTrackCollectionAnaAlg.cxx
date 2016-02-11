/*!
 * Title:   MCTrackCollectionAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg to put properties of collection of MCTracks in a tree.
 * 
 */

#include "lardata/MCBase/MCTrack.h"
#include "TTree.h"
#include <numeric>

#include "MCTrackCollectionAnaAlg.h"

void sim::MCTrackCollectionAnaAlg::SetOutputTree(TTree* tree, bool fill)
{
  fTree = tree;
  fFillTree = fill;

  fTree->Branch("mct_run",&fRun,"mct_run/i");
  fTree->Branch("mct_event",&fEvent,"mct_event/i");

  fTree->Branch("mct_nmctracks",&fNMCTracks,"mct_nmctracks/i");
  fTree->Branch("mct_dp",&fDParticle,"mct_dp/i");
  fTree->Branch("mct_dpfrac",&fDParticleFraction,"mct_dpfrac/F");

  fTree->Branch("mct_dporigin",&fDParticleOrigin,"mct_dporigin/I");
  fTree->Branch("mct_dppdg",&fDParticlePdgCode,"mct_dppdg/I");
  fTree->Branch("mct_dptrackid",&fDParticleTrackId,"mct_dptrackid/i");
  fTree->Branch("mct_dpmpdg",&fDParticleMotherPdgCode,"mct_dpmpdg/I");
  fTree->Branch("mct_dpmtrackid",&fDParticleMotherTrackId,"mct_dpmtrackid/i");
  fTree->Branch("mct_dpapdg",&fDParticleAncestorPdgCode,"mct_dpapdg/I");
  fTree->Branch("mct_dpatrackid",&fDParticleAncestorTrackId,"mct_dpatrackid/i");

  fTree->Branch("mct_dpstartx",&fDParticleStartX,"mct_dpstartx/F");
  fTree->Branch("mct_dpstarty",&fDParticleStartY,"mct_dpstarty/F");
  fTree->Branch("mct_dpstartz",&fDParticleStartZ,"mct_dpstartz/F");
  fTree->Branch("mct_dpstarte",&fDParticleStartE,"mct_dpstarte/F");
  fTree->Branch("mct_dpendx",&fDParticleEndX,"mct_dpendx/F");
  fTree->Branch("mct_dpendy",&fDParticleEndY,"mct_dpendy/F");
  fTree->Branch("mct_dpendz",&fDParticleEndZ,"mct_dpendz/F");
  fTree->Branch("mct_dpende",&fDParticleEndE,"mct_dpende/F");

  fTree->Branch("mct_coly",&fCollectionY,"mct_coly/F");
  fTree->Branch("mct_colz",&fCollectionZ,"mct_colz/F");
  fTree->Branch("mct_colx",&fCollectionX,"mct_colx/F");
  fTree->Branch("mct_colrmsy",&fCollectionRMSY,"mct_colrmsy/F");
  fTree->Branch("mct_colrmsz",&fCollectionRMSZ,"mct_colrmsz/F");
  fTree->Branch("mct_colrmsx",&fCollectionRMSX,"mct_colrmsx/F");
  fTree->Branch("mct_colE",&fCollectionEnergy,"mct_colE/F");

  fTree->Branch("mct_miny",&fMinY,"mct_miny/F");
  fTree->Branch("mct_maxy",&fMaxY,"mct_maxy/F");
  fTree->Branch("mct_minz",&fMinZ,"mct_minz/F");
  fTree->Branch("mct_maxz",&fMaxZ,"mct_maxz/F");
  fTree->Branch("mct_minx",&fMinX,"mct_minx/F");
  fTree->Branch("mct_maxx",&fMaxX,"mct_maxx/F");

}

void sim::MCTrackCollectionAnaAlg::FillTree(unsigned int run, unsigned int event,
					    const std::vector<sim::MCTrack>& mctVec)
{
  fRun = run;
  fEvent = event;

  fNMCTracks = mctVec.size();
  fCollectionEnergy = 0;

  fMinY = 99999;
  fMinZ = 99999;
  fMinX = 99999;
  fMaxY = -99999;
  fMaxZ = -99999;
  fMaxX = -99999;

  if(mctVec.size()==0){
    fTree->Fill();
    return;
  }
    

  std::vector<double> y_vals;
  std::vector<double> z_vals;
  std::vector<double> x_vals;
  std::vector<double> stepL_vals;

  std::vector<double> mct_length(mctVec.size(),0);
  int d_index = -1;
  for(size_t i_p=0; i_p<mctVec.size(); i_p++){


    fCollectionEnergy += mctVec[i_p].Start().E();

    //y_vals.reserve( y_vals.size() + mctVec[i_p].size()-1 );
    //z_vals.reserve( z_vals.size() + mctVec[i_p].size()-1 );
    //x_vals.reserve( x_vals.size() + mctVec[i_p].size()-1 );

    for(size_t i_s=1; i_s < mctVec[i_p].size(); i_s++){

      TVector3 const& vec1 = mctVec[i_p][i_s-1].Position().Vect();
      TVector3 const& vec2 = mctVec[i_p][i_s].Position().Vect();
      double stepL = (vec2-vec1).Mag();

      double thisy = 0.5*(vec1.Y() + vec2.Y());
      double thisz = 0.5*(vec1.Z() + vec2.Z());
      double thisx = 0.5*(vec1.X() + vec2.X());

      y_vals.push_back(thisy);
      z_vals.push_back(thisz);
      x_vals.push_back(thisx);
      stepL_vals.push_back(stepL);
      
      mct_length[i_p] += stepL;

      if(thisy > fMaxY) fMaxY = thisy;
      if(thisy < fMinY) fMinY = thisy;
      if(thisz > fMaxZ) fMaxZ = thisz;
      if(thisz < fMinZ) fMinZ = thisz;
      if(thisx > fMaxX) fMaxX = thisx;
      if(thisx < fMinX) fMinX = thisx;

    }//end loop over steps
  }//end loop over tracks

  double totalL = std::accumulate(mct_length.begin(),mct_length.end(),0.0);

  fDParticle = std::distance(mct_length.begin(),std::max_element(mct_length.begin(),mct_length.end()));
  fDParticleFraction = mct_length.at(fDParticle) / totalL;
  FillDominantParticleInfo(mctVec.at(fDParticle));

  double sumy=0,sumz=0,sumx=0;
  for(size_t i_step=0; i_step<stepL_vals.size(); i_step++){
    sumy += stepL_vals[i_step]*y_vals[i_step];
    sumz += stepL_vals[i_step]*z_vals[i_step];
    sumx += stepL_vals[i_step]*x_vals[i_step];
  }
  
  fCollectionY = sumy/totalL;
  fCollectionZ = sumz/totalL;
  fCollectionX = sumx/totalL;
  
  double sumy2=0,sumz2=0,sumx2=0;
  for(size_t i_step=0; i_step<stepL_vals.size(); i_step++){
    sumy2 += stepL_vals[i_step]*(y_vals[i_step]-fCollectionY)*(y_vals[i_step]-fCollectionY);
    sumz2 += stepL_vals[i_step]*(z_vals[i_step]-fCollectionZ)*(z_vals[i_step]-fCollectionZ);
    sumx2 += stepL_vals[i_step]*(x_vals[i_step]-fCollectionX)*(x_vals[i_step]-fCollectionX);
  }

  fCollectionRMSY = std::sqrt(sumy2/totalL);
  fCollectionRMSZ = std::sqrt(sumz2/totalL);
  fCollectionRMSX = std::sqrt(sumx2/totalL);

  if(fFillTree) fTree->Fill();
}

void sim::MCTrackCollectionAnaAlg::FillDominantParticleInfo(const sim::MCTrack& mctrack)
{
  size_t nsteps = mctrack.size();

  fDParticleOrigin = mctrack.Origin();
  fDParticlePdgCode = mctrack.PdgCode();
  fDParticleTrackId = mctrack.TrackID();
  fDParticleStartY = mctrack[0].Y();
  fDParticleStartZ = mctrack[0].Z();
  fDParticleStartX = mctrack[0].X();
  fDParticleStartE = mctrack[0].E();
  fDParticleEndY = mctrack[nsteps-1].Y();
  fDParticleEndZ = mctrack[nsteps-1].Z();
  fDParticleEndX = mctrack[nsteps-1].X();
  fDParticleEndE = mctrack[nsteps-1].E();
  fDParticleMotherPdgCode = mctrack.MotherPdgCode();
  fDParticleMotherTrackId = mctrack.MotherTrackID();
  fDParticleAncestorPdgCode = mctrack.AncestorPdgCode();
  fDParticleAncestorTrackId = mctrack.AncestorTrackID();
}
