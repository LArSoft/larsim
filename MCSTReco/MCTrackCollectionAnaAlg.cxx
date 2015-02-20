/*!
 * Title:   MCTrackCollectionAnaAlg Class
 * Author:  Wes Ketchum (wketchum@lanl.gov)
 *
 * Description:
 * Alg to put properties of collection of MCTracks in a tree.
 * 
 */

#include "MCBase/MCTrack.h"
#include "TTree.h"
#include <algorithm>

#include "MCTrackCollectionAnaAlg.h"

void sim::MCTrackCollectionAnaAlg::SetOutputTree(TTree* tree)
{
  fTree = tree;

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
  fTree->Branch("mct_dpstartp",&fDParticleStartP,"mct_dpstartp/F");
  fTree->Branch("mct_dpendx",&fDParticleEndX,"mct_dpendx/F");
  fTree->Branch("mct_dpendy",&fDParticleEndY,"mct_dpendy/F");
  fTree->Branch("mct_dpendz",&fDParticleEndZ,"mct_dpendz/F");
  fTree->Branch("mct_dpende",&fDParticleEndE,"mct_dpende/F");
  fTree->Branch("mct_dpendp",&fDParticleEndP,"mct_dpendp/F");

  fTree->Branch("mct_coly",&fCollectionY,"mct_coly/F"):
  fTree->Branch("mct_colz",&fCollectionZ,"mct_colz/F"):
  fTree->Branch("mct_colx",&fCollectionX,"mct_colx/F"):
  fTree->Branch("mct_colrmsy",&fCollectionRMSY,"mct_colrmsy/F"):
  fTree->Branch("mct_colrmsz",&fCollectionRMSZ,"mct_colrmsz/F"):
  fTree->Branch("mct_colrmsx",&fCollectionRMSX,"mct_colrmsx/F"):
  fTree->Branch("mct_colE",&fCollectionEnergy,"mct_colE/F"):

  fTree->Branch("mct_miny",&fMinY,"mct_miny/F"):
  fTree->Branch("mct_maxy",&fMaxY,"mct_maxy/F"):
  fTree->Branch("mct_minz",&fMinZ,"mct_minz/F"):
  fTree->Branch("mct_maxz",&fMaxZ,"mct_maxz/F"):
  fTree->Branch("mct_minx",&fMinX,"mct_minx/F"):
  fTree->Branch("mct_maxx",&fMaxX,"mct_maxx/F"):

}

void sim::MCTrackCollectionAnaAlg::FillTree(unsigned int run, unsigned int event,
					    const std::vector<sim::MCTrack>& mctVec)
{
  fRun = run;
  fEvent = event;

  fNMCTracks = mctVec.size();
  double totalE = 0;

  fMinY = 99999;
  fMinZ = 99999;
  fMinX = 99999;
  fMaxY = -99999;
  fMaxZ = -99999;
  fMaxX = -99999;

  std::vector<double> y_vals;
  std::vector<double> z_vals;
  std::vector<double> x_vals;

  std::vector<double> mct_length(mctVec.size(),0);
  int d_index = -1;
  for(size_t i_p=0; i_p<mctVec.size(); i_p++){


    totalE += mctVec[i_p].Start().E();

    y_vals.reserve( y_vals.size() + mctVec[i_p].size()-1 );
    z_vals.reserve( z_vals.size() + mctVec[i_p].size()-1 );
    x_vals.reserve( x_vals.size() + mctVec[i_p].size()-1 );

    for(size_t i_s=1; i_s < mctVec[i_p].size(); i_s++){

      TVector3 const& vec1 = mctVec[i_p][i_s-1].Position().Vect();
      TVector3 const& vec2 = mctVec[i_p][i_s].Position().Vect();
      double stepL = (vec2-vec1).Mag();

      double thisy = 0.5*(vec1.Y() + vec2.Y());
      double thisy = 0.5*(vec1.Z() + vec2.Z());
      double thisy = 0.5*(vec1.X() + vec2.X());

      y_vals.push_back(thisy*stepL);
      z_vals.push_back(thisz*stepL);
      x_vals.push_back(thisx*stepL);

      mct_length[i_p] += stepL;

      if(thisy > fMaxY) thisy = fMaxY;
      if(thisy < fMinY) thisy = fMinY;
      if(thisz > fMaxZ) thisz = fMaxZ;
      if(thisz < fMinZ) thisz = fMinZ;
      if(thisx > fMaxX) thisx = fMaxX;
      if(thisx < fMinX) thisx = fMinX;

    }//end loop over steps
  }//end loop over tracks

  double totalL = std::accumulate(mct_length.begin(),mct_length.end());

  fDParticle = std::distance(mct_length.begin(),std::max(mct_length.begin(),mct_length.end()));
  fDParticleFraction = mct_length.at(fDParticle) / totalL;
  FillDominantParticleInfo(mctVec.at(fDParticle));

  fCollectionY = std::accumulate(y_vals.begin(),y_vals.end())/totalL;
  fCollectionZ = std::accumulate(z_vals.begin(),z_vals.end())/totalL;
  fCollectionX = std::accumulate(x_vals.begin(),x_vals.end())/totalL;
  
  double sumy2=0;
  for(auto const& y : y_vals)
    sumy2 += (y-fCollectionY)*(y-fCollectionY);
  fCollectionRMSY = std::sqrt(sumy2/totalL);

  double sumz2=0;
  for(auto const& z : z_vals)
    sumz2 += (z-fCollectionZ)*(z-fCollectionZ);
  fCollectionRMSZ = std::sqrt(sumz2/totalL);

  double sumx2=0;
  for(auto const& x : x_vals)
    sumx2 += (x-fCollectionX)*(x-fCollectionX);
  fCollectionRMSX = std::sqrt(sumx2/totalL);
}

void sim::MCTrackCollectionAnaAlg::FillDominantParticleInfo(const sim::MCTrack& mctrack)
{
}
