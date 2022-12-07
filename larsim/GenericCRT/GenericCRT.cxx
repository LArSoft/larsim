/*!
 * Title:   GenericCRT Utility Class
 * Author:  Andrzej Szelc (andrzejs@fnal.gov)
 *
 * Description:
 * Class with Algorithms to convert sim::AuxDetHits to sim::AuxDetSimChannels
 *
 */

#include "GenericCRT.h"

#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Evaluator/Evaluator.h"

#include <algorithm> // std::find()
#include <utility>   // std::move()

sim::GenericCRTUtility::GenericCRTUtility(std::string const& energyUnitsScale,
                                          geo::AuxDetGeometryCore const& auxDetGeom)
  : fAuxDetGeom{&auxDetGeom}
{
  HepTool::Evaluator eval;
  eval.setStdMath();
  eval.setSystemOfUnits();

  const std::string scaleExpression = "MeV / " + energyUnitsScale;
  fEnergyUnitsScale = eval.evaluate(scaleExpression.c_str());

  if (eval.status() != 0) fEnergyUnitsScale = 1.;
}

sim::AuxDetIDE sim::GenericCRTUtility::toAuxDetIDE(sim::AuxDetHit const& adhit) const
{
  sim::AuxDetIDE outputIDE;

  outputIDE.trackID = adhit.GetTrackID();
  outputIDE.energyDeposited = adhit.GetEnergyDeposited() * fEnergyUnitsScale;
  outputIDE.entryX = adhit.GetEntryX();
  outputIDE.entryY = adhit.GetEntryY();
  outputIDE.entryZ = adhit.GetEntryZ();
  outputIDE.entryT = adhit.GetEntryT();
  outputIDE.exitX = adhit.GetExitX();
  outputIDE.exitY = adhit.GetExitY();
  outputIDE.exitZ = adhit.GetExitZ();
  outputIDE.exitT = adhit.GetExitT();
  outputIDE.exitMomentumX = adhit.GetExitMomentumX();
  outputIDE.exitMomentumY = adhit.GetExitMomentumY();
  outputIDE.exitMomentumZ = adhit.GetExitMomentumZ();

  return outputIDE;
}

std::vector<unsigned int> sim::GenericCRTUtility::GetAuxDetChannels(
  std::vector<sim::AuxDetHit> const& adhits) const
{
  std::vector<unsigned int> AuxDetChanNumber;
  AuxDetChanNumber.reserve(size(adhits));

  for (auto const& hit : adhits) {

    std::vector<unsigned int>::iterator Chanitr =
      std::find(AuxDetChanNumber.begin(), AuxDetChanNumber.end(), hit.GetID());

    if (Chanitr == AuxDetChanNumber.end()) { //If trackID is already in the map, update it
      //if channel ID is not in the set yet, add it
      AuxDetChanNumber.push_back(hit.GetID());
    }
  }

  return AuxDetChanNumber;
}

sim::AuxDetSimChannel sim::GenericCRTUtility::GetAuxDetSimChannelByNumber(
  std::vector<sim::AuxDetHit> const& adhits,
  unsigned int inputchannel) const
{

  //loop over sim::AuxDetHits and assign them to AuxDetSimChannels.
  std::vector<sim::AuxDetIDE> IDEvector;

  size_t ad_id_no = 9999;
  size_t ad_sen_id_no = 9999;

  for (auto const& auxDetHit : adhits) {

    double xcoordinate = (auxDetHit.GetEntryX() + auxDetHit.GetExitX()) / 2.0;
    double ycoordinate = (auxDetHit.GetEntryY() + auxDetHit.GetExitY()) / 2.0;
    double zcoordinate = (auxDetHit.GetEntryZ() + auxDetHit.GetExitZ()) / 2.0;
    geo::Point_t const worldPos{xcoordinate, ycoordinate, zcoordinate};

    if (auxDetHit.GetID() == inputchannel) // this is the channel we want.
    {
      // Find the IDs given the hit position
      fAuxDetGeom->FindAuxDetSensitiveAtPosition(worldPos, ad_id_no, ad_sen_id_no, 0.0001);

      mf::LogDebug("GenericCRTUtility")
        << "Found an AuxDetHit with ID " << auxDetHit.GetID() << " for AuxDet ID " << ad_id_no
        << " Sens ID " << ad_sen_id_no << std::endl;

      auto tempIDE = toAuxDetIDE(auxDetHit);

      std::vector<sim::AuxDetIDE>::iterator IDEitr =
        std::find(IDEvector.begin(), IDEvector.end(), tempIDE);

      if (IDEitr != IDEvector.end()) { //If trackID is already in the map, update it
        // Andrzej's note - following logic from AuxDetReadout in Legacy, but why are the
        // other paremeters getting overwritten like that?
        IDEitr->energyDeposited += tempIDE.energyDeposited;
        IDEitr->exitX = tempIDE.exitX;
        IDEitr->exitY = tempIDE.exitY;
        IDEitr->exitZ = tempIDE.exitZ;
        IDEitr->exitT = tempIDE.exitT;
        IDEitr->exitMomentumX = tempIDE.exitMomentumX;
        IDEitr->exitMomentumY = tempIDE.exitMomentumY;
        IDEitr->exitMomentumZ = tempIDE.exitMomentumZ;
      }
      else { //if trackID is not in the set yet, add it
        IDEvector.push_back(std::move(tempIDE));
      } //else

      // break;
    } // end if the AuxDetHit channel checks out.

  } // end main loop on AuxDetHit

  mf::LogDebug("GenericCRTUtility")
    << "Returning AuxDetSimChannel for ID " << ad_id_no << " " << ad_sen_id_no << ", with "
    << IDEvector.size() << " IDEs." << std::endl;

  //push back the AuxDetSimChannel Vector.
  //TODO check the last parameter values.
  return sim::AuxDetSimChannel(ad_id_no, std::move(IDEvector), ad_sen_id_no);
}

std::vector<sim::AuxDetSimChannel> sim::GenericCRTUtility::GetAuxDetSimChannels(
  std::vector<sim::AuxDetHit> const& adhits) const
{
  auto const auxDetChannels = GetAuxDetChannels(adhits);
  std::vector<sim::AuxDetSimChannel> auxDetVector;
  auxDetVector.reserve(size(auxDetChannels));

  for (auto const channelNum : auxDetChannels) {
    auxDetVector.push_back(GetAuxDetSimChannelByNumber(adhits, channelNum));
  }

  return auxDetVector;
}
