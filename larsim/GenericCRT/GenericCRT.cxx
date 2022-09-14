/*!
 * Title:   GenericCRT Utility Class
 * Author:  Andrzej Szelc (andrzejs@fnal.gov)
 *
 * Description:
 * Class with Algorithms to convert sim::AuxDetHits to sim::AuxDetSimChannels
 *
 */

#include "GenericCRT.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Evaluator/Evaluator.h"

#include <algorithm> // std::find()
#include <utility>   // std::move()

sim::GenericCRTUtility::GenericCRTUtility(const std::string energyUnitsScale)
{
  HepTool::Evaluator eval;
  eval.setStdMath();
  eval.setSystemOfUnits();

  const std::string scaleExpression = "MeV / " + energyUnitsScale;
  fEnergyUnitsScale = eval.evaluate(scaleExpression.c_str());

  if (eval.status() != 0) fEnergyUnitsScale = 1.;
}

sim::AuxDetIDE sim::GenericCRTUtility::toAuxDetIDE(const sim::AuxDetHit& InputHit) const
{
  sim::AuxDetIDE outputIDE;

  outputIDE.trackID = InputHit.GetTrackID();
  outputIDE.energyDeposited = InputHit.GetEnergyDeposited() * fEnergyUnitsScale;
  outputIDE.entryX = InputHit.GetEntryX();
  outputIDE.entryY = InputHit.GetEntryY();
  outputIDE.entryZ = InputHit.GetEntryZ();
  outputIDE.entryT = InputHit.GetEntryT();
  outputIDE.exitX = InputHit.GetExitX();
  outputIDE.exitY = InputHit.GetExitY();
  outputIDE.exitZ = InputHit.GetExitZ();
  outputIDE.exitT = InputHit.GetExitT();
  outputIDE.exitMomentumX = InputHit.GetExitMomentumX();
  outputIDE.exitMomentumY = InputHit.GetExitMomentumY();
  outputIDE.exitMomentumZ = InputHit.GetExitMomentumZ();

  return outputIDE;
}

std::vector<unsigned int> sim::GenericCRTUtility::GetAuxDetChannels(
  const std::vector<sim::AuxDetHit>& InputHitVector) const
{

  std::vector<unsigned int> AuxDetChanNumber;
  AuxDetChanNumber.reserve(size(InputHitVector));

  for (auto const& hit : InputHitVector) {

    std::vector<unsigned int>::iterator Chanitr =
      std::find(AuxDetChanNumber.begin(), AuxDetChanNumber.end(), hit.GetID());

    if (Chanitr == AuxDetChanNumber.end()) { //If trackID is already in the map, update it
      //if channel ID is not in the set yet, add it
      AuxDetChanNumber.push_back(hit.GetID());
    } //
  }

  return AuxDetChanNumber;
}

sim::AuxDetSimChannel sim::GenericCRTUtility::GetAuxDetSimChannelByNumber(
  const std::vector<sim::AuxDetHit>& InputHitVector,
  unsigned int inputchannel) const
{
  std::vector<sim::AuxDetIDE> IDEvector;
  //loop over sim::AuxDetHits and assign them to AuxDetSimChannels.

  size_t ad_id_no = 9999;
  size_t ad_sen_id_no = 9999;

  for (auto const& auxDetHit : InputHitVector) {

    double xcoordinate = (auxDetHit.GetEntryX() + auxDetHit.GetExitX()) / 2.0;
    double ycoordinate = (auxDetHit.GetEntryY() + auxDetHit.GetExitY()) / 2.0;
    double zcoordinate = (auxDetHit.GetEntryZ() + auxDetHit.GetExitZ()) / 2.0;
    geo::Point_t const worldPos{xcoordinate, ycoordinate, zcoordinate};

    if (auxDetHit.GetID() == inputchannel) // this is the channel we want.
    {
      // Find the IDs given the hit position
      fGeo->FindAuxDetSensitiveAtPosition(worldPos, ad_id_no, ad_sen_id_no, 0.0001);

      mf::LogDebug("GenericCRTUtility")
        << "Found an AuxDetHit with ID " << auxDetHit.GetID() << " for AuxDet ID " << ad_id_no
        << " Sens ID " << ad_sen_id_no << std::endl;

      auto tempIDE = toAuxDetIDE(auxDetHit);

      std::vector<sim::AuxDetIDE>::iterator IDEitr =
        std::find(IDEvector.begin(), IDEvector.end(), tempIDE);

      if (IDEitr != IDEvector.end()) { //If trackID is already in the map, update it
        //Andrzej's note - following logic from AuxDetReadout in Legacy, but why are the other paremeters getting overwritten like that?
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
  const std::vector<sim::AuxDetHit>& InputHitVector) const
{
  auto const auxDetChannels = GetAuxDetChannels(InputHitVector);
  std::vector<sim::AuxDetSimChannel> auxDetVector;
  auxDetVector.reserve(size(auxDetChannels));

  for (auto const channelNum : auxDetChannels) {
    auxDetVector.push_back(GetAuxDetSimChannelByNumber(InputHitVector, channelNum));
  }

  return auxDetVector;
}
