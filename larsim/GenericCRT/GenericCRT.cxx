/*!
 * Title:   GenericCRT Utility Class
 * Author:  Andrzej Szelc (andrzejs@fnal.gov)
 *
 * Description:
 * Class with Algorithms to convert sim::AuxDetHits to sim::AuxDetSimChannels
 *
*/


#include "GenericCRT.h"

#include <utility> // std::move()
#include <algorithm> // std::find()




// sim::GenericCRTUtility::GenericCRTUtility()
// {
//
//
// }



unsigned int sim::GenericCRTUtility::CopyAuxDetHittoAuxDetIDE(const sim::AuxDetHit &InputHit,sim::AuxDetIDE &outputIDE)
{


   outputIDE.trackID		    = InputHit.GetTrackID();
   outputIDE.energyDeposited	= InputHit.GetEnergyDeposited();
   outputIDE.entryX		= InputHit.GetEntryX();
   outputIDE.entryY		= InputHit.GetEntryY();
   outputIDE.entryZ		= InputHit.GetEntryZ();
   outputIDE.entryT		= InputHit.GetEntryT();
   outputIDE.exitX		= InputHit.GetExitX();
   outputIDE.exitY		= InputHit.GetExitY();
   outputIDE.exitZ		= InputHit.GetExitZ();
   outputIDE.exitT		= InputHit.GetExitT();
   outputIDE.exitMomentumX	= InputHit.GetExitMomentumX();
   outputIDE.exitMomentumY	= InputHit.GetExitMomentumY();
   outputIDE.exitMomentumZ	= InputHit.GetExitMomentumZ();


   return  InputHit.GetID();
}




unsigned int sim::GenericCRTUtility::GetNumberofAuxDetChannels(const std::vector<sim::AuxDetHit> &InputHitVector,std::vector<unsigned int> &AuxDetChanNumbers )
{

  AuxDetChanNumbers.reserve(InputHitVector.size());


  for(auto HitIter : InputHitVector)
  {

  std::vector<unsigned int>::iterator Chanitr
      = std::find(AuxDetChanNumbers.begin(), AuxDetChanNumbers.end(), HitIter.GetID());

   if(Chanitr == AuxDetChanNumbers.end()){ //If trackID is already in the map, update it
         //if channel ID is not in the set yet, add it
      AuxDetChanNumbers.push_back(HitIter.GetID());
        }//

  }

return AuxDetChanNumbers.size();

}



sim::AuxDetSimChannel const sim::GenericCRTUtility::GetAuxDetSimChannelByNumber( const std::vector<sim::AuxDetHit> &InputHitVector,unsigned int inputchannel)
{
        std::vector<sim::AuxDetIDE> IDEvector;
       IDEvector.resize(0);
        //loop over sim::AuxDetHits and assign them to AuxDetSimChannels.

        size_t ad_id_no = 9999;
        size_t ad_sen_id_no = 9999;

        for(auto AuxDetHitIter : InputHitVector)
        {

            unsigned int channel = AuxDetHitIter.GetID();
            double xcoordinate = (AuxDetHitIter.GetEntryX() + AuxDetHitIter.GetExitX())/2.0;
            double ycoordinate = (AuxDetHitIter.GetEntryY() + AuxDetHitIter.GetExitY())/2.0;
            double zcoordinate = (AuxDetHitIter.GetEntryZ() + AuxDetHitIter.GetExitZ())/2.0;
            double worldPos[3] = {xcoordinate,ycoordinate,zcoordinate};
            fGeo->FindAuxDetSensitiveAtPosition(worldPos, ad_id_no, ad_sen_id_no, 0.0001);
            std::cout << "AD,ADS = " << ad_id_no << "," << ad_sen_id_no << std::endl;
            if(channel == inputchannel)   // this is the channel we want.
            {
                sim::AuxDetIDE tempIDE;
                CopyAuxDetHittoAuxDetIDE(AuxDetHitIter,tempIDE);

                std::vector<sim::AuxDetIDE>::iterator IDEitr
                = std::find(IDEvector.begin(), IDEvector.end(), tempIDE);

                if(IDEitr != IDEvector.end()){ //If trackID is already in the map, update it
                    //Andrzej's note - following logic from AuxDetReadout in Legacy, but why are the other paremeters getting overwritten like that?
                    IDEitr->energyDeposited += tempIDE.energyDeposited;
                    IDEitr->exitX            = tempIDE.exitX;
                    IDEitr->exitY            = tempIDE.exitY;
                    IDEitr->exitZ            = tempIDE.exitZ;
                    IDEitr->exitT            = tempIDE.exitT;
                    IDEitr->exitMomentumX    = tempIDE.exitMomentumX;
                    IDEitr->exitMomentumY    = tempIDE.exitMomentumY;
                    IDEitr->exitMomentumZ    = tempIDE.exitMomentumZ;
                    }
                    else{  //if trackID is not in the set yet, add it
                        IDEvector.push_back(std::move(tempIDE));
                    }//else


            } // end if the AuxDetHit channel checks out.

        } // end main loop on AuxDetHit

        //push back the AuxDetSimChannel Vector.
        //TODO check the last parameter values.
       const sim::AuxDetSimChannel adsc=sim::AuxDetSimChannel(ad_id_no, std::move(IDEvector), ad_sen_id_no);

  return adsc;
}




void sim::GenericCRTUtility::FillAuxDetSimChannels(const std::vector<sim::AuxDetHit> &InputHitVector, std::vector<sim::AuxDetSimChannel> *AuxDetVector)
{

    std::vector<unsigned int> AuxDetChanNumbers;   // vector of channels that have AuxDetHits.

    //get number of NauxDetChannels and fill a map with their IDs.
    unsigned int NAuxDetChannels=GetNumberofAuxDetChannels(InputHitVector,AuxDetChanNumbers );
    AuxDetVector->reserve(NAuxDetChannels);


    //loop on AuxDetChannelNumbers
    for(unsigned iAuxDet=0;iAuxDet < AuxDetChanNumbers.size();iAuxDet++ )
    {
     const sim::AuxDetSimChannel adsc=GetAuxDetSimChannelByNumber(InputHitVector,AuxDetChanNumbers[iAuxDet]);
     AuxDetVector->push_back(std::move(adsc));

    }//end loop on AuxDetChanNumbers
}
