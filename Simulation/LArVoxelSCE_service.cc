///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelSCE_service.cc
/// \brief Implements offset in reconstructed position of ionization electrons due to SCE (Space Charge Effect)
///
/// \author  Michael Mooney (mooney@bnl.gov)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// This service provides tools to calculate the offset in reconstructed position of ioniziation electrons 
/// deposited within the TPC due to the space charge effect, which is utilized in the LArVoxelReadout class.  
/// The offsets in X (T), Y, and Z are calculated using an external program and it is expected that they will
/// be provided in a variety of different forms (parametric, matrix, interpolation engine, etc.) which can be
/// selected between via the ParameterSet.  For now the offsets are encoded in parametric form using a set of
/// TGraphs stored in a ROOT file.

/// This service is a temporary "dirty" hack and will eventually be superceded with a database solution located
/// elsewhere in the framework.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

#include "cetlib/exception.h"
#include "Simulation/LArVoxelSCE.h"

namespace sim {

  //----------------------------------------------------------------------------
  /// Constructor
  LArVoxelSCE::LArVoxelSCE(fhicl::ParameterSet const& pset, 
					 art::ActivityRegistry & /* reg */) 
  {
    this->reconfigure(pset);
  }

  //----------------------------------------------------------------------------
  /// Destructor
  LArVoxelSCE::~LArVoxelSCE()
  {
  }

  //----------------------------------------------------------------------------
  void LArVoxelSCE::reconfigure(fhicl::ParameterSet const& pset)
  {
    fEnableSCE = fLgpHandle->EnableSCE();

    if(fEnableSCE == true)
    {
      fRepresentationType = pset.get<std::string>("RepresentationType");
      fInputFilename = pset.get<std::string>("InputFilename");

      std::string fname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fInputFilename,fname);

      std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
      if(!infile->IsOpen()) throw art::Exception( art::errors::NotFound ) << "Could not find the space charge effect file " << fname << "!" << std::endl;

      if(fRepresentationType == "Parametric")
      {      
        for(int i = 0; i < 5; i++)
        {
          g1_x[i] = (TGraph*)infile->Get(Form("deltaX/g1_%d",i));
          g2_x[i] = (TGraph*)infile->Get(Form("deltaX/g2_%d",i));
          g3_x[i] = (TGraph*)infile->Get(Form("deltaX/g3_%d",i));   
          g4_x[i] = (TGraph*)infile->Get(Form("deltaX/g4_%d",i));
          g5_x[i] = (TGraph*)infile->Get(Form("deltaX/g5_%d",i));

          g1_y[i] = (TGraph*)infile->Get(Form("deltaY/g1_%d",i));
          g2_y[i] = (TGraph*)infile->Get(Form("deltaY/g2_%d",i));
          g3_y[i] = (TGraph*)infile->Get(Form("deltaY/g3_%d",i));   
          g4_y[i] = (TGraph*)infile->Get(Form("deltaY/g4_%d",i));
          g5_y[i] = (TGraph*)infile->Get(Form("deltaY/g5_%d",i));
          g6_y[i] = (TGraph*)infile->Get(Form("deltaY/g6_%d",i));

          g1_z[i] = (TGraph*)infile->Get(Form("deltaZ/g1_%d",i));
          g2_z[i] = (TGraph*)infile->Get(Form("deltaZ/g2_%d",i));
          g3_z[i] = (TGraph*)infile->Get(Form("deltaZ/g3_%d",i));   
          g4_z[i] = (TGraph*)infile->Get(Form("deltaZ/g4_%d",i));
        }

        g1_x[5] = (TGraph*)infile->Get("deltaX/g1_5");
        g2_x[5] = (TGraph*)infile->Get("deltaX/g2_5");
        g3_x[5] = (TGraph*)infile->Get("deltaX/g3_5");   
        g4_x[5] = (TGraph*)infile->Get("deltaX/g4_5");
        g5_x[5] = (TGraph*)infile->Get("deltaX/g5_5");

        g1_y[5] = (TGraph*)infile->Get("deltaY/g1_5");
        g2_y[5] = (TGraph*)infile->Get("deltaY/g2_5");
        g3_y[5] = (TGraph*)infile->Get("deltaY/g3_5");   
        g4_y[5] = (TGraph*)infile->Get("deltaY/g4_5");
        g5_y[5] = (TGraph*)infile->Get("deltaY/g5_5");
        g6_y[5] = (TGraph*)infile->Get("deltaY/g6_5");
        
        g1_x[6] = (TGraph*)infile->Get("deltaX/g1_6");
        g2_x[6] = (TGraph*)infile->Get("deltaX/g2_6");
        g3_x[6] = (TGraph*)infile->Get("deltaX/g3_6");   
        g4_x[6] = (TGraph*)infile->Get("deltaX/g4_6");
        g5_x[6] = (TGraph*)infile->Get("deltaX/g5_6");              
      }

      infile->Close();
    }

    return;
  }

  //----------------------------------------------------------------------------
  /// Primary working method of service that provides position offsets to be
  /// used in ionization electron drift
  std::vector<double> LArVoxelSCE::GetPosOffsets(double xVal, double yVal, double zVal) const
  {
    std::vector<double> thePosOffsets;

    if(fRepresentationType == "Parametric")
      thePosOffsets = GetPosOffsetsParametric(xVal,yVal,zVal);
    else
      thePosOffsets.resize(3,0.0);

    return thePosOffsets;
  }

  //----------------------------------------------------------------------------
  /// Provides position offsets using a parametric representation
  std::vector<double> LArVoxelSCE::GetPosOffsetsParametric(double xVal, double yVal, double zVal) const
  {
    std::vector<double> thePosOffsetsParametric;

    double xValNew = TransformX(xVal);
    double yValNew = TransformY(yVal);
    double zValNew = TransformZ(zVal);

    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"X"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Y"));
    thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Z"));

    return thePosOffsetsParametric;
  }

  //----------------------------------------------------------------------------
  /// Provides one position offset using a parametric representation, for a given
  /// axis
  double LArVoxelSCE::GetOnePosOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
  {      
    double parA[6][7];
    double parB[6];
    
    for(int j = 0; j < 6; j++)
    {
      for(int i = 0; i < 7; i++)
        parA[j][i] = 0.0;
    
      parB[j] = 0.0;
    }
    
    if(axis == "X")
    {
      for(int j = 0; j < 7; j++)
      {
        parA[0][j] = g1_x[j]->Eval(zValNew);
        parA[1][j] = g2_x[j]->Eval(zValNew);
        parA[2][j] = g3_x[j]->Eval(zValNew);
        parA[3][j] = g4_x[j]->Eval(zValNew);
        parA[4][j] = g5_x[j]->Eval(zValNew);
      }
    
      f1_x->SetParameters(parA[0]);
      f2_x->SetParameters(parA[1]);
      f3_x->SetParameters(parA[2]);
      f4_x->SetParameters(parA[3]);
      f5_x->SetParameters(parA[4]);
    }
    else if(axis == "Y")
    {
      for(int j = 0; j < 6; j++)
      {
        parA[0][j] = g1_y[j]->Eval(zValNew);
        parA[1][j] = g2_y[j]->Eval(zValNew);
        parA[2][j] = g3_y[j]->Eval(zValNew);
        parA[3][j] = g4_y[j]->Eval(zValNew);
        parA[4][j] = g5_y[j]->Eval(zValNew);
        parA[5][j] = g6_y[j]->Eval(zValNew);
      }
    
      f1_y->SetParameters(parA[0]);
      f2_y->SetParameters(parA[1]);
      f3_y->SetParameters(parA[2]);
      f4_y->SetParameters(parA[3]);
      f5_y->SetParameters(parA[4]);
      f6_y->SetParameters(parA[5]);
    }
    else if(axis == "Z")
    {
      for(int j = 0; j < 5; j++)
      {
        parA[0][j] = g1_z[j]->Eval(zValNew);
        parA[1][j] = g2_z[j]->Eval(zValNew);
        parA[2][j] = g3_z[j]->Eval(zValNew);
        parA[3][j] = g4_z[j]->Eval(zValNew);
      }
    
      f1_z->SetParameters(parA[0]);
      f2_z->SetParameters(parA[1]);
      f3_z->SetParameters(parA[2]);
      f4_z->SetParameters(parA[3]);
    }
    
    double aValNew;
    double bValNew;
    
    if(axis == "Y")
    {
      aValNew = xValNew;
      bValNew = yValNew;
    }
    else
    {
      aValNew = yValNew;
      bValNew = xValNew;
    }
    
    double offsetValNew = 0.0;
    if(axis == "X")
    {
      parB[0] = f1_x->Eval(aValNew-1.25);
      parB[1] = f2_x->Eval(aValNew-1.25);
      parB[2] = f3_x->Eval(aValNew-1.25);
      parB[3] = f4_x->Eval(aValNew-1.25);
      parB[4] = f5_x->Eval(aValNew-1.25);
    
      fFinal_x->SetParameters(parB);
      offsetValNew = 100.0*fFinal_x->Eval(bValNew-1.25);
    }
    else if(axis == "Y")
    {
      parB[0] = f1_y->Eval(aValNew-1.25);
      parB[1] = f2_y->Eval(aValNew-1.25);
      parB[2] = f3_y->Eval(aValNew-1.25);
      parB[3] = f4_y->Eval(aValNew-1.25);
      parB[4] = f5_y->Eval(aValNew-1.25);
      parB[5] = f6_y->Eval(aValNew-1.25);
    
      fFinal_y->SetParameters(parB);
      offsetValNew = 100.0*fFinal_y->Eval(bValNew-1.25);
    }
    else if(axis == "Z")
    {
      parB[0] = f1_z->Eval(aValNew-1.25);
      parB[1] = f2_z->Eval(aValNew-1.25);
      parB[2] = f3_z->Eval(aValNew-1.25);
      parB[3] = f4_z->Eval(aValNew-1.25);
    
      fFinal_z->SetParameters(parB);
      offsetValNew = 100.0*fFinal_z->Eval(bValNew-1.25);
    }
    
    return offsetValNew;
  }

  //----------------------------------------------------------------------------
  /// Transform X to SCE X coordinate:  [2.56,0.0] --> [0.0,2.50]
  double LArVoxelSCE::TransformX(double xVal) const
  {
    double xValNew;
    xValNew = 2.50 - (2.50/2.56)*(xVal/100.0);

    return xValNew;
  }

  //----------------------------------------------------------------------------
  /// Transform Y to SCE Y coordinate:  [-1.165,1.165] --> [0.0,2.50]
  double LArVoxelSCE::TransformY(double yVal) const
  {
    double yValNew;
    yValNew = (2.50/2.33)*((yVal/100.0)+1.165);

    return yValNew;
  }

  //----------------------------------------------------------------------------
  /// Transform Z to SCE Z coordinate:  [0.0,10.37] --> [0.0,10.0]
  double LArVoxelSCE::TransformZ(double zVal) const
  {
    double zValNew;
    zValNew = (10.0/10.37)*(zVal/100.0);

    return zValNew;
  }

} // namespace sim
 
namespace sim {

  DEFINE_ART_SERVICE(LArVoxelSCE)


} // namespace sim
