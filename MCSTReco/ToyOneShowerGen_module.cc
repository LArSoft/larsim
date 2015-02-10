////////////////////////////////////////////////////////////////////////
// Class:       ToyOneShowerGen
// Module Type: producer
// File:        ToyOneShowerGen_module.cc
//
// Generated at Mon Aug 11 14:14:59 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

#include "CLHEP/Random/RandFlat.h"

#include "TDatabasePDG.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TF1.h"

#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/sim.h"

class ToyOneShowerGen;

class ToyOneShowerGen : public art::EDProducer {
public:
  explicit ToyOneShowerGen(fhicl::ParameterSet const & p);
  virtual ~ToyOneShowerGen();

  void produce(art::Event & e) override;
  void beginRun(art::Run & run) override;
  std::vector<double> GetXYZDirection(double uz=0);
  std::vector<double> GetXYZPosition();

private:

  TF1 *fShapeMomentum;
  TF1 *fShapeTheta;
  CLHEP::RandFlat *fFlatRandom;
  size_t fNEvents;
  double fMass;
  double fTime;
  int    fPDGCode;
};


ToyOneShowerGen::ToyOneShowerGen(fhicl::ParameterSet const & p)
  : fShapeMomentum(nullptr), fShapeTheta(nullptr), fFlatRandom(nullptr)
{
  gRandom->SetSeed(0);
  TRandom3 r(0);
  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();

  // Time
  fTime = p.get<double>("Time");

  // PDGCode
  fPDGCode = p.get<int>("PDGCode");

  // Mass
  fMass = TDatabasePDG().GetParticle(fPDGCode)->Mass();
  if(fPDGCode != 22 && !fMass) throw std::exception();
  
  //
  // Momentum distribution
  //
  
  // Read-in formula
  fShapeMomentum = new TF1("fShapeMomentum",
			   (p.get<std::string>("MomentumShapeFormula")).c_str(),
			   p.get<double>("MomentumLowerBound"),
			   p.get<double>("MomentumUpperBound"));
  // Check formula
  if(!fShapeMomentum)

    throw cet::exception(__PRETTY_FUNCTION__) 
      << "Failed making momentum spectrum shape from provided formula!" << std::endl;

  // Read-in parameter values
  std::vector<double> shape_par_values = p.get<std::vector<double> >("MomentumShapeParameters");

  // Check parameter values
  if((int)(shape_par_values.size()) != fShapeMomentum->GetNpar())

    throw cet::exception(__PRETTY_FUNCTION__)
      << "Number of parameters provided to MomentumShapeFormula does not match with the formula!" << std::endl;

  // Set parameter values
  for(size_t i=0; i<shape_par_values.size(); ++i)

    fShapeMomentum->SetParameter(i,shape_par_values.at(i));

  //
  // Angular distribution
  //
  
  // Read-in formula
  fShapeTheta = new TF1("fShapeTheta",
			(p.get<std::string>("ThetaShapeFormula")).c_str(),
			p.get<double>("ThetaLowerBound"),
			p.get<double>("ThetaUpperBound"));
  // Check formula
  if(!fShapeTheta)

    throw cet::exception(__PRETTY_FUNCTION__) 
      << "Failed making momentum spectrum shape from provided formula!" << std::endl;

  // Read-in parameter values
  shape_par_values = p.get<std::vector<double> >("ThetaShapeParameters");

  // Check parameter values
  if((int)(shape_par_values.size()) != fShapeTheta->GetNpar())

    throw cet::exception(__PRETTY_FUNCTION__)
      << "Number of parameters provided to ThetaShapeFormula does not match with the formula!" << std::endl;

  // Set parameter values
  for(size_t i=0; i<shape_par_values.size(); ++i)

    fShapeTheta->SetParameter(i,shape_par_values.at(i));

  //
  // Random engine initialization
  //
  createEngine(sim::GetRandomNumberSeed());
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  fFlatRandom = new CLHEP::RandFlat(engine);
}

//------------------------------------------------------------------------------
void ToyOneShowerGen::beginRun(art::Run& run)
{
  // grab the geometry object to see what geometry we are using
  art::ServiceHandle<geo::Geometry> geo;

  std::unique_ptr<sumdata::RunData> runData(new sumdata::RunData(geo->DetectorName()));

  run.put(std::move(runData));

  return;
}

ToyOneShowerGen::~ToyOneShowerGen()
{
  delete fShapeTheta;
  delete fShapeMomentum;
  delete fFlatRandom;
}

std::vector<double> ToyOneShowerGen::GetXYZPosition() {

  std::vector<double> pos(3,0);

  //pos.at(0) = fFlatRandom->fire(170.,2390.);
  //pos.at(1) = fFlatRandom->fire(-995.,995.);
  //pos.at(2) = fFlatRandom->fire(170.,10190.);

  //pos.at(0) = fFlatRandom->fire(17.,239.);
  //pos.at(1) = fFlatRandom->fire(-99.5,99.5);
  //pos.at(2) = fFlatRandom->fire(17.,1019.);

  art::ServiceHandle<geo::Geometry> geo;

  pos.at(0) = fFlatRandom->fire(0.,2.*(geo->DetHalfWidth()));
  pos.at(1) = fFlatRandom->fire(-1.*(geo->DetHalfHeight()), geo->DetHalfHeight());
  pos.at(2) = fFlatRandom->fire(0.,geo->DetLength());

  return pos;
}

std::vector<double> ToyOneShowerGen::GetXYZDirection(double uz) {

  std::vector<double> udir(3,0);

  double phi = fFlatRandom->fire(0, 2 * 3.141592653589793238);

  udir[0] = TMath::Sin(TMath::ACos(uz)) * TMath::Cos(phi);
  udir[1] = TMath::Sin(TMath::ACos(uz)) * TMath::Sin(phi);
  udir[2] = uz;

  return udir;
}

void ToyOneShowerGen::produce(art::Event & e)
{
  std::unique_ptr< std::vector<simb::MCTruth> > mctArray(new std::vector<simb::MCTruth>);
  double mom = fShapeMomentum->GetRandom();
  double Uz   = TMath::Cos(fShapeTheta->GetRandom());

  std::vector<double> pos = GetXYZPosition();

  std::vector<double> dir = GetXYZDirection(Uz);

  simb::MCTruth truth;

  TLorentzVector pos_lorentz(pos.at(0), pos.at(1), pos.at(2), fTime);
  TLorentzVector mom_lorentz( dir.at(0) * mom,
			      dir.at(1) * mom,
			      dir.at(2) * mom,
			      sqrt(pow(mom,2)+pow(fMass,2)));

  simb::MCParticle part(0, fPDGCode, "primary", 0, fMass, 1);

  part.AddTrajectoryPoint(pos_lorentz, mom_lorentz);

  truth.Add(part);

  mctArray->push_back(truth);
  
  e.put(std::move(mctArray));

}

DEFINE_ART_MODULE(ToyOneShowerGen)
