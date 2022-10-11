////////////////////////////////////////////////////////////////////////
//
// SimWire class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
// mwang@fnal.gov (2/04/21)
//
// - Revised field response functions with quadratic and sinusoidal for
//   collection and induction planes, respectively
// - Updated electronics response with uBooNE spice-based model
// - Included more realistic noise models:
//   (a) modified uBooNE model -> "ModUBooNE" in fcl
//   (b) ArgoNeuT data driven model -> "ArgoNeuT" in fcl
// - Included following distributions for fluctuating magnitude of
//   noise frequency components:
//   (a) simple modified Poisson -> "SimplePoisson" in fcl
//   (b) weighted Poisson from ArgoNeuT DDN -> "WeightedPoisson" in fcl
// - Included choice for generating unique noise in each channel for
//   each event
// - Updated ConvoluteResponseFunctions() to do things in same fashion
//   as in a typical SignalShapingXXX service for a particular XXX
//   experiment
//
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include "TComplex.h"
#include "TMath.h"

// C++ includes
#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandFlat.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/Simulation/SimChannel.h"

// Detector simulation of raw signals on wires
namespace detsim {

  class SimWire : public art::EDProducer {
  public:
    explicit SimWire(fhicl::ParameterSet const& pset);

  private:
    void produce(art::Event& evt) override;
    void beginJob() override;

    void ConvoluteResponseFunctions(); ///< convolute electronics and field response

    void SetFieldResponse(); ///< response of wires to field
    void SetElectResponse(); ///< response of electronics

    void GenNoise(std::vector<float>& array, CLHEP::HepRandomEngine& engine);

    std::string fDriftEModuleLabel; ///< module making the ionization electrons
    raw::Compress_t fCompression;   ///< compression type to use

    int fNTicks;                          ///< number of ticks of the clock
    double fSampleRate;                   ///< sampling rate in ns
    unsigned int fNSamplesReadout;        ///< number of ADC readout samples in 1 readout frame
    double fCol3DCorrection;              ///< correction factor to account for 3D path of
                                          ///< electrons thru wires
    double fInd3DCorrection;              ///< correction factor to account for 3D path of
                                          ///< electrons thru wires
    unsigned int fNElectResp;             ///< number of entries from response to use
    double fInputFieldRespSamplingPeriod; ///< Sampling period in the input field response.
    double fShapeTimeConst;               ///< time constants for exponential shaping
    double fADCPerPCAtLowestASICGain;     ///< ADCs/pC at lowest gain setting of 4.7 mV/fC
    double fASICGainInMVPerFC;            ///< actual gain setting used in mV/fC
    int fNoiseNchToSim;                   ///< number of noise channels to generate
    std::string fNoiseModelChoice;        ///< choice for noise model
    std::string fNoiseFluctChoice;        ///< choice for noise freq component mag fluctuations

    std::vector<int> fFieldRespTOffset;     ///< field response time offset in ticks
    std::vector<int> fCalibRespTOffset;     ///< calib response time offset in ticks
    std::vector<float> fColFieldParams;     ///< collection plane field function parameterization
    std::vector<float> fIndFieldParams;     ///< induction plane field function parameterization
    std::vector<double> fColFieldResponse;  ///< response function for the field @ collection plane
    std::vector<double> fIndFieldResponse;  ///< response function for the field @ induction plane
    std::vector<TComplex> fColShape;        ///< response function for the field @ collection plane
    std::vector<TComplex> fIndShape;        ///< response function for the field @ induction plane
    std::vector<double> fElectResponse;     ///< response function for the electronics
    std::vector<std::vector<float>> fNoise; ///< noise on each channel for each time
    std::vector<double> fNoiseModelPar;     ///< noise model params
    std::vector<double> fNoiseFluctPar;     ///< Poisson noise fluctuations params

    TH1D* fIndFieldResp; ///< response function for the field @ induction plane
    TH1D* fColFieldResp; ///< response function for the field @ collection plane
    TH1D* fElectResp;    ///< response function for the electronics
    TH1D* fColTimeShape; ///< convoluted shape for field x electronics @ col plane
    TH1D* fIndTimeShape; ///< convoluted shape for field x electronics @ ind plane
    TH1D* fNoiseDist;    ///< distribution of noise counts
    TF1* fNoiseFluct;    ///< Poisson dist for fluctuations in magnitude of noise freq components

    CLHEP::HepRandomEngine& fEngine; ///< Random-number engine owned by art
  };                                 // class SimWire

}

namespace detsim {

  //-------------------------------------------------
  SimWire::SimWire(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fDriftEModuleLabel{pset.get<std::string>("DriftEModuleLabel")}
    , fCompression{pset.get<std::string>("CompressionType") == "Huffman" ? raw::kHuffman :
                                                                           raw::kNone}
    , fCol3DCorrection{pset.get<double>("Col3DCorrection")}
    , fInd3DCorrection{pset.get<double>("Ind3DCorrection")}
    , fInputFieldRespSamplingPeriod{pset.get<double>("InputFieldRespSamplingPeriod")}
    , fShapeTimeConst{pset.get<double>("ShapeTimeConst")}
    , fADCPerPCAtLowestASICGain{pset.get<double>("ADCPerPCAtLowestASICGain")}
    , fASICGainInMVPerFC{pset.get<double>("ASICGainInMVPerFC")}
    , fNoiseNchToSim{pset.get<int>("NoiseNchToSim")}
    , fNoiseModelChoice{pset.get<std::string>("NoiseModelChoice")}
    , fNoiseFluctChoice{pset.get<std::string>("NoiseFluctChoice")}
    , fFieldRespTOffset{pset.get<std::vector<int>>("FieldRespTOffset")}
    , fCalibRespTOffset{pset.get<std::vector<int>>("CalibRespTOffset")}
    , fColFieldParams{pset.get<std::vector<float>>("ColFieldParams")}
    , fIndFieldParams{pset.get<std::vector<float>>("IndFieldParams")}
    , fNoiseModelPar{pset.get<std::vector<double>>("NoiseModelPar")}
    , fNoiseFluctPar{pset.get<std::vector<double>>("NoiseFluctPar")}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService> {}->createEngine(*this, pset, "Seed"))
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    fSampleRate = sampling_rate(clockData);
    fNSamplesReadout = detProp.NumberTimeSamples();

    MF_LOG_WARNING("SimWire") << "SimWire is an example module that works for the "
                              << "MicroBooNE detector.  Each experiment should implement "
                              << "its own version of this module to simulate electronics "
                              << "response.";

    produces<std::vector<raw::RawDigit>>();
  }

  //-------------------------------------------------
  void SimWire::beginJob()
  {
    // ... get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    fNoiseDist = tfs->make<TH1D>("Noise", ";Noise (ADC);", 1000, -10., 10.);

    art::ServiceHandle<util::LArFFT const> fFFT;
    fNTicks = fFFT->FFTSize();

    // ... Poisson dist function for fluctuating magnitude of noise frequency component
    if (fNoiseFluctChoice == "SimplePoisson") {
      // .. simple modified Poisson with (x-1)! in denominator
      double params[1];
      fNoiseFluct = new TF1("_poisson", "[0]**(x) * exp(-[0]) / ROOT::Math::tgamma(x)", 0, 5.);
      params[0] = fNoiseFluctPar[0]; // Poisson mean
      fNoiseFluct->SetParameters(params);
    }
    else if (fNoiseFluctChoice == "WeightedPoisson") {
      // .. weighted Poisson in ArgoNeuT DDN model
      double params[3];
      fNoiseFluct = new TF1(
        "_poisson", "[0]*pow([1]/[2], x/[2])*exp(-[1]/[2])/ROOT::Math::tgamma(x/[2]+1.)", 0, 5.);
      params[0] = fNoiseFluctPar[0];
      params[1] = fNoiseFluctPar[1];
      params[2] = fNoiseFluctPar[2];
      fNoiseFluct->SetParameters(params);
    }
    else {
      throw cet::exception("SimWire::beginJob")
        << fNoiseFluctChoice << " is an unknown noise fluctuation choice" << std::endl;
    }

    // ... generate the noise in advance depending on value of fNoiseNchToSim:
    //     positive - generate N=fNoiseNchToSim channels & randomly pick from pool when adding to signal
    //     zero     - no noise
    //     negative - generate unique noise for each channel for each event
    if (fNoiseNchToSim > 0) {
      if (fNoiseNchToSim > 10000) {
        throw cet::exception("SimWire::beginJob")
          << fNoiseNchToSim << " noise channels requested exceeds 10000" << std::endl;
      }
      fNoise.resize(fNoiseNchToSim);
      for (unsigned int p = 0; p < fNoise.size(); ++p) {
        GenNoise(fNoise[p], fEngine);
        for (int i = 0; i < fNTicks; ++i) {
          fNoiseDist->Fill(fNoise[p][i]);
        }
      }
    }

    // ... set field response and electronics response, then convolute them
    SetFieldResponse();
    SetElectResponse();
    ConvoluteResponseFunctions();
  }

  //-------------------------------------------------
  void SimWire::produce(art::Event& evt)
  {

    art::ServiceHandle<geo::Geometry const> geo;

    // ... generate unique noise for each channel in each event
    if (fNoiseNchToSim < 0) {
      fNoise.clear();
      fNoise.resize(geo->Nchannels());
      for (unsigned int p = 0; p < geo->Nchannels(); ++p) {
        GenNoise(fNoise[p], fEngine);
      }
    }

    // ... make a vector of const sim::SimChannel* that has same number
    //     of entries as the number of channels in the detector
    //     and set the entries for the channels that have signal on them
    //     using the chanHandle
    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel, chanHandle);

    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for (size_t c = 0; c < chanHandle.size(); ++c) {
      channels[chanHandle[c]->Channel()] = chanHandle[c];
    }

    // ... make an unique_ptr of sim::SimDigits that allows ownership of the produced
    //     digits to be transferred to the art::Event after the put statement below
    auto digcol = std::make_unique<std::vector<raw::RawDigit>>();

    art::ServiceHandle<util::LArFFT> fFFT;

    // ... Add all channels
    CLHEP::RandFlat flat(fEngine);

    std::map<int, double>::iterator mapIter;
    for (unsigned int chan = 0; chan < geo->Nchannels(); chan++) {

      std::vector<short> adcvec(fNTicks, 0);
      std::vector<double> fChargeWork(fNTicks, 0.);

      if (channels[chan]) {

        // .. get the sim::SimChannel for this channel
        const sim::SimChannel* sc = channels[chan];

        // .. loop over the tdcs and grab the number of electrons for each
        for (int t = 0; t < fNTicks; ++t)
          fChargeWork[t] = sc->Charge(t);

        int time_offset = 0;

        // .. Convolve charge with appropriate response function
        if (geo->SignalType(chan) == geo::kInduction) {
          fFFT->Convolute(fChargeWork, fIndShape);
          time_offset = fFieldRespTOffset[1] + fCalibRespTOffset[1];
        }
        else {
          fFFT->Convolute(fChargeWork, fColShape);
          time_offset = fFieldRespTOffset[0] + fCalibRespTOffset[0];
        }

        // .. Apply field response offset
        std::vector<int> temp;
        if (time_offset <= 0) {
          temp.assign(fChargeWork.begin(), fChargeWork.begin() - time_offset);
          fChargeWork.erase(fChargeWork.begin(), fChargeWork.begin() - time_offset);
          fChargeWork.insert(fChargeWork.end(), temp.begin(), temp.end());
        }
        else {
          temp.assign(fChargeWork.end() - time_offset, fChargeWork.end());
          fChargeWork.erase(fChargeWork.end() - time_offset, fChargeWork.end());
          fChargeWork.insert(fChargeWork.begin(), temp.begin(), temp.end());
        }
      }

      // ... Add noise to signal depending on value of fNoiseNchToSim
      if (fNoiseNchToSim != 0) {
        int noisechan = chan;
        if (fNoiseNchToSim > 0) {
          noisechan = TMath::Nint(flat.fire() * (1. * (fNoise.size() - 1) + 0.1));
        }
        for (int i = 0; i < fNTicks; ++i) {
          adcvec[i] = (short)TMath::Nint(fNoise[noisechan][i] + fChargeWork[i]);
        }
      }
      else {
        for (int i = 0; i < fNTicks; ++i) {
          adcvec[i] = (short)TMath::Nint(fChargeWork[i]);
        }
      }

      adcvec.resize(fNSamplesReadout);

      // ... compress the adc vector using the desired compression scheme,
      //     if raw::kNone is selected nothing happens to adcvec
      //     This shrinks adcvec, if fCompression is not kNone.
      raw::Compress(adcvec, fCompression);

      // ... add this digit to the collection
      digcol->emplace_back(chan, fNTicks, move(adcvec), fCompression);

    } //end loop over channels

    evt.put(std::move(digcol));

    return;
  }

  //-------------------------------------------------
  void SimWire::ConvoluteResponseFunctions()
  {
    double tick, ticks, peak;

    std::vector<double> col(fNTicks, 0.);
    std::vector<double> ind(fNTicks, 0.);
    std::vector<TComplex> kern(fNTicks / 2 + 1);
    std::vector<double> delta(fNTicks, 0.);

    art::ServiceHandle<util::LArFFT> fFFT;

    // ... do collection plane
    fColShape.resize(fNTicks / 2 + 1);
    fFFT->DoFFT(fElectResponse, fColShape);

    fFFT->DoFFT(fColFieldResponse, kern);
    for (unsigned int i = 0; i < kern.size(); ++i)
      fColShape[i] *= kern[i];

    fFFT->DoInvFFT(fColShape, col);

    delta[0] = 1.0;
    peak = fFFT->PeakCorrelation(delta, col);
    tick = 0.;
    ticks = tick - peak;
    fFFT->ShiftData(fColShape, ticks);
    fFFT->DoInvFFT(fColShape, col);

    // ... do induction plane
    fIndShape.resize(fNTicks / 2 + 1);
    fFFT->DoFFT(fElectResponse, fIndShape);

    kern.clear();
    kern.resize(fNTicks / 2 + 1);
    fFFT->DoFFT(fIndFieldResponse, kern);
    for (unsigned int i = 0; i < kern.size(); ++i)
      fIndShape[i] *= kern[i];

    fFFT->DoInvFFT(fIndShape, ind);

    delta.resize(0);
    delta.resize(fNTicks, 0);
    delta[0] = 1.0;
    peak = fFFT->PeakCorrelation(delta, ind);
    tick = 0.;
    ticks = tick - peak;
    fFFT->ShiftData(fIndShape, ticks);
    fFFT->DoInvFFT(fIndShape, ind);

    // ... write the time-domain shapes out to a file
    art::ServiceHandle<art::TFileService const> tfs;
    fColTimeShape = tfs->make<TH1D>(
      "ConvolutedCollection", ";ticks; Electronics#timesCollection", fNTicks, 0, fNTicks);
    fIndTimeShape = tfs->make<TH1D>(
      "ConvolutedInduction", ";ticks; Electronics#timesInduction", fNTicks, 0, fNTicks);

    // ... check that you did the right thing
    for (unsigned int i = 0; i < ind.size(); ++i) {
      fColTimeShape->Fill(i, col[i]);
      fIndTimeShape->Fill(i, ind[i]);
    }

    fColTimeShape->Write();
    fIndTimeShape->Write();
  }

  //-------------------------------------------------
  void SimWire::GenNoise(std::vector<float>& noise, CLHEP::HepRandomEngine& engine)
  {
    CLHEP::RandFlat flat(engine);

    noise.clear();
    noise.resize(fNTicks, 0.);
    std::vector<TComplex> noiseFrequency(fNTicks / 2 + 1, 0.); // noise in frequency space

    double pval = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // .. width of frequencyBin in kHz
    double binWidth = 1.0 / (fNTicks * fSampleRate * 1.0e-6);

    for (int i = 0; i < fNTicks / 2 + 1; ++i) {

      double x = (i + 0.5) * binWidth;

      if (fNoiseModelChoice == "Legacy") {
        // ... Legacy exponential model kept here for reference:
        //     par[0]=NoiseFact, par[1]=NoiseWidth, par[2]=LowCutoff, par[3-7]=0
        //     example parameter values for fcl: NoiseModelPar:[ 1.32e-1,120,7.5,0,0,0,0,0 ]
        pval = fNoiseModelPar[0] * exp(-(double)i * binWidth / fNoiseModelPar[1]);
        double lofilter = 1.0 / (1.0 + exp(-(i - fNoiseModelPar[2] / binWidth) / 0.5));
        flat.fireArray(1, rnd, 0, 1);
        pval *= lofilter * (0.9 + 0.2 * rnd[0]);
      }
      else if (fNoiseModelChoice == "ModUBooNE") {
        // ... Modified uBooNE model with additive exp to account for low freq region:
        //     example parameter values for fcl: NoiseModelPar:[
        //                                         4450.,-530.,280.,110.,
        //                                         -0.85,18.,0.064,74. ]
        pval = fNoiseModelPar[0] * exp(-0.5 * pow((x - fNoiseModelPar[1]) / fNoiseModelPar[2], 2)) *
                 exp(-0.5 * pow(x / fNoiseModelPar[3], fNoiseModelPar[4])) +
               fNoiseModelPar[5] + exp(-fNoiseModelPar[6] * (x - fNoiseModelPar[7]));
        double randomizer = fNoiseFluct->GetRandom();
        pval = pval * randomizer / fNTicks;
      }
      else if (fNoiseModelChoice == "ArgoNeuT") {
        // ... ArgoNeuT data driven model:
        //     In fcl set parameters to: NoiseModelPar:[
        //                                 5000,-5.52058e2,2.81587e2,-5.66561e1,
        //                                 4.10817e1,1.76284e1,1e-1,5.97838e1 ]
        pval = fNoiseModelPar[0] * exp(-0.5 * pow((x - fNoiseModelPar[1]) / fNoiseModelPar[2], 2)) *
                 ((fNoiseModelPar[3] / (x + fNoiseModelPar[4])) + 1) +
               fNoiseModelPar[5] + exp(-fNoiseModelPar[6] * (x - fNoiseModelPar[7]));
        double randomizer = fNoiseFluct->GetRandom();
        pval = pval * randomizer / fNTicks;
      }
      else {
        throw cet::exception("SimWire::GenNoise")
          << fNoiseModelChoice << " is an unknown choice for the noise model" << std::endl;
      }

      flat.fireArray(1, rnd, 0, 1);
      phase = rnd[0] * 2. * TMath::Pi();

      TComplex tc(pval * cos(phase), pval * sin(phase));
      noiseFrequency[i] += tc;
    }

    // .. inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    // .. multiply each noise value by fNTicks as the InvFFT
    //    divides each bin by fNTicks assuming that a forward FFT
    //    has already been done.
    for (unsigned int i = 0; i < noise.size(); ++i)
      noise[i] *= 1. * fNTicks;
  }

  //-------------------------------------------------
  void SimWire::SetFieldResponse()
  {

    // ... Files to write the response functions to
    art::ServiceHandle<art::TFileService const> tfs;
    fIndFieldResp =
      tfs->make<TH1D>("InductionFieldResponse", ";t (ns);Induction Response", fNTicks, 0, fNTicks);
    fColFieldResp = tfs->make<TH1D>(
      "CollectionFieldResponse", ";t (ns);Collection Response", fNTicks, 0, fNTicks);

    fColFieldResponse.resize(fNTicks, 0.);
    fIndFieldResponse.resize(fNTicks, 0.);

    // ... First set response for collection plane
    int nbinc = fColFieldParams[0];

    double integral = 0.;
    for (int i = 1; i < nbinc; ++i) {
      fColFieldResponse[i] = i * i;
      integral += fColFieldResponse[i];
    }

    for (int i = 0; i < nbinc; ++i) {
      fColFieldResponse[i] *= fColFieldParams[1] / integral;
      fColFieldResp->Fill(i, fColFieldResponse[i]);
    }

    // ... Now set response for induction plane
    int nbini = fIndFieldParams[0];
    unsigned short lastbini = 2 * nbini;

    integral = 0;
    for (unsigned short i = 0; i < lastbini; ++i) {
      double ang = i * TMath::Pi() / nbini;
      fIndFieldResponse[i] = sin(ang);
      if (fIndFieldResponse[i] > 0) { integral += fIndFieldResponse[i]; }
      else {
        fIndFieldResponse[i] *=
          fIndFieldParams[2]; // scale the negative lobe by 10% (from ArgoNeuT)
      }
    }
    ++lastbini;

    for (unsigned short i = 0; i < lastbini; ++i) {
      fIndFieldResponse[i] *= fIndFieldParams[1] / integral;
      fIndFieldResp->Fill(i, fIndFieldResponse[i]);
    }

    // ... Save the field responses
    fColFieldResp->Write();
    fIndFieldResp->Write();
  }

  //-------------------------------------------------
  void SimWire::SetElectResponse()
  {
    fElectResponse.resize(fNTicks, 0.);
    std::vector<double> time(fNTicks, 0.);

    // ... Gain and shaping time variables from fcl file:
    double Ao = 1.0;
    double To = fShapeTimeConst; //peaking time

    // ... this is actually sampling time, in ns
    mf::LogInfo("SimWire::SetElectResponse")
      << "Check sampling intervals: " << fInputFieldRespSamplingPeriod << " ns"
      << "Check number of samples: " << fNTicks;

    // ... The following sets the microboone electronics response function in
    //     time-space. Function comes from BNL SPICE simulation of DUNE35t
    //     electronics. SPICE gives the electronics transfer function in
    //     frequency-space. The inverse laplace transform of that function
    //     (in time-space) was calculated in Mathematica and is what is being
    //     used below. Parameters Ao and To are cumulative gain/timing parameters
    //     from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain.
    //     They have been adjusted to make the SPICE simulation to match the
    //     actual electronics response. Default params are Ao=1.4, To=0.5us.
    double max = 0;

    for (size_t i = 0; i < fElectResponse.size(); ++i) {

      // ... convert time to microseconds, to match fElectResponse[i] definition
      time[i] = (1. * i) * fInputFieldRespSamplingPeriod * 1e-3;
      fElectResponse[i] =
        4.31054 * exp(-2.94809 * time[i] / To) * Ao -
        2.6202 * exp(-2.82833 * time[i] / To) * cos(1.19361 * time[i] / To) * Ao -
        2.6202 * exp(-2.82833 * time[i] / To) * cos(1.19361 * time[i] / To) *
          cos(2.38722 * time[i] / To) * Ao +
        0.464924 * exp(-2.40318 * time[i] / To) * cos(2.5928 * time[i] / To) * Ao +
        0.464924 * exp(-2.40318 * time[i] / To) * cos(2.5928 * time[i] / To) *
          cos(5.18561 * time[i] / To) * Ao +
        0.762456 * exp(-2.82833 * time[i] / To) * sin(1.19361 * time[i] / To) * Ao -
        0.762456 * exp(-2.82833 * time[i] / To) * cos(2.38722 * time[i] / To) *
          sin(1.19361 * time[i] / To) * Ao +
        0.762456 * exp(-2.82833 * time[i] / To) * cos(1.19361 * time[i] / To) *
          sin(2.38722 * time[i] / To) * Ao -
        2.6202 * exp(-2.82833 * time[i] / To) * sin(1.19361 * time[i] / To) *
          sin(2.38722 * time[i] / To) * Ao -
        0.327684 * exp(-2.40318 * time[i] / To) * sin(2.5928 * time[i] / To) * Ao +
        +0.327684 * exp(-2.40318 * time[i] / To) * cos(5.18561 * time[i] / To) *
          sin(2.5928 * time[i] / To) * Ao -
        0.327684 * exp(-2.40318 * time[i] / To) * cos(2.5928 * time[i] / To) *
          sin(5.18561 * time[i] / To) * Ao +
        0.464924 * exp(-2.40318 * time[i] / To) * sin(2.5928 * time[i] / To) *
          sin(5.18561 * time[i] / To) * Ao;

      if (fElectResponse[i] > max) max = fElectResponse[i];

    } // end loop over time buckets

    // ... "normalize" fElectResponse[i], before the convolution

    for (auto& element : fElectResponse) {
      element /= max;
      element *= fADCPerPCAtLowestASICGain * 1.60217657e-7;
      element *= fASICGainInMVPerFC / 4.7; // relative to lowest gain setting of 4.7 mV/fC
    }

    fNElectResp = fElectResponse.size();

    // ... write the response out to a file

    art::ServiceHandle<art::TFileService const> tfs;
    fElectResp = tfs->make<TH1D>(
      "ElectronicsResponse", ";t (ns);Electronics Response", fNElectResp, 0, fNElectResp);
    for (unsigned int i = 0; i < fNElectResp; ++i) {
      fElectResp->Fill(i, fElectResponse[i]);
    }

    fElectResp->Write();
  }

}

DEFINE_ART_MODULE(detsim::SimWire)
