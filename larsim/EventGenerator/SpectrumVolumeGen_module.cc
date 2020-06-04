/// \file  SpectrumVolumeGen_module.cc
/// \brief Generate anything in any volume of the geometry
///        Handy for quick radiological studies
/// \author  plasorak@FNAL.GOV
///          April 2020 PLasorak

#include "larsim/EventGenerator/BaseRadioGen.h"
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

namespace evgen {
  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes

  class SpectrumVolumeGen : public evgen::BaseRadioGen {
  public:
    explicit SpectrumVolumeGen(fhicl::ParameterSet const& pset);

  private:
    // This is called for each event.
    void produce_radio(art::Event& evt);
    int m_pdg;
    double m_mass;
    TH1D* m_spectrum;
  };

}

namespace evgen{

  SpectrumVolumeGen::SpectrumVolumeGen(fhicl::ParameterSet const& pset):
    BaseRadioGen(pset) {
    std::string isotope="";
    isotope = pset.get<std::string>("isotope");
    m_pdg = pset.get<int>("isotope");

    bool mass_specified = pset.get_if_present<double>("mass", m_mass);
    if (not mass_specified) {
      const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
      const TParticlePDG* definition = databasePDG->GetParticle(m_pdg);
      // Check that the particle is known to ROOT.  If not, this is
      // not a major error; Geant4 has an internal particle coding
      // scheme for nuclei that ROOT doesn't recognize.
      if (definition != 0){
        m_mass = definition->Mass();
      } else {
        throw cet::exception("SpectrumVolumeGen") << "Cannot find the particle pdg " << m_pdg;
      }
    }

    m_isotope.push_back(isotope);

    double min_p=0, max_p=0;
    std::vector<double> bins;
    bool minmax_mode = pset.get_if_present<double>("spectrum_p_min", min_p);
    if (minmax_mode) {
      max_p = pset.get<double>("spectrum_p_max");
    } else {
      bins = pset.get<std::vector<double>>("bins");
    }

    std::string function;
    std::vector<double> parameters;
    bool have_params;
    std::vector<double> spectrum;
    bool func_mode = pset.get_if_present<std::string>("function", function);

    if (func_mode) {
      have_params = pset.get_if_present<std::vector<double>>("parameters", parameters);
    } else {
      spectrum = pset.get<std::vector<double>>("spectrum");
    }

    if (func_mode and minmax_mode) {
      int nbins = pset.get<int>("nbins");
      TF1 the_func("spect_func",function.c_str(), min_p, max_p);
      if (have_params) the_func.SetParameters(parameters.data());

      m_spectrum = new TH1D("spectrum", ";p [GeV];PDF", nbins, min_p, max_p);

      for (int ibin=1; ibin<nbins+1; ++ibin) {
        double x = m_spectrum->GetBinCenter(ibin);
        double y = the_func.Eval(x);
        if (y<0) y=0;
        m_spectrum->SetBinContent(ibin, y);
      }

    } else if (func_mode and not minmax_mode) {
      if (bins.size()<=1) {
        throw cet::exception("SpectrumVolumeGen") << "You need to specify more than 1 bin edge to use the \"bins\" option.";
      }
      MF_LOG_WARNING("SpectrumVolumeGen") << "\033[31mYou are using a function and non constant binning...\n"
                                          << "Please please please check spectrum and momentum distribution are what you actually want!\033[0m\n";
      int nbins = bins.size()-1;

      TF1 the_func("spect_func",function.c_str(), *bins.begin(), *bins.rbegin());
      if (have_params) the_func.SetParameters(parameters.data());

      m_spectrum = new TH1D("spectrum", ";p [GeV];PDF", nbins, bins.data());

      for (int ibin=1; ibin<nbins+1; ++ibin) {
        double x = m_spectrum->GetBinCenter(ibin);
        double y = the_func.Eval(x);
        if (y<0) y=0;
        m_spectrum->SetBinContent(ibin, y);
      }

    } else if (not func_mode and minmax_mode) {
      int nbins = spectrum.size();

      m_spectrum = new TH1D("spectrum", ";p [GeV];PDF", nbins, min_p, max_p);
      for (int ibin=1; ibin<nbins+1; ++ibin) {
        m_spectrum->SetBinContent(ibin, spectrum.at(ibin-1));
      }

    } else {
      int nbins = spectrum.size();
      if ((unsigned)nbins != spectrum.size())
        throw cet::exception("SpectrumVolumeGen") << "spectrum and bins don't have coherent size!";

      m_spectrum = new TH1D("spectrum", ";p [GeV];PDF", nbins, bins.data());

      for (int ibin=1; ibin<nbins+1; ++ibin) {
        m_spectrum->SetBinContent(ibin, spectrum.at(ibin-1));
      }
    }

    m_spectrum->Scale(1. / m_spectrum->Integral());

    art::ServiceHandle<art::TFileService> tfs;
    TH1D* sp = tfs->make<TH1D>(*m_spectrum);
    (void)sp;

  }

  //____________________________________________________________________________
  void SpectrumVolumeGen::produce_radio(art::Event& evt) {

    //unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    int track_id=-1;
    const std::string primary_str("primary");

    int n_decay = GetNDecays();

    for (int iDecay=0; iDecay<n_decay; ++iDecay) {
      TLorentzVector position;

      if (GetGoodPositionTime(position)) {
        simb::MCParticle part(track_id, m_pdg, primary_str);

        TLorentzVector momentum = dirCalc(m_spectrum->GetRandom(), m_mass);
        part.AddTrajectoryPoint(position, momentum);

        truth.Add(part);
        FillHistos(part);

        track_id--;
      } // GetGoodPosition
    } // idecay


    MF_LOG_DEBUG("SpectrumVolumeGen") << truth;
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  }

//Calculate an arbitrary direction with a given magnitude p

}//end namespace evgen

DEFINE_ART_MODULE(evgen::SpectrumVolumeGen)
