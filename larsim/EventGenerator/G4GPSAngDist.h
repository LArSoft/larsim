/**
 * @author      : Daniele Guffanti (University and INFN Milano-Bicocca)
 * @file        : G4GPSAngDist.h
 * @created     : Friday Oct 16, 2026 15:34:16 CEST
 * @brief Direction-generation tool wrapping Geant4's GPS angular distribution.
 *
 * Configures a G4SPSAngDistribution from FHiCL and samples directions using
 * the CLHEP random engine supplied by the caller.
 */

#ifndef LARSIM_G4GPSANGDIST_H

#define LARSIM_G4GPSANGDIST_H

#include "DirectionGenerator.h"
#include "art/Utilities/ToolConfigTable.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geant4/G4RandomTools.hh"
#include "Geant4/G4SPSAngDistribution.hh"
#include "Geant4/G4SPSPosDistribution.hh"
#include "Geant4/G4SPSRandomGenerator.hh"

#include "TFile.h"
#include "TH1.h"

namespace evgen {
  namespace {
    /**
     * @brief Temporarily installs a CLHEP engine as the Geant4 random engine.
     *
     * Restores the previous engine when the guard is destroyed.
     */
    class G4EngineGuard {
    public:
      explicit G4EngineGuard(CLHEP::HepRandomEngine& engine) : fPrevious{G4Random::getTheEngine()}
      {
        G4Random::setTheEngine(&engine);
      }
      ~G4EngineGuard() { G4Random::setTheEngine(fPrevious); }
      G4EngineGuard(G4EngineGuard const&) = delete;
      G4EngineGuard& operator=(G4EngineGuard const&) = delete;

    private:
      CLHEP::HepRandomEngine* fPrevious = {};
    };
  }

  /**
   * @brief Samples particle directions using a Geant4 GPS angular distribution.
   *
   * Supported distribution types are `iso`, `cos`, `planar`, `beam1d`,
   * `beam2d`, and `user`. Angular limits and rotation axes are configured
   * through FHiCL. For `user`, theta and/or phi histograms can be supplied
   * as arrays or read from ROOT histograms.
   *
   * @note Histogram bin contents represent probability per interval, not
   * probability density. The tool removes weight outside the configured
   * angular range, splitting bins at the range boundaries if necessary.
   */
  class G4GPSAngDist : public IDirectionGenerator {
  public:
    /// FHiCL parameters for the angular distribution.
    struct Config {
      fhicl::Atom<std::string> type{
        fhicl::Name{"distType"},
        fhicl::Comment{"Type of angular distribution to use in G4GPSAngDistribution\n"
                       "Options: iso [default], cos, planar, beam1d, beam2d, user"},
        "iso"};
      fhicl::Sequence<double, 3> rot1{
        fhicl::Name{"rot1"},
        fhicl::Comment{
          "Defines the first (x' direction) rotation vector AR1 [default (1,0,0)]\n"
          "for the angular distribution and is not necessarily a unit vector.\n"
          "Used with /gps/ang/rot2 to compute the angular distribution rotation matrix."},
        {1, 0, 0}};
      fhicl::Sequence<double, 3> rot2{
        fhicl::Name{"rot2"},
        fhicl::Comment{
          "Defines the second rotation vector AR2 [default (0,1,0)]\n"
          "for the angular distribution and is not necessarily a unit vector.\n"
          "Used with /gps/ang/rot1 to compute the angular distribution rotation matrix."},
        {0, 1, 0}};
      fhicl::Atom<double> minTheta{
        fhicl::Name{"minTheta"},
        fhicl::Comment{
          "Sets a minimum value [default 0 rad] for the θ distribution. Units are rad."},
        0.};
      fhicl::Atom<double> maxTheta{
        fhicl::Name{"maxTheta"},
        fhicl::Comment{
          "Sets a maximum value [default π rad] for the θ distribution. Units are rad."},
        CLHEP::pi};
      fhicl::Atom<double> minPhi{
        fhicl::Name{"minPhi"},
        fhicl::Comment{
          "Sets a minimum value [default 0 rad] for the φ distribution. Units are rad."},
        0.};
      fhicl::Atom<double> maxPhi{
        fhicl::Name{"maxPhi"},
        fhicl::Comment{
          "Sets a maximum value [default 2π rad] for the φ distribution. Units are rad."},
        2 * CLHEP::pi};
      fhicl::Atom<double> sigma_r{
        fhicl::Name{"sigma_r"},
        fhicl::Comment{"Sets the standard deviation [default 0 rad] of beam directional profile in "
                       "radial. Units are rad."},
        0.};
      fhicl::Atom<double> sigma_x{fhicl::Name{"sigma_x"},
                                  fhicl::Comment{"Sets the standard deviation [default 0 rad] of "
                                                 "beam directional profile in x. Units are rad."},
                                  0.};
      fhicl::Atom<double> sigma_y{fhicl::Name{"sigma_y"},
                                  fhicl::Comment{"Sets the standard deviation [default 0 rad] of "
                                                 "beam directional profile in y. Units are rad."},
                                  0.};
      fhicl::Atom<int> verbosity{
        fhicl::Name{"verbosity"},
        fhicl::Comment{"Sets the verbosity level [default 0] for the G4GPSAngDist tool."},
        0};

      /**
         * @brief Configuration for a user-defined angular histogram.
         *
         * Supply either `x` and `y`, or `file` and `name`. For arrays, `x`
         * contains the lower edge followed by upper bin edges, in radians;
         * `y` contains interval weights and must start with zero. ROOT files
         * are searched for in `FW_SEARCH_PATH`.
         */
      struct UserHistConfig {
        fhicl::OptionalSequence<double> x{
          fhicl::Name("x"),
          fhicl::Comment("Points: lower edge first, then upper bin edges, ascending [rad]")};
        fhicl::OptionalSequence<double> y{
          fhicl::Name("y"),
          fhicl::Comment("Probability per interval; y[0] must be 0; same length as x")};
        fhicl::OptionalAtom<std::string> file{
          fhicl::Name("file"),
          fhicl::Comment("ROOT file, searched in FW_SEARCH_PATH")};
        fhicl::OptionalAtom<std::string> name{
          fhicl::Name("name"),
          fhicl::Comment("TH1 name; bin contents are the probability per bin")};
      };

      fhicl::OptionalTable<UserHistConfig> thetaHist{
        fhicl::Name("thetaHist"),
        fhicl::Comment("User-defined theta distribution (for \"user\" type); must provide either x "
                       "and y or file and name")};
      fhicl::OptionalTable<UserHistConfig> phiHist{
        fhicl::Name("phiHist"),
        fhicl::Comment("User-defined phi distribution (for \"user\" type); must provide either x "
                       "and y or file and name")};
    };

    /// Histogram points as (angle in radians, interval weight) pairs.
    using HistData_t = std::vector<std::pair<double, double>>;

    using Parameters = art::ToolConfigTable<Config>;

    /**
       * @brief Configures the Geant4 angular distribution.
       * @param config Validated FHiCL tool configuration.
       * @throws art::Exception If configuration values or histograms are invalid.
       */
    explicit G4GPSAngDist(Parameters const& config);

    /**
       * @brief Samples a direction using the supplied random engine.
       * @param engine CLHEP engine temporarily installed for the Geant4 call.
       * @return Sampled direction vector.
       */
    geo::Vector_t generate(CLHEP::HepRandomEngine& engine) override;

  protected:
    /// Writes the configured distribution type, angular limits, and axes.
    void dumpSettingsImpl(std::ostream& out) const override;

  private:
    geo::Vector_t fAngRef1 = {};
    geo::Vector_t fAngRef2 = {};

    std::unique_ptr<G4SPSAngDistribution> fAngDist;
    std::unique_ptr<G4SPSRandomGenerator> fRandGen;
    std::unique_ptr<G4SPSPosDistribution> fPosDist;

    /**
       * @brief Reads histogram points from FHiCL arrays or a ROOT histogram.
       * @param config Histogram source configuration.
       * @param label Histogram name used in diagnostic messages.
       * @return Histogram edges and interval weights.
       */
    static HistData_t readHist(const Config::UserHistConfig& config, const std::string& label);

    /// Copies ROOT histogram bin edges and contents into histogram points.
    static HistData_t fromTH1(const TH1& h);

    /**
       * @brief Validates and restricts histogram weights to an angular range.
       * @param pts Histogram points, modified in place.
       * @param label Histogram name used in diagnostic messages.
       * @param min Lower angular limit in radians.
       * @param max Upper angular limit in radians.
       * @throws art::Exception If the histogram fails validation.
       */
    static void validate_user_hist(HistData_t& pts,
                                   const std::string& label,
                                   double min,
                                   double max);
  };

  G4GPSAngDist::G4GPSAngDist(Parameters const& config)
    : fAngDist{std::make_unique<G4SPSAngDistribution>()}
    , fRandGen{std::make_unique<G4SPSRandomGenerator>()}
    , fPosDist{std::make_unique<G4SPSPosDistribution>()}
  {
    // Set verbosity
    fAngDist->SetVerbosity(config().verbosity());

    // Disable user-defined angular distribution and rotation axes
    fAngDist->SetUserWRTSurface(false);
    fAngDist->SetUseUserAngAxis(false);

    // Validate angular ranges
    if (config().minTheta() < 0. || config().maxTheta() > CLHEP::pi ||
        config().minTheta() >= config().maxTheta()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid theta range [" << config().minTheta() << ", "
        << config().maxTheta() << "] rad. Must satisfy 0 <= minTheta < maxTheta <= pi.";
    }
    if (config().minPhi() < 0. || config().maxPhi() > 2 * CLHEP::pi ||
        config().minPhi() >= config().maxPhi()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid phi range [" << config().minPhi() << ", " << config().maxPhi()
        << "] rad. Must satisfy 0 <= minPhi < maxPhi <= 2pi.";
    }
    fAngDist->SetMinTheta(config().minTheta());
    fAngDist->SetMaxTheta(config().maxTheta());
    fAngDist->SetMinPhi(config().minPhi());
    fAngDist->SetMaxPhi(config().maxPhi());

    // Validate type
    static const std::array<std::string, 6> validTypes{
      "iso", "cos", "planar", "beam1d", "beam2d", "user"};
    const std::string type = config().type();
    if (std::find(validTypes.begin(), validTypes.end(), type) == validTypes.end()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid 'type' \"" << type << "\".\n"
        << "Valid options: iso, cos, planar, beam1d, beam2d, user.";
    }
    fAngDist->SetAngDistType(type);

    // Validate rotation vectors
    // Rotation vectors must be non-zero
    geo::Vector_t v1{config().rot1()[0], config().rot1()[1], config().rot1()[2]};
    geo::Vector_t v2{config().rot2()[0], config().rot2()[1], config().rot2()[2]};
    if (v1.Mag2() < 1e-8 || v2.Mag2() < 1e-8) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: 'rot1' and 'rot2' must be non-zero vectors.";
    }
    // Rotation vectors must not be parallel
    if (std::abs(v1.Cross(v2).Mag2()) < 1e-8) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: 'rot1' and 'rot2' must not be parallel vectors.";
    }
    fAngRef1.SetCoordinates(config().rot1()[0], config().rot1()[1], config().rot1()[2]);
    fAngDist->DefineAngRefAxes(
      "angref1", G4ThreeVector(config().rot1()[0], config().rot1()[1], config().rot1()[2]));
    fAngRef2.SetCoordinates(config().rot2()[0], config().rot2()[1], config().rot2()[2]);
    fAngDist->DefineAngRefAxes(
      "angref2", G4ThreeVector(config().rot2()[0], config().rot2()[1], config().rot2()[2]));

    fAngDist->SetBeamSigmaInAngR(config().sigma_r());
    fAngDist->SetBeamSigmaInAngX(config().sigma_x());
    fAngDist->SetBeamSigmaInAngY(config().sigma_y());

    // Set the user histograms if provided
    if (type == "user") {
      mf::LogInfo("G4GPSAngDist") << "Using user-defined angular distribution"
                                  << "thetaHist present: " << config().thetaHist().has_value()
                                  << ", phiHist present: " << config().phiHist().has_value();

      Config::UserHistConfig th_config, phi_config;
      bool const hasTheta = config().thetaHist(th_config);
      bool const hasPhi = config().phiHist(phi_config);
      if (!hasTheta && !hasPhi)
        throw art::Exception(art::errors::Configuration)
          << "distType user needs thetaHist and/or phiHist\n";
      if (hasTheta) {
        auto pts = readHist(th_config, "thetaHist");
        validate_user_hist(pts, "thetaHist", config().minTheta(), config().maxTheta());
        for (auto const& [x, w] : pts)
          fAngDist->UserDefAngTheta(G4ThreeVector{x, w, 0.});
      }
      if (hasPhi) {
        auto pts = readHist(phi_config, "phiHist");
        validate_user_hist(pts, "phiHist", config().minPhi(), config().maxPhi());
        for (auto const& [x, w] : pts)
          fAngDist->UserDefAngPhi(G4ThreeVector{x, w, 0.});
      }
    }

    // Set the random generator for pos and ang distributions to use the same engine
    fPosDist->SetBiasRndm(fRandGen.get());
    fAngDist->SetBiasRndm(fRandGen.get());
    fAngDist->SetPosDistribution(fPosDist.get());
  }

  geo::Vector_t G4GPSAngDist::generate(CLHEP::HepRandomEngine& engine)
  {
    G4EngineGuard guard{engine};
    G4ParticleMomentum dir = fAngDist->GenerateOne();
    return {dir.x(), dir.y(), dir.z()};
  }

  G4GPSAngDist::HistData_t G4GPSAngDist::readHist(const Config::UserHistConfig& config,
                                                  const std::string& label)
  {
    std::vector<double> xs, ys;
    std::string file, name;
    bool const hasXY = config.x(xs) && config.y(ys);
    bool const hasRoot = config.file(file) && config.name(name);
    if (hasXY == hasRoot)
      throw art::Exception(art::errors::Configuration)
        << label << ": give either x and y, or file and name\n";
    if (hasXY) {
      if (xs.size() != ys.size())
        throw art::Exception(art::errors::Configuration) << label << ": x and y differ in length\n";
      G4GPSAngDist::HistData_t pts;
      for (std::size_t i = 0; i < xs.size(); ++i)
        pts.emplace_back(xs[i], ys[i]);
      return pts;
    }
    cet::search_path search_path{"FW_SEARCH_PATH"};
    std::string path;
    if (!search_path.find_file(file, path))
      throw art::Exception(art::errors::Configuration)
        << label << ": cannot find " << file << " in FW_SEARCH_PATH\n";
    std::unique_ptr<TFile> f{TFile::Open(path.c_str(), "READ")};
    auto* h = f ? f->Get<TH1>(name.c_str()) : nullptr;
    if (!h)
      throw art::Exception(art::errors::Configuration)
        << label << ": histogram " << name << " not found in " << path << '\n';
    return fromTH1(*h); // values copied; the file can close afterwards
  }

  G4GPSAngDist::HistData_t G4GPSAngDist::fromTH1(const TH1& h)
  {
    G4GPSAngDist::HistData_t pts;
    pts.emplace_back(h.GetXaxis()->GetBinLowEdge(1), 0.);
    for (int i = 1; i <= h.GetNbinsX(); ++i)
      pts.emplace_back(h.GetXaxis()->GetBinUpEdge(i), h.GetBinContent(i));
    return pts;
  }

  void G4GPSAngDist::validate_user_hist(G4GPSAngDist::HistData_t& data,
                                        std::string const& label,
                                        double min,
                                        double max)
  {
    auto bad = [&](char const* why) {
      throw art::Exception(art::errors::Configuration) << label << ": " << why << '\n';
    };
    if (data.size() < 2 || data.size() > 1024) bad("need 2..1024 points");
    if (data.front().second != 0.) bad("y[0] must be 0 (it would be a point mass)");

    G4GPSAngDist::HistData_t adjusted_data;
    adjusted_data.reserve(data.size() + 2);
    adjusted_data.push_back(data.front());

    bool support = false;
    bool removedWeight = false;

    for (std::size_t i = 1; i < data.size(); ++i) {
      const double& lower = data[i - 1].first;
      const double& upper = data[i].first;
      const double& weight = data[i].second;

      if (!std::isfinite(upper) || !(upper > lower))
        bad("x values must be finite and strictly increasing");
      if (!std::isfinite(weight) || weight < 0.) bad("y values must be finite and non-negative");

      double left = lower;
      auto append = [&](double right) {
        const bool is_inside = (left >= min && right <= max);
        const double fraction = is_inside ? (right - left) / (upper - lower) : 0.;
        adjusted_data.emplace_back(right, weight * fraction);
        if (is_inside && fraction > 0.) support = true;
        if (!is_inside && weight > 0.) removedWeight = true;
        left = right;
      };

      if (min > lower && min < upper) append(min);
      if (max > lower && max < upper) append(max);
      append(upper);
    }

    if (!support) bad("no support in the specified range");
    if (adjusted_data.size() < 2) bad("not enough points in the specified range");
    if (adjusted_data.size() > 1024) bad("too many points in the specified range");

    data.swap(adjusted_data);

    if (removedWeight) {
      mf::LogWarning("G4GPSAngDist")
        << label << ": set histogram weight outside [" << min << ", " << max << "] to zero.";
    }
  }

  void G4GPSAngDist::dumpSettingsImpl(std::ostream& out) const
  {
    out << "G4GPSAngDist settings:\n"
        << "  type: " << fAngDist->GetDistType() << "\n"
        << "  minTheta: " << fAngDist->GetMinTheta() << " rad\n"
        << "  maxTheta: " << fAngDist->GetMaxTheta() << " rad\n"
        << "  minPhi: " << fAngDist->GetMinPhi() << " rad\n"
        << "  maxPhi: " << fAngDist->GetMaxPhi() << " rad\n"
        << "  rot1: (" << fAngRef1.x() << ", " << fAngRef1.y() << ", " << fAngRef1.z() << ")\n"
        << "  rot2: (" << fAngRef2.x() << ", " << fAngRef2.y() << ", " << fAngRef2.z() << ")\n";
  }
} // namespace evgen

#endif /* end of include guard LARSIM_G4GPSANGDIST_H */
