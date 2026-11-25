#define BOOST_TEST_MODULE MergeTrackIDs
#include "boost/test/unit_test.hpp"

#include "lardataobj/Simulation/sim.h"
#include "larsim/MergeSimSources/MergeSimSources.h"

#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace {
  using Ancestry = std::map<int, std::set<int>>;

  Ancestry mergeAncestry(Ancestry const& input, int offset)
  {
    sim::MergeSimSourcesUtility merger{{offset}};
    std::vector<sim::ParticleAncestryMap> output;
    merger.MergeParticleAncestryMaps(output, sim::ParticleAncestryMap{input}, 0);
    BOOST_REQUIRE_EQUAL(output.size(), 1U);
    return output.front().GetMap();
  }

  sim::SimEnergyDeposit mergeDeposit(int assigned, int original, int offset)
  {
    sim::MergeSimSourcesUtility merger{{offset}};
    std::vector<sim::SimEnergyDeposit> const input{
      {3, 4, 0.5, 0.25, {1., 2., 3.}, {4., 5., 6.}, 7., 8., assigned, 11, original}};
    std::vector<sim::SimEnergyDeposit> output;
    merger.MergeSimEnergyDeposits(output, input, 0);
    BOOST_REQUIRE_EQUAL(output.size(), 1U);
    return output.front();
  }
}

BOOST_AUTO_TEST_CASE(AncestrySentinelsAndSignedIDs)
{
  int const missing = sim::NoParticleId;
  Ancestry const input{{missing, {2, missing}}, {7, {3, -4}}, {-8, {5}}};
  for (int const offset : {0, 10000000, 20000000}) {
    Ancestry const expected{{missing, {2 + offset, missing}},
                            {7 + offset, {3 + offset, -4 - offset}},
                            {-8 - offset, {5 + offset}}};
    BOOST_CHECK(mergeAncestry(input, offset) == expected);
  }
}

BOOST_AUTO_TEST_CASE(DepositSentinelsAndSignedIDs)
{
  int const missing = sim::NoParticleId;
  for (int const offset : {0, 10000000, 20000000}) {
    auto const noAssigned = mergeDeposit(missing, 12, offset);
    BOOST_CHECK_EQUAL(noAssigned.TrackID(), missing);
    BOOST_CHECK_EQUAL(noAssigned.OrigTrackID(), 12 + offset);

    auto const noOriginal = mergeDeposit(-7, missing, offset);
    BOOST_CHECK_EQUAL(noOriginal.TrackID(), -7 - offset);
    BOOST_CHECK_EQUAL(noOriginal.OrigTrackID(), missing);

    auto const bothMissing = mergeDeposit(missing, missing, offset);
    BOOST_CHECK_EQUAL(bothMissing.TrackID(), missing);
    BOOST_CHECK_EQUAL(bothMissing.OrigTrackID(), missing);
    BOOST_CHECK_EQUAL(bothMissing.Energy(), 0.25);
    BOOST_CHECK_EQUAL(bothMissing.NumElectrons(), 4);
    BOOST_CHECK_EQUAL(bothMissing.T0(), 7.);
    BOOST_CHECK_EQUAL(bothMissing.T1(), 8.);
  }
}

BOOST_AUTO_TEST_CASE(RejectOverflowAndSentinelCollisions)
{
  int const high = std::numeric_limits<int>::max();
  int const low = std::numeric_limits<int>::min() + 1;
  BOOST_CHECK(mergeAncestry({{high - 1, {low + 1}}}, 1) == Ancestry({{high, {low}}}));
  BOOST_CHECK_THROW(mergeAncestry({{high, {1}}}, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeAncestry({{1, {high}}}, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeAncestry({{low, {1}}}, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeAncestry({{1, {low}}}, 2), std::overflow_error);
  BOOST_CHECK_THROW(mergeDeposit(high, 1, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeDeposit(1, high, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeDeposit(low, 1, 1), std::overflow_error);
  BOOST_CHECK_THROW(mergeDeposit(1, low, 2), std::overflow_error);
}
