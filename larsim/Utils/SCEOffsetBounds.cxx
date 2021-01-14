/**
 * @file larsim/Utils/SCEOffsetBounds.cxx
 *
 * @brief Implementation of SCE Offset Bounds Checking method(s)
 *
 * @author Tom Junk (trj@fnal.gov)
 *
 * $Log: $
 */

// LArSoft
#include "larsim/Utils/SCEOffsetBounds.h"

bool larsim::Utils::SCE::out_of_bounds(geo::Vector_t const &offset)
{
  constexpr double limit{1.e9};
  return (std::abs(offset.X()) > limit ||
          std::abs(offset.Y()) > limit ||
          std::abs(offset.Z()) > limit);
}
