/**
 * @file larsim/Utils/SCEOffsetBounds.h
 *
 * @brief Utility function for testing if Space Charge offsets are out of bounds
 *
 * @author Tom Junk (trj@fnal.gov)
 *
 * $log: $
 */
#ifndef LARSIMSCEOFFSETBOUNDS_H_SEEN
#define LARSIMSCEOFFSETBOUNDS_H_SEEN

//LArSoft
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

namespace larsim {
  namespace Utils {
    namespace SCE {
      bool out_of_bounds(geo::Vector_t const& offset);
    }
  }
}

#endif // #ifndef LARSIMSCEOFFSETBOUNDS_H_SEEN
