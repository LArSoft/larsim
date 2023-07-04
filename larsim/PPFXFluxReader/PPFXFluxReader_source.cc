#include "PPFXFluxReader.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"

namespace fluxr {
  typedef art::Source<PPFXFluxReader> PPFXFluxReaderSource;
}

DEFINE_ART_INPUT_SOURCE(fluxr::PPFXFluxReaderSource)
