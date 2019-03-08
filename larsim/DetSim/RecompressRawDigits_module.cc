/**
 * @file   RecompressRawDigits_module.cc
 * @brief  Module saving raw::RawDigit with a different compression algorithm.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 14, 2016
 *
 * Provides:
 * 
 * * `raw::RecompressRawDigits` module
 * 
 */

// LArSoft libraries
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h" // raw::Compress(), raw::Uncompress()
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::Compress_t

// framework libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" // art::ValidHandle
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <algorithm> // std::transform()
#include <memory> // std::make_unique()
#include <cctype> // std::toupper()


namespace raw {
  
  /**
   * @brief Writes the input `raw::RawDigit` with a different compression
   * 
   * This module writes a collection of `raw::RawDigit` in a new data product,
   * using the specified compression mode.
   * 
   * 
   * Input
   * ------
   * 
   * A single collection of `raw::RawDigit` objects.
   * 
   * 
   * Output
   * -------
   * 
   * A single collection of `raw::RawDigit` objects.
   * 
   * 
   * Configuration
   * --------------
   * 
   * * *rawDigitLabel* (input tag, _mandatory_): the input tag for the original
   *   `raw::RawDigit` collection
   * * *compressionType* (string, _mandatory_): the compression mode code
   *   (use names as in the enumerator `raw::Compress_t`)
   * * *instanceName* (string, optional): if specified, the output collection is
   *   saved with the specified product instance name (by default, none is used)
   * 
   */
  class RecompressRawDigits: public art::EDProducer {
      public:
    
    struct Config {
      
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> rawDigitLabel{
        Name("rawDigitLabel"),
        Comment("input tag for the original raw::RawDigit collection")
        };
      
      fhicl::Atom<std::string> compressionType{
        Name("compressionType"),
        Comment("compression mode code (e.g. \"kHuffman\")")
        };
      
      fhicl::Atom<std::string> instanceName{
        Name("instanceName"),
        Comment("instance name for the output data product (none by default)"),
        std::string() // default
        };
      
    }; // Config
    
    
    using Parameters = art::EDProducer::Table<Config>;
    
    
    /// Constructor; see module documentation for configuration directions
    explicit RecompressRawDigits(Parameters const& config)
      : EDProducer{config}
      , fRawDigitLabel(config().rawDigitLabel())
      , fCompressionType(parseCompressionType(config().compressionType()))
      , fInstanceName(config().instanceName())
      {
        produces<std::vector<raw::RawDigit>>(fInstanceName);
      }
    
    
    virtual void produce(art::Event& event) override;
    
    
    /**
     * @brief Returns a `RawDigit` with its waveform compressed in `newFormat`
     * @param digit `RawDigit` to be recompressed
     * @param newFormat compression format for the returned `RawDigit`
     * @param force run even if input and target compression formats match
     * @return `RawDigit` with the waveform in a new compression format
     * 
     * A copy of `RawDigit` is always returned, with the waveform stored in the
     * `newFormat` compression format.
     * If the input `RawDigit` is already in that format, the waveform is simply
     * copied, unless `force` is specified `true`, in which case data is
     * uncompressed and compressed back.
     * 
     */
    static raw::RawDigit recompress
      (raw::RawDigit const& digit, Compress_t newFormat, bool force = false);

    
    /**
     * @brief Returns the compression mode corresponding to the specified string
     * @param spec
     * @return the compression mode spec describes
     * @throw art::Exception (type: art::errors::Configuration) if invalid spec
     * 
     * This method returns the compression type described by the specification.
     * Specification is not case sensitive and can omit the trailing letter `k`.
     * Therefore, for example to specify `kNone` all the following
     * specifications are valid: `"kNone"`, `"kNONE"`, `"none"`.
     * 
     */
    static raw::Compress_t parseCompressionType(std::string spec);
    
    
      private:
    art::InputTag fRawDigitLabel; ///< tag of input raw digit collection
    raw::Compress_t fCompressionType; ///< type of compression to be applied
    std::string fInstanceName; ///< instance name of output data product
    
    
  }; // class RecompressRawDigits
  
  
} // namespace raw




//------------------------------------------------------------------------------
//---  implementation
//------------------------------------------------------------------------------
namespace {
  std::string toupper(std::string s) {
    // act directly on s (which is a copy of the original);
    // explicit cast is needed to pick which std::toupper() to use
    std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) std::toupper);
    return s;
  } // raw::RecompressRawDigits::toupper()
} // local namespace


//------------------------------------------------------------------------------
void raw::RecompressRawDigits::produce(art::Event& event) {
  
  //
  // prepare input and output container
  //
  auto oldRawDigitHandle
    = event.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  
  auto newRawDigits = std::make_unique<std::vector<raw::RawDigit>>();
  
  //
  // create a new compressed raw digit for each existing one
  //
  std::set<raw::Compress_t> formats;
  for (raw::RawDigit const& oldDigits: *oldRawDigitHandle) {
    
    // keep track of all the formats
    formats.insert(oldDigits.Compression());
    
    // recompress the digit, and put it into the future data product
    newRawDigits->emplace_back(recompress(oldDigits, fCompressionType));
    
  } // for oldDigits
  
  //
  // store the new digits into the event
  //
  event.put(std::move(newRawDigits), fInstanceName);
  
} // raw::RecompressRawDigits::produce()



//------------------------------------------------------------------------------
raw::RawDigit raw::RecompressRawDigits::recompress
  (raw::RawDigit const& digit, Compress_t newFormat, bool force /* = false */)
{
  
  if ((newFormat == digit.Compression()) && !force) 
    return digit; // return a copy
  
  // uncompress the data
  raw::RawDigit::ADCvector_t ADCs(digit.Samples());
  raw::Uncompress(digit.ADCs(), ADCs, digit.Compression());    
  
  // compress the data (it happens in place)
  raw::Compress(ADCs, newFormat);
  
  // we force a copy in a buffer sized just big enough to host the data;
  // in case the buffer is already of the right size (as it happens when
  // newFormat specifies no compression), then no copy happens
  ADCs.shrink_to_fit();
  
  // create a new digit with the new buffer
  return raw::RawDigit{
    digit.Channel(),
    digit.Samples(),
    std::move(ADCs),
    newFormat
    };
  
} // raw::RecompressRawDigits::recompress()


//------------------------------------------------------------------------------
raw::Compress_t raw::RecompressRawDigits::parseCompressionType(std::string spec)
{
  std::string SPEC(toupper(spec));
  
  if ((SPEC == "NONE"           ) || (SPEC == "KNONE"           )) return raw::kNone           ;
  if ((SPEC == "HUFFMAN"        ) || (SPEC == "KHUFFMAN"        )) return raw::kHuffman        ;
  if ((SPEC == "ZEROSUPPRESSION") || (SPEC == "KZEROSUPPRESSION")) return raw::kZeroSuppression;
  if ((SPEC == "ZEROHUFFMAN"    ) || (SPEC == "KZEROHUFFMAN"    )) return raw::kZeroHuffman    ;
  if ((SPEC == "DYNAMICDEC"     ) || (SPEC == "KDYNAMICDEC"     )) return raw::kDynamicDec     ;
  
  // if a valid one triggers this error, it needs to be added above
  throw art::Exception(art::errors::Configuration)
    << "Unrecognized compression type: '" << spec << "'\n";
  
} // raw::RecompressRawDigits::parseCompressionType()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(raw::RecompressRawDigits)



//------------------------------------------------------------------------------
