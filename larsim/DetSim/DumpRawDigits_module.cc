/**
 * @file   DumpRawDigits_module.cc
 * @brief  Dumps on screen the content of the raw digits
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 14th, 2015
 */

// C//C++ standard libraries
#include <string>
#include <algorithm> // std::min(), std::copy_n()
#include <ios> // std::fixed
#include <iomanip> // std::setprecision(), std::setw()
#include <memory> // std::unique_ptr<>

// support libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// LArSoft includes
#include "larcore/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/SimpleTypesAndConstants/geo_types.h" // geo::View_t
#include "lardata/RawData/raw.h" // raw::Uncompress()
#include "lardata/RawData/RawDigit.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

namespace {
  
  // FIXME there should be an utility like this in lardata/Utilities soon
  /**
   * @brief Records the minimum and maximum of a value
   * @param T type of the values to record
   */
  template <typename T>
  class ExtremaAccumulator {
      public:
    using Data_t = T; ///< type of the extrema to be accumulated
    
    /// Default constructor: nothing accumulated yet
    ExtremaAccumulator(): minimum(), maximum(), n(0) {}
    
    ///@{
    ///@name Accessors
    Data_t samples() const { return n; }
    Data_t min() const { return hasStats()? minimum: T(); }
    Data_t max() const { return hasStats()? maximum: T(); }
    ///@}
    
    /// Returns whether this object has accumulated any statistics
    bool hasStats() const { return samples() > 0; }
    
    /// Returns whether all the accumulated objects have the same value
    bool isConstant() const { return !hasStats() || (min() == max()); }
    
    /**
     * @brief Adds a single value to the sample
     * @param value the value to be recorded
     * @return this object
     */
    ExtremaAccumulator& add(const Data_t& value)
      {
        if (n == 0) minimum = maximum = value;
        else {
          if (minimum > value) minimum = value;
          else if (maximum < value) maximum = value;
        }
        ++n;
        return *this;
      } // add(Data_t)
    
    /**
     * @brief Adds a sequence of values to the sample
     * @param ITER type of iterator
     * @param begin (constant) iterator pointing to the first element
     * @param end (constant) iterator pointing after the last element
     * @return this object
     */
    template <typename ITER>
    ExtremaAccumulator& add(ITER begin, ITER end)
      {
        for (ITER iter = begin; iter != end; ++iter) add(*iter);
        return *this;
      } // add(ITER)
    
      private:
    Data_t minimum, maximum; ///< extrema recorded so far
    unsigned int n; ///< number of recorded values
  }; // class ExtremaAccumulator<>
  
} // local namespace


namespace detsim {

  /**
   * @brief Prints the content of all the raw digits on screen
   *
   * This analyser prints the content of all the raw digits into the
   * LogInfo/LogVerbatim stream.
   * 
   * <b>Configuration parameters</b>
   * 
   * - <b>DetSimModuleLabel</b> (string, default: "daq"): label of the
   *   producer used to create the raw::RawDigits collection
   * - <b>OutputCategory</b> (string, default: "DumpDigits"): the category used
   *   for the output (useful for filtering)
   * - <b>DigitsPerLine</b> (integer, default: 20): the dump of digits and ticks
   *   will put this many of them for each line
   * - <b>IgnoreFilters</b> (boolean, default: false): if true, channel filters
   *   will be ignored; by default, only raw digits on channels that are not bad
   *   are printed out
   * - <b>Pedestal</b> (integer, default: 0): digit values are written relative
   *   to this number
   *
   */
  class DumpRawDigits: public art::EDAnalyzer {
      public:
    
    /// Default constructor
    explicit DumpRawDigits(fhicl::ParameterSet const& pset); 
    
    /// Does the printing
    void analyze (const art::Event& evt);
    
      private:
    
    std::string fDetSimModuleLabel; ///< name of module that produced the digits
    std::string fOutputCategory; ///< category for LogInfo output
    unsigned int fDigitsPerLine; ///< ticks/digits per line in the output
    short fPedestal; ///< ADC pedestal, will be subtracted from digits
    bool bIgnoreFilters; ///< use all the wires, don't filter them
    
  }; // class DumpRawDigits
  
} // namespace detsim


namespace detsim {

  //-------------------------------------------------
  DumpRawDigits::DumpRawDigits(fhicl::ParameterSet const& pset) 
    : EDAnalyzer         (pset)
    , fDetSimModuleLabel (pset.get<std::string>("DetSimModuleLabel", "daq"))
    , fOutputCategory    (pset.get<std::string>("OutputCategory", "DumpDigits"))
    , fDigitsPerLine     (pset.get<unsigned int>("DigitsPerLine", 20))
    , fPedestal          (pset.get<short>("Pedestal", 0))
    , bIgnoreFilters     (pset.get<bool>("IgnoreFilters", false))
    {}


  //-------------------------------------------------
  void DumpRawDigits::analyze(const art::Event& evt) {
    
    // fetch the data to be dumped on screen
    art::ValidHandle<std::vector<raw::RawDigit>> Digits
      = evt.getValidHandle<std::vector<raw::RawDigit>>(fDetSimModuleLabel);
    
    // channel filter: create one only if requested
    lariov::ChannelStatusProvider const* channelStatus = bIgnoreFilters
      ? nullptr
      : art::ServiceHandle<lariov::ChannelStatusService>()->GetProviderPtr();
    
    mf::LogInfo(fOutputCategory)
      << "The event contains " << Digits->size() << " raw digits";
    if (fPedestal != 0) {
      mf::LogVerbatim(fOutputCategory) << "A pedestal of " << fPedestal
        << " will be subtracted from all raw digits";
    } // if pedestal
    
    // a portable version of code dumping a wire is in
    // uboonecode/uboone/CalWireROI_module.cc
    
    std::vector<float> DigitBuffer(fDigitsPerLine), LastBuffer;
    for (const raw::RawDigit& digits: *Digits) {
      // uncompress the digits
      raw::RawDigit::ADCvector_t ADCs(digits.Samples());
      raw::Uncompress(digits.ADCs(), ADCs, digits.Compression());
      
      // print a header for the raw digits
      { // limit the scope of out:
        mf::LogVerbatim out(fOutputCategory);
        out << "  #" << digits.Channel() << ":";
        if (channelStatus && channelStatus->IsBad(digits.Channel())) {
          out << " bad channel";
          continue;
        }
        out << " " << ADCs.size() << " time ticks";
        if (digits.Samples() != ADCs.size())
          out << " [!!! EXPECTED " << digits.Samples() << "] ";
        out
          << " (" << digits.NADC() << " after compression); compression type: ";
        switch (digits.Compression()) {
          case raw::kNone:            out << "no compression"; break;
          case raw::kHuffman:         out << "Huffman encoding" ; break;
          case raw::kZeroSuppression: out << "zero suppression"; break;
          case raw::kZeroHuffman:     out << "zero suppression + Huffman encoding";
                                      break;
          case raw::kDynamicDec:      out << "dynamic decimation"; break;
          default:
            out << "unknown (#" << ((int) digits.Compression()) << ")"; break;
        } // switch
      } // block
      
      // print the content of the channel
      if (fDigitsPerLine > 0) {
        unsigned int repeat_count = 0; // additional lines like the last one
        unsigned int index = 0;
        LastBuffer.clear();
        ExtremaAccumulator<float> Extrema;
        mf::LogVerbatim(fOutputCategory) << "  content of the channel ("
          << fDigitsPerLine << " ticks per line):";
        auto iTick = ADCs.cbegin(), tend = ADCs.cend(); // const iterators
        while (iTick != tend) {
          // the next line will show at most fDigitsPerLine ticks
          unsigned int line_size
            = std::min(fDigitsPerLine, (unsigned int) ADCs.size() - index);
          if (line_size == 0) break; // no more ticks
          
          // fill the new buffer (iTick will move forward)
          DigitBuffer.resize(line_size);
          std::vector<float>::iterator iBuf = DigitBuffer.begin(),
            bend = DigitBuffer.end();
          while ((iBuf != bend) && (iTick != tend))
            Extrema.add(*(iBuf++) = *(iTick++) - fPedestal);
          index += line_size;
          
          // if the new buffer is the same as the old one, just mark it
          if (DigitBuffer == LastBuffer) {
            repeat_count += 1;
            continue;
          }
          
          // if there are previous repeats, write that on screen
          // before the new, different line
          if (repeat_count > 0) {
            mf::LogVerbatim(fOutputCategory)
              << "      [ ... repeated " << repeat_count << " more times ]";
            repeat_count = 0;
          }
          
          // dump the new line of ticks
          mf::LogVerbatim line_out(fOutputCategory);
          line_out << "   ";
          for (float digit: DigitBuffer)
            line_out << " " << std::setw(4) << digit;
          
          // quick way to assign DigitBuffer to LastBuffer
          // (we don't care we lose the former)
          std::swap(LastBuffer, DigitBuffer);
          
        } // while
        if (repeat_count > 0) {
          mf::LogVerbatim(fOutputCategory)
            << "      [ ... repeated " << repeat_count << " more times to the end]";
        }
        if (!Extrema.isConstant()) {
          mf::LogVerbatim(fOutputCategory)
            << "    range of " << Extrema.samples()
            << " samples: [" << Extrema.min() << ";" << Extrema.max() << "]";
        }
      } // if dumping the ticks
      
    } // for wire
    
  } // DumpRawDigits::analyze()

  DEFINE_ART_MODULE(DumpRawDigits)

} // namespace detsim
