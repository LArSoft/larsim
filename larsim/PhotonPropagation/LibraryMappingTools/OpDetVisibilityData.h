/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/OpDetVisibilityData.h
 * @brief  Data container for use with photon visibility libraries.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 29, 2019
 * 
 * This is a header-only library.
 * 
 */

#ifndef LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_OPDETVISIBILITYDATA_H
#define LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_OPDETVISIBILITYDATA_H

// LArSoft libraries
#include "lardataalg/Utilities/MappedContainer.h"

// C++ standard libraries
#include <iterator> // std::size()
#include <utility> // std::forward()
#include <type_traits> // std::enable_if_t<>


namespace phot {
  
  /**
   * @brief Trait for determining the validity of library data in a container.
   * @tparam Cont type of container being tested
   * @tparam Enable dummy type for `std::enable_if_t` hooks
   * 
   * This class is expected to be specialized to provide specific behavior for
   * custom classes.
   * The required interface consists of the single static method:
   * * `static bool isValid(Cont const&)` returning whether `cont` is valid
   * 
   * 
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T, std::size_t N>
   * struct MyArray;
   * 
   * template <typename T, std::size_t N>
   * struct LibraryDataValidatorStruct<MyArray<T, N>, void> {
   *   static bool isValid(Cont const&) { return N > 0U; }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * or
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct MyVector;
   * 
   * template <typename T>
   * struct LibraryDataValidatorStruct
   *   <T, std::enable_if_t<util::is_instance_of<MyVector, T>>>
   * {
   *   static bool isValid(Cont const& cont) { return !cont.isInvalid(); }
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Cont, typename Enable = void>
  struct LibraryDataValidatorStruct;


  /**
   * @brief Returns the validity of content of library data in the container.
   * @tparam Cont type of the container being tested
   * @param cont a reference to the container being tested
   * @return whether the content of library data in `cont` is valid.
   * 
   * Customization of this behaviour should be through specialization of
   * `phot::LibraryDataValidatorStruct`.
   * 
   */
  template <typename Cont>
  bool isValidLibraryData(Cont&& cont);
  
  
  /**
   * @brief A container for photon visibility mapping data.
   * 
   * This container is designed to replace the use of plain pointers to library
   * data, adding an intermediate mapping layer.
   * 
   * A proper mapping allows to use symmetries of the detector to reduce the
   * redundant size of the library, at the cost of some run time.
   * 
   * This container is used as a `util::MappedContainer` object, with the
   * addition of being aware of whether the contained data is valid or not.
   * The validity is estimated by considering the original data that is being
   * mapped, that is the data from the library. If the data is invalid or empty,
   * as it may be for a null pointer or an empty vector, the container data is
   * also considered to be invalid.
   * 
   * 
   * Implementation notes
   * ---------------------
   * 
   * In the current implementation, the concept of data being invalid is
   * implemented by forcing the nominal size of the container to be zero.
   * 
   * 
   * Known limitations
   * ------------------
   * 
   * Currently the mapping of the optical detectors is forced to be a static
   * collection of detector IDs. If more flexibility is needed, the design of
   * the container internals need to be updated to allow run time polymorphism.
   * 
   * Given the type of validity detection described above, it takes extra steps
   * to provide a default mapping for areas not covered by the library. While
   * with plain `util::MappedContainer` that can be achieved by a empty data
   * container and a mapping returning `InvalidIndex` for all optical detectors,
   * here this would just yield to an invalid container. Instead, together with
   * that mapping, some fake data needs to be supplied to
   * `phot::OpDetVisibilityData`, which will be basically ignored and may be of
   * any size larger than zero (or more precisely, 
   * `phot::isValidLibraryData()` must return `false`).
   * 
   */
  template <typename Cont, typename Mapping>
  class OpDetVisibilityData
    : public util::MappedContainer<Cont, Mapping>
  {
    using ContainerBase_t = util::MappedContainer<Cont, Mapping>;
    
      public:
    
    // import types from base class
    using typename ContainerBase_t::DataContainer_t;
    using typename ContainerBase_t::Mapping_t;
    using typename ContainerBase_t::size_type;
    using typename ContainerBase_t::value_type;
    
      private:
    
    /// Determines whether the current content should be considered valid.
    static size_type effectiveSize
      (DataContainer_t const& cont, size_type allegedSize);
    
    /// Determines whether the current content should be considered valid.
    static size_type effectiveSize
      (DataContainer_t const& cont, Mapping_t const& mapping);
    
      public:
    
    
    // --- BEGIN Constructors --------------------------------------------------
    
    /// @name Validity information
    /// @{
    /// Default constructor: data will be invalid.
    OpDetVisibilityData() = default;
    
    /**
     * @brief Constructor: acquires data, mapping and a default value.
     * @param cont container with the data to be mapped
     * @param mapping the mapping to be used
     * @param size the size of the container after mapping
     * @param defValue value to be used as default
     * @see `MappedContainer(DataContainer_t const&, Mapping_t const&, size_type, value_type)`
     * 
     * This works as the equivalent `util::MappedContainer` constructor, except
     * that if the data in `cont` is deemed invalid the size is overridden to
     * zero.
     */
    OpDetVisibilityData(
      DataContainer_t const& cont,
      Mapping_t const& mapping,
      size_type size,
      value_type defValue
      )
      : ContainerBase_t(cont, mapping, effectiveSize(cont, size), defValue)
      {}
    
    /**
     * @brief Constructor: acquires data and mapping.
     * @param cont container with the data to be mapped
     * @param mapping the mapping to be used
     * @param size the size of the container after mapping
     * @see `MappedContainer(DataContainer_t const&, Mapping_t const&, size_type)`
     * 
     * This works as the equivalent `util::MappedContainer` constructor, except
     * that if the data in `cont` is deemed invalid the size is overridden to
     * zero.
     */
    OpDetVisibilityData
      (DataContainer_t const& cont, Mapping_t const& mapping, size_type size)
      : ContainerBase_t(cont, mapping, effectiveSize(cont, size))
      {}
    
    /**
     * @brief Constructor: acquires data and mapping.
     * @param cont container with the data to be mapped
     * @param mapping the mapping to be used
     * @see `util::MappedContainer(DataContainer_t const&, Mapping_t const&, size_type)`
     * 
     * The size of the container is declared to be the minimal one
     * (see `minimal_size()`), unless data of the library `cont` is invalid,
     * in which case the full container is invalid.
     * 
     * The default value is assigned as in
     * `util::MappedContainer(DataContainer_t const&, Mapping_t const&, size_type)`.
     */
    OpDetVisibilityData(DataContainer_t const& cont, Mapping_t const& mapping)
      : ContainerBase_t(cont, mapping, effectiveSize(cont, mapping))
      {}
    
    
    // --- END Constructors ----------------------------------------------------
    
    
    // --- BEGIN Validity information ------------------------------------------
    /// @name Validity information
    /// @{
    
    /// Returns whether the contained data is valid or not.
    bool isValid() const { return !ContainerBase_t::empty(); }
    
    /// Returns whether the container has valid data.
    operator bool() const { return isValid(); }
    
    /// Returns whether the container has no valid data.
    bool operator!() const { return !isValid(); }
    
    /// @}
    // --- END Validity information --------------------------------------------
    
    
  }; // OpDetVisibilityData
  
  
} // namespace phot


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------

namespace phot {
  
  namespace details {
  
    //--------------------------------------------------------------------------
    template <typename Cont>
    auto generic_size(Cont&& cont)
      { using std::size; return size(std::forward<Cont>(cont)); }
    
    
    //--------------------------------------------------------------------------
    // layered implementation of `LibraryDataValidatorStructStandardImpl`:
    
    // This is the last layer: if compiler points here,
    // the type `Cont` is not being supported yet:
    // - std::size() (or better, `generic_size()`)
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_size;
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_size
      <Cont, std::enable_if_t<util::always_true_v<decltype(
        generic_size(std::declval<Cont>())
      )>>
      >
    {
      static bool isValid(Cont const& cont) { return generic_size(cont) > 0U; }
    }; // struct LibraryDataValidatorStructStandardImpl_size
    
    
    // - Cont::empty()
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_empty
      : LibraryDataValidatorStructStandardImpl_size<Cont>
      {};
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_empty
      <Cont, std::enable_if_t<util::always_true_v<decltype(
        std::declval<Cont const>().empty()
      )>>
      >
    {
      static bool isValid(Cont const& cont) { return !cont.empty(); }
    }; // struct LibraryDataValidatorStructStandardImpl_empty
    
    
    // - Cont::is_valid()
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_is_valid
      : LibraryDataValidatorStructStandardImpl_empty<Cont>
      {};
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_is_valid
      <Cont, std::enable_if_t<util::always_true_v<decltype(
        std::declval<Cont const>().is_valid()
      )>>
      >
    {
      static bool isValid(Cont const& cont) { return cont.is_valid(); }
    }; // struct LibraryDataValidatorStructStandardImpl_is_valid
    
    
    // - Cont::isValid()
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_isValid
      : LibraryDataValidatorStructStandardImpl_is_valid<Cont>
      {};
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_isValid
      <Cont, std::enable_if_t<util::always_true_v<decltype(
        std::declval<Cont const>().isValid()
      )>>
      >
    {
      static bool isValid(Cont const& cont) { return cont.isValid(); }
    }; // struct LibraryDataValidatorStructStandardImpl_isValid
    
    
    // - std::unique_ptr<> (using `enable_if` to catch all qualifiers)
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_unique_ptr
      : LibraryDataValidatorStructStandardImpl_isValid<Cont>
      {};
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_unique_ptr
      <Cont, std::enable_if_t<util::is_unique_ptr_v<Cont>>>
    {
      static bool isValid(Cont const& cont) { return bool(cont); }
    }; // struct LibraryDataValidatorStructStandardImpl_unique_ptr
    
    
    // C pointer types:
    template <typename T>
    struct LibraryDataValidatorStructStandardImpl_pointer_branch;
    
    // - T*
    template <typename T>
    struct LibraryDataValidatorStructStandardImpl_pointer_branch<T*>
      { static bool isValid(T* ptr) { return bool(ptr); } };
    
    // - T[]
    template <typename T>
    struct LibraryDataValidatorStructStandardImpl_pointer_branch<T[]>
      : LibraryDataValidatorStructStandardImpl_size<T*> {};
    
    // - T[N]
    template <typename T, std::size_t N>
    struct LibraryDataValidatorStructStandardImpl_pointer_branch<T[N]>
      { static bool isValid(T (&) [N]) { return N > 0U; } };
    
    // - entry point
    template <typename Cont, typename = void>
    struct LibraryDataValidatorStructStandardImpl_pointer
      : public LibraryDataValidatorStructStandardImpl_unique_ptr<Cont>
    {};
    
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl_pointer
      <Cont, std::enable_if_t<std::is_pointer_v<std::decay_t<Cont>>>>
      : LibraryDataValidatorStructStandardImpl_pointer_branch
        <std::decay_t<Cont>>
    {};
    
    
    // - entry point
    template <typename Cont>
    struct LibraryDataValidatorStructStandardImpl
      : public LibraryDataValidatorStructStandardImpl_pointer<Cont>
    {};
    
    
    //--------------------------------------------------------------------------
    
  } // namespace details


  //----------------------------------------------------------------------------
  // we pick a standard implementation of our own, and let users add here
  template <typename Cont, typename /* = void */>
  struct LibraryDataValidatorStruct
    : public details::LibraryDataValidatorStructStandardImpl
     <std::remove_reference_t<Cont>>
  {};

    
  //----------------------------------------------------------------------------
  template <typename Cont>
  bool isValidLibraryData(Cont&& cont) {
    return
      LibraryDataValidatorStruct<std::remove_reference_t<Cont>>::isValid(cont);
  } // isValidLibraryData()


  //----------------------------------------------------------------------------
  template <typename Cont, typename Mapping>
  auto OpDetVisibilityData<Cont, Mapping>::effectiveSize 
    (DataContainer_t const& cont, size_type allegedSize) -> size_type 
    { return isValidLibraryData(cont)? allegedSize: 0U; }


  template <typename Cont, typename Mapping>
  auto OpDetVisibilityData<Cont, Mapping>::effectiveSize 
    (DataContainer_t const& cont, Mapping_t const& mapping) -> size_type
  { 
    return isValidLibraryData(cont)
      ? ContainerBase_t::minimal_size(cont, mapping)
      : 0U
      ;
  } // OpDetVisibilityData::effectiveSize(DataContainer_t, Mapping_t)


  //----------------------------------------------------------------------------
  
} // namespace phot


#endif // LARSIM_PHOTONPROPAGATION_LIBRARYMAPPINGTOOLS_OPDETVISIBILITYDATA_H
