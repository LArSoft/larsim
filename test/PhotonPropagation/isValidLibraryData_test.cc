/**
 * @file    isValidLibraryData_test.cc
 * @brief   Unit test for `phot::isValidLibraryData()`.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 29, 2019
 * @see     `larsim/PhotonPropagation/LibraryMappingTools/OpDetVisibilityData.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ( isValidLibraryData_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK()

// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/OpDetVisibilityData.h"
#include "larcorealg/CoreUtils/MetaUtils.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>
#include <functional> // std::reference_wrapper
#include <type_traits>


//------------------------------------------------------------------------------
class SizedContainer {
    protected:
  std::size_t fSize = 0U;
  
    public:
  SizedContainer(std::size_t size = 0U): fSize(size) {}
  
}; // SizedContainer

struct ContainerWith_size: public SizedContainer {
  using SizedContainer::SizedContainer;
  
  std::size_t size() const { return fSize; }
  
}; // ContainerWith_size

struct ContainerWith_empty: public SizedContainer {
  using SizedContainer::SizedContainer;
  
  bool empty() const { return fSize == 0U; }
  
}; // ContainerWith_empty

struct ContainerWith_isValid: public SizedContainer {
  using SizedContainer::SizedContainer;
  
  bool isValid() const { return fSize > 0U; }
  
}; // ContainerWith_isValid

struct ContainerWith_is_valid: public SizedContainer {
  using SizedContainer::SizedContainer;
  
  bool is_valid() const { return fSize > 0U; }
  
}; // ContainerWith_is_valid


//------------------------------------------------------------------------------
void isValidLibraryData_test() {
  
  std::vector<float> valid_vector(4U, 1.0);
  std::vector<float> invalid_vector;
  
  std::array<float, 4U> valid_STLarray;
  // all STL arrays are valid to us
  
  float valid_array[4U] = { 1.0, 2.0, 3.0, 4.0 };
  // all C arrays are valid to us
  
  float const* valid_cptr = valid_array;
  float const* invalid_cptr = nullptr;
  
  float const* const valid_cptrc = valid_array;
  float const* const invalid_cptrc = nullptr;
  
  std::unique_ptr<float> valid_uptr { new float[4U] };
  std::unique_ptr<float> invalid_uptr;
  
  ContainerWith_size       valid_with_size    { 4U };
  ContainerWith_size     invalid_with_size    ;
  
  ContainerWith_empty      valid_with_empty   { 4U };
  ContainerWith_empty    invalid_with_empty   ;
  
  ContainerWith_is_valid   valid_with_is_valid{ 4U };
  ContainerWith_is_valid invalid_with_is_valid;
  
  ContainerWith_isValid    valid_with_isValid { 4U };
  ContainerWith_isValid  invalid_with_isValid ;
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_vector));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_vector));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_STLarray));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_array));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_cptr));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_cptr));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_cptrc));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_cptrc));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_uptr));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_uptr));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_with_size    ));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_with_size    ));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_with_empty   ));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_with_empty   ));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_with_is_valid));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_with_is_valid));
  
  BOOST_CHECK( phot::isValidLibraryData(  valid_with_isValid ));
  BOOST_CHECK(!phot::isValidLibraryData(invalid_with_isValid ));
  
} // isValidLibraryData_test()


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(isValidLibraryData_TestCase) {
  isValidLibraryData_test();
} // BOOST_AUTO_TEST_CASE(isValidLibraryData_TestCase)

//------------------------------------------------------------------------------
