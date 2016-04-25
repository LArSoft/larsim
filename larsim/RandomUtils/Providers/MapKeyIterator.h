/**
 * @file   MapKeyIterator.h
 * @brief  Provides a key iterator from a map-type object
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 18th, 2015
 *
 * The classes provided here are templates.
 */
#ifndef LARSIM_RANDOMUTILS_PROVIDERS_MAPKEYITERATOR_H
#define LARSIM_RANDOMUTILS_PROVIDERS_MAPKEYITERATOR_H 1

// C++ standard libraries
#include <iterator>
#include <tuple>
#include <type_traits>


namespace sim {
  namespace SeedServiceHelper {
    
    
    /** ************************************************************************
     * @brief Wraps a map-like object iterator to iterate on its first element
     * @tparam Iter the iterator type
     *
     * This object can wrap a map-like iterator type to iterate through its key
     * in a for loop.
     * 
     * The Iter::value_type must support a std::get<0>() call, that is what's
     * used to extract the key.
     *
     * Example of usage:
     *     
     *     std::map<int, double> my_map;
     *     // ... fill the map ...
     *     // print all the kets
     *     for (MapKeyIterator<map_id_t::const_iterator> it = map_id.cbegin();
     *       it != map_id.cend(); ++it)
     *       std::cout << "Key: " << *it << std::endl;
     *     
     * Although the non-constant iterators are supported, this is of little use
     * since typically the key is declared constant anyway (in the example
     * above, the key is `const int`).
     */
    template <typename Iter>
    class MapKeyIterator: public Iter {
        public:
      using map_iterator_t = Iter;
      using iterator_traits_t = std::iterator_traits<map_iterator_t>;
      
      // iterator traits
      using difference_type = typename iterator_traits_t::difference_type;
      using value_type
        = typename std::tuple_element<0, typename iterator_traits_t::value_type>::type;
      using pointer = typename std::add_pointer<value_type>::type;
      using reference = typename std::add_lvalue_reference<value_type>::type;
      using iterator_category = typename iterator_traits_t::iterator_category;
      
      /// Default constructor is really default
      MapKeyIterator() = default;
      
      // accept quietly all the other constructery C++ provides
      
      /// Initialize from a original iterator type
      MapKeyIterator(map_iterator_t const& from):
        map_iterator_t(from) {}
      
      /// Return the key of the original iterator value
      reference operator*() { return get_key(get_value()); }
      
      /// Return a pointer to the key of the original iterator value
      pointer operator->() { return &get_key(get_value()); }
      
        protected:
      using pair_reference = typename iterator_traits_t::reference;
      
      pair_reference get_value() const { return map_iterator_t::operator*(); }
      
      static reference get_key(pair_reference value) 
        { return std::get<0>(value); }
      
    }; // MapKeyIterator<>
    
    
    
    /**
     * @brief Provides iterators for std::begin() and std::end()
     * @tparam Map type of the map to be wrapped
     *
     * This object can wrap a map-like type to iterate through its key in a
     * range-based for loop.
     * 
     * Example of usage:
     *     
     *     std::map<int, double> my_map;
     *     // ... fill the map ...
     *     // print all the kets
     *     for (int const& key: MapKeyConstIteratorBox<std::map<int, double>>(my_map))
     *       std::cout << "Key: " << key << std::endl;
     *     
     * Only constant iteration on keys is provided.
     */
    template <typename Map>
    class MapKeyConstIteratorBox {
        public:
      using map_t = Map;
      
      using const_iterator = MapKeyIterator<typename map_t::const_iterator>;
      
      MapKeyConstIteratorBox(map_t const& map):
        beginKey(map.cbegin()), endKey(map.cend())
        {}
      
      const_iterator begin()  { return beginKey; }
      const_iterator end()    { return endKey; }
      const_iterator cbegin() { return beginKey; }
      const_iterator cend()   { return endKey; }
      
        protected:
      const_iterator beginKey; ///< iterator pointing to the begin key
      const_iterator endKey;   ///< iterator pointing to the endkey
      
    }; // class MapKeyConstIteratorBox<>
  
  } // namespace SeedServiceHelper
} // namespace sim

#endif // LARSIM_RANDOMUTILS_PROVIDERS_MAPKEYITERATOR_H
