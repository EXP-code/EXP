#include <dirent.h>

#include <algorithm>
#include <exception>
#include <memory>
#include <vector>
#include <utility>

namespace scandirpp {
  typedef std::unique_ptr<struct dirent> EntryPtr;
  typedef std::vector<EntryPtr> EntryPtrVector;

  class ScandirException : public std::system_error {
  public:
    ScandirException(const std::string& dir) :
      std::system_error(errno, std::system_category(), std::string("Error scanning directory \"" + dir + "\"").c_str())
    {}
  };

  inline EntryPtrVector get_entry_ptrs(const std::string& dir) {
    struct dirent** ptrs;
    int n = scandir(dir.c_str(), &ptrs, nullptr, alphasort);
    if (n < 0) {
      throw ScandirException(dir);
    }
    return EntryPtrVector(ptrs, ptrs+n);
  }

  template<class T>
  bool default_filter(const T&) { return true; }

  inline struct dirent extract_entry(const struct dirent& entry) {
    return {entry};
  }

  inline std::string extract_name(const struct dirent& entry) {
    return entry.d_name;
  }

  inline ino_t extract_ino(const struct dirent& entry) {
    return entry.d_ino;
  }

  template<class OutputIterator,
           class ValueExtractor,
           class ValueFilter = bool (*)(const decltype(std::declval<ValueExtractor>()(std::declval<const struct dirent>()))&),
           class EntryFilter = bool (*)(const struct dirent&)>
  inline void get_values(const std::string& dir,
                         OutputIterator result,
                         ValueExtractor value_extractor,
                         ValueFilter value_filter = default_filter,
                         EntryFilter entry_filter = default_filter)
  {
    const EntryPtrVector& entry_ptrs = get_entry_ptrs(dir);

    auto copy_matching = [&result, &value_extractor, &entry_filter, &value_filter](const EntryPtr& entry_ptr) {
      struct dirent& entry = *(entry_ptr.get());
      if (entry_filter(entry)) {
        auto value = value_extractor(entry);
        if (value_filter(value)) {
          *result = value;
          result++;
        }
      }
    };

    std::for_each(entry_ptrs.begin(), entry_ptrs.end(), copy_matching);
  }

  template<class ValueExtractor,
           class ValueType   = decltype(std::declval<ValueExtractor>()(std::declval<const struct dirent>())),
           class ValueFilter = bool (*)(const ValueType&),
           class EntryFilter = bool (*)(const struct dirent&)>
  inline std::vector<ValueType> get_vector(const std::string& dir,
                                           ValueExtractor&& value_extractor,
                                           ValueFilter&& value_filter = default_filter,
                                           EntryFilter&& entry_filter = default_filter)
  {
    std::vector<ValueType> result;
    get_values(dir, std::back_inserter(result), std::forward<ValueExtractor>(value_extractor),
               std::forward<ValueFilter>(value_filter), std::forward<EntryFilter>(entry_filter));
    return result;
  }

  template<class ValueFilter = bool (*)(const struct dirent&), class EntryFilter = bool (*)(const struct dirent&)>
  inline std::vector<struct dirent> get_entries(const std::string& dir,
                                                ValueFilter&& value_filter = default_filter,
                                                EntryFilter&& entry_filter = default_filter)
  {
    return get_vector(dir, extract_entry, std::forward<ValueFilter>(value_filter), std::forward<EntryFilter>(entry_filter));
  }

  template<class ValueFilter = bool (*)(const std::string&), class EntryFilter = bool (*)(const struct dirent&)>
  inline std::vector<std::string> get_names(const std::string& dir,
                                            ValueFilter&& value_filter = default_filter,
                                            EntryFilter&& entry_filter = default_filter)
  {
    return get_vector(dir, extract_name, std::forward<ValueFilter>(value_filter), std::forward<EntryFilter>(entry_filter));
  }

  template<class ValueFilter = bool (*)(const ino_t&), class EntryFilter = bool (*)(const struct dirent&)>
  inline std::vector<ino_t> get_inos(const std::string& dir,
                                     ValueFilter&& value_filter = default_filter,
                                     EntryFilter&& entry_filter = default_filter)
  {
    return get_vector(dir, extract_ino, std::forward<ValueFilter>(value_filter), std::forward<EntryFilter>(entry_filter));
  }

}
