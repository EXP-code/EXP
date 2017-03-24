// Compile: g++ -g -o jt json_test.cc
// Run: ./jt

#include <iostream>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string/predicate.hpp>

// Custom translator for bool (only supports std::string)
struct BoolTranslator
{
    typedef std::string internal_type;
    typedef bool        external_type;

    // Converts a string to bool
    boost::optional<external_type> get_value(const internal_type& str)
    {
        if (!str.empty())
        {
            using boost::algorithm::iequals;

            if (iequals(str, "true") || iequals(str, "yes") || str == "1")
                return boost::optional<external_type>(true);
            else
                return boost::optional<external_type>(false);
        }
        else
            return boost::optional<external_type>(boost::none);
    }

    // Converts a bool to string
    boost::optional<internal_type> put_value(const external_type& b)
    {
        return boost::optional<internal_type>(b ? "true" : "false");
    }
};

/*  Specialize translator_between so that it uses our custom translator for
    bool value types. Specialization must be in boost::property_tree
    namespace. */
namespace boost {
namespace property_tree {

template<typename Ch, typename Traits, typename Alloc> 
struct translator_between<std::basic_string< Ch, Traits, Alloc >, bool>
{
    typedef BoolTranslator type;
};

} // namespace property_tree
} // namespace boost

int main()
{
    boost::property_tree::ptree pt;

    read_json("test.json", pt);
    double f = pt.get<double>("value");
    int  i = pt.get<int>("number");
    bool b = pt.get<bool>("enabled");
    bool g = pt.get<bool>("good");
    std::string s = pt.get<std::string>("name");
    std::cout << "f=" << f << std::endl
	      << "i=" << i << std::endl
	      << "b=" << std::boolalpha << b << std::endl
	      << "g=" << std::boolalpha << g << std::endl
	      << "s=" << s << std::endl;
}
