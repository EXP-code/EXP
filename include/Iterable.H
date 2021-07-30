#include <iterator>

class myit
{
private:
    int value_;
    class intholder
    {
        int value_;
    public:
        intholder(int value): value_(value) {}
        int operator*() { return value_; }
    };
public:
    // Previously provided by std::iterator - see update below
    typedef int                     value_type;
    typedef std::ptrdiff_t          difference_type;
    typedef int*                    pointer;
    typedef int&                    reference;
    typedef std::input_iterator_tag iterator_category;

    explicit myit(int value) : value_(value) {}
    int operator*() const { return value_; }
    bool operator==(const myit& other) const { return value_ == other.value_; }
    bool operator!=(const myit& other) const { return !(*this == other); }
    intholder operator++(int)
    {
        intholder ret(value_);
        ++*this;
        return ret;
    }
    myit& operator++()
    {
        ++value_;
        return *this;
    }
};
