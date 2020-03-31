#include <algorithm>
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <tuple>
#include <map>

// Wrapper class for std::map to keep a finite size
//
template<typename key_t, typename value_t>
class Cache : public std::map<key_t, value_t>
{
public:
  Cache(int N=1) : max_size(N) {}

  void resize(int N) { max_size = N; }

  void add( key_t key, value_t value )
  {
    (*this)[key] = value;
    if (this->size() > max_size) {
      this->erase(--(*this).end());
    }
  }

private:
  int max_size;
};


// Class for representing a point
// Coordinate_type must be a numeric type
template<typename coordinate_type, size_t dimensions>
class point
{
public:
  // Null initializer
  point(std::array<coordinate_type, dimensions> c) : coords_(c) {}
  
  // List constructor
  point(std::initializer_list<coordinate_type> list)
  {
    size_t n = std::min(dimensions, list.size());
    std::copy_n(list.begin(), n, coords_.begin());
  }

  // Copy constructor
  point(const point& p)
  {
    for (size_t n=0; n<dimensions; n++) coords_[n] = p.coords_[n];
  }

  /**
   * Returns the coordinate in the given dimension.
   *
   * @param index dimension index (zero based)
   * @return coordinate in the given dimension
   */
  coordinate_type get(size_t index) const
  {
    return coords_[index];
  }
  
  /**
   * Returns the distance squared from this point to another
   * point.
   *
   * @param pt another point
   * @return distance squared from this point to the other point
   */
  double distance(const point& pt) const
  {
    double dist = 0;
    for (size_t i = 0; i < dimensions; ++i)
      {
	double d = get(i) - pt.get(i);
	dist += d * d;
      }
    return dist;
  }
private:
  std::array<coordinate_type, dimensions> coords_;
};

// For iostream printing of points
//
template<typename coordinate_type, size_t dimensions>
std::ostream& operator<<(std::ostream& out,
			 const point<coordinate_type, dimensions>& pt)
{
  out << '(';
  for (size_t i = 0; i < dimensions; ++i) {
    if (i > 0) out << ", ";
    out << pt.get(i);
  }
  out << ')';
  return out;
}

// k-d tree implementation
//
template<typename coordinate_type, size_t dimensions>
class kdtree
{
public:
  typedef point<coordinate_type, dimensions> point_type;

private:
  struct node
  {
    node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr)
    {
    }
    
    coordinate_type get(size_t index) const
    {
      return point_.get(index);
    }
    
    double distance(const point_type& pt) const
    {
      return point_.distance(pt);
    }
    point_type point_;
    node* left_;
    node* right_;
  };

  node* root_;

  Cache<double, node*> best_;

  size_t visited_;
  std::vector<node> nodes_;
  
  struct node_cmp
  {
    node_cmp(size_t index) : index_(index)
    {
    }

    bool operator()(const node& n1, const node& n2) const
    {
      return n1.point_.get(index_) < n2.point_.get(index_);
    }
    size_t index_;
  };
  
  node* make_tree(size_t begin, size_t end, size_t index)
  {
    if (end <= begin) return nullptr;
    size_t n = begin + (end - begin)/2;
    std::nth_element(&nodes_[begin], &nodes_[n], &nodes_[end], node_cmp(index));
    index = (index + 1) % dimensions;
    nodes_[n].left_ = make_tree(begin, n, index);
    nodes_[n].right_ = make_tree(n + 1, end, index);
    return &nodes_[n];
  }
  
  void nearestN(node* root, const point_type& point, size_t index, int N)
  {
    if (root == nullptr) return;

    ++visited_;

    double d = root->distance(point);
    if (best_.size()==0 || d < best_.rbegin()->first) {
      best_.add(d, root);
    }

    // This is only correct is the test point is never in the data set . . .
    // if (best_.begin()->first == 0) return;

    double dx = root->get(index) - point.get(index);
    index = (index + 1) % dimensions;
    nearestN(dx > 0 ? root->left_ : root->right_, point, index, N);

    if (dx * dx >= best_.rbegin()->first) return;
    nearestN(dx > 0 ? root->right_ : root->left_, point, index, N);
  }
  
public:
  //@{
  // Copy constructor2
  kdtree(const kdtree&) = delete;
  kdtree& operator=(const kdtree&) = delete;
  //@}
  
  /**
   * Constructor taking a pair of iterators. Adds each
   * point in the range [begin, end) to the tree.
   *
   * @param begin start of range
   * @param end end of range
   */
  template<typename iterator>
  kdtree(iterator begin, iterator end)
  {
    best_.resize(1);
    visited_ = 0;
    nodes_.reserve(std::distance(begin, end));
    for (auto i = begin; i != end; ++i)
      nodes_.emplace_back(*i);
    root_ = make_tree(0, nodes_.size(), 0);
  }
  
  /**
   * Constructor taking a function object that generates
   * points. The function object will be called n times
   * to populate the tree.
   *
   * @param f function that returns a point
   * @param n number of points to add
   */
  template<typename func>
  kdtree(func&& f, size_t n)
  {
    best_.resize(1);
    visited_ = 0;
    nodes_.reserve(n);
    for (size_t i = 0; i < n; ++i)
      nodes_.emplace_back(f());
    root_ = make_tree(0, nodes_.size(), 0);
  }
  
  /**
   * Returns true if the tree is empty, false otherwise.
   */
  bool empty() const
  {
    return nodes_.empty();
  }
  
  /**
   * Returns the number of nodes visited by the last call
   * to nearest().
   */
  size_t visited() const
  {
    return visited_;
  }
  
  /**
   * Returns the distance between the input point and return value
   * from the last call to nearest().
   */
  double distance() const
  {
    return std::sqrt(best_.begin()->first);
  }
  
  /**
   * Finds the nearest point in the tree to the given point.
   * It is not valid to call this function if the tree is empty.
   *
   * @param pt a point
   * @param the nearest point in the tree to the given point
   */
  std::tuple<point_type, point_type, double>
  nearestN(const point_type& pt, int N)
  {
    if (root_ == nullptr) throw std::logic_error("tree is empty");
    best_.clear();
    best_.resize(N);
    visited_ = 0;
    nearestN(root_, pt, 0, N);
    return {best_.begin()->second->point_,
	    (best_.begin()++)->second->point_,
	    std::sqrt(best_.rbegin()->first)};
  }

  std::vector<double> getDist()
  {
    std::vector<double> ret;
    for (auto v : best_) ret.push_back(v.first);
    return ret;
  }

};

