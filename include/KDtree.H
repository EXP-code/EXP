#include <algorithm>
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <tuple>
#include <map>

namespace KDtree
{
  //! Wrapper class for std::map to keep a finite size
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


  /** Class for representing a point
      Coordinate_type must be a numeric type
      Field (weight) is a double
  */
  template<typename coordinate_type, size_t dimensions>
  class point
  {
  public:
    //! Null initializer
    point() {}
    
    //! Constructor
    point(std::array<coordinate_type, dimensions> c, double m=1,
	  unsigned long indx=0) : coords_(c), mass_(m), indx_(indx) {}
  
    point(std::array<coordinate_type, dimensions> c,
	  std::array<coordinate_type, dimensions> v,
	  double m=1, unsigned long indx=0) : coords_(c), vels_(v),
					      mass_(m), indx_(indx) {}
  
    //! List constructor
    point(std::initializer_list<coordinate_type> list, double m=1,
	  unsigned long indx=0) : mass_(m), indx_(indx)
    {
      size_t n = std::min(dimensions, list.size());
      std::copy_n(list.begin(), n, coords_.begin());
    }
    
    point(std::initializer_list<coordinate_type> list,
	  std::initializer_list<coordinate_type> vlst,
	  double m=1,
	  unsigned long indx=0) : mass_(m), indx_(indx)
    {
      size_t n = std::min(dimensions, list.size());
      std::copy_n(list.begin(), n, coords_.begin());
      std::copy_n(vlst.begin(), n, vels_.begin());
    }

    //! Copy constructor
    point(const point& p)
    {
      for (size_t n=0; n<dimensions; n++) coords_[n] = p.coords_[n];
      for (size_t n=0; n<dimensions; n++) vels_  [n] = p.vels_  [n];
      mass_ = p.mass_;
      indx_ = p.indx_;
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
    
    coordinate_type vel(size_t index) const
    {
      return vels_[index];
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

    /**
     * Returns the squared velocity between this point to another
     * point.
     *
     * @param pt another point
     * @return squared speed between this point and the other point
     */
    double speed(const point& pt) const
    {
      double rvel = 0;
      for (size_t i = 0; i < dimensions; ++i)
	{
	  double d = vel(i) - pt.vel(i);
	  rvel += d * d;
	}
      return rvel;
    }
    
    double mass() const
    {
      return mass_;
    }
    
    unsigned long indx() const
    {
      return indx_;
    }

  private:
    std::array<coordinate_type, dimensions> coords_;
    std::array<coordinate_type, dimensions> vels_;
    double mass_;
    unsigned long indx_;
  };

  //! For iostream printing of points
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

  //! k-d tree implementation
  template<typename coordinate_type, size_t dimensions>
  class kdtree
  {
  public:

    using point_type = point<coordinate_type, dimensions>;

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
    
    using cache = Cache<double, node*>;

  private:

    node* root_;

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
      nodes_[n].left_  = make_tree(begin, n, index);
      nodes_[n].right_ = make_tree(n + 1, end, index);
      return &nodes_[n];
    }
  
    void nearestN(node* root, const point_type& point, size_t index, int N,
		  cache& best_)
    {
      if (root == nullptr) return;

      ++visited_;		// For diagnostic info only

      double d = root->distance(point);
      if (best_.size()<N || d < best_.rbegin()->first) {
	best_.add(d, root);
      }
      
      // This is only correct if the test point is never in the data set . . .
      // if (best_.begin()->first == 0) return;
      
      double dx = root->get(index) - point.get(index);
      index = (index + 1) % dimensions;
      
      nearestN(dx > 0 ? root->left_  : root->right_, point, index, N, best_);
      
      if (best_.size()>=N and dx * dx >= best_.rbegin()->first) return;
      nearestN(dx > 0 ? root->right_ : root->left_,  point, index, N, best_);
    }
  
    void accum(node* root, std::vector<std::vector<unsigned long>>& ret, unsigned cur)
    {
      ret[cur].push_back(root->point_.indx());
      if (root->left_)  accum(root->left_,  ret, cur);
      if (root->right_) accum(root->right_, ret, cur);
    }

    void walk(node* root, std::vector<std::vector<unsigned long>>& ret,
	      unsigned int cur, unsigned int lev, unsigned int level)
    {
      if (lev == level) accum(root, ret, cur);
      else {
	ret[cur].push_back(root->point_.indx());
	walk(root->left_,  ret, 2*cur + 0, lev + 1, level);
	walk(root->right_, ret, 2*cur + 1, lev + 1, level);
      }
    }


  public:

    //@{
    //! Copy constructor2
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
    double distance(cache& best_) const
    {
      return std::sqrt(best_.begin()->first);
    }
    
    /**
     * Finds the nearest N points in the tree to the given point.  It is
     * not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @param N is the number of nearest points
     *
     * Returns: tuple of the first points, summed weight, and the radius of the Nth
     * point
     */
    std::tuple<point_type, double, double, cache>
    nearestN(const point_type& pt, int N)
    {
      if (root_ == nullptr) throw std::logic_error("tree is empty");
      cache best_(N);
      visited_ = 0;
      nearestN(root_, pt, 0, N, best_);

      double wgt = 0.0;		// Sum weights
      for (auto b : best_) wgt += b.second->point_.mass();

#if __GNUC__ > 6 or defined(__clang__)
      return {best_.begin()->second->point_, wgt, std::sqrt(best_.rbegin()->first), best_};
#else
      std::tuple<point_type, double, double, cache> ret;
      std::get<0>(ret) = best_.begin()->second->point_;
      std::get<1>(ret) = wgt;
      std::get<2>(ret) = std::sqrt(best_.rbegin()->first);
      std::get<3>(ret) = best_;
      return ret;
#endif
    }

    /**
     * Finds the nearest N points in the tree to the given point.  It is
     * not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @param N is the number of nearest points
     
     * Returns: tuple of the nearest point list and the radius of the
     * Nth point
     */
    std::tuple<std::vector<point_type>, double, cache>
    nearestList(const point_type& pt, int N)
    {
      if (root_ == nullptr) throw std::logic_error("tree is empty");
      cache best_(N);
      visited_ = 0;
      nearestN(root_, pt, 0, N, best_);

      std::vector<point_type> pts; // The returned point list
      for (auto b : best_) pts.push_back(b.second->point_);

#if __GNUC__ > 6
      return {pts, std::sqrt(best_.rbegin()->first), best_};
#else
      std::tuple<std::vector<point_type>, double> ret;
      std::get<0>(ret) = pts;
      std::get<1>(ret) = std::sqrt(best_.rbegin()->first);
      std::get<2>(ret) = best_;
      return ret;
#endif
    }

    std::vector<double> getDist(cache& best_)
    {
      std::vector<double> ret;
      for (auto v : best_) ret.push_back(v.first);
      return ret;
    }
    
    std::vector<std::vector<unsigned long>> getPartition(unsigned level)
    {
      int partitions = (1 << level);
      std::vector<std::vector<unsigned long>> ret(partitions);
      walk(root_, ret, 0, 0, level);
      return ret;
    }

  };

}
// END: namespace KDtree
