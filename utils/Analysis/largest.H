#ifndef _largest_H
#define _largest_H

#include <queue>

// Encapsulate a priority queue for maximum value list
//
template <class U>
class largestValues
{
public:

  //! Construtor: size is the desired size of the heap
  largestValues(int size) : maxSize(size) {}

  //! Add a new value set
  void operator()(U data)
  {
    if (pq.size()<maxSize) {
      pq.push(data);
      return;
    }

    if (data.val > pq.top().val) {
      pq.pop();		 // Get rid of the root
      pq.push(data); // Priority queue will automatically restructure
    }
  }

  //! Get the current list of values
  std::vector<U> getValues()
  {
    std::vector<U> ret;		// Make a return vector
    auto cq = pq;		// Copy the queue
    while(cq.size()!=0) {	// Keep going until empty
      ret.push_back(cq.top());	// Push back the root element
      cq.pop();			// Remove the root element
    }
    return ret;			// Return the vector
  }

  //! Empty the queue
  void clear() {
    while(pq.size()!=0) pq.pop();
  }

private:

  std::priority_queue<U, std::vector<U>, std::greater<U>> pq;
  int maxSize;

};

#endif
