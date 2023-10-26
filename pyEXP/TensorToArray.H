#ifndef _TENSOR_TO_ARRAY_H_
#define _TENSOR_TO_ARRAY_H_

//! Helper function that maps the Eigen::Tensor<T, 3> into an numpy.ndarray
template <typename T>
py::array_t<T> make_ndarray(Eigen::Tensor<T, 3>& mat)
{
  // Get the tensor dimenions
  auto dims = mat.dimensions();

  // Check rank
  if (dims.size() != 3) {
    std::ostringstream sout;
    sout << "make_ndarray: tensor rank must be 3, found " << dims.size();
    throw std::runtime_error(sout.str());
  }
  
  // Make the memory mapping
  return py::array_t<T>
    (
     // shape
     {dims[0], dims[1], dims[2]},
     // C-style contiguous strides for double
     {sizeof(T), dims[0]*sizeof(T), dims[0]*dims[1]*sizeof(T)},
     // the data pointer
     mat.data()
     );
}

//! Helper function that maps the Eigen::Tensor<T, 4> into an numpy.ndarray
template <typename T>
py::array_t<T> make_ndarray4(Eigen::Tensor<T, 4>& mat)
{
  // Get the tensor dimenions
  auto dims = mat.dimensions();

  // Check rank
  if (dims.size() != 4) {
    std::ostringstream sout;
    sout << "make_ndarray: tensor rank must be 4, found " << dims.size();
    throw std::runtime_error(sout.str());
  }
  
  // Make the memory mapping
  return py::array_t<T>
    (
     // shape
     {dims[0], dims[1], dims[2], dims[3]},
     // C-style contiguous strides for double
     {sizeof(T), dims[0]*sizeof(T), dims[0]*dims[1]*sizeof(T), dims[0]*dims[1]*dims[2]*sizeof(T)},
     // the data pointer
     mat.data()
     );
}

#endif
