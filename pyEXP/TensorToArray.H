#ifndef _TENSOR_TO_ARRAY_H_
#define _TENSOR_TO_ARRAY_H_

//! Helper function that maps the Eigen::Tensor<T, 3> into an numpy.ndarray
template <typename T>
py::array_t<T> make_ndarray3(Eigen::Tensor<T, 3>& mat)
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

template <typename T>
Eigen::Tensor<T, 3> make_tensor3(py::array_t<T>& in)
{
  // Request a buffer descriptor from Python
  py::buffer_info buffer_info = in.request();

  // Extract the data 
  T *data = static_cast<T *>(buffer_info.ptr);

  // Get the data shape
  std::vector<ssize_t> shape = buffer_info.shape;

  // Check rank
  if (shape.size() != 3) {
    std::ostringstream sout;
    sout << "make_tensor3: tensor rank must be 3, found "
	 << shape.size();
    throw std::runtime_error(sout.str());
  }

  // Reorder the data to satisfy the col-major Eigen::Tensor ordering
  //
  Eigen::Tensor<T, 3> tensor(shape[0], shape[1], shape[2]);
  for (int i=0, c=0; i < shape[0]; i++) {
    for (int j=0; j < shape[1]; j++) {
      for (int k=0; k < shape[2]; k++, c++) {
	tensor(i, j, k) = data[c];
      }
    }
  }

  return tensor;
}

template <typename T>
Eigen::Tensor<T, 4> make_tensor4(py::array_t<T> array)
{
  // Request a buffer descriptor from Python
  py::buffer_info buffer_info = array.request();

  // Get the array dimenions
  T *data = static_cast<T *>(buffer_info.ptr);
  std::vector<ssize_t> shape = buffer_info.shape;

  // Check rank
  if (shape.size() != 4) {
    std::ostringstream sout;
    sout << "make_tensor4: tensor rank must be 4, found "
	 << shape.size();
    throw std::runtime_error(sout.str());
  }

  // Build result tensor with col-major ordering
  Eigen::Tensor<T, 4> tensor(shape[0], shape[1], shape[2], shape[3]);
  for (int i=0, l=0; i < shape[0]; i++) {
    for (int j=0; j < shape[1]; j++) {
      for (int k=0; k < shape[2]; k++) {
	for (int l=0; l < shape[3]; k++) {
	  tensor(i, j, k, l) = data[l++];
	}
      }
    }
  }

  return tensor;
}


#endif
