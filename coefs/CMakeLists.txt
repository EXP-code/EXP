set(TESTH5_SRC testh5.cc)

include_directories($<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/>
  ${HighFive_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
add_executable(testh5 ${TESTH5_SRC})
target_link_libraries(testh5 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
