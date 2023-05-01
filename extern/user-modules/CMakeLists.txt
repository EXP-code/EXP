# Add all subdirectories

# Get all the files including directories from this level into the
# V_GLOB vector

file(GLOB FVEC LIST_DIRECTORIES true "*")

# Run through items in the FVEC vector and add directories for
# compilation

foreach(item ${FVEC})
  if(IS_DIRECTORY ${item})
    add_subdirectory(${item})
  endif()
endforeach()
      
