add_subdirectory(c)

if(SWIG_PYTHON)
    add_subdirectory(python)
endif()

if(SWIG_MATLAB)
    add_subdirectory(matlab)
endif()
