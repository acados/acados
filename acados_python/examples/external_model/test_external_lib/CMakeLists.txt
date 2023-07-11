CMAKE_MINIMUM_REQUIRED ( VERSION 3.8 )
PROJECT("ExternalOdeCasADi" CXX C)

# Set folder to look path to
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/modules)

# Set a default build type to 'Release' if none was specified
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
ENDIF()

find_package(Casadi REQUIRED)
link_directories(${filename}
  ${CMAKE_INSTALL_PREFIX}/lib
)


add_library(external_ode_casadi SHARED src/external_ode_casadi.cpp)

target_link_libraries(external_ode_casadi
    ${CASADI_LIBRARIES}
  )
    
set_target_properties(external_ode_casadi PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_include_directories(external_ode_casadi PUBLIC
    include
    ${CMAKE_INSTALL_PREFIX}/include
)

# Write down all files to include here
set(source_files
    test_my_lib.cpp
)

# Setup each project
foreach(src_file ${source_files})
    get_filename_component(filename ${src_file} NAME_WE)

    ADD_EXECUTABLE (${filename}
        ${filename}.cpp
    )
    target_include_directories(${filename} PUBLIC
        ${CASADI_INCLUDE_DIR}
    )

    target_link_libraries(${filename}
        external_ode_casadi
        ${CASADI_LIBRARIES}
      )
endforeach()


install(TARGETS external_ode_casadi DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
