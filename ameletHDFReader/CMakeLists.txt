PROJECT(ameletHDFReader)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#-------------------------------------------------------------------------------
# Found all dependaces.
# ------------------------------------------------------------------------------
FIND_PACKAGE(HDF5 REQUIRED)

FIND_PACKAGE(AmeletHDF REQUIRED)
LINK_DIRECTORIES(${AMELETHDF_LIBRARY_DIRS})

set(VTKAH5_LIBRARIES_DIRS
  ${VTK_LIBRARY_DIRS}
  ${AMELETHDF_LIBRARY_DIRS}
  ${HDF5_LIBRARY_DIR})

IF(WIN32)
   set(LIBAMELETHDFC amelethdfc.lib)
ELSE(WIN32)
   set(LIBAMELETHDFC libamelethdfc.a)
ENDIF(WIN32)

set(VTKAH5_LIBRARIES
  amelethdf_readers
  common
  ${LIBAMELETHDFC}
  hdf5 hdf5_hl
  ${HDF5_LIBRARIES})


set(VTKAH5_LIBRARIES_INCLUDE_DIRS
  ${PARAVIEW_INCLUDE_DIRS}
  ${PARAVIEW_GUI_INCLUDE_DIRS}
  ${VTKCLIENTSERVER_INCLUDE_DIR}
  ${PVSERVERCOMMON_INCLUDE_DIR}
  ${AMELETHDF_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIR})

# INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
  ${VTKAH5_LIBRARIES_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/ahdfreaders
  ${CMAKE_CURRENT_SOURCE_DIR}/common
  )

#-------------------------------------------------------------------------------
# Configure the Paraview plugin
#-------------------------------------------------------------------------------
# This is just a simple variable that lists the source files.
SET(SRCS
  vtkAmeletHDFReader.cxx
  )

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
ADD_PARAVIEW_PLUGIN(ameletHDFReader "0.2"
  SERVER_MANAGER_SOURCES ${SRCS}
  SERVER_MANAGER_XML ameletHDFReaderSM.xml
#  GUI_RESOURCE_FILES ameletHDFReaderGUI.xml
  )


#-------------------------------------------------------------
# Configure compiler
#-------------------------------------------------------------
if ( "${CMAKE_C_COMPILER}" MATCHES "clang" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" )
  set ( CMAKE_COMPILER_IS_CLANG 1 )
endif()

IF(CMAKE_COMPILER_IS_GNUCXX)
  set(CXX_COMPILE_FLAGS           "-Wall -Wextra -fPIC")
  set(CXX_DEBUG_COMPILE_FLAGS     "-g3 -D_DEBUG")
  set(CXX_RELEASE_COMPILE_FLAGS   "-O3 -g0 -DNDEBUG") # -Werror

  set(CXX_COVERAGE_COMPILE_FLAGS  "-fprofile-arcs -ftest-coverage")
  set(CXX_COVERAGE_LINK_FLAGS     "-lgcov -fprofile-arcs -ftest-coverage")

ELSEIF (MSVC10) # check is VC on win
  set(CXX_COMPILE_FLAGS           "/nologo /W1 /Oy- /Gm- /EHsc /GS /fp:precise /Zc:wchar_t- /Zc:forScope /Gd")
  set(CXX_DEBUG_COMPILE_FLAGS     "/Zi /WX- /Od /MDd /RTC1 /Oy-")
  set(CXX_RELEASE_COMPILE_FLAGS   "/D_WIN32 /O2 /MD ")

ELSEIF (CMAKE_COMPILER_IS_CLANG)
  set(CXX_COMPILE_FLAGS           "-Wall -Wextra -fPIC")
  set(CXX_DEBUG_COMPILE_FLAGS     "-g3 -D_DEBUG")
  set(CXX_RELEASE_COMPILE_FLAGS   "-O3 -g0 -DNDEBUG") # -Werror

ELSE ()
  MESSAGE( STATUS "Unsupported compiler defaut configuration is used!")
ENDIF ()

set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} ${CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG
  "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_DEBUG_COMPILE_FLAGS} ${CXX_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_RELEASE_COMPILE_FLAGS}")

SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG
  "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${CXX_COVERAGE_LINK_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS_DEBUG
  "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${CXX_COVERAGE_LINK_FLAGS}")


#-------------------------------------------------------------------------------
# Build the library and Paraview plugin
# ------------------------------------------------------------------------------
# We need to add library dependencies for the plugin.
ADD_SUBDIRECTORY(ahdfreaders)
ADD_SUBDIRECTORY(common)

TARGET_LINK_LIBRARIES(ameletHDFReader ${VTKAH5_LIBRARIES})


#-------------------------------------------------------------------------------
# Build Tools
# ------------------------------------------------------------------------------
# find_package(VTK REQUIRED)
# Define here the extra dependance for tools
set(VTKAH5_EXTRA_TOOLS_LIBRARIES_DIRS)
set(VTKAH5_EXTRA_TOOLS_LIBRARIES ${VTK_LIBRARIES})
set(VTKAH5_EXTRA_TOOLS_LIBRARIES_INCLUDE)

ADD_SUBDIRECTORY(tools)


#-------------------------------------------------------------------------------
# Build tests
# ------------------------------------------------------------------------------
ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

# Define here the extra dependance for tests
set(VTKAH5_EXTRA_TESTS_LIBRARIES_DIRS ${Boost_LIBRARY_DIRS})
set(VTKAH5_EXTRA_TESTS_LIBRARIES ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${VTK_LIBRARIES})
set(VTKAH5_EXTRA_TESTS_LIBRARIES_INCLUDE ${Boost_INCLUDE_DIRS})

ADD_SUBDIRECTORY(tests)
