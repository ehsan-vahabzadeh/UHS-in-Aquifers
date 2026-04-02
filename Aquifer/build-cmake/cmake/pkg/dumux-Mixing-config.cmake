if(NOT dumux-Mixing_FOUND)
# Whether this module is installed or not
set(dumux-Mixing_INSTALLED ON)

# Settings specific to the module

# Package initialization

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was dumux-Mixing-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#report other information
set_and_check(dumux-Mixing_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dumux-Mixing_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set(dumux-Mixing_CMAKE_CONFIG_VERSION "2.10")
set(dumux-Mixing_CXX_FLAGS "")
set(dumux-Mixing_CXX_FLAGS_DEBUG "-g")
set(dumux-Mixing_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dumux-Mixing_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dumux-Mixing_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dumux-Mixing_DEPENDS "dumux")
set(dumux-Mixing_SUGGESTS "")
set(dumux-Mixing_MODULE_PATH "${PACKAGE_PREFIX_DIR}/share/dune/cmake/modules")
set(dumux-Mixing_PYTHON_WHEELHOUSE "${PACKAGE_PREFIX_DIR}/share/dune/wheelhouse")
set(dumux-Mixing_LIBRARIES "")
set(dumux-Mixing_HASPYTHON 0)
set(dumux-Mixing_PYTHONREQUIRES "")

# Resolve dune dependencies
include(CMakeFindDependencyMacro)
macro(find_and_check_dune_dependency module version)
  find_dependency(${module})
  list(PREPEND CMAKE_MODULE_PATH "${dune-common_MODULE_PATH}")
  include(DuneModuleDependencies)
  list(POP_FRONT CMAKE_MODULE_PATH)
  if(dune-common_VERSION VERSION_GREATER_EQUAL "2.10")
    dune_check_module_version(${module} QUIET REQUIRED VERSION "${version}")
  endif()
endmacro()

find_and_check_dune_dependency(dumux " ")

# Set up DUNE_LIBS, DUNE_FOUND_DEPENDENCIES, DUNE_*_FOUND, and HAVE_* variables
if(dumux-Mixing_LIBRARIES)
  message(STATUS "Setting dumux-Mixing_LIBRARIES=${dumux-Mixing_LIBRARIES}")
  list(PREPEND DUNE_LIBS ${dumux-Mixing_LIBRARIES})
endif()
list(APPEND DUNE_FOUND_DEPENDENCIES dumux-Mixing)
set(DUNE_dumux-Mixing_FOUND TRUE)
set(HAVE_DUMUX_MIXING TRUE)

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


# If this file is found in a super build that includes dumux-Mixing, the 
# `dumux-Mixing-targets.cmake`-file has not yet been generated. This variable
# determines whether the configuration of dumux-Mixing has been completed.
get_property(dumux-Mixing_IN_CONFIG_MODE GLOBAL PROPERTY dumux-Mixing_LIBRARIES DEFINED)

#import the target
if(dumux-Mixing_LIBRARIES AND NOT dumux-Mixing_IN_CONFIG_MODE)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dumux-Mixing-targets.cmake")
endif()

endif()
