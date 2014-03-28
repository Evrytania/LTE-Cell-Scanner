# For OpenCL LTE Scanner (https://github.com/JiaoXianjun/LTE-Cell-Scanner)
# Jiao Xianjun (putaoshu@gmail.com)

# - Find OpenCL
# Find the native OpenCL includes and library
# This module defines
#  OPENCL_INCLUDE_DIR, where to find cl.h, etc.
#  OPENCL_LIBRARIES, the libraries needed to use OpenCL.
#  OPENCL_FOUND, If false, do not try to use OpenCL.

FIND_PATH(OPENCL_INCLUDE_DIR CL/cl.h
  /usr/include
  /opt/AMDAPP/include/
  /usr/local/include
  /usr/local/cuda/include/
  NO_DEFAULT_PATH
)

FIND_LIBRARY(OPENCL_LIBRARY libOpenCL.so
  /usr/lib64
  /opt/AMDAPP/lib/x86_64
  /usr/lib
  /usr/local/lib
  NO_DEFAULT_PATH
)

IF (OPENCL_LIBRARY AND OPENCL_INCLUDE_DIR)
  SET(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
  SET(OPENCL_FOUND "YES")
ELSE (OPENCL_LIBRARY AND OPENCL_INCLUDE_DIR)
  SET(OPENCL_FOUND "NO")
  MESSAGE(STATUS "OPENCL not found.")
ENDIF (OPENCL_LIBRARY AND OPENCL_INCLUDE_DIR)

IF (OPENCL_FOUND)
  IF (NOT OPENCL_FIND_QUIETLY)
    MESSAGE(STATUS "Found OpenCL: ${OPENCL_LIBRARIES}")
  ENDIF (NOT OPENCL_FIND_QUIETLY)
ELSE (OPENCL_FOUND)
  IF (OPENCL_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find OpenCL library")
  ENDIF (OPENCL_FIND_REQUIRED)
ENDIF (OPENCL_FOUND)

# Deprecated declarations.
GET_FILENAME_COMPONENT (NATIVE_OPENCL_LIB_PATH ${OPENCL_LIBRARY} PATH)

MARK_AS_ADVANCED(
  OPENCL_LIBRARY
  OPENCL_INCLUDE_DIR
)

