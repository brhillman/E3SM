/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONFIG_HPP
#define HOMMEXX_CONFIG_HPP

#ifdef HOMMEXX_CONFIG_IS_CMAKE
# include "Hommexx_config.h"
# ifdef HAVE_CONFIG_H
#  include "config.h.c"
# endif
#elif !defined(HOMMEXX_VECTOR_SIZE) 
// Establish a good candidate vector size for eam builds
# ifdef HOMMEXX_ENABLE_GPU
#  define HOMMEXX_VECTOR_SIZE 1
# else
#  define HOMMEXX_VECTOR_SIZE 8
# endif
#endif

#if ! defined HOMMEXX_CUDA_SPACE && ! defined HOMMEXX_OPENMP_SPACE && ! defined HOMMEXX_THREADS_SPACE && ! defined HOMMEXX_SERIAL_SPACE && ! defined HOMMEXX_HIP_SPACE && ! defined HOMMEXX_SYCL_SPACE
# define HOMMEXX_DEFAULT_SPACE
#endif

#ifndef HOMMEXX_MPI_ON_DEVICE
# define HOMMEXX_MPI_ON_DEVICE 1
#endif

#include <Kokkos_Core.hpp>

#ifdef HOMMEXX_ENABLE_GPU 
# ifndef HOMMEXX_CUDA_MIN_WARP_PER_TEAM
#  define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 8
# endif
# ifndef HOMMEXX_CUDA_MAX_WARP_PER_TEAM
#  define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 16
# endif
#elif !defined(HOMMEXX_CONFIG_IS_CMAKE)
# define HOMMEXX_CUDA_MIN_WARP_PER_TEAM 1
# define HOMMEXX_CUDA_MAX_WARP_PER_TEAM 1
#endif

#if defined KOKKOS_COMPILER_GNU
// See https://github.com/kokkos/kokkos-kernels/issues/129
# define ConstExceptGnu
#else
# define ConstExceptGnu const
#endif

#endif // HOMMEXX_CONFIG_HPP
