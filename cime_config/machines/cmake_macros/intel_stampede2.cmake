string(APPEND CMAKE_C_FLAGS " -xCORE-AVX2")
string(APPEND CPPDEFS " -DLINUX")
if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
string(APPEND CPPDEFS " -DARCH_MIC_KNL")
string(APPEND CMAKE_Fortran_FLAGS " -fp-model consistent -fimf-use-svml")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -qno-opt-dynamic-align")
string(APPEND CMAKE_Fortran_FLAGS " -xCORE-AVX2")
string(APPEND CMAKE_EXE_LINKER_FLAGS " -lpthread")