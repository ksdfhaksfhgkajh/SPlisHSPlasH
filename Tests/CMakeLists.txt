set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(Kernel)
add_subdirectory(ReadWriteState)
if (FOUND_SIMD)
	add_subdirectory(AVX)
endif()
