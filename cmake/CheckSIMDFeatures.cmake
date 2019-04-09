include(CheckCXXSourceRuns)

function(check_marchnative result)
set(CMAKE_REQUIRED_FLAGS "-march=native")

check_cxx_source_runs("
	int main()
	{
		return 0;
	}" FIND_MARCHNATIVE)
set(${result}  ${FIND_MARCHNATIVE} PARENT_SCOPE)
endfunction()


function(check_avx2 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-mavx2")
elseif(NOT CMAKE_CL_64)
	set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
endif()

check_cxx_source_runs("
	#include <immintrin.h>
	int main()
	{
		__m256i a = _mm256_set_epi32 (-1, 2, -3, 4, -1, 2, -3, 4);
		__m256i result = _mm256_abs_epi32 (a);
		return 0;
	}" FIND_AVX_20)

set(${result} ${FIND_AVX_20} PARENT_SCOPE)
endfunction()


function(check_avx1 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-mavx")
elseif(NOT CMAKE_CL_64)
	set(CMAKE_REQUIRED_FLAGS "/arch:AVX")
endif()

check_cxx_source_runs("
	#include <immintrin.h>
	int main()
	{
		__m256 a = _mm256_set_ps (-1.0f, 2.0f, -3.0f, 4.0f, -1.0f, 2.0f, -3.0f, 4.0f);
		__m256 b = _mm256_set_ps (1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
		__m256 result = _mm256_add_ps (a, b);
		return 0;
	}" FIND_AVX_10)

set(${result} ${FIND_AVX_10} PARENT_SCOPE)
endfunction()



function(check_sse42 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-msse4.2")
endif()

check_cxx_source_runs("
	#include <emmintrin.h>
	#include <nmmintrin.h>
	int main ()
	{
		long long a[2] = {  1, 2 };
		long long b[2] = { -1, 3 };
		long long c[2];
		__m128i va = _mm_loadu_si128 ((__m128i*)a);
		__m128i vb = _mm_loadu_si128 ((__m128i*)b);
		__m128i vc = _mm_cmpgt_epi64 (va, vb);

		_mm_storeu_si128 ((__m128i*)c, vc);
		if (c[0] == -1LL && c[1] == 0LL)
			return (0);
		else
			return (1);
	}"
	HAVE_SSE4_2_EXTENSIONS)

set(${result} ${HAVE_SSE4_2_EXTENSIONS} PARENT_SCOPE)
endfunction()



function(check_sse41 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-msse4.1")
endif()

check_cxx_source_runs("
	#include <smmintrin.h>
	int main ()
	{
		__m128 a, b;
		float vals[4] = {1, 2, 3, 4};
		const int mask = 123;
		a = _mm_loadu_ps (vals);
		b = a;
		b = _mm_dp_ps (a, a, mask);
		_mm_storeu_ps (vals,b);
		return (0);
	}"
	HAVE_SSE4_1_EXTENSIONS)

set(${result} ${HAVE_SSE4_1_EXTENSIONS} PARENT_SCOPE)
endfunction()



function(check_sse3 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-msse3")
endif()

check_cxx_source_runs("
	#include <pmmintrin.h>
	int main ()
	{
		volatile __m128d a, b;
		double vals[2] = {0};
		a = _mm_loadu_pd (vals);
		b = _mm_hadd_pd (a,a);
		_mm_storeu_pd (vals, b);
		return (0);
	}"
	HAVE_SSE3_EXTENSIONS)

set(${result} ${HAVE_SSE3_EXTENSIONS} PARENT_SCOPE)
endfunction()



function(check_sse2 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-msse2")
elseif(NOT CMAKE_CL_64)
	set(CMAKE_REQUIRED_FLAGS "/arch:SSE2")
endif()

check_cxx_source_runs("
	#include <emmintrin.h>
	int main ()
	{
		__m128d a, b;
		double vals[2] = {0};
		a = _mm_loadu_pd (vals);
		b = _mm_add_pd (a,a);
		_mm_storeu_pd (vals,b);
		return (0);
	}"
	HAVE_SSE2_EXTENSIONS)

set(${result} ${HAVE_SSE2_EXTENSIONS} PARENT_SCOPE)
endfunction()



function(check_sse1 result)
set(CMAKE_REQUIRED_FLAGS)

if(NOT MSVC)
	set(CMAKE_REQUIRED_FLAGS "-msse")
elseif(NOT CMAKE_CL_64)
	set(CMAKE_REQUIRED_FLAGS "/arch:SSE")
endif()

check_cxx_source_runs("
	#include <xmmintrin.h>
	int main ()
	{
		volatile __m128 a, b;
		float vals[4] = {0};
		a = _mm_loadu_ps (vals);
		b = a;
		b = _mm_add_ps (a,b);
		_mm_storeu_ps (vals,b);
		return (0);
	}"
	HAVE_SSE_EXTENSIONS)

set(${result} ${HAVE_SSE_EXTENSIONS} PARENT_SCOPE)
endfunction()


macro(CGOGN_CHECK_FOR_SSE)
	set(CGOGN_SSE_FLAGS)
	set(CGOGN_SSE_DEFINITIONS)
	set(CMAKE_REQUIRED_FLAGS)

	# check_marchnative(HAVE_MARCH_NATIVE)
	check_avx2(HAVE_AVX2_EXTENSIONS)
	check_avx1(HAVE_AVX1_EXTENSIONS)
	check_sse42(HAVE_SSE4_2_EXTENSIONS)
	check_sse41(HAVE_SSE4_1_EXTENSIONS)
	check_sse3(HAVE_SSE3_EXTENSIONS)
	check_sse2(HAVE_SSE2_EXTENSIONS)
	check_sse1(HAVE_SSE_EXTENSIONS)
	if(NOT MSVC)
		# if (HAVE_MARCH_NATIVE)
		# 	message(STATUS "-march=native flag support detected.")
		# 	set(CGOGN_SSE_FLAGS "${CGOGN_SSE_FLAGS} -march=native")
		# endif()

		if (HAVE_AVX2_EXTENSIONS)
			message(STATUS "avx2 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -mavx2 -mfpmath=sse)
		elseif (HAVE_AVX1_EXTENSIONS)
			message(STATUS "avx support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -mavx -mfpmath=sse)
		elseif(HAVE_SSE4_2_EXTENSIONS)
			message(STATUS "sse4.2 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -msse4.2 -mfpmath=sse)
		elseif(HAVE_SSE4_1_EXTENSIONS)
			message(STATUS "sse4.1 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -msse4.1 -mfpmath=sse)
		elseif(HAVE_SSSE3_EXTENSIONS)
			message(STATUS "sse3 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -mssse3 -mfpmath=sse)
		elseif(HAVE_SSE2_EXTENSIONS)
			message(STATUS "sse2 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -msse2 -mfpmath=sse)
		elseif(HAVE_SSE_EXTENSIONS)
			message(STATUS "sse support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} -msse -mfpmath=sse)
		else()
			# Setting -ffloat-store to alleviate 32bit vs 64bit discrepancies on non-SSE platforms.
			set(CGOGN_SSE_FLAGS -ffloat-store)
		endif()
	else(NOT MSVC)
		if (HAVE_AVX2_EXTENSIONS)
			message(STATUS "avx2 support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} /arch:AVX2)
		elseif (HAVE_AVX1_EXTENSIONS)
			message(STATUS "avx support detected.")
			set(CGOGN_SSE_FLAGS ${CGOGN_SSE_FLAGS} /arch:AVX)
		endif()
	endif()
endmacro()
