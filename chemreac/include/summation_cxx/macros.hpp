#pragma once // -*- eval: (read-only-mode); -*-
#line 2 "/work/summation_cxx/macros.hpp"
#ifndef SMMTNCXX_RESTRICT
  #if defined(__GNUC__)
    #define SMMTNCXX_RESTRICT __restrict__
  #elif defined(_MSC_VER) && _MSC_VER >= 1400
    #define SMMTNCXX_RESTRICT __restrict
  // #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
  //   #define SMMTNCXX_RESTRICT restrict
  #else
    #define SMMTNCXX_RESTRICT
  #endif
#endif

#ifndef SMMTNCXX_PREFER_INLINE
  #if defined(__GNUC__)
    #define SMMTNCXX_PREFER_INLINE __attribute__((flatten))
  #elif defined(_MSC_VER) && _MSC_VER >= 1400
    #define SMMTNCXX_PREFER_INLINE __forceinline
  #else
    #define SMMTNCXX_PREFER_INLINE
  #endif
#endif

// Math macros to support e.g. __float128 without std lib support:
#ifndef SMMTNCXX_ABS
#define SMMTNCXX_ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef SMMTNCXX_NEUMAIER_BRANCH
// see test/bench.cpp
#define SMMTNCXX_NEUMAIER_BRANCH 1
#endif

// #define SMMTNCXX_FRIEND_OPERATORS(CLS, RET, TYP)
//     friend RET operator+(const CLS&, const TYP&);
//     friend RET operator+(const TYP&, const CLS&);

// #define SMMTNCXX_OPERATORS_IMPL(CLS, RET, TYP)
//     friend RET operator+(const CLS&, const TYP&);
//     friend RET operator+(const TYP&, const CLS&);
