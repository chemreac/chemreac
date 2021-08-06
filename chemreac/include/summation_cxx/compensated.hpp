#pragma once // -*- eval: (read-only-mode); -*-
#line 2 "/work/summation_cxx/compensated.hpp"
#ifdef __FAST_MATH__
#error fast math enabled (/fp:fast, -ffast-math), this would negate compensation.
#endif
#include "summation_cxx/macros.hpp"
#include <cstddef> // std::size_t

namespace summation_cxx {
    enum class Compensation { NONE, KAHAN, NEUMAIER };

    namespace /* anonymous */ {
        template<typename T>
        SMMTNCXX_PREFER_INLINE void accum_kahan_destructive(
            T& SMMTNCXX_RESTRICT accu,
            T& SMMTNCXX_RESTRICT carry,
            T& SMMTNCXX_RESTRICT elem)
        {
            elem -= carry;
            const T tmp = accu + elem;
            carry = T{tmp - accu} - elem;
            accu = tmp;
        }
        template<typename T>
        SMMTNCXX_PREFER_INLINE void accum_kahan(
            T& SMMTNCXX_RESTRICT accu,
            T& SMMTNCXX_RESTRICT carry,
            const T& SMMTNCXX_RESTRICT elem)
        {
            T y = elem;
            accum_kahan_destructive(accu, carry, y);
        }

        template<typename T>
        SMMTNCXX_PREFER_INLINE void accum_neumaier(
            T& SMMTNCXX_RESTRICT acm,
            T& SMMTNCXX_RESTRICT carry,
            const T& SMMTNCXX_RESTRICT elem)
        {
            const T tmp = acm + elem;
#if SMMTNCXX_NEUMAIER_BRANCH == 1
            if (SMMTNCXX_ABS(tmp) > SMMTNCXX_ABS(elem)) {
              carry += T{acm - tmp} + elem;
            } else {
                carry += T{elem - tmp} + acm;
            }
#else
            T cases[2] = {T{elem - tmp} + acm, T{acm - tmp} + elem};
            carry += cases[SMMTNCXX_ABS(tmp) > SMMTNCXX_ABS(elem)];
#endif
            acm = tmp;
        }
    }
}
