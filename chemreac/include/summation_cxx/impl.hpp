#pragma once
#include <summation_cxx/compensated.hpp>
#include <array>
#include <cassert>


namespace summation_cxx {
    template <typename T, summation_cxx::Compensation Scheme> struct Accumulator;
    template <typename T, summation_cxx::Compensation Scheme> struct AccuView;
    namespace detail {
        template <typename T, Compensation scheme, typename Derived> struct Operators;
    }
}
namespace summation_cxx::detail {
    template <typename T, Compensation scheme, typename Derived>
    struct Operators {
        typedef T underlying_type;
        typedef Accumulator<T, scheme> accumulator_type;
        typedef AccuView<T, scheme> view_type;

#define ACCUM(cv_qual) static_cast<cv_qual Derived *>(this)->accum()
#define CARRY(cv_qual) static_cast<cv_qual Derived *>(this)->carry()
        template <typename U>
        U to() const {
            if constexpr (Derived::compensation_scheme == Compensation::KAHAN) {
                return ACCUM(const);
            } else if constexpr (Derived::compensation_scheme == Compensation::NEUMAIER) {
                if constexpr (sizeof(T) > sizeof(U)) {
                    return ACCUM(const) + CARRY(const);
                } else {
                    return static_cast<U>(ACCUM(const)) + static_cast<U>(CARRY(const));
                }
            } else {
                assert(false);
            }
        }
        Derived& operator+=(T arg) {
          if constexpr (scheme == Compensation::KAHAN) {
            accum_kahan_destructive(ACCUM(), CARRY(), arg);
          } else if constexpr (scheme == Compensation::NEUMAIER) {
            accum_neumaier(ACCUM(), CARRY(), arg);
          } else {
            assert(false);
          }
          return *(static_cast<Derived *>(this));
        }
        Derived& operator-=(T arg) {
            Derived& self = *(static_cast<Derived *>(this));
            self += -arg;
            return self;
        }

        // template<std::enable_if<scheme == Compensation::KAHAN>>
        // Derived& operator+=(T&& arg) {
        //     accum_kahan_destructive(ACCUM(), CARRY(), arg);
        //     return *this;
        // }

        // template<std::enable_if<scheme == Compensation::NEUMAIER>>
        // Derived& operator+=(T&& arg) {
        //     accum_neumaier(ACCUM(), CARRY(), arg);
        //     return *this;
        // }
        // template<std::enable_if<scheme == Compensation::KAHAN>>
        // Derived& operator+=(const T& arg) {
        //     accum_kahan(ACCUM(), CARRY(), arg);
        //     return *this;
        // }

        // template<std::enable_if<scheme == Compensation::NEUMAIER>>
        // Derived& operator+=(const T& arg) {
        //     accum_neumaier(ACCUM(), CARRY(), arg);
        //     return *this;
        // }

        void operator=(const T& arg) {
            ACCUM() = arg;
            CARRY() = 0;
        }
        void operator/=(const T& arg) {
            ACCUM() /= arg;
            CARRY() /= arg;
        }
        void operator*=(const T& arg) {
            ACCUM() *= arg;
            CARRY() *= arg;
        }
        void operator+=(const accumulator_type& other) {
            *this += other.accum();
            CARRY() += other.carry();
        }
        void operator-=(const accumulator_type& other) {
            *this -= other.accum();
            CARRY() -= other.carry();
        }
        accumulator_type operator*(const T &arg) const {
            return accumulator_type(ACCUM(const)*arg, CARRY(const)*arg);
        }
        accumulator_type operator*(const accumulator_type& other) const {
            return accumulator_type(ACCUM(const)*other.accum(),
                                    CARRY(const)*other.accum() + ACCUM(const)*other.carry() + CARRY(const)*other.carry());
        }
        accumulator_type operator/(const accumulator_type& other) const {
            const T denom = other.template to<T>();
            accumulator_type result {ACCUM()/denom, CARRY()/denom};
        }
        accumulator_type operator+(const accumulator_type& other) const {
            return accumulator_type(ACCUM(const)+other.accum(), CARRY(const)+other.carry());
        }
        accumulator_type operator+() const {
            return accumulator_type(ACCUM(const), CARRY(const));
        }
        accumulator_type operator-() const {
            return accumulator_type(-ACCUM(const), -CARRY(const));
        }
    };
    template<typename T, typename Derived>
    typename Derived::accumulator_type operator*(const T& factor_a, const Derived& factor_b)
    {
        return factor_b*factor_a; // multiplication is commutative
    }

#undef ACCUM
#undef CARRY


}
