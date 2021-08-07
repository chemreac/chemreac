#pragma once
#include <array>
#include <cmath>
#include <cstring>
#include "summation_cxx/compensated.hpp"
#include <summation_cxx/impl.hpp>

namespace summation_cxx {

    template <typename T, Compensation scheme>
    struct Accumulator : public detail::Operators<T, scheme, Accumulator<T, scheme>> {
        static constexpr Compensation compensation_scheme { scheme };
    private:
        std::array<T, 2> data {};
    public:
        T& accum() { return data.data()[0]; }
        T& carry() { return data.data()[1]; }
        const T& accum() const { return data.data()[0]; }
        const T& carry() const { return data.data()[1]; }
    public:
        Accumulator() = default;
        Accumulator(T accum) {
            data[0] = accum;
        }
        explicit Accumulator(T accum, T carry) {
            data[0] = accum;
            data[1] = carry;
        }
        void clear() { data.clear(); }

        template <int n>
        static constexpr T sum(const std::array<T, n>& arr) {
            Accumulator ta {};
#if defined(__clang__)
#pragma unroll 16
#elif defined(__GNUC__)
#pragma GCC unroll 16
#endif
            for (const auto& e: arr) {
                ta += e;
            }
            return to<T>(ta);
        }
    };
    template<typename T> using AccumulatorKahan = Accumulator<T, Compensation::KAHAN>;
    template<typename T> using AccumulatorNeumaier = Accumulator<T, Compensation::NEUMAIER>;

    template<typename T, Compensation scheme>
    T pow(const Accumulator<T, scheme>& base, T exponent) {
        return std::pow(base.accum(), exponent);
    }
    template<typename T, Compensation scheme>
    T pow(const Accumulator<T, scheme>& base, int exponent) {
        return std::pow(base.accum(), static_cast<T>(exponent));
    }

}
