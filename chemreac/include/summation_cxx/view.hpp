#pragma once // -*- eval: (read-only-mode); -*-
#line 2 "/work/summation_cxx/view.hpp"
#include <summation_cxx/impl.hpp>

namespace summation_cxx {

    template <typename T, Compensation scheme>
    struct AccuView : //public summation_cxx::detail::AccuBase<T>,
                      public detail::Operators<T, scheme, AccuView<T, scheme>> {
        static constexpr Compensation compensation_scheme { scheme };
    private:
        T * ptr;
    public:
        T& accum() { return ptr[0]; }
        T& carry() { return ptr[1]; }
        const T& accum() const { return ptr[0]; }
        const T& carry() const { return ptr[1]; }
    public:
        AccuView() = delete;
        AccuView(T * data) : ptr(data) {
            assert(data);
        }
        Accumulator<T, scheme> deepcopy() {
            return Accumulator<T, scheme>{this->accum(), this->carry()};
        }
    };


    template<typename T> using AccuViewKahan = AccuView<T, Compensation::KAHAN>;
    template<typename T> using AccuViewNeumaier = AccuView<T, Compensation::NEUMAIER>;
}
