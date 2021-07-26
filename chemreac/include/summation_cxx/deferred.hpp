// -*- eval: (read-only-mode); -*-
#pragma once
#include <cstddef>
#include <cstring> // std::memset
#include <memory> // std::make_unique
#include <summation_cxx/compensated.hpp>

namespace summation_cxx {
    namespace detail {
        template <typename T, typename Base>
        struct RefBase {
        private:
            T * ptr;
        protected:
            T& accum() { return ptr[0]; }
            T& carry() { return ptr[1]; }
        public:
            RefBase(T * data) : ptr(data) {}
            void operator=(T arg) {
                this->accum() = arg;
                this->carry() = 0;
            }
            void operator*=(T arg) {
                this->accum() *= arg;
                this->carry() *= arg;
            }
            void operator/=(T arg) {
                this->accum() /= arg;
                this->carry() /= arg;
            }
        };
    }
    template <typename T>
    struct KahanRef : detail::RefBase<T, KahanRef<T>> {
        void operator+=(T arg) {
            accum_kahan_destructive(this->accum(), this->carry(), arg);
        }
    };

    template <typename T>
    struct NeumaierRef : detail::RefBase<T, NeumaierRef<T>> {
        using detail::RefBase<T, NeumaierRef<T>>::operator=;
        void operator+=(T arg) {
            accum_neumaier(this->accum(), this->carry(), arg);
        }
    };

    namespace detail {
        template <typename T, typename Ref, typename Base>
        struct RangedAccumulatorBase {
        protected:
            T * tgt;
            std::unique_ptr<T[]> storage;
            std::size_t sz;
        public:
            RangedAccumulatorBase(std::size_t sz) : storage(std::make_unique<T[]>(sz*2)), sz(sz) {}
            void init(T * target) {
                tgt = target;
                std::memset(storage.get(), 0x00, sizeof(T)*sz*2);
            }
            Ref operator[](std::size_t idx) { return Ref{&storage[idx * 2]}; }
        };
    }
    template <typename T>
    struct RangedAccumulatorKahan : detail::RangedAccumulatorBase<T, KahanRef<T>, RangedAccumulatorKahan<T>> {
        using detail::RangedAccumulatorBase<T, KahanRef<T>, RangedAccumulatorKahan<T>>::RangedAccumulatorBase;
        void commit() {
            for (std::size_t i=0; i<this->sz; ++i) {
                this->tgt[i] = this->storage[i*2];
            }
        }
    };
    template <typename T>
    struct RangedAccumulatorNeumaier : detail::RangedAccumulatorBase<T, NeumaierRef<T>, RangedAccumulatorNeumaier<T>> {
        using detail::RangedAccumulatorBase<T, NeumaierRef<T>, RangedAccumulatorNeumaier<T>>::RangedAccumulatorBase;
        void commit() {
            for (std::size_t i=0; i<this->sz; ++i) {
                this->tgt[i] = this->storage[i*2] + this->storage[i*2 + 1];
            }
        }
    };

    namespace detail {
        template <typename T, typename Base>
        struct ScopedAccumulatorBase {
        private:
            T accum {0}, carry {0};
            T* tgt {nullptr};
        public:
            ScopedAccumulatorBase(T *target) : tgt(target) {}
            void clear() { accum = T{0}; carry = T{0}; }
            void operator*=(T arg) {
                accum *= arg;
                carry *= arg;
            }
            void operator/=(T arg) {
                accum /= arg;
                carry /= arg;
            }
        };
    }
    template <typename T>
    struct ScopedAccumulatorKahan : detail::ScopedAccumulatorBase<T, ScopedAccumulatorKahan<T>>{
        void operator+=(T arg) {
            accum_kahan_destructive(this->accum, this->carry, arg);
        }
        ~ScopedAccumulatorKahan(){
            *(this->tgt) = this->accum;
        }
    };
    template <typename T>
    struct ScopedAccumulatorNeumaier : detail::ScopedAccumulatorBase<T, ScopedAccumulatorNeumaier<T>>{
        void operator+=(T arg) {
            accum_neumaier(this->accum, this->carry, arg);
        }
        ~ScopedAccumulatorNeumaier(){
            *(this->tgt) = this->accum + this->carry;
        }
    };
}
