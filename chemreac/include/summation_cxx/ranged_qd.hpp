#pragma once
#include <cstring> // std::memset
#include <memory> // make_unique
#include <qd/dd_real.h>

/// patch qd library with unary plus:
dd_real operator+(const dd_real& arg)
{
    return arg;
}

namespace summation_cxx {
static_assert(sizeof(dd_real) == 2 * sizeof(double));
struct RangedAccumulatorDD {
    typedef double underlying_type;
    typedef double target_type;
    typedef dd_real accumulator_type;
    typedef dd_real& view_type;
    typedef const dd_real& const_view_type;

private:
    target_type* tgt {};
    std::unique_ptr<dd_real[]> storage {};
    std::size_t sz {};
    bool cumulative {};

public:
    RangedAccumulatorDD() = default;
    explicit RangedAccumulatorDD(std::size_t sz)
        : storage(std::make_unique<dd_real[]>(sz))
        , sz(sz)
    {
    }
    void init(target_type* target, bool cumulative = false)
    {
        tgt = target;
        this->cumulative = cumulative;
        if (sz > 0 /* UB to call memset over zero bytes. */) {
            // doing this only makes sense if commit() is not always called.
            // std::fill(storage.get(), storage.get()+sz, 0);
            std::memset(reinterpret_cast<void*>(storage.get()), 0x00, sizeof(dd_real) * sz);
        }
    }
    view_type operator[](std::size_t idx)
    {
        return storage[idx];
    }
    const_view_type operator[](std::size_t idx) const
    {
        return storage[idx];
    }
    void commit() const
    {
#define SXX_LOOP(OP)                            \
    for (std::size_t i = 0; i < sz; ++i) {           \
        this->tgt[i] OP to_double(this->storage[i]); \
    }
        if (cumulative) {
            SXX_LOOP(+=)
        } else {
            SXX_LOOP(=)
        }
#undef SXX_LOOP
    }
};
}
