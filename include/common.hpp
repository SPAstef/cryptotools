#pragma once

#include <array>
#include <bit>
#include <cinttypes>
#include <limits>

namespace cryptotools
{
    template<size_t bits>
    using uint = unsigned _BitInt(bits);

    template<size_t bits>
    constexpr size_t bitsof(uint<bits>)
    {
        return bits;
    }

    template<typename T = size_t>
    constexpr T exp(size_t bits)
    {
        return T{1} << bits;
    }

    template<size_t bits>
    constexpr size_t hamming_weight(uint<bits> x)
    {
        if constexpr (bits <= 64)
            return std::popcount(static_cast<uint64_t>(x));
        else
        {
            size_t w = 0;

            do
                w += std::popcount(static_cast<uint64_t>(x));
            while (x >>= 64);

            return w;
        }
    }

    template<size_t bits>
    constexpr bool hamming_parity(uint<bits> x)
    {
        return hamming_weight(x) & 1;
    }

    template<size_t n, typename T>
    constexpr std::array<T, n> make_array(T val)
    {
        return []<size_t... Is>(T val, std::index_sequence<Is...>) -> std::array<T, n>
        { return {{(static_cast<void>(Is), val)...}}; }(val, std::make_index_sequence<n>());
    }

} // namespace cryptotools
