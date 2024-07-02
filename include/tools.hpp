#pragma once

#include "common.hpp"
#include <algorithm>
#include <iomanip>
#include <ostream>

namespace cryptotools
{
    template<size_t bits>
    using Sbox = std::array<uint<bits>, exp(bits)>;

    using ssize_t = ptrdiff_t;

    template<typename T, size_t r, size_t c = r>
    class Table : public std::array<T, r * c>
    {
    private:
        using super = std::array<T, r * c>;

    public:
        using super::super;
        using super::operator[];

        Table(const std::array<T, r * c> &arr) : super{arr} {}
        Table(std::array<T, r * c> &&arr) : super{arr} {}

        T &operator[](size_t i, size_t j) noexcept { return (*this)[i * c + j]; }

        const T &operator[](size_t i, size_t j) const noexcept { return (*this)[i * c + j]; }

        constexpr size_t rows() const { return r; }
        constexpr size_t cols() const { return c; }

        friend std::ostream &operator<<(std::ostream &os, const Table &tab)
        {
            size_t w = std::ceil(std::log10(std::ranges::max(tab))) + (std::ranges::min(tab) < 0);

            if constexpr (std::is_floating_point_v<T>)
                w += 5;

            for (size_t i = 0; i < tab.rows(); ++i)
            {
                os << "[ ";
                for (size_t j = 0; j < tab.cols(); ++j)
                    os << std::setw(w) << std::setprecision(4)  << tab[i, j] << ' ';
                os << "]\n";
            }

            return os;
        }
    };


    template<typename T, size_t r, size_t c = r>
    class Matrix : public Table<T, r, c>
    {
    private:
        using super = Table<T, r, c>;

    public:
        using super::super;
        using super::operator[];
    };

    template<size_t bits>
    constexpr Table<size_t, exp(bits)> difference_distribution_table(const Sbox<bits> &sbox)
    {
        static constexpr size_t N = exp(bits);

        Table<size_t, N> ddt{};

        for (uint<bits> x{};;)
        {
            for (uint<bits> y{};;)
            {
                ++ddt[x ^ y, sbox[x] ^ sbox[y]];
                if (!++y)
                    break;
            }

            if (!++x)
                break;
        }

        return ddt;
    }

    template<size_t bits>
    constexpr size_t differential_uniformity(const Sbox<bits> &sbox)
    {
        auto ddt{difference_distribution_table(sbox)};

        ddt[0] = 0;

        return std::ranges::max(ddt);
    }

    template<size_t bits>
    constexpr size_t differential_branch_number(const Sbox<bits> &sbox)
    {
        size_t b = std::numeric_limits<size_t>::max();

        for (uint<bits> x{};;)
        {
            for (uint<bits> y{x + 1}; y != x; ++y)
            {
                size_t w = hamming_weight(x ^ y) + hamming_weight(sbox[x] ^ sbox[y]);

                if (w < b)
                    b = w;

                if (++y == x)
                    break;
            }
            if (!++x)
                break;
        }

        return b;
    }

    template<size_t bits>
    constexpr Table<ssize_t, exp(bits)> linear_approximation_table(const Sbox<bits> &sbox)
    {
        static constexpr size_t N = exp(bits);

        Table<ssize_t, N> lat;

        lat.fill(-ssize_t{N >> 1});

        for (uint<bits> x{};;)
        {
            for (uint<bits> mx{};;)
            {
                for (uint<bits> my{};;)
                {
                    lat[mx, my] += !hamming_parity(x & mx ^ sbox[x] & my);
                    if (!++my)
                        break;
                }

                if (!++mx)
                    break;
            }

            if (!++x)
                break;
        }

        return lat;
    }

    template<size_t bits>
    constexpr ssize_t linear_bias(const Sbox<bits> &sbox)
    {
        auto lat{linear_approximation_table(sbox)};

        lat[0] = 0;
        auto [min, max] = std::ranges::minmax(lat);

        return std::abs(min) > std::abs(max) ? min : max;
    }

    template<size_t bits>
    constexpr Table<double, exp(bits)> linear_bias_table(const Sbox<bits> &sbox)
    {
        auto lat{linear_approximation_table(sbox)};
        Table<double, exp(bits)> lbt;

        std::ranges::transform(lat, lbt.begin(),
                               [](ssize_t x) { return double{x} / (exp(bits)); });

        return lbt;
    }
} // namespace cryptotools
