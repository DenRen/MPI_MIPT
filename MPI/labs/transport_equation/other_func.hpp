#pragma once

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <random>
#include <array>

#include "print_lib.hpp"

template <typename T>
bool operator != (const std::vector <T>& lhs,
                  const std::vector <T>& rhs)
{
    const std::size_t size = lhs.size ();
    if (size != rhs.size ()) {
        return true;
    }

    for (std::size_t i = 0; i < size; ++i) {
        if (lhs[i] != rhs[i]) {
            return true;
        }
    }

    return false;
} // bool operator != (const std::vector <T>& lhs, const std::vector <T>& rhs)

template <typename T>
bool operator == (const std::vector <T>& lhs,
                  const std::vector <T>& rhs)
{
    const std::size_t size = lhs.size ();
    if (size != rhs.size ()) {
        return false;
    }

    for (std::size_t i = 0; i < size; ++i) {
        if (lhs[i] != rhs[i]) {
            return false;
        }
    }

    return true;
} // bool operator == (const std::vector <T>& lhs, const std::vector <T>& rhs)

template <typename T, typename Rand>
std::vector <T>
getRandFillVector (std::size_t size,
                   Rand& rand)
{
    std::vector <T> vec (size);
    for (auto& item : vec) {
        item = rand ();
    }

    return vec; // RVO
} // getRandFillVector (std::size_t size, Rand& rand)

template <typename T, typename Rand>
std::vector <T>
getUniqRandFillVector (std::size_t size,
                       T module,
                       Rand& rand)
{
    std::set <T> set;
    while (set.size () != size) {
        set.insert (rand () % module);
    }

    std::vector <T> res;
    res.reserve (size);
    for (auto&& item : set) {
        res.emplace_back (std::move (item));
    }

    return res; // RVO
} // getRandFillVector (std::size_t size, Rand& rand)s

template <typename T, typename Rand>
std::vector <T>
getRandFillVector (std::size_t size,
                   Rand& rand,
                   T module)
{
    std::vector <T> vec (size);
    for (auto& item : vec) {
        item = rand () % module;
    }

    return vec; // RVO
} // getRandFillVector (std::size_t size, Rand& rand, T module)

namespace seclib {
    class RandomGenerator {
        std::mt19937 rand_;
        std::mt19937_64 rand_64_;

        RandomGenerator (std::random_device::result_type seed) :
            rand_ (seed),
            rand_64_ (seed)
        {}

        template <typename F, typename... Args>
        std::vector <std::invoke_result_t <F, Args...>>
        getFilledVector (std::size_t size,
                         F&& func, Args&&... args)
        {
            std::vector <std::invoke_result_t <F, Args...>> vec (size);

            for (std::size_t i = 0; i < size; ++i) {
                vec[i] = func (std::forward <Args> (args)...);
            }

            return vec;
        }

        template <typename F, typename... Args>
        std::vector <std::invoke_result_t <F, Args...>>
        getUniqueFilledVector (std::size_t size,
                               F&& func, Args&&... args)
        {
            using T = std::invoke_result_t <F, Args...>;

            std::set <T> unqie_elems;
            while (unqie_elems.size () != size) {
                T value = std::invoke (std::forward <F> (func),
                                       std::forward <Args> (args)...);
                unqie_elems.emplace (std::move (value));
            }

            std::vector <T> vec;
            vec.reserve (size);
            std::move (std::begin (unqie_elems), std::end (unqie_elems),
                       std::back_inserter (vec));

            std::shuffle (std::begin (vec), std::end (vec), rand_);

            return vec;
        }

    public:
        RandomGenerator () :
            RandomGenerator (std::random_device {}())
        {}

        template <typename T>
        T
        get_rand_val () {
            if constexpr (sizeof (T) <= 4) {
                return rand_ ();
            } else {
                return rand_64_ ();
            }
        }

        template <typename T>
        T
        get_rand_val (T module) {
            return get_rand_val <T> () % module;
        }

        template <typename T>
        T
        get_rand_val (T min,
                      T max) {
            return min + get_rand_val <T> (max + 1 - min);
        }

        template <typename T>
        std::vector <T>
        get_vector (std::size_t size)
        {
            return getFilledVector (size, [this] () {
                return get_rand_val <T> ();
            });
        }

        template <typename T>
        std::vector <T>
        get_vector (std::size_t size,
                   T module)
        {
            return getFilledVector (size, [this, module] () {
                return get_rand_val <T> (module);
            });
        }

        template <typename T>
        std::vector <T>
        get_vector_uniq (std::size_t size)
        {
            return getUniqueFilledVector (size, [this] () {
                return get_rand_val <T> ();
            });
        }

        template <typename T>
        std::vector <T>
        get_vector_uniq (std::size_t size,
                         T module)
        {
            return getUniqueFilledVector (size, [this, module] () {
                return get_rand_val <T> (module);
            });
        }

        std::string
        get_string (std::size_t len,
                    char min,
                    char max)
        {
            if (min > max) {
                throw std::invalid_argument ("min > max");
            }

            std::string str;
            str.reserve (len);

            auto it = str.begin ();
            while (len >= 8) {
                uint64_t value = get_rand_val <uint64_t> ();

                for (int i = 0; i < 8; ++i) {
                    str += static_cast <char> (min + value % (max - min + 1));
                    value >>= 8;
                }

                len -= 8;
            }

            uint64_t value = get_rand_val <uint64_t> ();

            for (int i = 0; i < len; ++i) {
                str += static_cast <char> (value % (max - min + 1) + min);
                value >>= 8;
            }


            return str;
        }
    }; // class RandomGenerator

    class BitEnumerator {
        typedef uint64_t type;
        type mask = 0;

    public:
        BitEnumerator () = default;
        BitEnumerator (type mask) :
            mask (mask)
        {}

        BitEnumerator&
        operator = (type mask) {
            this->mask = mask;
            return *this;
        }

        bool
        get (int n) const {
            if (n >= sizeof (mask) * 8) {
                throw std::invalid_argument ("overbound request");
            }
            return (mask >> n) & 1;
        }

        void
        operator ++ () {
            ++mask;
        }

        void
        operator ++ (int) {
            ++(*this);
        }
    }; // class BitEnumerator
} // namespace seclib