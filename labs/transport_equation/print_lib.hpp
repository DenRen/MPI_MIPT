#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <iterator>

template <typename T>
std::ostream&
operator << (std::ostream& os,
             const std::vector <T>& vec)
{
    const std::size_t size = vec.size ();
    if (size == 0) {
        return os;
    }

    for (std::size_t i = 0; i + 1 < size; ++i) {
        os << vec[i] << " ";
    }

    return os << vec[size - 1];
}

template <typename F, typename S>
std::ostream&
operator << (std::ostream& os,
             const std::pair <F, S>& pair)
{
    return os << pair.first << ": " << pair.second;
}

template <typename K, typename V>
std::ostream&
operator << (std::ostream& os,
             const std::map <K, V>& map)
{
    os << "{";

    const std::size_t size = map.size ();
    auto iter = map.cbegin ();
    for (std::size_t i = 0; i + 1 < size; ++i) {
        os << *iter++ << ", ";
    }

    return  os << *iter << "}";
}

template <typename... Args>
std::ostream&
print_all (std::ostream& os,
           const Args&... args)
{
    return (os << ... << args);
}

namespace detail {

template <typename T, typename SepT = std::string>
class AddSeparator {
    const T& ref_;
    const SepT& sep_;

public:
    AddSeparator (const T& ref, SepT sep = " "):
        ref_ (ref),
        sep_ (sep)
    {}

    std::ostream& print (std::ostream& os) const {
        return os << sep_ << ref_;
    }
};

} // namespace detail

template <typename T>
std::ostream&
operator << (std::ostream& os,
             const detail::AddSeparator <T>& val)
{
    return val.print (os);
}

template <typename Arg, typename... Args>
std::ostream&
print_all_sep (std::string sep,
               std::ostream& os,
               const Arg& arg,
               const Args&... args)
{
    return ((os << arg) << ... << detail::AddSeparator (args, sep));
}

template <typename Arg, typename... Args>
std::ostream&
print_all_sep (std::ostream& os,
               const Arg& arg,
               const Args&... args)
{
    return print_all_sep (" ", os, arg, args...);
}

template <typename Arg, typename... Args>
std::ostream&
print_all_sep (const Arg& arg,
               const Args&... args)
{
    return print_all_sep (std::cout, arg, args...);
}

namespace detail {

template <int N, typename... Args>
struct PrintTuple {
    template <typename OStreamT>
    void
    print (OStreamT& os,
           const std::tuple <Args...>& tuple)
    {
        PrintTuple <N-1, Args...> a;
        a.print (os, tuple);
        os << ", " << std::get <N-1> (tuple);
    }
};

template <typename... Args>
struct PrintTuple <1, Args...> {
    template <typename OStreamT>
    void
    print (OStreamT& os,
           const std::tuple <Args...>& tuple)
    {
        os << std::get <0> (tuple);
    }
};

}

template <typename... Args>
std::ostream&
operator << (std::ostream& os,
             const std::tuple <Args...>& tuple)
{
    detail::PrintTuple <sizeof... (Args), Args...> _print_tuple;

    os << "{";
    _print_tuple.print (os, tuple);
    os << "}";

    return os;
}

template <typename T>
std::ostream&
operator << (std::ostream& os,
             const std::set <T>& set)
{
    auto size = set.size ();
    auto iter = std::cbegin (set);

    os <<  "{";
    if  (size != 0) {    
        while (size-- != 1) {
            os << *iter++ << ", ";
        }
    }

    return os  << *iter <<  "}";
}

template <typename T, std::size_t N>
std::ostream&
operator << (std::ostream& os, const std::array <T, N>& arr) {
    if constexpr (N != 0) {
        os << arr[0];
        for (std::size_t i = 1; i < N; ++i) {
            os << ", " << arr[i];
        }
    }

    return os;
}