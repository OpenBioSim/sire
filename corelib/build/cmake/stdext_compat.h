#pragma once
// stdext::make_checked_array_iterator and make_unchecked_array_iterator were
// removed in MSVC 14.38 (VS 2022 17.8). Qt5's qcompilerdetection.h maps
// QT_MAKE_CHECKED_ARRAY_ITERATOR and QT_MAKE_UNCHECKED_ARRAY_ITERATOR to these
// functions unconditionally on MSVC, so provide no-op stubs.
#if defined(__cplusplus) && defined(_MSC_VER) && _MSC_VER >= 1938
#include <cstddef>
namespace stdext
{
    template <typename T>
    T *make_checked_array_iterator(T *ptr, std::size_t, std::size_t index = 0)
    {
        return ptr + index;
    }

    template <typename T>
    T *make_unchecked_array_iterator(T *ptr)
    {
        return ptr;
    }
}
#endif
