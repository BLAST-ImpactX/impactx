#ifndef SYNERGIA_UTILS_KOKKOS_TYPES_H
#define SYNERGIA_UTILS_KOKKOS_TYPES_H

#include <Kokkos_Core.hpp>

namespace kt {
  template <class T, size_t N>
  struct array_type {
    T data[N] = {};

    KOKKOS_INLINE_FUNCTION
    constexpr size_t
    size() const
    {
      return N;
    }

    KOKKOS_INLINE_FUNCTION
    T&
    operator[](size_t i)
    {
      return data[i];
    }

    KOKKOS_INLINE_FUNCTION
    T const&
    operator[](size_t i) const
    {
      return data[i];
    }

    KOKKOS_INLINE_FUNCTION // initialize data[] to 0
      void
      init()
    {
      for (int i = 0; i < N; i++) data[i] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    T*
    begin()
    {
      return data;
    }

    KOKKOS_INLINE_FUNCTION
    T const*
    begin() const
    {
      return data;
    }

    KOKKOS_INLINE_FUNCTION
    T*
    end()
    {
      return data + N;
    }

    KOKKOS_INLINE_FUNCTION
    T const*
    end() const
    {
      return data + N;
    }

    KOKKOS_INLINE_FUNCTION
    array_type&
    operator+=(const array_type& src)
    {
      for (int i = 0; i < N; i++) data[i] += src.data[i];
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    void
    operator+=(const volatile array_type& src) volatile
    {
      for (int i = 0; i < N; i++) data[i] += src.data[i];
    }
  };

  // alias
  template <class T, int N>
  using arr_t = array_type<T, N>;

  template <typename T>
  KOKKOS_INLINE_FUNCTION T
  qpow(T x, int i)
  {
    T retval{1};
    while (i--) { retval = retval * x; }
    return retval;
  }

  template <class T, int N>
  struct SumArray {
  public:
    // Required
    typedef SumArray reducer;
    typedef array_type<T, N> value_type;
    typedef Kokkos::View<value_type*, Kokkos::MemoryUnmanaged> result_view_type;

  private:
    value_type& value;

  public:
    KOKKOS_INLINE_FUNCTION
    SumArray(value_type& value_) : value(value_) {}

    // Required
    KOKKOS_INLINE_FUNCTION
    void
    join(value_type& dest, const value_type& src) const
    {
      dest += src;
    }

    KOKKOS_INLINE_FUNCTION
    void
    join(volatile value_type& dest, const volatile value_type& src) const
    {
      dest += src;
    }

    KOKKOS_INLINE_FUNCTION
    void
    init(value_type& val) const
    {
      val.init();
    }

    KOKKOS_INLINE_FUNCTION
    value_type&
    reference() const
    {
      return value;
    }

    KOKKOS_INLINE_FUNCTION
    result_view_type
    view() const
    {
      return result_view_type(&value);
    }

    KOKKOS_INLINE_FUNCTION
    bool
    references_scalar() const
    {
      return true;
    }
  };
}

#endif
