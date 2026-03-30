// Generic SIMD Vector
#ifndef GSVECTOR_H_
#define GSVECTOR_H_

#include <Kokkos_Core.hpp>
#include <synergia/foundation/trigon_traits.h>

// no simd when build for CUDA
#ifdef SYNERGIA_ENABLE_CUDA
#undef GSV_SSE
#undef GSV_AVX
#undef GSV_V4D
#undef GSV_MIC
#endif

// helper
namespace detail {
  template <class T, class E = void>
  struct VectorHelper {
    KOKKOS_INLINE_FUNCTION
    static constexpr int
    size()
    {
      return 1;
    }

    KOKKOS_INLINE_FUNCTION
    static T
    ld(const double* p)
    {
      return *p;
    }

    KOKKOS_INLINE_FUNCTION
    static void
    st(double* p, const T& v)
    {
      *p = v;
    }
  };
} // namespace detail

// expression class
template <typename E, class T>
struct VecExpr {
  typedef T vec_t;

  KOKKOS_INLINE_FUNCTION
  vec_t
  cal() const
  {
    return static_cast<E const&>(*this).cal();
  }

  KOKKOS_INLINE_FUNCTION
  operator E&() { return static_cast<E&>(*this); }

  KOKKOS_INLINE_FUNCTION
  operator E const&() const { return reinterpret_cast<const E&>(*this); }
  //{ return static_cast<const E &>(*this); }
};

// the vector wrapper base class
template <class T>
struct GSVec : public VecExpr<GSVec<T>, T> {
  using data_t = T;

  T data;

  KOKKOS_INLINE_FUNCTION
  static constexpr int
  size()
  {
    return detail::VectorHelper<T>::size();
  }

  template <typename U = T>
  KOKKOS_INLINE_FUNCTION
  GSVec(const T* t, typename std::enable_if<is_trigon<U>::value>::type* = 0)
    : data(*t)
  {}

  KOKKOS_INLINE_FUNCTION
  GSVec(const double d) : data(d) {}

  KOKKOS_INLINE_FUNCTION
  GSVec(const double* p) : data(detail::VectorHelper<T>::ld(p)) {}

  KOKKOS_INLINE_FUNCTION
  void
  load(const double* p)
  {
    data = detail::VectorHelper<T>::ld(p);
  }

  KOKKOS_INLINE_FUNCTION
  void
  store(double* p) const
  {
    detail::VectorHelper<T>::st(p, data);
  }

  template <typename U = T>
  KOKKOS_INLINE_FUNCTION void
  store(T* p, typename std::enable_if<is_trigon<U>::value>::type* = 0) const
  {
    *p = data;
  }

  KOKKOS_INLINE_FUNCTION
  T&
  cal()
  {
    return data;
  }

  KOKKOS_INLINE_FUNCTION
  T
  cal() const
  {
    return data;
  }

  template <typename E>
  KOKKOS_INLINE_FUNCTION
  GSVec(VecExpr<E, T> const& vec)
  {
    E const& v = vec;
    data = v.cal();
  }

#if 0
  KOKKOS_INLINE_FUNCTION
    operator double() const
    { return data; }
#endif
};

template <class T>
KOKKOS_INLINE_FUNCTION bool
operator==(GSVec<T> const& lhs, double rhs)
{
  return lhs.data == rhs;
}

template <class T>
KOKKOS_INLINE_FUNCTION bool
operator<(GSVec<T> const& lhs, double rhs)
{
  return lhs.data < rhs;
}

template <class T>
KOKKOS_INLINE_FUNCTION bool
operator>(GSVec<T> const& lhs, double rhs)
{
  return lhs.data > rhs;
}

template <class T>
KOKKOS_INLINE_FUNCTION bool
operator<=(GSVec<T> const& lhs, double rhs)
{
  return lhs.data <= rhs;
}

template <class T>
KOKKOS_INLINE_FUNCTION bool
operator>=(GSVec<T> const& lhs, double rhs)
{
  return lhs.data >= rhs;
}

// expression classes
template <typename E1, typename E2, class T, class E = void>
struct VecAdd : public VecExpr<VecAdd<E1, E2, T>, T> {
  E1 const& _u;
  E2 const& _v;

  KOKKOS_INLINE_FUNCTION
  VecAdd(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v) : _u(u), _v(v) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecAdd<E1, E2, T>, T>::vec_t
  cal() const
  {
    return _u.cal() + _v.cal();
  }
};

template <typename E1, typename E2, class T>
struct VecSub : public VecExpr<VecSub<E1, E2, T>, T> {
  E1 const& _u;
  E2 const& _v;

  KOKKOS_INLINE_FUNCTION
  VecSub(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v) : _u(u), _v(v) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecSub<E1, E2, T>, T>::vec_t
  cal() const
  {
    return _u.cal() - _v.cal();
  }
};

template <typename E1, typename E2, class T>
struct VecMul : public VecExpr<VecMul<E1, E2, T>, T> {
  E1 const& _u;
  E2 const& _v;

  KOKKOS_INLINE_FUNCTION
  VecMul(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v) : _u(u), _v(v) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecMul<E1, E2, T>, T>::vec_t
  cal() const
  {
    return _u.cal() * _v.cal();
  }
};

template <typename E1, typename E2, class T>
struct VecDiv : public VecExpr<VecDiv<E1, E2, T>, T> {
  E1 const& _u;
  E2 const& _v;

  KOKKOS_INLINE_FUNCTION
  VecDiv(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v) : _u(u), _v(v) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecDiv<E1, E2, T>, T>::vec_t
  cal() const
  {
    return _u.cal() / _v.cal();
  }
};

template <typename E, class T>
struct VecNeg : public VecExpr<VecNeg<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecNeg(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecNeg<E, T>, T>::vec_t
  cal() const
  {
    return -(_u.cal());
  }
};

template <typename E, class T>
struct VecSqrt : public VecExpr<VecSqrt<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecSqrt(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecSqrt<E, T>, T>::vec_t
  cal() const
  {
    return sqrt(_u.cal());
  }
};

template <typename E, class T>
struct GSVecLog : public VecExpr<GSVecLog<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  GSVecLog(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<GSVecLog<E, T>, T>::vec_t
  cal() const
  {
    return log(_u.cal());
  }
};

template <typename E, class T>
struct GSVecExp : public VecExpr<GSVecExp<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  GSVecExp(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<GSVecExp<E, T>, T>::vec_t
  cal() const
  {
    return exp(_u.cal());
  }
};

template <typename E, class T>
struct VecSin : public VecExpr<VecSin<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecSin(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecSin<E, T>, T>::vec_t
  cal() const
  {
    return sin(_u.cal());
  }
};

template <typename E, class T>
struct VecCos : public VecExpr<VecCos<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecCos(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecCos<E, T>, T>::vec_t
  cal() const
  {
    return cos(_u.cal());
  }
};

template <typename E, class T>
struct VecAsin : public VecExpr<VecAsin<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecAsin(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecAsin<E, T>, T>::vec_t
  cal() const
  {
    return asin(_u.cal());
  }
};

template <typename E, class T>
struct VecAtan : public VecExpr<VecAtan<E, T>, T> {
  E const& _u;

  KOKKOS_INLINE_FUNCTION
  VecAtan(VecExpr<E, T> const& u) : _u(u) {}

  KOKKOS_INLINE_FUNCTION
  typename VecExpr<VecAtan<E, T>, T>::vec_t
  cal() const
  {
    return atan(_u.cal());
  }
};

// overload operators
//

template <typename E1, typename E2, class T>
KOKKOS_INLINE_FUNCTION VecAdd<E1, E2, T> const
operator+(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v)
{
  return VecAdd<E1, E2, T>(u, v);
}

template <typename E1, typename E2, class T>
KOKKOS_INLINE_FUNCTION VecSub<E1, E2, T> const
operator-(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v)
{
  return VecSub<E1, E2, T>(u, v);
}

template <typename E1, typename E2, class T>
KOKKOS_INLINE_FUNCTION VecMul<E1, E2, T> const
operator*(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v)
{
  return VecMul<E1, E2, T>(u, v);
}

template <typename E1, typename E2, class T>
KOKKOS_INLINE_FUNCTION VecDiv<E1, E2, T> const
operator/(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v)
{
  return VecDiv<E1, E2, T>(u, v);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecNeg<E, T> const
operator-(VecExpr<E, T> const& u)
{
  return VecNeg<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecSqrt<E, T> const
sqrt(VecExpr<E, T> const& u)
{
  return VecSqrt<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION GSVecLog<E, T> const
log(VecExpr<E, T> const& u)
{
  return GSVecLog<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION GSVecExp<E, T> const
exp(VecExpr<E, T> const& u)
{
  return GSVecExp<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecSin<E, T> const
sin(VecExpr<E, T> const& u)
{
  return VecSin<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecCos<E, T> const
cos(VecExpr<E, T> const& u)
{
  return VecCos<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecAsin<E, T> const
asin(VecExpr<E, T> const& u)
{
  return VecAsin<E, T>(u);
}

template <typename E, class T>
KOKKOS_INLINE_FUNCTION VecAtan<E, T> const
atan(VecExpr<E, T> const& u)
{
  return VecAtan<E, T>(u);
}

// specialization for different platforms

class Vec2d;
class Vec4d;
class Vec8d;
class vector4double;

// headers
#if defined(GSV_SSE) || defined(GSV_AVX) || defined(GSV_AVX512)

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wattributes"
#endif

#include "vectorclass/vectorclass.h"
#include "vectorclass/vectormath_exp.h"
#include "vectorclass/vectormath_trig.h"

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#elif defined(GSV_QPX)
#include <mass_simd.h>
#endif

namespace detail {
  // specialization of helper class
  template <class T>
  struct VectorHelper<
    T,
    typename std::enable_if<std::is_same<T, Vec2d>::value>::type> {
    KOKKOS_INLINE_FUNCTION
    static constexpr int
    size()
    {
      return 2;
    }

    KOKKOS_INLINE_FUNCTION
    static T
    ld(const double* p)
    {
      T t;
      t.load_a(p);
      return t;
    }

    KOKKOS_INLINE_FUNCTION
    static void
    st(double* p, const T& v)
    {
      v.store_a(p);
    }
  };

  template <class T>
  struct VectorHelper<
    T,
    typename std::enable_if<std::is_same<T, Vec4d>::value>::type> {
    KOKKOS_INLINE_FUNCTION
    static constexpr int
    size()
    {
      return 4;
    }

    KOKKOS_INLINE_FUNCTION
    static T
    ld(const double* p)
    {
      T t;
      t.load_a(p);
      return t;
    }

    KOKKOS_INLINE_FUNCTION
    static void
    st(double* p, const T& v)
    {
      v.store_a(p);
    }
  };

  template <class T>
  struct VectorHelper<
    T,
    typename std::enable_if<std::is_same<T, Vec8d>::value>::type> {
    KOKKOS_INLINE_FUNCTION
    static constexpr int
    size()
    {
      return 8;
    }

    KOKKOS_INLINE_FUNCTION
    static T
    ld(const double* p)
    {
      T t;
      t.load_a(p);
      return t;
    }

    KOKKOS_INLINE_FUNCTION
    static void
    st(double* p, const T& v)
    {
      v.store_a(p);
    }
  };
} // namespace detail

// operations
template <typename E1, typename E2, class T>
struct VecAdd<
  E1,
  E2,
  T,
  typename std::enable_if<std::is_same<T, vector4double>::value>::type>
  : public VecExpr<VecAdd<E1, E2, T>, T> {
  E1 const& _u;
  E2 const& _v;
  VecAdd(VecExpr<E1, T> const& u, VecExpr<E2, T> const& v) : _u(u), _v(v) {}
  typename VecExpr<VecAdd<E1, E2, T>, T>::vec_t
  cal() const
  {
    return vec_add(_u.cal(), _v.cal());
  }
};

// stream operator
template <class T>
inline std::ostream&
operator<<(std::ostream& out, GSVec<T> const& v)
{
  out << "(" << v.data << ")";
  return out;
}

template <class T>
inline std::enable_if_t<std::is_same<T, Vec2d>::value, std::ostream&>
operator<<(std::ostream& out, T const& v)
{
  out << v[0] << ", " << v[1];
  return out;
}

template <class T>
inline std::enable_if_t<std::is_same<T, Vec4d>::value, std::ostream&>
operator<<(std::ostream& out, T const& v)
{
  out << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3];
  return out;
}

template <class T>
inline std::enable_if_t<std::is_same<T, Vec8d>::value, std::ostream&>
operator<<(std::ostream& out, T const& v)
{
  out << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << ", " << v[4]
      << ", " << v[5] << ", " << v[6] << ", " << v[7];
  return out;
}

// define the GSVector type
#if defined(GSV_SSE)
typedef GSVec<Vec2d> GSVector;
#elif defined(GSV_AVX)
typedef GSVec<Vec4d> GSVector;
#elif defined(GSV_AVX512)
typedef GSVec<Vec8d> GSVector;
#elif defined(GSV_QPX)
typedef GSVec<vector4double> GSVector;
#else
typedef GSVec<double> GSVector;
#endif

#endif
