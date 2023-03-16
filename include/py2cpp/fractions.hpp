#pragma once

/** @file include/fractions.hpp
 *  This is a C++ Library header.
 */

// #include <boost/operators.hpp>
// #include <cmath>
#include <numeric>
#include <type_traits>
#include <utility>

// #include "common_concepts.h"

#if __cpp_constexpr >= 201304
#define CONSTEXPR14 constexpr
#else
#define CONSTEXPR14 inline
#endif

namespace fun {

/**
 * @brief absolute (unsigned)
 *
 * @tparam T
 * @param[in] a
 */
template <typename T>
CONSTEXPR14 auto abs(const T &a) ->
    typename std::enable_if<std::is_unsigned<T>::value, T>::type {
  return a;
}

/**
 * @brief absolute (signed)
 *
 * @tparam T
 * @param[in] a
 */
template <typename T>
CONSTEXPR14 auto abs(const T &a) ->
    typename std::enable_if<!std::is_unsigned<T>::value, T>::type {
  return (a < T(0)) ? -a : a;
}

/**
 * @brief Greatest common divider
 *
 * @tparam _Mn
 * @param[in] __m
 * @param[in] __n
 * @return _Mn
 */
template <typename _Mn>
CONSTEXPR14 auto gcd_recur(const _Mn &__m, const _Mn &__n) -> _Mn {
  if (__n == 0) {
    return abs(__m);
  }
  return gcd_recur(__n, __m % __n);
}

/**
 * @brief Greatest common divider
 *
 * @tparam _Mn
 * @param[in] __m
 * @param[in] __n
 * @return _Mn
 */
template <typename _Mn>
CONSTEXPR14 auto gcd(const _Mn &__m, const _Mn &__n) -> _Mn {
  if (__m == 0) {
    return abs(__n);
  }
  return gcd_recur(__m, __n);
}

/**
 * @brief Least common multiple
 *
 * @tparam _Mn
 * @param[in] __m
 * @param[in] __n
 * @return _Mn
 */
template <typename _Mn>
CONSTEXPR14 auto lcm(const _Mn &__m, const _Mn &__n) -> _Mn {
  if (__m == 0 || __n == 0) {
    return 0;
  }
  return (abs(__m) / gcd(__m, __n)) * abs(__n);
}

/**
 * @brief Fraction
 *
 * @tparam Z
 */
template <typename Z> struct Fraction {
  Z _num;
  Z _den;

  /**
   * @brief Construct a new Fraction object
   *
   * @param[in] num
   * @param[in] den
   */
  CONSTEXPR14 Fraction(Z num, Z den)
      : _num{std::move(num)}, _den{std::move(den)} {
    this->normalize();
  }

  /**
   * @brief normalize to a canonical form
   *
   * denominator is always non-negative and co-prime with numerator
   */
  CONSTEXPR14 auto normalize() -> Z {
    this->normalize1();
    return this->normalize2();
  }

  /**
   * @brief normalize to a canonical form
   *
   * denominator is always non-negative
   */
  CONSTEXPR14 void normalize1() {
    if (this->_den < Z(0)) {
      this->_num = -this->_num;
      this->_den = -this->_den;
    }
  }

  /**
   * @brief normalize to a canonical form
   *
   * denominator is always co-prime with numerator
   */
  CONSTEXPR14 auto normalize2() -> Z {
    Z common = gcd(this->_num, this->_den);
    if (common == Z(1) || common == Z(0)) {
      return common;
    }
    this->_num /= common;
    this->_den /= common;
    return common;
  }

  /**
   * @brief Construct a new Fraction object
   *
   * @param[in] num
   */
  CONSTEXPR14 explicit Fraction(Z &&num) : _num{std::move(num)}, _den(Z(1)) {}

  /**
   * @brief Construct a new Fraction object
   *
   * @param[in] num
   */
  CONSTEXPR14 explicit Fraction(const Z &num) : _num{num}, _den(1) {}

  /**
   * @brief Construct a new Fraction object
   *
   * @param[in] num
   */
  CONSTEXPR14 Fraction() : _num(0), _den(1) {}

  /**
   * @brief
   *
   * @return const Z&
   */
  CONSTEXPR14 auto num() const noexcept -> const Z & { return _num; }

  /**
   * @brief
   *
   * @return const Z&
   */
  CONSTEXPR14 auto den() const noexcept -> const Z & { return _den; }

  /**
   * @brief cross product
   *
   * @param[in] rhs
   * @return Z
   */
  CONSTEXPR14 auto cross(const Fraction &rhs) const -> Z {
    return this->_num * rhs._den - this->_den * rhs._num;
  }

  /** @name Comparison operators
   *  ==, !=, <, >, <=, >= etc.
   */
  ///@{

  /**
   * @brief Equal to
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator==(const Fraction &lhs, const Z &rhs)
      -> bool {
    if (lhs._den == Z(1) || rhs == Z(0)) {
      return lhs._num == rhs;
    }
    auto lhs2{lhs};
    auto rhs2{rhs};
    std::swap(lhs2._den, rhs2);
    lhs2.normalize2();
    return lhs2._num < lhs2._den * rhs2;
  }

  /**
   * @brief Less than
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator<(const Fraction &lhs, const Z &rhs) -> bool {
    if (lhs._den == Z(1) || rhs == Z(0)) {
      return lhs._num < rhs;
    }
    auto lhs2{lhs};
    auto rhs2{rhs};
    std::swap(lhs2._den, rhs2);
    lhs2.normalize2();
    return lhs2._num < lhs2._den * rhs2;
  }

  /**
   * @brief Less than
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator<(const Z &lhs, const Fraction &rhs) -> bool {
    if (rhs._den == Z(1) || lhs == Z(0)) {
      return lhs < rhs._num;
    }
    auto lhs2{lhs};
    auto rhs2{rhs};
    std::swap(rhs2._den, lhs2);
    rhs2.normalize2();
    return rhs2._den * lhs2 < rhs2._num;
  }

  /**
   * @brief Equal to
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator==(const Z &lhs, const Fraction &rhs)
      -> bool {
    return rhs == lhs;
  }

  /**
   * @brief Equal to
   *
   * @param[in] rhs
   * @return true
   * @return false
   */

  /**
   * @brief Equal to
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator==(const Fraction &lhs, const Fraction &rhs)
      -> bool {
    if (lhs._den == rhs._den) {
      return lhs._num == rhs._num;
    }
    auto lhs2{lhs};
    auto rhs2{rhs};
    std::swap(lhs2._den, rhs2._num);
    lhs2.normalize2();
    rhs2.normalize2();
    return lhs2._num * rhs2._den == lhs2._den * rhs2._num;
  }

  /**
   * @brief Less than
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator<(const Fraction &lhs, const Fraction &rhs)
      -> bool {
    if (lhs._den == rhs._den) {
      return lhs._num < rhs._num;
    }
    auto lhs2{lhs};
    auto rhs2{rhs};
    std::swap(lhs2._den, rhs2._num);
    lhs2.normalize2();
    rhs2.normalize2();
    return lhs2._num * rhs2._den < lhs2._den * rhs2._num;
  }

  /**
   * @brief
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator!=(const Fraction &rhs) const -> bool {
    return !(*this == rhs);
  }

  /**
   * @brief Greater than
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator>(const Fraction &rhs) const -> bool {
    return rhs < *this;
  }

  /**
   * @brief Greater than or euqal to
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator>=(const Fraction &rhs) const -> bool {
    return !(*this < rhs);
  }

  /**
   * @brief Less than or equal to
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator<=(const Fraction &rhs) const -> bool {
    return !(rhs < *this);
  }

  /**
   * @brief Greater than
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator>(const Z &rhs) const -> bool { return rhs < *this; }

  /**
   * @brief Less than or equal to
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator<=(const Z &rhs) const -> bool {
    return !(rhs < *this);
  }

  /**
   * @brief Greater than or equal to
   *
   * @param[in] rhs
   * @return true
   * @return false
   */
  CONSTEXPR14 auto operator>=(const Z &rhs) const -> bool {
    return !(*this < rhs);
  }

  /**
   * @brief Greater than
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator>(const Z &lhs, const Fraction &rhs) -> bool {
    return rhs < lhs;
  }

  /**
   * @brief Less than or equal to
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator<=(const Z &lhs, const Fraction &rhs)
      -> bool {
    return !(rhs < lhs);
  }

  /**
   * @brief Greater than or euqal to
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return true
   * @return false
   */
  friend CONSTEXPR14 auto operator>=(const Z &lhs, const Fraction &rhs)
      -> bool {
    return !(lhs < rhs);
  }

  ///@}

  /**
   * @brief reciprocal
   *
   */
  CONSTEXPR14 void reciprocal() {
    std::swap(this->_num, this->_den);
    this->normalize1();
  }

  /**
   * @brief multiply and assign
   *
   * @param[in] rhs
   * @return Fraction&
   */
  CONSTEXPR14 auto operator*=(Fraction rhs) -> Fraction & {
    std::swap(this->_num, rhs._num);
    this->normalize2();
    rhs.normalize2();
    this->_num *= rhs._num;
    this->_den *= rhs._den;
    return *this;
  }

  /**
   * @brief multiply
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator*(Fraction lhs, const Fraction &rhs)
      -> Fraction {
    return lhs *= rhs;
  }

  /**
   * @brief multiply and assign
   *
   * @param[in] rhs
   * @return Fraction&
   */
  CONSTEXPR14 auto operator*=(Z rhs) -> Fraction & {
    std::swap(this->_num, rhs);
    this->normalize2();
    this->_num *= rhs;
    return *this;
  }

  /**
   * @brief multiply
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator*(Fraction lhs, const Z &rhs) -> Fraction {
    return lhs *= rhs;
  }

  /**
   * @brief multiply
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator*(const Z &lhs, Fraction rhs) -> Fraction {
    return rhs *= lhs;
  }

  /**
   * @brief divide and assign
   *
   * @param[in] rhs
   * @return Fraction&
   */
  CONSTEXPR14 auto operator/=(Fraction rhs) -> Fraction & {
    std::swap(this->_den, rhs._num);
    this->normalize();
    rhs.normalize2();
    this->_num *= rhs._den;
    this->_den *= rhs._num;
    return *this;
  }

  /**
   * @brief divide
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator/(Fraction lhs, const Fraction &rhs)
      -> Fraction {
    return lhs /= rhs;
  }

  /**
   * @brief divide and assign
   *
   * @param[in] rhs
   * @return Fraction&
   */
  CONSTEXPR14 auto operator/=(Z rhs) -> Fraction & {
    std::swap(this->_den, rhs);
    this->normalize();
    this->_den *= rhs;
    return *this;
  }

  /**
   * @brief divide
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator/(Fraction lhs, const Z &rhs) -> Fraction {
    return lhs /= rhs;
  }

  /**
   * @brief divide
   *
   * @param[in] lhs
   * @param[in] rhs
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator/(const Z &lhs, Fraction rhs) -> Fraction {
    rhs.reciprocal();
    return rhs *= lhs;
  }

  /**
   * @brief Negate
   *
   * @return Fraction
   */
  CONSTEXPR14 auto operator-() const -> Fraction {
    auto res = Fraction(*this);
    res._num = -res._num;
    return res;
  }

  /**
   * @brief Add
   *
   * @param[in] rhs
   * @return Fraction
   */
  CONSTEXPR14 auto operator+(const Fraction &rhs) const -> Fraction {
    if (this->_den == rhs._den) {
      return Fraction(this->_num + rhs._num, this->_den);
    }
    const auto common = gcd(this->_den, rhs._den);
    if (common == Z(0)) {
      return Fraction(rhs._den * this->_num + this->_den * rhs._num, Z(0));
    }
    const auto l = this->_den / common;
    const auto r = rhs._den / common;
    auto d = this->_den * r;
    auto n = r * this->_num + l * rhs._num;
    return Fraction(std::move(n), std::move(d));
  }

  /**
   * @brief Subtract
   *
   * @param[in] frac
   * @return Fraction
   */
  CONSTEXPR14 auto operator-(const Fraction &frac) const -> Fraction {
    return *this + (-frac);
  }

  /**
   * @brief Add
   *
   * @param[in] frac
   * @param[in] i
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator+(Fraction frac, const Z &i) -> Fraction {
    return frac += i;
  }

  /**
   * @brief Add
   *
   * @param[in] i
   * @param[in] frac
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator+(const Z &i, Fraction frac) -> Fraction {
    return frac += i;
  }

  /**
   * @brief
   *
   * @param[in] i
   * @return Fraction
   */
  CONSTEXPR14 auto operator-(const Z &i) const -> Fraction {
    return *this + (-i);
  }

  /**
   * @brief
   *
   * @param[in] rhs
   * @return Fraction
   */
  CONSTEXPR14 auto operator+=(const Fraction &rhs) -> Fraction & {
    return *this -= (-rhs);
  }

  /**
   * @brief
   *
   * @param[in] rhs
   * @return Fraction
   */
  CONSTEXPR14 auto operator-=(const Fraction &rhs) -> Fraction & {
    if (this->_den == rhs._den) {
      this->_num -= rhs._num;
      this->normalize2();
      return *this;
    }

    auto other{rhs};
    std::swap(this->_den, other._num);
    auto common_n = this->normalize2();
    auto common_d = other.normalize2();
    std::swap(this->_den, other._num);
    this->_num = this->cross(other);
    this->_den *= other._den;
    std::swap(this->_den, common_d);
    this->normalize2();
    this->_num *= common_n;
    this->_den *= common_d;
    this->normalize2();
    return *this;
  }

  /**
   * @brief
   *
   * @param[in] i
   * @return Fraction
   */
  CONSTEXPR14 auto operator+=(const Z &i) -> Fraction & {
    return *this -= (-i);
  }

  /**
   * @brief
   *
   * @param[in] rhs
   * @return Fraction
   */
  CONSTEXPR14 auto operator-=(const Z &rhs) -> Fraction & {
    if (this->_den == Z(1)) {
      this->_num -= rhs;
      return *this;
    }

    auto other{rhs};
    std::swap(this->_den, other);
    auto common_n = this->normalize2();
    std::swap(this->_den, other);
    this->_num -= other * this->_den;
    this->_num *= common_n;
    this->normalize2();
    return *this;
  }

  /**
   * @brief
   *
   * @param[in] c
   * @param[in] frac
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator-(const Z &c, const Fraction &frac)
      -> Fraction {
    return c + (-frac);
  }

  /**
   * @brief
   *
   * @param[in] c
   * @param[in] frac
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator+(int &&c, const Fraction &frac) -> Fraction {
    return frac + Z(c);
  }

  /**
   * @brief
   *
   * @param[in] c
   * @param[in] frac
   * @return Fraction
   */
  friend CONSTEXPR14 auto operator-(int &&c, const Fraction &frac) -> Fraction {
    return (-frac) + Z(c);
  }

  /**
   * @brief
   *
   * @param[in] c
   * @param[in] frac
   * @return Fraction<Z>
   */
  friend CONSTEXPR14 auto operator*(int &&c, const Fraction &frac) -> Fraction {
    return frac * Z(c);
  }

  /**
   * @brief
   *
   * @tparam _Stream
   * @tparam Z
   * @param[in] os
   * @param[in] frac
   * @return _Stream&
   */
  template <typename _Stream>
  friend auto operator<<(_Stream &os, const Fraction &frac) -> _Stream & {
    os << "(" << frac.num() << "/" << frac.den() << ")";
    return os;
  }
};

// For template deduction
// typename{Z} Fraction(const Z &, const Z &) noexcept -> Fraction<Z>;

} // namespace fun
