#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <utility>
#include <vector>

struct big_integer {
private:
  static constexpr std::size_t GROUP_SIZE = 9;
  static constexpr std::int32_t FROM_RADIX = 10;
  static constexpr std::int32_t BASE_POWER = 32; // assuming BASE = 2 ** BASE_POWER
  static constexpr std::int64_t BASE = 1ll << BASE_POWER;
  static constexpr std::int64_t INT64_T_MIN = std::numeric_limits<std::int64_t>::min();

  static constexpr std::array<std::int64_t, GROUP_SIZE + 1> RADIX_POWERS = [] {
    std::array<std::int64_t, GROUP_SIZE + 1> arr{};
    arr[0] = 1;
    for (std::size_t i = 1; i <= GROUP_SIZE; ++i) {
      arr[i] = FROM_RADIX * arr[i - 1];
    }
    return arr;
  }();

private:
  bool negative;
  std::vector<std::uint32_t> magnitude;

  void strip_zeros() noexcept;
  bool is_zero() const noexcept;
  bool compare_magnitudes(const big_integer& other) const;
  big_integer& bitwise_op_base(const big_integer& other,
                               const std::function<std::uint32_t(std::uint32_t, std::uint32_t)>& op);

  static std::uint32_t inverse_if_necessary(std::uint32_t val, bool necessary, bool negate) noexcept;

  static std::function<std::int64_t(std::int64_t, std::int64_t)> get_additive_op(bool sub);
  static std::pair<std::int64_t, std::int64_t> additive_mod(
      std::int64_t val, const std::function<std::int64_t(std::int64_t, std::int64_t)>& op);
  big_integer& proceed_additive(const big_integer& rhs, bool subtract);
  big_integer& mul_eq_short(std::uint32_t b);
  big_integer& divide_eq_short(std::uint32_t b);
  big_integer& add_eq_short(std::uint32_t b);
  big_integer& sub_eq_short(std::uint32_t b);
  big_integer mul_short(std::uint32_t b) const;
  friend std::pair<big_integer, std::uint32_t> divmod(const big_integer& a, std::uint32_t b);
  friend std::pair<big_integer, big_integer> divmod(const big_integer& a, const big_integer& b);

public:
  big_integer() noexcept;
  big_integer(const big_integer& other);
  big_integer(int a);
  big_integer(unsigned int a);
  big_integer(long a);
  big_integer(unsigned long a);
  big_integer(long long a);
  big_integer(unsigned long long a);
  explicit big_integer(const std::string& str);
  ~big_integer() noexcept;
  big_integer& operator=(const big_integer& other);

  void swap(big_integer& other) noexcept;
  big_integer& operator+=(const big_integer& rhs);
  big_integer& operator-=(const big_integer& rhs);
  big_integer& operator*=(const big_integer& rhs);
  big_integer& operator/=(const big_integer& rhs);
  big_integer& operator%=(const big_integer& rhs);

  big_integer& operator&=(const big_integer& rhs);
  big_integer& operator|=(const big_integer& rhs);
  big_integer& operator^=(const big_integer& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend big_integer operator+(const big_integer& a, const big_integer& b);
  friend big_integer operator-(const big_integer& a, const big_integer& b);
  friend big_integer operator*(const big_integer& a, const big_integer& b);
  friend big_integer operator/(const big_integer& a, const big_integer& b);
  friend big_integer operator%(const big_integer& a, const big_integer& b);

  friend big_integer operator&(const big_integer& a, const big_integer& b);
  friend big_integer operator|(const big_integer& a, const big_integer& b);
  friend big_integer operator^(const big_integer& a, const big_integer& b);

  friend big_integer operator<<(const big_integer& a, int b);
  friend big_integer operator>>(const big_integer& a, int b);

  friend bool operator==(const big_integer& a, const big_integer& b);
  friend bool operator!=(const big_integer& a, const big_integer& b);
  friend bool operator<(const big_integer& a, const big_integer& b);
  friend bool operator>(const big_integer& a, const big_integer& b);
  friend bool operator<=(const big_integer& a, const big_integer& b);
  friend bool operator>=(const big_integer& a, const big_integer& b);

  friend std::string to_string(const big_integer& a);
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, const big_integer& b);
big_integer operator*(const big_integer& a, const big_integer& b);
big_integer operator/(const big_integer& a, const big_integer& b);
big_integer operator%(const big_integer& a, const big_integer& b);

big_integer operator&(const big_integer& a, const big_integer& b);
big_integer operator|(const big_integer& a, const big_integer& b);
big_integer operator^(const big_integer& a, const big_integer& b);

big_integer operator<<(const big_integer& a, int b);
big_integer operator>>(const big_integer& a, int b);

bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
bool operator<(const big_integer& a, const big_integer& b);
bool operator>(const big_integer& a, const big_integer& b);
bool operator<=(const big_integer& a, const big_integer& b);
bool operator>=(const big_integer& a, const big_integer& b);

std::string to_string(const big_integer& a);
std::ostream& operator<<(std::ostream& out, const big_integer& a);
