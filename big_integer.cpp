#include "big_integer.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>

big_integer::big_integer() noexcept : negative(false) {}

big_integer::big_integer(const big_integer& other) = default;

big_integer::big_integer(int a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(unsigned int a) : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(long a) : big_integer(static_cast<long long>(a)) {}

big_integer::big_integer(unsigned long a) : big_integer(static_cast<unsigned long long>(a)) {}

big_integer::big_integer(long long a)
    : big_integer(static_cast<std::uint64_t>((a >= 0 || a == big_integer::INT64_T_MIN) ? a : -a)) {
  negative = a < 0;
}

big_integer::big_integer(unsigned long long a) : negative(false) {
  do {
    magnitude.push_back(a % BASE);
    a >>= BASE_POWER;
  } while (a);
}

big_integer::big_integer(const std::string& str) : big_integer() {
  if (str.empty()) {
    throw std::invalid_argument("Unexpected empty string in ctor");
  } else if (str[0] == '-' && str.size() == 1) {
    throw std::invalid_argument("Unable to construct integer from '-'");
  }

  magnitude.push_back(0);
  std::size_t actual_group_size = 0;
  for (std::size_t i = str[0] == '-' || str[0] == '+'; i < str.size(); i += actual_group_size) {
    if (str[i] < '0' || str[i] > '9') {
      std::string message = "Expected digit, found ";
      message.push_back(str[i]);
      throw std::invalid_argument(message);
    }
    int value = std::stoi(str.substr(i, big_integer::GROUP_SIZE), &actual_group_size, big_integer::FROM_RADIX);
    mul_eq_short(big_integer::RADIX_POWERS[actual_group_size]);
    add_eq_short(value);
  }
  if (str[0] == '-') {
    negative = true;
  }
}

big_integer::~big_integer() noexcept = default;

big_integer& big_integer::operator=(const big_integer& other) {
  if (&other != this) {
    big_integer(other).swap(*this);
  }
  return *this;
}

std::uint32_t big_integer::inverse_if_necessary(std::uint32_t val, bool necessary, bool negate) noexcept {
  if (necessary) {
    val = negate ? (-val) : (~val);
  }
  return val;
}

void big_integer::swap(big_integer& other) noexcept {
  std::swap(negative, other.negative);
  std::swap(magnitude, other.magnitude);
}

bool big_integer::is_zero() const noexcept {
  return magnitude.empty() || (magnitude.size() == 1 && magnitude[0] == 0);
}

void big_integer::strip_zeros() noexcept {
  while (magnitude.size() > 1 && magnitude.back() == 0) {
    magnitude.pop_back();
  }
  if (is_zero()) {
    negative = false;
  }
}

std::pair<big_integer, big_integer> divmod(const big_integer& a, const big_integer& b) {
  // Using BasecaseDivRem approach

  if (a.magnitude.size() < b.magnitude.size()) {
    return {0, a};
  }

  big_integer divisor(b);
  big_integer dividend(a);
  divisor.negative = false;
  dividend.negative = false;

  // normalization
  std::uint32_t factor = 1;
  if (divisor.magnitude.back() < big_integer::BASE / 2) {
    factor = big_integer::BASE / (divisor.magnitude.back() + 1);
    divisor.mul_eq_short(factor);
    dividend.mul_eq_short(factor);
  }

  big_integer quot;
  std::size_t quot_size = dividend.magnitude.size() - divisor.magnitude.size() + 1;
  quot.magnitude.assign(quot_size, 0);
  quot.negative = a.negative ^ b.negative;

  std::uint64_t d = divisor.magnitude.back();
  std::size_t n = divisor.magnitude.size();

  divisor <<= big_integer::BASE_POWER * (quot_size - 1);
  if (dividend >= divisor) {
    quot.magnitude[quot_size - 1] = 1;
    dividend -= divisor;
  }

  for (std::size_t j = quot_size - 1; j > 0; --j) {
    std::uint64_t q = 0;
    if (n + j - 1 < dividend.magnitude.size()) {
      q += static_cast<std::uint64_t>(dividend.magnitude[n + j - 1]) * big_integer::BASE;
    }
    if (n + j - 2 < dividend.magnitude.size()) {
      q += dividend.magnitude[n + j - 2];
    }
    q /= d;
    quot.magnitude[j - 1] = std::min(q, static_cast<std::uint64_t>(big_integer::BASE - 1));

    divisor >>= big_integer::BASE_POWER;
    dividend -= divisor.mul_short(quot.magnitude[j - 1]);
    while (dividend < 0) {
      quot.magnitude[j - 1]--;
      dividend += divisor;
    }
  }
  dividend.divide_eq_short(factor);
  dividend.negative = a.negative;
  quot.strip_zeros();
  dividend.strip_zeros();
  return {quot, dividend};
}

std::pair<big_integer, std::uint32_t> divmod(const big_integer& a, std::uint32_t b) {
  std::uint64_t carry = 0;
  big_integer quot;
  quot.magnitude.assign(a.magnitude.size(), 0);
  quot.negative = a.negative;

  for (std::size_t i = a.magnitude.size(); i > 0; --i) {
    std::uint64_t current = a.magnitude[i - 1] + carry * big_integer::BASE;
    quot.magnitude[i - 1] = current / b;
    carry = current % b;
  }
  quot.strip_zeros();
  return {quot, carry};
}

std::function<std::int64_t(std::int64_t, std::int64_t)> big_integer::get_additive_op(bool sub) {
  std::function<std::int64_t(std::int64_t, std::int64_t)> op;
  if (sub) {
    op = std::minus();
  } else {
    op = std::plus();
  }
  return op;
}

std::pair<std::int64_t, std::int64_t> big_integer::additive_mod(
    std::int64_t val, const std::function<std::int64_t(std::int64_t, std::int64_t)>& op) {
  std::pair<std::int64_t, std::int64_t> res = {val, 0};
  if (val < 0 || val >= big_integer::BASE) {
    res.first = op(res.first, big_integer::BASE);
    res.second = op(0, 1);
  }
  return res;
}

big_integer& big_integer::add_eq_short(std::uint32_t b) {
  if (magnitude.empty() || (negative && magnitude.size() == 1 && magnitude[0] <= b)) {
    if (magnitude.empty()) {
      magnitude.push_back(b);
    } else {
      magnitude[0] = b - magnitude[0];
    }
    negative = false;
    return *this;
  }
  std::function<std::int64_t(std::int64_t, std::int64_t)> op = get_additive_op(negative);
  std::pair<std::int64_t, std::int64_t> now = {0, op(0, b)};
  magnitude.resize(magnitude.size() + 1);
  for (std::size_t i = 0; now.second; ++i) {
    now = additive_mod(magnitude[i] + now.second, op);
    magnitude[i] = now.first;
  }
  strip_zeros();
  return *this;
}

big_integer& big_integer::sub_eq_short(std::uint32_t b) {
  negative ^= 1;
  add_eq_short(b);
  negative ^= 1;
  return *this;
}

big_integer& big_integer::divide_eq_short(std::uint32_t b) {
  std::pair<big_integer, std::uint32_t> res = divmod(*this, b);
  swap(res.first);
  return *this;
}

big_integer& big_integer::mul_eq_short(std::uint32_t b) {
  std::uint64_t carry = 0;
  std::size_t i = 0;

  magnitude.resize(magnitude.size() + 1);
  while (i < magnitude.size()) {
    std::uint64_t mul = static_cast<std::uint64_t>(magnitude[i]) * b + carry;
    magnitude[i] = mul % big_integer::BASE;
    carry = mul / big_integer::BASE;
    ++i;
  }

  strip_zeros();
  return *this;
}

big_integer big_integer::mul_short(std::uint32_t b) const {
  big_integer tmp(*this);
  tmp.mul_eq_short(b);
  return tmp;
}

big_integer& big_integer::proceed_additive(const big_integer& rhs, bool subtract) {
  const bool subtracting = negative != (rhs.negative ^ subtract);

  const bool less = compare_magnitudes(rhs); // abs(*this) < abs(rhs)
  const big_integer& greater = less ? rhs : *this;
  const big_integer& lower = less ? *this : rhs;
  std::function<std::int64_t(std::int64_t, std::int64_t)> op = get_additive_op(subtracting);
  std::size_t greater_size = greater.magnitude.size();
  magnitude.resize(greater_size + 1);
  std::pair<std::int64_t, std::int64_t> now = {0, 0};
  for (std::size_t i = 0; i < greater_size; ++i) {
    std::int64_t sum = now.second + greater.magnitude[i];
    if (i < lower.magnitude.size()) {
      sum = op(sum, lower.magnitude[i]);
    }
    now = additive_mod(sum, op);
    magnitude[i] = now.first;
  }
  if (now.second) {
    magnitude[greater_size] = now.second;
  }
  if (subtracting && less) {
    // |a| < |b|; either (a < 0, b > 0, (a + b) > 0) or (a > 0, b < 0, (a + b) < 0)
    negative = !negative;
  }
  strip_zeros();
  return *this;
}

big_integer& big_integer::operator+=(const big_integer& rhs) {
  return proceed_additive(rhs, false);
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  return proceed_additive(rhs, true);
}

big_integer& big_integer::operator*=(const big_integer& rhs) {
  std::size_t sz = magnitude.size();
  magnitude.resize(magnitude.size() + rhs.magnitude.size() + 1);
  for (std::size_t i = sz; i > 0; --i) {
    std::uint64_t block = magnitude[i - 1];
    std::uint32_t carry = 0;
    magnitude[i - 1] = 0;
    for (std::size_t j = 0; j < rhs.magnitude.size() || carry; ++j) {
      // (2^32 - 1) * (2^32 - 1) + (2^32 - 1) + (2^32 - 1) = ull::max()
      std::uint64_t mul = block * (j < rhs.magnitude.size() ? rhs.magnitude[j] : 0) + magnitude[i + j - 1] + carry;
      magnitude[i + j - 1] = mul % big_integer::BASE;
      carry = mul / big_integer::BASE;
    }
  }
  negative ^= rhs.negative;
  strip_zeros();
  return *this;
}

big_integer& big_integer::operator/=(const big_integer& rhs) {
  std::pair<big_integer, big_integer> res = divmod(*this, rhs);
  swap(res.first);
  return *this;
}

big_integer& big_integer::operator%=(const big_integer& rhs) {
  std::pair<big_integer, big_integer> res = divmod(*this, rhs);
  swap(res.second);
  return *this;
}

big_integer& big_integer::bitwise_op_base(const big_integer& other,
                                          const std::function<std::uint32_t(std::uint32_t, std::uint32_t)>& op) {
  /*
   * -x = ~x + 1, so in the first block we're manipulating with -block,
   * in the next blocks there are two cases
   * 1) if all previous block were 0, so ~0 = 111.111 in binary ~0 + 1 = 0 and the "carry" to another block,
   *    thus we have to manipulate with -block as well
   * 2) in other case we're manipulating just with ~block
   */
  bool new_negativity = op(negative, other.negative);
  bool full_inverse_this = true;
  bool full_inverse_other = true;
  bool full_inverse_res = true;

  magnitude.resize(std::max(magnitude.size(), other.magnitude.size()));
  for (std::size_t i = 0; i < std::max(magnitude.size(), other.magnitude.size()); ++i) {
    std::uint32_t other_block = i < other.magnitude.size() ? other.magnitude[i] : 0;
    std::uint32_t res = inverse_if_necessary(op(inverse_if_necessary(magnitude[i], negative, full_inverse_this),
                                                inverse_if_necessary(other_block, other.negative, full_inverse_other)),
                                             new_negativity, full_inverse_res);
    if (magnitude[i] != 0) {
      full_inverse_this = false;
    }
    if (other.magnitude[i] != 0) {
      full_inverse_other = false;
    }

    magnitude[i] = res;
    if (magnitude[i] != 0) {
      full_inverse_res = false;
    }
  }
  negative = new_negativity;
  strip_zeros();
  return *this;
}

big_integer& big_integer::operator&=(const big_integer& rhs) {
  return bitwise_op_base(rhs, std::bit_and());
}

big_integer& big_integer::operator|=(const big_integer& rhs) {
  return bitwise_op_base(rhs, std::bit_or());
}

big_integer& big_integer::operator^=(const big_integer& rhs) {
  return bitwise_op_base(rhs, std::bit_xor());
}

big_integer& big_integer::operator<<=(int rhs) {
  /*
   * Since we have Base = 2 ** 32, shifting by 32 * x bits leads to adding x zero blocks in the start of the magnitude
   * So we do this rhs / 32 times
   * Than we're multiplying our value by 2 ** (rhs % 32) (i.e. by small (int) number)
   */

  magnitude.insert(magnitude.begin(), rhs / BASE_POWER, 0);
  return mul_eq_short(1u << static_cast<std::uint32_t>(rhs % BASE_POWER));
}

big_integer& big_integer::operator>>=(int rhs) {
  /*
   * Same as <<, but with dividing
   * N.B. Dividing rounds to 0, when ">>" to -inf, so we the result after diving may be equals to (real_answer + 1)
   * Thus if value is negative, we're checking the (value & (1 << rhs)) bit of the number,
   * where value is presented as two's complement number.
   * -x = ~x + 1, and there is two cases
   * 1) when we want to look at -block
   *    it's when it is a first block or all previous block equal to 0
   *    (i.e. ~0 = 11..11 in binary representation, thus if we add 1 to it, it will give "carry" to the next block),
   * 2) otherwise, we want to look at ~block
   */

  bool full_inverse = true;
  for (std::size_t i = 0; !is_zero() && i < (rhs / BASE_POWER); ++i) {
    if (magnitude[0] != 0) {
      full_inverse = false;
    }
  }
  magnitude.erase(magnitude.begin(), magnitude.begin() + (rhs / BASE_POWER));
  rhs %= BASE_POWER;
  // looking at the lowest bit as if we were shifting in two's complement
  bool odd = inverse_if_necessary(is_zero() ? 0 : magnitude[0], negative, full_inverse) & (1u << rhs);
  divide_eq_short(1u << rhs);
  // if that bit doesn't match our lowest bit, it means that division rounded wrong, and we should correct the result
  if (odd != ((is_zero() ? 0 : magnitude[0]) & 1)) {
    sub_eq_short(1);
  }
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  big_integer res(*this);
  if (!is_zero()) {
    res.negative = !negative;
  }
  return res;
}

big_integer big_integer::operator~() const {
  big_integer res(*this);
  res.negative = !res.negative;
  res.sub_eq_short(1);
  return res;
}

big_integer& big_integer::operator++() {
  return add_eq_short(1);
}

big_integer big_integer::operator++(int) {
  big_integer tmp = *this;
  ++(*this);
  return tmp;
}

big_integer& big_integer::operator--() {
  return sub_eq_short(1);
}

big_integer big_integer::operator--(int) {
  big_integer tmp = *this;
  --(*this);
  return tmp;
}

big_integer operator+(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp += b;
  return tmp;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp -= b;
  return tmp;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp *= b;
  return tmp;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  return divmod(a, b).first;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return divmod(a, b).second;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp &= b;
  return tmp;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp |= b;
  return tmp;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  big_integer tmp(a);
  tmp ^= b;
  return tmp;
}

big_integer operator<<(const big_integer& a, int b) {
  big_integer tmp(a);
  tmp <<= b;
  return tmp;
}

big_integer operator>>(const big_integer& a, int b) {
  big_integer tmp(a);
  tmp >>= b;
  return tmp;
}

bool big_integer::compare_magnitudes(const big_integer& other) const {
  if (magnitude.size() != other.magnitude.size()) {
    return magnitude.size() < other.magnitude.size();
  }
  return std::lexicographical_compare(magnitude.rbegin(), magnitude.rend(), other.magnitude.rbegin(),
                                      other.magnitude.rend());
}

bool operator==(const big_integer& a, const big_integer& b) {
  return (a.is_zero() && b.is_zero()) || (a.negative == b.negative && a.magnitude == b.magnitude);
}

bool operator!=(const big_integer& a, const big_integer& b) {
  return !(a == b);
}

bool operator<(const big_integer& a, const big_integer& b) {
  return a != b && a.negative ^ (a.negative == b.negative && a.compare_magnitudes(b));
}

bool operator>(const big_integer& a, const big_integer& b) {
  return b < a;
}

bool operator<=(const big_integer& a, const big_integer& b) {
  return !(a > b);
}

bool operator>=(const big_integer& a, const big_integer& b) {
  return !(a < b);
}

std::string to_string(const big_integer& a) {
  std::string res;

  big_integer to_write(a);
  do {
    std::pair<big_integer, std::uint32_t> divmod_res =
        divmod(to_write, big_integer::RADIX_POWERS[big_integer::GROUP_SIZE]);
    std::string string_mod = std::to_string(divmod_res.second);
    std::size_t ensure_length = (divmod_res.first.is_zero() ? 0 : big_integer::GROUP_SIZE - string_mod.size());
    res.append(string_mod.rbegin(), string_mod.rend());
    res.append(ensure_length, '0');
    to_write.swap(divmod_res.first);
  } while (!to_write.is_zero());

  if (a.negative && !a.is_zero()) {
    res.push_back('-');
  }

  std::reverse(res.begin(), res.end());
  return res;
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}
