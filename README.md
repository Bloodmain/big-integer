# Bigint

This repository implements the class of the big signed integer. 

This class has:
* (non-allocating) Default constructor, initializing the number with 0.
* Copy constructor.
* Constructors from the integer types.
* Explicit constructor from `std::string`
* Assignment operators.
* Comparison operators.
* Arithmetic operators (add, sub, mul, div, unary minus and plus).
* Increment and decrement.
* Bitwise operators (and, or, xor, not).
* Bitwise shifts.
* `std::string to_string(big_integer const&)`

## Some realisation info
* Mul and div have O(nm) time complexity
* Digits of the number are represented with 32 bits numbers. 
* Division uses BasecaseDivRem approach from ["Modern Computer Arithmetic"](https://members.loria.fr/PZimmermann/mca/mca-0.5.pdf)