# Floating point arithmetic

an IEEE 754 conforming software implementation of basic floating point arithmetic written in pure C
## Features

- support for addition, subtraction, multiplication and division
- calculations can be carried out with half and full precision
- prints the result in the same format as the ```%a``` format specifier

## Build

```sh
make all
```

## Usage

```
float <precision (h | f)> <rounding> <first number (in hex IEEE form)> [<operation> <second number (in hex IEEE form)>]
```
