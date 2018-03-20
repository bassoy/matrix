# C++ Matrix Template Library

The C++ Matrix Library is header-only library with two header files `matrix.h`, `timer.h` and a `test.cpp` file.
The `matrix.h` header file provides a matrix template class along with expression templates and generic lambdas of the C++14 standard.

## Requirements

- compiler supporting the C++17 standard (uses `if constexpr`)
- tested with clang 4.0.1 and gcc 7.2.0

## Template class `iosb::matrix`

- has template parameters `<T,M,N,S>` where `T` is the type, `S`, is the storage format, `M` and `N` the number of rows and columns of the matrix
- stores elements either in column-major or row-major storage format

## Template class `iosb::timer`

- encapsulates `std::chrono::high_resolution_clock` for convenient runtime measurements.
- has one template parameter `D` which can be e.g. of type `iosb::nanoseconds` or `iosb::seconds`
- simply wrap your expression with the static methods `tic()` and `toc()` to measure runtime.



## Overloaded operators

- elementwise arithmetic operators `+`,`-`,`*`,`/` with scalars and matrices
- elementwise assignment operators `+=`,`-=`,`*=`,`/=` with scalars and matrices
- matrix-multiplication operator `|`
- matrix-transposition operator `!`

## Testing and profiling

- uses openmp if `-fopenmp` flag is provided
- compile with `clang++ -Wall -Werror -std=c++1z -O3 -fopenmp -march=native test.cpp -o test` or
- compile with `g++ -Wall -Werror -std=c++1z -fopenmp -march=native -O3 test.cpp -o test`
- outputs the runtime of multiple matrix expressions 
- for smaller dimensions the `test` programm produces a `check.m` script that can be run in Octave or Matlab to verify the numerical results
