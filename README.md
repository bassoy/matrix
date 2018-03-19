# C++ Matrix Template Library

The [C++ Matrix Library] is header-only library with one header file `matrix.h` and a `test.cpp` file.
It provides a matrix template class along with expression templates and generic lambdas of the C++14 standard.

## Requirements

- compiler supporting the C++14 standard
- tested with clang 4.0.1 and gcc 7.2.0

## Template class `iosb::matrix`

- template parameters `<T,M,N>` where `T` is the type, `M` and `N` the number of rows and columns of the matrix
- stores elements in column-major storage format (this could be generalized)


## Overloaded operators

- elementwise arithmetic operators `+`,`-`,`*`,`/` with scalars and matrices
- elementwise assignment operators `+=`,`-=`,`*=`,`/=` with scalars and matrices
- matrix-multiplication operator `|`
- matrix-transposition operator `!`

## Testing and profiling

- uses openmp if `-fopenmp` flag is provided
- compile with `clang++ -Wall -Werror -std=c++14 -O3 -fopenmp -march=native test.cpp -o test` or
- compile with `g++ -Wall -Werror -std=c++14 -fopenmp -march=native -O3 test.cpp -o test`
- outputs the runtime of multiple matrix expressions 
- for smaller dimensions the `test` programm produces a `check.m` script that can be run in Octave or Matlab to verify the numerical results
