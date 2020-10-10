![Python package](https://github.com/timleslie/pyblas/workflows/Unit%20tests/badge.svg)
![test-black](https://github.com/timleslie/pyblas/workflows/Linting%20-%20Black/badge.svg)

# 🔢 PyBLAS

PyBLAS is a python port of the [netlib reference BLAS implementation](http://www.netlib.org/blas/).

## Usage

```shell
pip install numpy pyblas
```

```python
import numpy as np
from pyblas.level1 import dswap

x = np.array([1.2, 2.3, 3.4])
y = np.array([5.6, 7.8, 9.0])
N = len(x)  # The length of the vectors x and y
incx = 1  # The index spacing of the vector x
incy = 1  # The index spacing of the vector y

# Swap the values of the vectors x and y
dswap(N, x, incx, y, incy)
print(x, y)
```

For more details on how to use the PyBLAS library, please consult our [docs](/docs)

## What is BLAS

The **Basic Linear Algebra Subprograms** (BLAS) are a collection of functions which form the basis of many modern numerical computing packages, including numpy, scipy, and matlab.
They provide functions for performing basic calculations on vectors and matrices, which form the basis for more complex calculations such as solving systems of linear equations.

## Who is this for?

 * Maths and computer science students who are learning about linear algebra and want to see how to implement simple operations.
 * Algorithm developers who want to prototype their calculations in a high level language with the same APIs they will use in C or Fortran.
 * Data scientists who want to better understand what is going on under the hood of the algorithms they use.
 * Educators who want an easy to use BLAS implementation when teaching numerical computing courses.

## Performance

Optimal performance is *not* a goal of the PyBLAS project.
The project aim is to have a simple and readable implementation of the BLAS standards.
As such, we often forego optimization opporunities in favour of simplicity.

The project matches overall algorithmic complexity with the reference implementation for all functions.

## Accuracy

The project aims to match the numerical accuracy of the reference BLAS implementation.
A future goal is to have benchmarks which can be run against the system BLAS libraries to verify the numerical accuract of the functions

### Level 1

| single (s)           | double (d)    | single complex (c) | double complex (z) |                                           |
| :---:                | :---:         | :---:          | :---:          | :---:                                             |
| scopy                | dcopy         | ccopy          | zcopy          | y := x                                            |
| sswap                | dswap         | cswap          | zswap          | x, y := y                                         |
| sscal                | dscal         | cscal (csscal) | zscal (zdscal) | x := a*x                                          |
| saxpy                | daxpy         | caxpy          | zaxpy          | y := a*x + y                                      |
|                      |               | scabs1         | dcabs1         | `=> |Re(x_i)| + |Im(x_i)|`                        |
| sasum                | dasum         |                |                | `=> sum(|x_i|)`                                   |
|                      |               | scasum         | dzasum         | `=> sum(|Re(x_i)| + |Im(x_i)})`                   |
| sdot (sdsdot, dsdot) | ddot          | cdotu          | zdotu          | `=>  <x, y>`                                      |
|                      |               | cdotc          | zdotc          | `=>  <x^H, y>`                                    |
| snrm2                | dnrm2         |                |                | `=> sqrt(<x, x>`)                                 |
|                      |               | scnrm2         | dznrm2         | `=> sqrt(<x^H, x>`)                               |
| srot                 | drot          | csrot          | zdrot          | `[x_i, y_i] := [c s; -s c] * [x_i, y_i]`          |
| srotg                | drotg         | crotg          | zrotg          | `c := a/r, s:= b/r, a :=r, b := "z"`              |
| srotm                | drotm         |                |                | `[x_i, y_i] := [h_1 h_2; h_3, h_4] * [x_i, y_i] ` |
| srotmg               | drotmg        |                |                |                                                   |

### Level 2

| single (s)           | double (d)    | single complex (c) | double complex (z) |                                           |
| :---:                | :---:         | :---:          | :---:          | :---:                                             |
| ssyr   | dsyr   |        |       | `A := a*x*x^T + A` (sym)                  |
|        |        | cher   | zher  | `A := a*x*x^H + A` (her)                  |
| sspr   | dspr   |        |       | `A := a*x*x^T + A` (sym-packaged)         |
|        |        | chpr   | zhpr  | `A := a*x*x^H + A` (her-packaged)         |
| sger   | dger   | cgeru  | zgeru | `A := a*x*y^T + A`                        |
|        |        | cgerc  | zgerc | `A := a*x*y^H + A`                        |
| ssyr2  | dsyr2  |        |       | `A := a*x*y^T + a*y*x^T + A` (sym)        |
|        |        | cher2  | zher2 | `A := a*x*y^H + a^+*y*x^H + A` (her)      |
| sspr2  | dspr2  |        |       | `A := a*x*y^T + a*y*x^T + A` (sym-packed) |
|        |        | chpr2  | zhpr2 | `A := a*x*y^H + a*y*x^H + A` (her-packed) |

| | |  |
| strsv  | dtrsv  | ctrsv  | ztrsv | `x := A^-1*b` or `x := A^[TH]^-1*b` (tri) |
| stbsv  | dtbsv  | ctbsv  | ztbsv | `x := A^-1*b` or `x := A^[TH]^-1*b` (band) |
| stpsv  | dtpsv  | ctpsv  | ztpsv | `x := A^-1*b` or `x := A^[TH]^-1*b` (tri-packed) |
| | | |
| stpmv  | dtpmv  | ctpmv  | ztpmv | `x := A^[1TH]*x` (sym-packed)   |
| stbmv  | dtbmv  | ctbmv  | ztbmv | `x := A^[1TH]*x` (tri-band) |
| strmv  | dtrmv  | ctrmv  | ztrmv | `x := A^[1TH]*x` (tri) |
| sgemv  | dgemv  | cgemv  | zgemv | `y := a*A^[1TH]*x + b*y` |
|        |        | chemv  | zhemv | `y := a*A*x + b*y` (herm) |
| sgbmv  | dgbmv  | cgbmv  | zgbmv | `y := a*A^[1TH]*x + b*y` (band) |
|        |        | chbmv  | zhbmv | `y := a*A*x + b*y` (band-herm)   |
| ssymv  | dsymv  |                | `y := a*A*x + b*y` (sym)        |
| sspmv  | dspmv  |                | `y := a*A*x + b*y` (sym-packed) |
|        |        | chpmv  | zhpmv | `y := a*A*x + b*y` (her-packed) |
| ssbmv  | dsbmv  |        |       | `y := a*A*x + b*y` (sym-band)   |
| isamax | idamax |        |     | | ` => argmax(|x_i|)` |
|        |        | icamax |izamax | ` => argmax(|Re(x_i)| + |Im(x_i)|)` |


### Level 3
