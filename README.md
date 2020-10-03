# pyblas
A python port of the netlib reference BLAS implementation


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
