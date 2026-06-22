# sml-matrix

[![CI](https://github.com/sjqtentacles/sml-matrix/actions/workflows/ci.yml/badge.svg)](https://github.com/sjqtentacles/sml-matrix/actions/workflows/ci.yml)

Dense real linear algebra for general `m x n` matrices in Standard ML. Where a
graphics library fixes matrices at 2/3/4 dimensions, `sml-matrix` works with
arbitrary shapes: construction, arithmetic, LU decomposition with partial
pivoting, determinant, linear solves, inverse, QR decomposition, Cholesky,
SVD, symmetric eigendecomposition, least squares, pseudoinverse, numerical
rank, matrix norms, and condition number.

## Storage

A matrix is a dense `real array` in **row-major** order (entry `(i, j)` at flat
index `i * cols + j`) carried alongside its row and column counts. The type is
opaque; build values with `fromRows`, `make`, `zeros`, or `identity`.

## Portability

Pure Standard ML using only the Basis library -- no FFI, no threads. Verified
on **MLton** and **Poly/ML**.

## Building and testing

```sh
make test        # build + run the suite under MLton (default)
make test-poly   # run the suite under Poly/ML
make all-tests   # run under both
make clean
```

## Installing with smlpkg

```sh
smlpkg add github.com/sjqtentacles/sml-matrix
smlpkg sync
```

Then reference the library basis from your own `.mlb`:

```
lib/github.com/sjqtentacles/sml-matrix/sml-matrix.mlb
```

For Poly/ML, `use` the `matrix.sig` and `matrix.sml` sources in order.

## Usage

```sml
val a = Matrix.fromRows [[1.0, 2.0], [3.0, 4.0]]
val b = Matrix.fromRows [[5.0, 6.0], [7.0, 8.0]]

val c   = Matrix.mul (a, b)            (* matrix product            *)
val d   = Matrix.det a                 (* ~2.0                      *)
val x   = Matrix.solve a [5.0, 11.0]   (* solve A x = b             *)
val ai  = Matrix.inv a                 (* inverse                   *)
val {q, r} = Matrix.qr a               (* Q orthonormal, R upper    *)

(* SPD Cholesky, A = L * L^T *)
val spd = Matrix.fromRows [[4.0, 12.0, ~16.0], [12.0, 37.0, ~43.0], [~16.0, ~43.0, 98.0]]
val l   = Matrix.cholesky spd

(* Singular value decomposition: A = u * diag(s) * vt *)
val {u, s, vt} = Matrix.svd a

(* Symmetric eigenvalues (ascending) and orthonormal eigenvectors *)
val {values, vectors} = Matrix.eigSym (Matrix.fromRows [[2.0, 1.0], [1.0, 2.0]])

(* Least squares for an overdetermined system, and the pseudoinverse *)
val coeffs = Matrix.lstsq design b      (* min ||A x - b||_2 via QR *)
val ap     = Matrix.pinv a              (* Moore-Penrose pseudoinverse *)

(* Norms, condition number, rank, trace *)
val n2 = Matrix.norm2 a                 (* spectral norm *)
val k  = Matrix.cond a                  (* 2-norm condition number *)
val rk = Matrix.rank a                  (* numerical rank *)
val tr = Matrix.trace a                 (* sum of the diagonal *)
```

Dimension errors raise `Matrix.Dim`; singular/non-SPD systems raise
`Matrix.Singular`.

## API summary

| Function | Description |
| --- | --- |
| `fromRows : real list list -> t` | Build from rows (equal lengths). |
| `make : int * int -> real -> t` | `r x c` filled with a constant. |
| `zeros : int * int -> t` | `r x c` of zeros. |
| `identity : int -> t` | `n x n` identity. |
| `rows`/`cols : t -> int` | Dimensions. |
| `sub : t * int * int -> real` | Element `(i, j)` (zero-indexed). |
| `update : t * int * int * real -> unit` | Set element in place. |
| `toRows : t -> real list list` | Rows as lists. |
| `add`/`sub' : t * t -> t` | Elementwise sum / difference. |
| `scale : real -> t -> t` | Scalar multiply. |
| `mul : t * t -> t` | Matrix product (dimension-checked). |
| `transpose : t -> t` | Transpose. |
| `lu : t -> {l, u, p, sign}` | LU with partial pivoting (`P*A = L*U`). |
| `det : t -> real` | Determinant via LU. |
| `solve : t -> real list -> real list` | Solve `A x = b`. |
| `inv : t -> t` | Inverse. |
| `qr : t -> {q, r}` | QR (Gram-Schmidt), `Q*R = A`. |
| `cholesky : t -> t` | Lower factor `L` of SPD `A = L*L^T`. |
| `svd : t -> {u, s, vt}` | Singular values/vectors, `A = u*diag(s)*vt`. |
| `eigSym : t -> {values, vectors}` | Symmetric eigenpairs (ascending values). |
| `lstsq : t -> real list -> real list` | Least-squares solution via QR. |
| `pinv : t -> t` | Moore-Penrose pseudoinverse via SVD. |
| `rank : t -> int` | Numerical rank from the SVD. |
| `trace : t -> real` | Sum of the diagonal. |
| `norm1`/`normInf`/`normFro`/`norm2 : t -> real` | 1-, inf-, Frobenius, spectral norms. |
| `cond : t -> real` | 2-norm condition number. |

## License

MIT. See [LICENSE](LICENSE).
