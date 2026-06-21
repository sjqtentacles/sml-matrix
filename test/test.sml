(* Tests for sml-matrix. All real comparisons use an absolute tolerance via
   `near`/`checkNear` so the printed ok/FAIL lines are identical across MLton
   and Poly/ML (no raw Real.toString in labels). *)

structure MatrixTests =
struct
  open Harness
  structure M = Matrix

  val eps = 1E~9

  fun near (a, b) = Real.abs (a - b) <= eps

  fun checkNear name (expected, actual) =
    check name (near (expected, actual))

  (* Every entry of two equal-shaped matrices agrees within eps. *)
  fun matNear (a, b) =
    M.rows a = M.rows b andalso M.cols a = M.cols b andalso
    let
      val ok = ref true
      val () =
        List.app
          (fn i =>
             List.app
               (fn j =>
                  if near (M.sub (a, i, j), M.sub (b, i, j)) then ()
                  else ok := false)
               (List.tabulate (M.cols a, fn j => j)))
          (List.tabulate (M.rows a, fn i => i))
    in
      !ok
    end

  fun checkMat name (expected, actual) = check name (matNear (expected, actual))

  fun run () =
    let
      val () = section "construction and shape"
      val a = M.fromRows [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
      val () = checkInt "rows" (2, M.rows a)
      val () = checkInt "cols" (3, M.cols a)
      val () = checkNear "sub (0,0)" (1.0, M.sub (a, 0, 0))
      val () = checkNear "sub (1,2)" (6.0, M.sub (a, 1, 2))
      val () = checkMat "toRows roundtrip"
        (a, M.fromRows (M.toRows a))

      val () = checkRaises "fromRows ragged"
        (fn () => M.fromRows [[1.0, 2.0], [3.0]])
      val () = checkRaises "fromRows empty" (fn () => M.fromRows [])
      val () = checkRaises "sub out of range" (fn () => M.sub (a, 2, 0))

      val z = M.zeros (2, 2)
      val () = checkNear "zeros entry" (0.0, M.sub (z, 1, 1))
      val f = M.make (2, 3) 7.0
      val () = checkNear "make entry" (7.0, M.sub (f, 0, 1))
      val i3 = M.identity 3
      val () = checkNear "identity diag" (1.0, M.sub (i3, 1, 1))
      val () = checkNear "identity off-diag" (0.0, M.sub (i3, 0, 2))

      val () = section "update"
      val u = M.zeros (2, 2)
      val () = M.update (u, 0, 1, 9.0)
      val () = checkNear "update sets" (9.0, M.sub (u, 0, 1))

      val () = section "add / sub / scale"
      val b = M.fromRows [[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]
      val () = checkMat "add"
        (M.fromRows [[7.0, 7.0, 7.0], [7.0, 7.0, 7.0]], M.add (a, b))
      val () = checkMat "sub'"
        (M.fromRows [[~5.0, ~3.0, ~1.0], [1.0, 3.0, 5.0]], M.sub' (a, b))
      val () = checkMat "scale"
        (M.fromRows [[2.0, 4.0, 6.0], [8.0, 10.0, 12.0]], M.scale 2.0 a)
      val () = checkRaises "add mismatch"
        (fn () => M.add (a, M.zeros (3, 2)))

      val () = section "transpose"
      val () = checkMat "transpose shape"
        (M.fromRows [[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]], M.transpose a)
      val () = checkMat "transpose twice = original"
        (a, M.transpose (M.transpose a))

      val () = section "mul"
      (* Hand-computed 2x3 * 3x2 product.
         [[1 2 3],[4 5 6]] * [[7 8],[9 10],[11 12]]
         = [[58 64],[139 154]] *)
      val p23 = M.fromRows [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
      val q32 = M.fromRows [[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]]
      val () = checkMat "2x3 * 3x2"
        (M.fromRows [[58.0, 64.0], [139.0, 154.0]], M.mul (p23, q32))
      val () = checkRaises "mul dim mismatch"
        (fn () => M.mul (p23, M.zeros (2, 2)))

      val () = section "identity is multiplicative identity"
      val g = M.fromRows [[2.0, ~1.0, 0.0], [3.0, 5.0, 1.0]]
      val () = checkMat "A * I = A" (g, M.mul (g, M.identity 3))
      val () = checkMat "I * A = A" (M.mul (M.identity 2, g), g)

      val () = section "lu / det"
      (* Triangular det = product of diagonal. *)
      val tri = M.fromRows
        [[2.0, 7.0, 1.0], [0.0, 3.0, 5.0], [0.0, 0.0, 4.0]]
      val () = checkNear "det triangular" (24.0, M.det tri)
      (* Known 3x3 det. det = ~306.0 for this classic example. *)
      val d3 = M.fromRows
        [[6.0, 1.0, 1.0], [4.0, ~2.0, 5.0], [2.0, 8.0, 7.0]]
      val () = checkNear "det 3x3" (~306.0, M.det d3)
      val () = checkNear "det identity" (1.0, M.det (M.identity 4))
      (* LU reconstructs: P*A = L*U, i.e. L*U has rows of A permuted by p. *)
      val {l, u, p, sign = _} = M.lu d3
      val lu' = M.mul (l, u)
      val () = checkBool "P*A = L*U"
        (true,
         List.all
           (fn i =>
              List.all
                (fn j => near (M.sub (lu', i, j),
                              M.sub (d3, Array.sub (p, i), j)))
                (List.tabulate (3, fn j => j)))
           (List.tabulate (3, fn i => i)))
      val () = checkRaises "det non-square" (fn () => M.det (M.zeros (2, 3)))

      val () = section "solve / inv"
      (* Hand-checked system: A x = b with
         A = [[2 1 ~1],[~3 ~1 2],[~2 1 2]], b = [8 ~11 ~3]
         solution x = [2 3 ~1]. *)
      val sa = M.fromRows
        [[2.0, 1.0, ~1.0], [~3.0, ~1.0, 2.0], [~2.0, 1.0, 2.0]]
      val sx = M.solve sa [8.0, ~11.0, ~3.0]
      val () = checkBool "solve 3x3"
        (true,
         ListPair.all near (sx, [2.0, 3.0, ~1.0])
         andalso List.length sx = 3)
      (* A * inv(A) ~= I. *)
      val inva = M.inv sa
      val () = checkMat "A * inv(A) = I" (M.identity 3, M.mul (sa, inva))
      val () = checkMat "inv(A) * A = I" (M.identity 3, M.mul (inva, sa))
      (* Singular matrix raises. *)
      val singular = M.fromRows
        [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 1.0, 1.0]]
      val () = checkRaises "solve singular"
        (fn () => M.solve singular [1.0, 2.0, 3.0])
      val () = checkRaises "inv singular" (fn () => M.inv singular)
      val () = checkRaises "solve rhs length"
        (fn () => M.solve sa [1.0, 2.0])

      val () = section "qr"
      val qa = M.fromRows
        [[12.0, ~51.0, 4.0], [6.0, 167.0, ~68.0], [~4.0, 24.0, ~41.0]]
      val {q, r} = M.qr qa
      (* Q*R ~= A. *)
      val () = checkMat "Q*R = A" (qa, M.mul (q, r))
      (* R is upper-triangular: below-diagonal entries ~ 0. *)
      val () = checkBool "R upper-triangular"
        (true,
         List.all
           (fn i =>
              List.all
                (fn j => j >= i orelse near (M.sub (r, i, j), 0.0))
                (List.tabulate (3, fn j => j)))
           (List.tabulate (3, fn i => i)))
      (* Q's columns orthonormal: Q^T Q ~= I. *)
      val () = checkMat "Q^T Q = I"
        (M.identity 3, M.mul (M.transpose q, q))
      (* Tall matrix (m > n). *)
      val tall = M.fromRows
        [[1.0, 0.0], [1.0, 1.0], [1.0, 2.0], [1.0, 3.0]]
      val {q = qt, r = rt} = M.qr tall
      val () = checkMat "tall Q*R = A" (tall, M.mul (qt, rt))
      val () = checkMat "tall Q^T Q = I"
        (M.identity 2, M.mul (M.transpose qt, qt))
      val () = checkRaises "qr wide raises" (fn () => M.qr (M.zeros (2, 3)))
    in
      ()
    end
end
