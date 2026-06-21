(* matrix.sml

   Dense real linear algebra over general m x n matrices. See matrix.sig for
   the contract. Storage is a flat `real array` in row-major order plus the
   row/column counts: entry (i, j) lives at flat index i * c + j. *)

structure Matrix :> MATRIX =
struct
  exception Dim of string
  exception Singular

  val eps = 1E~9

  (* A matrix is its row count, column count, and a row-major data array. *)
  type t = {r : int, c : int, a : real array}

  fun rows ({r, ...} : t) = r
  fun cols ({c, ...} : t) = c

  fun idx (c, i, j) = i * c + j

  fun checkIndex ({r, c, ...} : t, i, j) =
    if i < 0 orelse i >= r orelse j < 0 orelse j >= c then
      raise Dim ("index (" ^ Int.toString i ^ ", " ^ Int.toString j ^
                 ") out of range for " ^ Int.toString r ^ "x" ^ Int.toString c)
    else ()

  fun sub (m as {c, a, ...} : t, i, j) =
    (checkIndex (m, i, j); Array.sub (a, idx (c, i, j)))

  fun update (m as {c, a, ...} : t, i, j, x) =
    (checkIndex (m, i, j); Array.update (a, idx (c, i, j), x))

  (* --- construction --- *)

  fun make (r, c) x =
    if r <= 0 orelse c <= 0 then
      raise Dim ("make: non-positive dimension " ^
                 Int.toString r ^ "x" ^ Int.toString c)
    else {r = r, c = c, a = Array.array (r * c, x)}

  fun zeros (r, c) = make (r, c) 0.0

  fun identity n =
    let
      val m = zeros (n, n)
      val {a, ...} = m
      fun loop i = if i >= n then () else (Array.update (a, idx (n, i, i), 1.0); loop (i + 1))
    in
      loop 0; m
    end

  fun fromRows rs =
    let
      val r = List.length rs
      val () = if r = 0 then raise Dim "fromRows: empty matrix" else ()
      val c = List.length (List.hd rs)
      val () = if c = 0 then raise Dim "fromRows: empty row" else ()
      val () =
        List.app
          (fn row =>
             if List.length row <> c then
               raise Dim "fromRows: ragged rows"
             else ())
          rs
      val a = Array.array (r * c, 0.0)
      val _ =
        List.foldl
          (fn (row, i) =>
             (ignore (List.foldl
                (fn (x, j) => (Array.update (a, idx (c, i, j), x); j + 1))
                0 row);
              i + 1))
          0 rs
    in
      {r = r, c = c, a = a}
    end

  fun toRows ({r, c, a} : t) =
    List.tabulate
      (r, fn i => List.tabulate (c, fn j => Array.sub (a, idx (c, i, j))))

  (* --- arithmetic --- *)

  fun sameShape (x : t, y : t) = rows x = rows y andalso cols x = cols y

  fun zipWith name f (x as {r, c, a} : t, y : t) =
    if not (sameShape (x, y)) then
      raise Dim (name ^ ": shape mismatch " ^
                 Int.toString (rows x) ^ "x" ^ Int.toString (cols x) ^ " vs " ^
                 Int.toString (rows y) ^ "x" ^ Int.toString (cols y))
    else
      let
        val {a = b, ...} = y
        val out = Array.array (r * c, 0.0)
        fun loop k =
          if k >= r * c then ()
          else (Array.update (out, k, f (Array.sub (a, k), Array.sub (b, k)));
                loop (k + 1))
      in
        loop 0; {r = r, c = c, a = out}
      end

  fun add (x, y) = zipWith "add" (op +) (x, y)
  fun sub' (x, y) = zipWith "sub'" (op -) (x, y)

  fun scale s ({r, c, a} : t) =
    let
      val out = Array.array (r * c, 0.0)
      fun loop k =
        if k >= r * c then ()
        else (Array.update (out, k, s * Array.sub (a, k)); loop (k + 1))
    in
      loop 0; {r = r, c = c, a = out}
    end

  fun transpose ({r, c, a} : t) =
    let
      val out = Array.array (r * c, 0.0)
      fun loop (i, j) =
        if i >= r then ()
        else if j >= c then loop (i + 1, 0)
        else (Array.update (out, idx (r, j, i), Array.sub (a, idx (c, i, j)));
              loop (i, j + 1))
    in
      loop (0, 0); {r = c, c = r, a = out}
    end

  fun mul (x as {r = rx, c = cx, a = ax} : t, y as {r = ry, c = cy, a = ay} : t) =
    if cx <> ry then
      raise Dim ("mul: inner dimensions disagree " ^
                 Int.toString rx ^ "x" ^ Int.toString cx ^ " * " ^
                 Int.toString ry ^ "x" ^ Int.toString cy)
    else
      let
        val out = Array.array (rx * cy, 0.0)
        fun cell (i, j) =
          let
            fun acc (k, s) =
              if k >= cx then s
              else acc (k + 1,
                        s + Array.sub (ax, idx (cx, i, k))
                            * Array.sub (ay, idx (cy, k, j)))
          in
            acc (0, 0.0)
          end
        fun loop (i, j) =
          if i >= rx then ()
          else if j >= cy then loop (i + 1, 0)
          else (Array.update (out, idx (cy, i, j), cell (i, j)); loop (i, j + 1))
      in
        loop (0, 0); {r = rx, c = cy, a = out}
      end

  (* --- LU with partial pivoting --- *)

  fun requireSquare name (m : t) =
    if rows m <> cols m then
      raise Dim (name ^ ": matrix is not square (" ^
                 Int.toString (rows m) ^ "x" ^ Int.toString (cols m) ^ ")")
    else rows m

  (* Returns (lu, perm, sign): `lu` is a fresh n x n array holding L (below the
     diagonal, unit diagonal implied) and U (on and above the diagonal) packed
     together; `perm` maps result row -> original row; `sign` is the
     permutation parity. Does not raise on singular input; a zero pivot simply
     yields a zero on U's diagonal. *)
  fun factor (m : t) =
    let
      val n = requireSquare "lu" m
      val {a = src, ...} = m
      val a = Array.tabulate (n * n, fn k => Array.sub (src, k))
      val perm = Array.tabulate (n, fn i => i)
      val sign = ref 1.0
      fun get (i, j) = Array.sub (a, idx (n, i, j))
      fun set (i, j, x) = Array.update (a, idx (n, i, j), x)
      fun swapRows (p, q) =
        if p = q then ()
        else
          let
            val () =
              List.app
                (fn j =>
                   let val t = get (p, j) in set (p, j, get (q, j)); set (q, j, t) end)
                (List.tabulate (n, fn j => j))
            val tp = Array.sub (perm, p)
          in
            Array.update (perm, p, Array.sub (perm, q));
            Array.update (perm, q, tp);
            sign := ~ (!sign)
          end
      fun pivot k =
        let
          fun best (i, bi, bv) =
            if i >= n then bi
            else
              let val v = Real.abs (get (i, k))
              in if v > bv then best (i + 1, i, v) else best (i + 1, bi, bv) end
        in
          best (k + 1, k, Real.abs (get (k, k)))
        end
      fun eliminate k =
        if k >= n then ()
        else
          let
            val () = swapRows (k, pivot k)
            val pivv = get (k, k)
          in
            if Real.abs pivv <= eps then eliminate (k + 1)
            else
              let
                fun rowloop i =
                  if i >= n then ()
                  else
                    let
                      val f = get (i, k) / pivv
                      val () = set (i, k, f)
                      fun colloop j =
                        if j >= n then ()
                        else (set (i, j, get (i, j) - f * get (k, j));
                              colloop (j + 1))
                    in
                      colloop (k + 1); rowloop (i + 1)
                    end
              in
                rowloop (k + 1); eliminate (k + 1)
              end
          end
    in
      eliminate 0;
      {n = n, a = a, perm = perm, sign = !sign}
    end

  fun lu (m : t) =
    let
      val {n, a, perm, sign} = factor m
      val L = identity n
      val U = zeros (n, n)
      val () =
        List.app
          (fn i =>
             List.app
               (fn j =>
                  let val v = Array.sub (a, idx (n, i, j))
                  in
                    if j < i then update (L, i, j, v)
                    else update (U, i, j, v)
                  end)
               (List.tabulate (n, fn j => j)))
          (List.tabulate (n, fn i => i))
    in
      {l = L, u = U, p = perm, sign = sign}
    end

  fun det (m : t) =
    let
      val {n, a, sign, ...} = factor m
      fun prod (i, acc) =
        if i >= n then acc else prod (i + 1, acc * Array.sub (a, idx (n, i, i)))
    in
      sign * prod (0, 1.0)
    end

  (* Solve using a precomputed factorisation against a single right-hand side. *)
  fun solveFactored ({n, a, perm, ...}, b : real array) =
    let
      val () =
        List.app
          (fn i =>
             if Real.abs (Array.sub (a, idx (n, i, i))) <= eps then
               raise Singular
             else ())
          (List.tabulate (n, fn i => i))
      (* Permute b: y[i] starts as b[perm[i]]. *)
      val y = Array.tabulate (n, fn i => Array.sub (b, Array.sub (perm, i)))
      (* Forward substitution: L y = Pb, unit lower diagonal. *)
      fun fwd i =
        if i >= n then ()
        else
          let
            fun acc (j, s) =
              if j >= i then s
              else acc (j + 1, s + Array.sub (a, idx (n, i, j)) * Array.sub (y, j))
          in
            Array.update (y, i, Array.sub (y, i) - acc (0, 0.0)); fwd (i + 1)
          end
      (* Back substitution: U x = y. *)
      val x = Array.array (n, 0.0)
      fun back i =
        if i < 0 then ()
        else
          let
            fun acc (j, s) =
              if j >= n then s
              else acc (j + 1, s + Array.sub (a, idx (n, i, j)) * Array.sub (x, j))
          in
            Array.update (x, i,
              (Array.sub (y, i) - acc (i + 1, 0.0)) / Array.sub (a, idx (n, i, i)));
            back (i - 1)
          end
    in
      fwd 0; back (n - 1); x
    end

  fun solve (m : t) b =
    let
      val n = requireSquare "solve" m
      val () =
        if List.length b <> n then
          raise Dim ("solve: rhs length " ^ Int.toString (List.length b) ^
                     " <> matrix size " ^ Int.toString n)
        else ()
      val f = factor m
      val bv = Array.fromList b
      val x = solveFactored (f, bv)
    in
      List.tabulate (n, fn i => Array.sub (x, i))
    end

  fun inv (m : t) =
    let
      val n = requireSquare "inv" m
      val f = factor m
      val out = zeros (n, n)
      val () =
        List.app
          (fn j =>
             let
               val e = Array.tabulate (n, fn i => if i = j then 1.0 else 0.0)
               val col = solveFactored (f, e)
             in
               List.app (fn i => update (out, i, j, Array.sub (col, i)))
                 (List.tabulate (n, fn i => i))
             end)
          (List.tabulate (n, fn j => j))
    in
      out
    end

  (* --- QR via classical Gram-Schmidt with one reorthogonalisation pass --- *)

  fun qr (m : t) =
    let
      val mm = rows m
      val nn = cols m
      val () =
        if mm < nn then
          raise Dim ("qr: needs rows >= cols, got " ^
                     Int.toString mm ^ "x" ^ Int.toString nn)
        else ()
      val Q = zeros (mm, nn)
      val R = zeros (nn, nn)
      fun colDot (Acol, Bcol) =
        let
          fun acc (i, s) =
            if i >= mm then s
            else acc (i + 1, s + sub (Q, i, Acol) * sub (Q, i, Bcol))
        in
          acc (0, 0.0)
        end
      (* Copy column j of A into column j of Q to start. *)
      fun initCol j =
        List.app (fn i => update (Q, i, j, sub (m, i, j)))
          (List.tabulate (mm, fn i => i))
      (* Subtract projection onto already-orthonormal column k, accumulating
         the coefficient into R[k,j]. *)
      fun project (j, k) =
        let
          val d = colDot (k, j)
          val () = update (R, k, j, sub (R, k, j) + d)
        in
          List.app
            (fn i => update (Q, i, j, sub (Q, i, j) - d * sub (Q, i, k)))
            (List.tabulate (mm, fn i => i))
        end
      fun norm j =
        let
          fun acc (i, s) =
            if i >= mm then s else acc (i + 1, s + sub (Q, i, j) * sub (Q, i, j))
        in
          Math.sqrt (acc (0, 0.0))
        end
      fun gs j =
        if j >= nn then ()
        else
          let
            val () = initCol j
            (* Two passes of orthogonalisation against columns 0..j-1. *)
            val cols' = List.tabulate (j, fn k => k)
            val () = List.app (fn k => project (j, k)) cols'
            val () = List.app (fn k => project (j, k)) cols'
            val nj = norm j
            val () = if nj <= eps then raise Singular else ()
            val () = update (R, j, j, nj)
            val () =
              List.app (fn i => update (Q, i, j, sub (Q, i, j) / nj))
                (List.tabulate (mm, fn i => i))
          in
            gs (j + 1)
          end
    in
      gs 0; {q = Q, r = R}
    end
end
