(* demo.sml - dense linear algebra on small fixed integer matrices. Every
   value shown here is an exact integer (products, transpose, determinant and
   trace of integer matrices), so entries are printed with Real.round and no
   real ever reaches Real.toString. Deterministic across runs and compilers. *)

structure M = Matrix

(* All demo matrices are integer-valued, so every result is an exact integer
   real; render it as an int to sidestep cross-compiler float formatting. *)
fun ri (r : real) = Int.toString (Real.round r)

fun showMat name m =
  let
    val rs = M.toRows m
  in
    print (name ^ " (" ^ Int.toString (M.rows m) ^ "x" ^ Int.toString (M.cols m) ^ "):\n");
    List.app
      (fn row =>
         print ("  [ " ^ String.concatWith " " (List.map ri row) ^ " ]\n"))
      rs
  end

val a = M.fromRows [[1.0, 2.0], [3.0, 4.0]]
val b = M.fromRows [[5.0, 6.0], [7.0, 8.0]]

val () = showMat "A" a
val () = showMat "B" b
val () = showMat "A * B" (M.mul (a, b))
val () = showMat "transpose A" (M.transpose a)
val () = print ("det A   = " ^ ri (M.det a) ^ "\n")
val () = print ("trace A = " ^ ri (M.trace a) ^ "\n")
