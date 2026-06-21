fun runAllSuites () =
  ( Harness.reset ()
  ; MatrixTests.run ()
  ; Harness.run () )

fun main () =
  OS.Process.exit
    (if runAllSuites () then OS.Process.success else OS.Process.failure)
