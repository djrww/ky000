fn main {
  let led = @lib.Ledger::new()
  @lib.num_register(led)
  @lib.perm_register(led)
  @lib.group_register(led)
  @lib.matrix_register(led)
  @lib.poly_register(led)
  @lib.alt_register(led)
  @lib.lie_register(led)
  @lib.sporadic_register(led)
  @lib.sheaf_register(led)
  @lib.hopf_register(led)
  @lib.twolocal_register(led)
  @lib.enrich_register(led)
  println("╔══════════════════════════════════════════════════════════════╗")
  println("║  有限單群 × 作用：公理 A1–A5 與義務自證稽核（MoonBit）        ║")
  println("╚══════════════════════════════════════════════════════════════╝")
  println("")
  println("公理：")
  for ax in @lib.all_axioms() {
    println("  " + ax.tag() + " — " + ax.statement())
  }
  println("")
  println("義務自證逐條：")
  let (p, t) = led.audit(true)
  println("")
  println("通過 " + p.to_string() + "/" + t.to_string())
  println("演算法條數（去重）=" + led.algorithm_count().to_string())
  println("算式分布：")
  for pr in led.calculus_histogram() {
    println("  " + pr.0 + " : " + pr.1.to_string())
  }
}
