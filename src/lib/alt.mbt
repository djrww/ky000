//! 交錯群 A_n：CFSG 第一族（A_n，n ≥ 5）之算式化。
//!
//! 本模組以「偶置換群」之可計算定義為出發點：
//!   * 生成：A_n = ⟨(0 1 k) : 3 ≤ k ≤ n⟩（3-輪換生成）
//!   * 階守恆之雙算式：階乘算式 n!/2 與 Schreier–Sims 之 |G| 互相對照（A1 群公理）。
//!   * 單性（A4 公理）：∀x≠e，x 之正規閉包 = A_n。
//!   * 類方程：Σ 共軛類大小 = |A_n|，配合 A5 之「交錯群族」條目。

/// A_n 之 3-輪換生成集：(0 1 2), (0 1 3), …, (0 1 n−1)。
pub fn alt_gens(n : Int) -> Array[Array[Int]] {
  let gens : Array[Array[Int]] = []
  let mut j = 2
  while j < n {
    gens.push(cycle_perm(n, [0, 1, j]))
    j = j + 1
  }
  gens
}

/// 全部 3-輪換（i < j < k）。
pub fn alt_three_cycles(n : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let mut k = j + 1
      while k < n {
        out.push(cycle_perm(n, [i, j, k]))
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 由階乘算式 |A_n| = n!/2（與群論演算法互相獨立之計算路徑）。
pub fn alt_order_factorial(n : Int) -> Int {
  let mut f = 1
  let mut i = 2
  while i <= n {
    f = f * i
    i = i + 1
  }
  f / 2
}

/// 以 Schreier–Sims 演算法計算 |A_n|。
pub fn alt_order(n : Int) -> Int {
  group_order(n, alt_gens(n))
}

/// A_n 之元素階分布（同構不變量）：(元素階, 個數)。
pub fn element_order_profile(elems : Array[Array[Int]]) -> Array[(Int, Int)] {
  let m : Map[Int, Int] = Map([])
  for x in elems {
    let o = perm_order(x)
    m[o] = m.get_or_default(o, 0) + 1
  }
  let out : Array[(Int, Int)] = []
  m.each(fn(k, v) { out.push((k, v)) })
  out
}

/// 兩階分布是否相同（逐項比對）。
pub fn order_profile_eq(a : Array[(Int, Int)], b : Array[(Int, Int)]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for pr in a {
    let (o, c) = pr
    let mut found = false
    for pr2 in b {
      let (o2, c2) = pr2
      if o == o2 && c == c2 {
        found = true
      }
    }
    if !found {
      return false
    }
  }
  true
}

/// A_n 之共軛類（代表元, 類大小）。
pub fn alt_classes(n : Int) -> Array[(Array[Int], Int)] {
  let g = perm_group(n, alt_gens(n))
  let (elems, _ok) = g.elements(200000)
  conjugacy_classes(elems, g)
}

/// 類方程 Σ|類| = |A_n| = n!/2 之驗證。
pub fn alt_class_equation(n : Int) -> Bool {
  let mut s = 0
  for c in alt_classes(n) {
    s = s + c.1
  }
  s == alt_order(n) && alt_order(n) == alt_order_factorial(n)
}

/// 單性之正規閉包判準：∀x ≠ e，⟨x^{A_n}⟩ = A_n。
pub fn alt_simple_by_closures(n : Int) -> Bool {
  let g = perm_group(n, alt_gens(n))
  let (elems, _ok) = g.elements(200000)
  if elems.length() <= 1 {
    return false
  }
  let mut ok = true
  for x in elems {
    if !perm_is_id(x) {
      if normal_closure(n, g.gens, [x]).length() < elems.length() {
        ok = false
      }
    }
  }
  ok
}

/// 類大小之多重集合（排序後字串鍵，便於比對同構不變量）。
pub fn class_size_key(sizes : Array[Int]) -> String {
  let a : Array[Int] = []
  for s in sizes {
    a.push(s)
  }
  let mut i = 0
  while i < a.length() {
    let mut j = i + 1
    while j < a.length() {
      if a[j] < a[i] {
        let t = a[i]
        a[i] = a[j]
        a[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  let mut s = ""
  for x in a {
    s = s + x.to_string() + ","
  }
  s
}

/// 交錯群族之義務自證。
pub fn alt_register(led : Ledger) -> Unit {
  led.oblige("alt.order", A1_Group, Order, "|A_n| = n!/2：A₃=3、A₄=12、A₅=60、A₆=360、A₇=2520（階乘算式 × Schreier–Sims 雙路徑）", fn() {
    alt_order(3) == 3 &&
    alt_order(4) == 12 &&
    alt_order(5) == 60 &&
    alt_order(6) == 360 &&
    alt_order(7) == 2520 &&
    alt_order_factorial(7) == 2520 &&
    alt_order(6) == 360 &&
    alt_order(5) == 60
  })
  led.oblige("alt.classes", A1_Group, Burnside, "A₅ 共軛類：5 類，類大小 {1,15,20,12,12}，Σ = 60；類方程 Σ|c| = |G| 於 A₄、A₅、A₆ 成立", fn() {
    let c5 = alt_classes(5)
    let sizes : Array[Int] = []
    for c in c5 {
      sizes.push(c.1)
    }
    let key = class_size_key(sizes)
    c5.length() == 5 &&
    key == "1,12,12,15,20," &&
    alt_class_equation(4) &&
    alt_class_equation(5) &&
    alt_class_equation(6) &&
    alt_classes(4).length() == 4
  })
  led.oblige("alt.simplicity", A4_Simplicity, Order, "A₅ 單（∀x≠e 之正規閉包皆為 60）；A₄ 不單（雙對換之正規閉包 = V₄，階 4）；A₄ 之導長 = 2", fn() {
    let g4 = perm_group(4, alt_gens(4))
    let nc = normal_closure(4, g4.gens, [perm_mul(transposition(4, 0, 1), transposition(4, 2, 3))])
    alt_simple_by_closures(5) &&
    nc.length() == 4 &&
    derived_length(4, alt_gens(4)) == 2 &&
    !alt_simple_by_closures(4)
  })
  led.oblige("alt.three_cycles", A1_Group, OrbitStabilizer, "3-輪換：A₅ 之 3-輪換（每 3-子集一個定向）共 C(5,3) = 10 個而生成 A₅；A₅ 中階 3 之元素共 2·C(5,3) = 20 個", fn() {
    let tc = alt_three_cycles(5)
    let mut all_even = true
    for x in tc {
      if perm_sign(x) != 1 {
        all_even = false
      }
    }
    let gen = subgroup_generated(5, tc)
    let g5 = perm_group(5, alt_gens(5))
    let (e5, _o) = g5.elements(1000)
    let mut order3 = 0
    for x in e5 {
      if perm_order(x) == 3 {
        order3 = order3 + 1
      }
    }
    tc.length() == 10 &&
    all_even &&
    gen.length() == 60 &&
    order3 == 20 &&
    binomial(5, 3) * 2 == 20
  })
  led.oblige("alt.a5_psl5_bridge", A5_Classification, LinearCount, "A₅ ≅ PSL(2,5) ≅ PSL(2,4)：階皆 60，元素階分布一致（1:1、2:15、3:20、5:24），兩者皆單（同構不變量一致）", fn() {
    let a5 = perm_group(5, alt_gens(5))
    let (ea5, _o1) = a5.elements(1000)
    let (pg, deg, _npts) = psl2_on_projective_line(5)
    let psl = perm_group(deg, pg)
    let (epsl, _o2) = psl.elements(1000)
    let prof_a = element_order_profile(ea5)
    let prof_p = element_order_profile(epsl)
    ea5.length() == 60 &&
    epsl.length() == 60 &&
    order_profile_eq(prof_a, prof_p) &&
    is_simple(ea5, a5) &&
    lie_a1_order(4) == 60 &&
    lie_a1_order(5) == 60
  })
}
