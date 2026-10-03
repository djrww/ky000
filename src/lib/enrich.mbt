//! 富化作用（enriched action）：Maschke／Reynolds 平均算子、特徵標內積、
//! Frobenius 互反與不變環（Molien 級數）之互相對照。
//!
//! 「富化」= 在作用之上再賦予線性（表示）結構：G → GL(V)。
//!   * Reynolds 算子 π = (1/|G|) Σ_g ρ(g)（p ∤ |G| 時存在）：π² = π、im π = V^G、
//!     rank π = dim V^G；由此給出 Maschke 定理之構造性證明（不變子空間必有互補）。
//!   * 特徵標內積 ⟨χ,ψ⟩ = (1/|G|) Σ_c |c| χ(c)ψ(c)（有理數精確運算）：既約特徵標正交。
//!   * Frobenius 互反：⟨Ind_H^G 1, χ⟩ = ⟨1, Res χ⟩ = χ 中平凡表示之重數。
//!   * Molien 級數 ↔ 對稱冪不變子空間維數 ↔ 整數分拆數（三重獨立算式）。

/// 由 n 點置換之生成元建立置換矩陣（每個生成元一個 n×n 矩陣）。
pub fn perm_rep_matrices(gens : Array[Array[Int]], n : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  for g in gens {
    let m = Array::make(n * n, 0)
    let mut i = 0
    while i < n {
      m[g[i] * n + i] = 1
      i = i + 1
    }
    out.push(m)
  }
  out
}

/// 展平之置換矩陣列（給定型別 g 之所有元素之矩陣）。
pub fn perm_matrices_of_group(g : PermGroup) -> Array[Array[Int]] {
  let (elems, _ok) = g.elements(100000)
  let out : Array[Array[Int]] = []
  for x in elems {
    let m = Array::make(g.n * g.n, 0)
    let mut i = 0
    while i < g.n {
      m[x[i] * g.n + i] = 1
      i = i + 1
    }
    out.push(m)
  }
  out
}

/// Reynolds 算子 π = (1/|G|) Σ ρ(g)（須 gcd(|G| mod p, p) = 1）。
pub fn reynolds_operator(mats : Array[Array[Int]], dim : Int, p : Int) -> Array[Int] {
  let mut acc = Array::make(dim * dim, 0)
  for a in mats {
    acc = mat_add(acc, a, p)
  }
  mat_scale(acc, mod_inv(mats.length() % p, p), p)
}

/// 矩陣之跡（mod p）。
pub fn mat_trace(m : Array[Int], n : Int, p : Int) -> Int {
  let mut s = 0
  let mut i = 0
  while i < n {
    s = s + m[i * n + i]
    i = i + 1
  }
  mod_pos(s, p)
}

/// 有理數特徵標內積：(1/|G|) Σ_c |c| χ(c) ψ(c)。
pub fn char_inner_product(
  sizes : Array[Int],
  chi : Array[Int],
  psi : Array[Int],
  order : Int
) -> Rat {
  let mut acc = rat(0, 1)
  let mut i = 0
  while i < sizes.length() {
    acc = rat_add(acc, rat(sizes[i] * chi[i] * psi[i], 1))
    i = i + 1
  }
  rat_div(acc, rat(order, 1))
}

/// 分拆數 p_k(n)：n 之分拆中每部分 ≤ k（動態規劃）。
pub fn partitions_with_parts_le(n : Int, k : Int) -> Int {
  let dp = Array::make(n + 1, 0)
  dp[0] = 1
  let mut part = 1
  while part <= k {
    let mut s = part
    while s <= n {
      dp[s] = dp[s] + dp[s - part]
      s = s + 1
    }
    part = part + 1
  }
  dp[n]
}

/// 富化作用之義務自證。
pub fn enrich_register(led : Ledger) -> Unit {
  led.oblige("enrich.reynolds_projection", A2_Action, Enriched, "Reynolds 算子（3 維置換表示，S₃ 於 F₅）：π² = π、π 與作用交換、rank π = dim V^G = 1（常數函數）、π 固定不變向量", fn() {
    let s3 = symmetric_group(3)
    let mats = perm_matrices_of_group(s3)
    let pi = reynolds_operator(mats, 3, 5)
    let pi2 = mat_mul(pi, pi, 3, 5)
    let mut commutes = true
    for a in mats {
      if !mat_eq(mat_mul(pi, a, 3, 5), mat_mul(a, pi, 3, 5)) {
        commutes = false
      }
    }
    let fixdim = joint_fixed_space(mats, 3, 5).length()
    let ones : Array[Int] = [1, 1, 1]
    let _ = ones
    let v : Array[Int] = [1, 1, 1]
    let pv = mat_vec(pi, v, 3, 5)
    fixdim == 1 &&
    mat_eq(pi2, pi) &&
    commutes &&
    matrix_rank(pi, 3, 3, 5) == 1 &&
    pv[0] == 1 && pv[1] == 1 && pv[2] == 1 &&
    mats.length() == 6
  })
  led.oblige("enrich.maschke_complement", A2_Action, Enriched, "Maschke 定理之構造性驗證：G = ⟨diag(1,−1)⟩ 作用於 F₅²（p ∤ |G|），不變直線 W = ⟨e₁⟩ 之 π-互補 ker π = ⟨e₂⟩ 亦不變，且 V = W ⊕ ker π", fn() {
    let rho : Array[Int] = [1, 0, 0, mod_pos(-1, 5)]
    let id = mat_id(2, 5)
    let pi0 = mat_add(id, rho, 5)
    let pi = mat_scale(pi0, mod_inv(2, 5), 5)
    let e1 : Array[Int] = [1, 0]
    let e2 : Array[Int] = [0, 1]
    let pe2 = mat_vec(pi, e2, 2, 5)
    let re2 = mat_vec(rho, e2, 2, 5)
    let pre2 = mat_vec(pi, re2, 2, 5)
    mat_eq(mat_mul(pi, pi, 2, 5), pi) &&
    mat_eq(mat_vec(pi, e1, 2, 5), e1) &&
    pe2[0] == 0 && pe2[1] == 0 &&
    re2[0] == 0 &&
    pre2[0] == 0 && pre2[1] == 0 &&
    matrix_rank(pi, 2, 2, 5) == 1 &&
    mat_trace(pi, 2, 5) == 1
  })
  led.oblige("enrich.char_orthogonality", A1_Group, Enriched, "特徵標正交（有理數精確）：S₃ 之平凡、符號、標準（2 維）特徵標兩兩內積 = δ；Σ d² = 1+1+4 = 6 = |S₃|；A₅ 之 5 個次數 1,3,3,4,5 滿足 Σ d² = 60", fn() {
    let sizes : Array[Int] = [1, 3, 2]
    let triv : Array[Int] = [1, 1, 1]
    let sign : Array[Int] = [1, -1, 1]
    let std : Array[Int] = [2, 0, -1]
    let n_tt = char_inner_product(sizes, triv, triv, 6)
    let n_ss = char_inner_product(sizes, sign, sign, 6)
    let n_dd = char_inner_product(sizes, std, std, 6)
    let c_ts = char_inner_product(sizes, triv, sign, 6)
    let c_td = char_inner_product(sizes, triv, std, 6)
    let c_sd = char_inner_product(sizes, sign, std, 6)
    let deg_a5 : Array[Int] = [1, 3, 3, 4, 5]
    let mut sum = 0
    for d in deg_a5 {
      sum = sum + d * d
    }
    rat_eq(n_tt, rat(1, 1)) &&
    rat_eq(n_ss, rat(1, 1)) &&
    rat_eq(n_dd, rat(1, 1)) &&
    rat_eq(c_ts, rat(0, 1)) &&
    rat_eq(c_td, rat(0, 1)) &&
    rat_eq(c_sd, rat(0, 1)) &&
    sum == 60 &&
    alt_classes(5).length() == 5
  })
  led.oblige("enrich.frobenius", A3_Composition, Enriched, "Frobenius 互反：Ind_{S₂}^{S₃} 1 之特徵標 χ_perm = χ_triv + χ_std（3 維），⟨χ_perm,χ_triv⟩ = 1、⟨χ_perm,χ_std⟩ = 1、⟨χ_perm,χ_perm⟩ = 2，且 [S₃:S₂] = 3", fn() {
    let sizes : Array[Int] = [1, 3, 2]
    let perm : Array[Int] = [3, 1, 0]
    let triv : Array[Int] = [1, 1, 1]
    let std : Array[Int] = [2, 0, -1]
    let a = char_inner_product(sizes, perm, triv, 6)
    let b = char_inner_product(sizes, perm, std, 6)
    let c = char_inner_product(sizes, perm, perm, 6)
    let mut decomposed = true
    let mut i = 0
    while i < 3 {
      if perm[i] != triv[i] + std[i] {
        decomposed = false
      }
      i = i + 1
    }
    rat_eq(a, rat(1, 1)) &&
    rat_eq(b, rat(1, 1)) &&
    rat_eq(c, rat(2, 1)) &&
    decomposed &&
    group_order(3, symmetric_group(3).gens) / group_order(2, symmetric_group(2).gens) == 3
  })
  led.oblige("enrich.invariants_molien", A2_Action, Molien, "不變子空間之三重算式一致（S₃ 於 3 變數，d ≤ 6）：對稱冪不變維數 = Molien 級數係數 = 部分 ≤ 3 之整數分拆數（1,1,2,3,4,5,7）", fn() {
    let s3 = symmetric_group(3)
    let gens = s3.gens
    let genmats = perm_rep_matrices(gens, 3)
    let (series, ok) = molien_series(gens, 3, 7)
    if !ok {
      return false
    }
    let mut all_ok = true
    let mut d = 0
    while d <= 6 {
      let inv = invariant_dim_linalg(genmats, 3, d, 5)
      let part = partitions_with_parts_le(d, 3)
      let coef = series[d]
      if inv != part {
        all_ok = false
      }
      if !rat_eq(coef, rat(part, 1)) {
        all_ok = false
      }
      d = d + 1
    }
    all_ok &&
    partitions_with_parts_le(6, 3) == 7 &&
    partitions_with_parts_le(4, 3) == 4
  })
}
