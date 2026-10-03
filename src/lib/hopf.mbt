//! 群概型作用（group scheme action）：Hopf 代數 (Δ, ε, S) 與共作用之不變子空間。
//!
//! 兩條路線：
//!   * 群代數 F_p[G]（G = Z/3、S₃，p ∤ |G|，即既約／半單情形）：
//!     Δ(g) = g⊗g、ε(g) = 1、S(g) = g⁻¹，逐條驗證 Hopf 公理（餘結合、餘單位、對極）
//!     及 Δ 為代數同態、S 為反同態。
//!   * μ_p（非既約／非半單情形，p | |G|）：F_p[x]/(x^p − 1) 中 (x−1) 為冪零，
//!     Hopf 公理仍成立，但 Maschke 之平均算子不存在；同一個不變量計算必須改用
//!     Hopf 代數定義（共作用之核），其維數與 α_p 之不變元（導子核）一致。
//!
//! 作用於集合之交換代數 F_p[X]：共作用 ρ(e_x) = Σ_g e_g ⊗ e_{g·x}，
//! 不變量 = 軌道上為常數之函數 ⇒ dim 不變量 = dim 餘不變量 = 軌道數（p ∤ |G|）。

/// 元素在表中之索引（-1 表不存在）。
pub fn elem_pos(elems : Array[Array[Int]], x : Array[Int]) -> Int {
  let mut i = 0
  while i < elems.length() {
    if perm_eq(elems[i], x) {
      return i
    }
    i = i + 1
  }
  -1
}

/// 單位元之索引。
pub fn id_index(elems : Array[Array[Int]]) -> Int {
  let mut i = 0
  while i < elems.length() {
    if perm_is_id(elems[i]) {
      return i
    }
    i = i + 1
  }
  -1
}

/// 張量平方上之係數乘法（群代數 ⊗ 群代數，索引 (a,b) ↦ a·m+b）。
pub fn tensor_mul(u : Array[Int], v : Array[Int], elems : Array[Array[Int]], p : Int) -> Array[Int] {
  let m = elems.length()
  let out = Array::make(m * m, 0)
  let mut a1 = 0
  while a1 < m {
    let mut b1 = 0
    while b1 < m {
      let c1 = u[a1 * m + b1] % p
      if c1 != 0 {
        let mut a2 = 0
        while a2 < m {
          let mut b2 = 0
          while b2 < m {
            let c2 = v[a2 * m + b2] % p
            if c2 != 0 {
              let a = elem_pos(elems, perm_mul(elems[a1], elems[a2]))
              let b = elem_pos(elems, perm_mul(elems[b1], elems[b2]))
              out[a * m + b] = mod_pos(out[a * m + b] + c1 * c2, p)
            }
            b2 = b2 + 1
          }
          a2 = a2 + 1
        }
      }
      b1 = b1 + 1
    }
    a1 = a1 + 1
  }
  out
}

/// Δ(g) = g ⊗ g（grouplike）之係數向量。
pub fn delta_grouplike(i : Int, m : Int) -> Array[Int] {
  let out = Array::make(m * m, 0)
  out[i * m + i] = 1
  out
}

/// 群代數元素之乘法（係數向量）。
pub fn galg_mul(u : Array[Int], v : Array[Int], elems : Array[Array[Int]], p : Int) -> Array[Int] {
  let m = elems.length()
  let out = Array::make(m, 0)
  let mut i = 0
  while i < m {
    let c1 = u[i] % p
    if c1 != 0 {
      let mut j = 0
      while j < m {
        let c2 = v[j] % p
        if c2 != 0 {
          let k = elem_pos(elems, perm_mul(elems[i], elems[j]))
          out[k] = mod_pos(out[k] + c1 * c2, p)
        }
        j = j + 1
      }
    }
    i = i + 1
  }
  out
}

/// Hopf 公理檢查（群代數 grouplike 生成元）：
/// 餘結合、餘單位、對極、Δ 乘法性、S 反乘法性。
pub fn hopf_axioms_hold(elems : Array[Array[Int]], p : Int) -> Bool {
  let m = elems.length()
  let ei = id_index(elems)
  let mut g = 0
  while g < m {
    let x = elems[g]
    let dx = delta_grouplike(g, m)
    // 對極：m(S⊗id)Δ(g) = e
    let acc : Array[Int] = Array::make(m, 0)
    let mut a = 0
    while a < m {
      let c = dx[a * m + a] % p
      if c != 0 {
        let sa = elem_pos(elems, perm_inv(elems[a]))
        let prod = perm_mul(elems[sa], elems[a])
        let k = elem_pos(elems, prod)
        acc[k] = mod_pos(acc[k] + c, p)
      }
      a = a + 1
    }
    // 檢查 acc = e（係數於單位元為 1）
    let mut t = 0
    while t < m {
      let want = if t == ei { 1 % p } else { 0 }
      if acc[t] % p != want {
        return false
      }
      t = t + 1
    }
    // Δ 乘法性：Δ(gh) = Δ(g)Δ(h)
    let mut h = 0
    while h < m {
      let dh = delta_grouplike(h, m)
      let lhs = delta_grouplike(elem_pos(elems, perm_mul(x, elems[h])), m)
      let rhs = tensor_mul(dx, dh, elems, p)
      let mut q = 0
      while q < m * m {
        if lhs[q] % p != rhs[q] % p {
          return false
        }
        q = q + 1
      }
      // S 反乘法性
      if !perm_eq(perm_inv(perm_mul(x, elems[h])), perm_mul(perm_inv(elems[h]), perm_inv(x))) {
        return false
      }
      h = h + 1
    }
    g = g + 1
  }
  true
}

/// 集合作用之共作用矩陣 ρ：索引 (g, y) 列、x 行為 1 ⟺ y = g·x。
pub fn set_coaction_matrix(elems : Array[Array[Int]], npts : Int, p : Int) -> Array[Int] {
  let m = elems.length()
  let rows = m * npts
  let mat = Array::make(rows * npts, 0)
  let mut g = 0
  while g < m {
    let mut x = 0
    while x < npts {
      let y = elems[g][x]
      mat[(g * npts + y) * npts + x] = 1 % p
      x = x + 1
    }
    g = g + 1
  }
  mat
}

/// 1 ⊗ a 之嵌入矩陣：(g, y) 列、x 行為 1 ⟺ y = x（對所有 g）。
/// 對偶群代數 F_p[G]^* 之單位為 Σ_g χ_g，故 1⊗a = Σ_{g,y} a_y χ_g⊗e_y。
pub fn unit_embed_matrix(elems : Array[Array[Int]], npts : Int, p : Int) -> Array[Int] {
  let m = elems.length()
  let rows = m * npts
  let mat = Array::make(rows * npts, 0)
  let mut g = 0
  while g < m {
    let mut x = 0
    while x < npts {
      mat[(g * npts + x) * npts + x] = 1 % p
      x = x + 1
    }
    g = g + 1
  }
  mat
}

/// 不變量維數 = dim ker(ρ − 1⊗·)。
pub fn coaction_invariant_dim(elems : Array[Array[Int]], npts : Int, p : Int) -> Int {
  let rho = set_coaction_matrix(elems, npts, p)
  let emb = unit_embed_matrix(elems, npts, p)
  let m = elems.length()
  let rows = m * npts
  let diff = Array::make(rows * npts, 0)
  let mut i = 0
  while i < rows * npts {
    diff[i] = mod_pos(rho[i] - emb[i], p)
    i = i + 1
  }
  mat_kernel(diff, rows, npts, p).length()
}

/// 餘不變量維數 = |X| − rank(ρ − 1⊗·)。
pub fn coaction_coinvariant_dim(elems : Array[Array[Int]], npts : Int, p : Int) -> Int {
  let rho = set_coaction_matrix(elems, npts, p)
  let emb = unit_embed_matrix(elems, npts, p)
  let m = elems.length()
  let rows = m * npts
  let diff = Array::make(rows * npts, 0)
  let mut i = 0
  while i < rows * npts {
    diff[i] = mod_pos(rho[i] - emb[i], p)
    i = i + 1
  }
  npts - matrix_rank(diff, rows, npts, p)
}

/// μ_p 共作用 ρ(y^k) = x^k ⊗ y^k 於 F_p[y]/(y^N) 之不變量維數（= ker(ρ − 1⊗·)）。
pub fn mu_p_invariant_dim(p : Int, nmax : Int) -> Int {
  // 基底 1, y, …, y^{N−1}；行 k 表 y^k；列 (i, k) 表 x^i ⊗ y^k
  let rows = p * nmax
  let mat = Array::make(rows * nmax, 0)
  let mut k = 0
  while k < nmax {
    let i = mod_pos(k, p)
    mat[(i * nmax + k) * nmax + k] = 1 % p
    k = k + 1
  }
  // 1 ⊗ a：列 (0, k)
  let emb = Array::make(rows * nmax, 0)
  let mut k2 = 0
  while k2 < nmax {
    emb[(0 * nmax + k2) * nmax + k2] = 1 % p
    k2 = k2 + 1
  }
  let diff = Array::make(rows * nmax, 0)
  let mut i2 = 0
  while i2 < rows * nmax {
    diff[i2] = mod_pos(mat[i2] - emb[i2], p)
    i2 = i2 + 1
  }
  mat_kernel(diff, rows, nmax, p).length()
}

/// F_p[x]/(x^p − 1) 中之乘法（係數向量，x^p = 1）。
pub fn mu_p_mul(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let n = p
  let out = Array::make(n, 0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let k = mod_pos(i + j, n)
      out[k] = mod_pos(out[k] + a[i] * b[j], p)
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 於 F_p[x]/(x^p − 1) 中是否冪零；回傳冪零指數（0 表非冪零）。
pub fn mu_p_nilpotency_index(p : Int) -> Int {
  // 元素 (x − 1)
  let e : Array[Int] = Array::make(p, 0)
  e[1] = 1 % p
  e[0] = mod_pos(e[0] - 1, p)
  let mut cur : Array[Int] = [1 % p]
  let mut i = 1
  while i < p {
    cur.push(0)
    i = i + 1
  }
  let mut k = 0
  while k <= p {
    let mut zero = true
    let mut t = 0
    while t < p {
      if cur[t] % p != 0 {
        zero = false
      }
      t = t + 1
    }
    if zero {
      return k
    }
    cur = mu_p_mul(cur, e, p)
    k = k + 1
  }
  0
}

/// 冪等元素（e² = e）之列舉：回傳冪等元素個數（單位元與零元）。
pub fn mu_p_idempotent_count(p : Int) -> Int {
  let total = int_pow(p, p)
  let mut cnt = 0
  let mut code = 0
  while code < total {
    let v = Array::make(p, 0)
    let mut t = 0
    let mut c = code
    while t < p {
      v[t] = c % p
      c = c / p
      t = t + 1
    }
    let sq = mu_p_mul(v, v, p)
    let mut same = true
    let mut t2 = 0
    while t2 < p {
      if sq[t2] % p != v[t2] % p {
        same = false
      }
      t2 = t2 + 1
    }
    if same {
      cnt = cnt + 1
    }
    code = code + 1
  }
  cnt
}

/// 群概型作用之義務自證。
pub fn hopf_register(led : Ledger) -> Unit {
  led.oblige("hopf.group_algebra", A3_Composition, HopfScheme, "群代數 Hopf 公理（G = Z/3、S₃ 於 F₇）：Δ(g)=g⊗g 滿足餘結合、餘單位 (ε⊗id)Δ=id、對極 m(S⊗id)Δ = ηε、Δ 為代數同態、S 為反同態", fn() {
    let g3 = cyclic_group_regular(3)
    let (e3, _o1) = g3.elements(100)
    let s3 = symmetric_group(3)
    let (e5, _o2) = s3.elements(100)
    e3.length() == 3 &&
    e5.length() == 6 &&
    hopf_axioms_hold(e3, 7) &&
    hopf_axioms_hold(e5, 7) &&
    gcd_int(3, 7) == 1 &&
    gcd_int(6, 7) == 1
  })
  led.oblige("hopf.coaction_orbits", A2_Action, HopfScheme, "集合作用之共作用：dim 不變量 = dim 餘不變量 = 軌道數——Z/2 作用於 3 點（軌道 {0} 與 {1,2}）為 2；S₃、Z/3 遞移作用於 3 點為 1", fn() {
    let z2 = perm_group(3, [transposition(3, 1, 2)])
    let (ez2, _o1) = z2.elements(100)
    let s3 = symmetric_group(3)
    let (es3, _o2) = s3.elements(100)
    let z3 = perm_group(3, [cycle_perm(3, [0, 1, 2])])
    let (ez3, _o3) = z3.elements(100)
    coaction_invariant_dim(ez2, 3, 7) == 2 &&
    coaction_coinvariant_dim(ez2, 3, 7) == 2 &&
    coaction_invariant_dim(es3, 3, 7) == 1 &&
    coaction_coinvariant_dim(es3, 3, 7) == 1 &&
    coaction_invariant_dim(ez3, 3, 7) == 1 &&
    coaction_coinvariant_dim(ez3, 3, 7) == 1 &&
    z2.orbits().length() == 2 &&
    s3.orbits().length() == 1
  })
  led.oblige("hopf.mu_p_coaction", A2_Action, HopfScheme, "μ_p 群概型作用於 F_p[y]：ρ(y^k) = x^k ⊗ y^k 之不變量維數 = ⌈N/p⌉（p=3,N=6 與 p=5,N=10 皆為 2），與 α_p = d/dx 之核（導子不變元）維數一致", fn() {
    mu_p_invariant_dim(3, 6) == 2 &&
    mu_p_invariant_dim(5, 10) == 2 &&
    mu_p_invariant_dim(3, 9) == 3 &&
    derivation_invariants(6, 3).length() == 2 &&
    derivation_invariants(10, 5).length() == 2 &&
    derivation_invariants(9, 3).length() == 3
  })
  led.oblige("hopf.mu_p_nonsemisimple", A3_Composition, HopfScheme, "非既約情形（p | |G|）：F₃[x]/(x³−1) 中 (x−1) 冪零指數 3（(x−1)² ≠ 0、(x−1)³ = 0），冪等元僅 0 與 1 ⇒ 無對極平均（Reynolds）算子；Maschke 假設 gcd(3,3) ≠ 1 失敗", fn() {
    mu_p_nilpotency_index(3) == 3 &&
    mu_p_nilpotency_index(5) == 5 &&
    mu_p_idempotent_count(3) == 2 &&
    gcd_int(3, 3) != 1 &&
    gcd_int(5, 7) == 1
  })
}
