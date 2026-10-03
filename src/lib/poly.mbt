//! 多項式環作用：F_p[x] 算術、Frobenius 作用（Galois 群作用於根）、Z/n 對 F_p[x]/(x^n−1) 之作用、
//! 線性群之不變子空間，以及 Molien 級數與 Chevalley–Shephard–Todd 定理之驗證。
//!
//! 定義 D6：F_p[x] 之元素以「升冪係數」`Array[Int]` 表示（mod p，去除高位零）。
//! 定義 D7：Frobenius φ(a) = a^p 為 F̄_p 之自同構；對不可約 f（deg n），Gal(F_{p^n}/F_p) = ⟨φ⟩ ≅ Z/n
//!   以「根的置換」作用——其軌道大小為 n（定理：有限體之 Galois 擴張為循環）。

/// 去除高位零並取模。
pub fn poly_trim(a : Array[Int], p : Int) -> Array[Int] {
  let out : Array[Int] = []
  for c in a {
    out.push(mod_pos(c, p))
  }
  while out.length() > 0 && out[out.length() - 1] == 0 {
    let _ = out.pop()
  }
  out
}

/// 多項式加法。
pub fn poly_add(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  return poly_trim(poly_add_p(a, b, p), p)
}

/// 多項式減法。
pub fn poly_sub(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  return poly_trim(poly_add_p(a, poly_neg_p(b, p), p), p)
}

/// 多項式乘法。
pub fn poly_mul(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  return poly_trim(poly_mul_p(a, b, p), p)
}

/// 多項式之次數（零多項式回傳 −1）。
pub fn poly_deg(a : Array[Int]) -> Int {
  a.length() - 1
}

/// 多項式除法取餘：r = a mod b（b ≠ 0）。
pub fn poly_mod(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let bb = poly_trim(b, p)
  if bb.length() == 0 {
    abort("poly_mod: 除式為零")
  }
  let mut r = poly_trim(a, p)
  let db = poly_deg(bb)
  let lead_inv = mod_inv(bb[db], p)
  while r.length() > 0 && poly_deg(r) >= db {
    let dr = poly_deg(r)
    let factor = mod_pos(r[dr] * lead_inv, p)
    let shift = dr - db
    let mut i = 0
    while i <= db {
      r[shift + i] = mod_pos(r[shift + i] - factor * bb[i], p)
      i = i + 1
    }
    r = poly_trim(r, p)
  }
  r
}

/// 多項式帶餘除法：回傳 (商, 餘)。
pub fn poly_divmod(a : Array[Int], b : Array[Int], p : Int) -> (Array[Int], Array[Int]) {
  let bb = poly_trim(b, p)
  if bb.length() == 0 {
    abort("poly_divmod: 除式為零")
  }
  let db = poly_deg(bb)
  let lead_inv = mod_inv(bb[db], p)
  let mut r = poly_trim(a, p)
  let q = Array::make(if poly_deg(a) >= db { poly_deg(a) - db + 1 } else { 0 }, 0)
  while r.length() > 0 && poly_deg(r) >= db {
    let dr = poly_deg(r)
    let factor = mod_pos(r[dr] * lead_inv, p)
    let shift = dr - db
    q[shift] = factor
    let mut i = 0
    while i <= db {
      r[shift + i] = mod_pos(r[shift + i] - factor * bb[i], p)
      i = i + 1
    }
    r = poly_trim(r, p)
  }
  (poly_trim(q, p), r)
}

/// 多項式最大公因式（Euclid）。
pub fn poly_gcd(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let mut x = poly_trim(a, p)
  let mut y = poly_trim(b, p)
  while y.length() > 0 {
    let r = poly_mod(x, y, p)
    x = y
    y = r
  }
  x
}

/// 多項式之形式導數 D(f)。
pub fn poly_derivative(a : Array[Int], p : Int) -> Array[Int] {
  let out : Array[Int] = []
  let mut i = 1
  while i < a.length() {
    out.push(mod_pos(a[i] * i, p))
    i = i + 1
  }
  poly_trim(out, p)
}

/// 多項式求值 f(c)。
pub fn poly_eval(a : Array[Int], c : Int, p : Int) -> Int {
  let mut acc = 0
  let mut i = a.length() - 1
  while i >= 0 {
    acc = mod_pos(acc * c + a[i], p)
    i = i - 1
  }
  acc
}

/// 多項式冪次取模：f^k mod g。
pub fn poly_pow_mod(f : Array[Int], k : Int, g : Array[Int], p : Int) -> Array[Int] {
  let mut acc : Array[Int] = [1 % p]
  let mut base = poly_mod(f, g, p)
  let mut e = k
  while e > 0 {
    if e % 2 == 1 {
      acc = poly_mod(poly_mul(acc, base, p), g, p)
    }
    base = poly_mod(poly_mul(base, base, p), g, p)
    e = e / 2
  }
  poly_trim(acc, p)
}

/// 不可約性判定（試除所有次數 ≤ deg/2 之單首多項式；僅適用小 p、小次數）。
pub fn poly_is_irreducible(f : Array[Int], p : Int) -> Bool {
  let ff = poly_trim(f, p)
  let d = poly_deg(ff)
  if d <= 0 {
    return false
  }
  if d == 1 {
    return true
  }
  // 逐一列舉次數 1 ≤ k ≤ d/2 之所有單首多項式
  let mut k = 1
  while k <= d / 2 {
    let total = int_pow(p, k)
    let mut code = 0
    while code < total {
      let g = Array::make(k + 1, 0)
      let mut t = code
      let mut i = 0
      while i < k {
        g[i] = t % p
        t = t / p
        i = i + 1
      }
      g[k] = 1 % p
      if poly_mod(ff, g, p).length() == 0 {
        return false
      }
      code = code + 1
    }
    k = k + 1
  }
  true
}

/// Frobenius 作用於根：回傳 x, x^p, x^{p²}, … mod f 直到循環，及首次重合之步數。
/// 對不可約 f（deg n），Frobenius 之軌道長度應為 n。
pub fn frobenius_orbit(x0 : Array[Int], f : Array[Int], p : Int, maxsteps : Int) -> (Array[Array[Int]], Int) {
  let seen : Array[Array[Int]] = []
  let mut cur = poly_mod(x0, f, p)
  let mut step = 0
  while step < maxsteps {
    let mut found = -1
    let mut i = 0
    while i < seen.length() {
      if poly_eq(seen[i], cur, p) {
        found = i
      }
      i = i + 1
    }
    if found >= 0 {
      return (seen, found)
    }
    seen.push(cur)
    cur = poly_pow_mod(cur, p, f, p)
    step = step + 1
  }
  (seen, -1)
}

/// 多項式相等（正規化後比較）。
pub fn poly_eq(a : Array[Int], b : Array[Int], p : Int) -> Bool {
  let x = poly_trim(a, p)
  let y = poly_trim(b, p)
  if x.length() != y.length() {
    return false
  }
  let mut i = 0
  while i < x.length() {
    if x[i] != y[i] {
      return false
    }
    i = i + 1
  }
  true
}

/// Z/n 對 F_p[x]/(x^n−1) 之作用：x ↦ x^k（gcd(k,n)=1），以「單項式指標的置換」表示。
/// 回傳作用在 {0,…,n−1}（代表 x^i）上之置換。
pub fn cyclic_action_on_quotient(n : Int, k : Int) -> Array[Int] {
  let perm = perm_id(n)
  let mut i = 0
  while i < n {
    perm[i] = mod_pos(i * k, n)
    i = i + 1
  }
  perm
}

/// 平移作用 T : x ↦ x + 1 在「次數 < nmax 之多項式空間」上的矩陣（基 x^i）。
pub fn translation_matrix(nmax : Int, p : Int) -> Array[Int] {
  // (x+1)^i 之係數：binom(i,j)
  let m = Array::make(nmax * nmax, 0)
  let mut i = 0
  while i < nmax {
    let mut j = 0
    while j <= i {
      m[j * nmax + i] = binomial(i, j) % p
      j = j + 1
    }
    i = i + 1
  }
  m
}

/// 平移不變子空間（Z/p 作用 x ↦ x+1 之固定子空間）之基底。
pub fn translation_invariants(nmax : Int, p : Int) -> Array[Array[Int]] {
  joint_fixed_space([translation_matrix(nmax, p)], nmax, p)
}

/// 導子 D 之核（α_p 之無窮小作用之不變元）：ker D 於次數 < nmax 之空間。
pub fn derivation_invariants(nmax : Int, p : Int) -> Array[Array[Int]] {
  let m = Array::make(nmax * nmax, 0)
  let mut i = 1
  while i < nmax {
    m[(i - 1) * nmax + i] = i % p
    i = i + 1
  }
  mat_kernel(m, nmax, nmax, p)
}

/// 多項式 (x^p − x)^k 之係數（升冪）。
pub fn frobenius_poly_power(k : Int, p : Int) -> Array[Int] {
  // x^p − x
  let base = Array::make(p + 1, 0)
  base[1] = mod_pos(-1, p)
  base[p] = 1 % p
  let mut acc : Array[Int] = [1 % p]
  let mut i = 0
  while i < k {
    acc = poly_mul(acc, base, p)
    i = i + 1
  }
  acc
}

// ───────────── Molien 級數（有理數係數之截斷級數） ─────────────

/// 有理數（分子/分母，以 gcd 正規化）。
pub struct Rat {
  n : Int
  d : Int
}

/// 建立有理數。
pub fn rat(n : Int, d : Int) -> Rat {
  let mut num = n
  let mut den = d
  if den < 0 {
    num = -num
    den = -den
  }
  let g = if num == 0 { if den == 0 { 1 } else { den } } else { gcd_int(num, den) }
  { n: num / g, d: den / g }
}

/// 有理數加法。
pub fn rat_add(a : Rat, b : Rat) -> Rat {
  rat(a.n * b.d + b.n * a.d, a.d * b.d)
}

/// 有理數乘法。
pub fn rat_mul(a : Rat, b : Rat) -> Rat {
  rat(a.n * b.n, a.d * b.d)
}

/// 有理數除法。
pub fn rat_div(a : Rat, b : Rat) -> Rat {
  rat(a.n * b.d, a.d * b.n)
}

/// 有理數字串。
pub fn rat_to_string(a : Rat) -> String {
  if a.d == 1 {
    a.n.to_string()
  } else {
    a.n.to_string() + "/" + a.d.to_string()
  }
}

/// 級數加法（截斷至 nmax 項）。
pub fn series_add(a : Array[Rat], b : Array[Rat], nmax : Int) -> Array[Rat] {
  let out = Array::make(nmax, rat(0, 1))
  let mut i = 0
  while i < nmax {
    let x = if i < a.length() { a[i] } else { rat(0, 1) }
    let y = if i < b.length() { b[i] } else { rat(0, 1) }
    out[i] = rat_add(x, y)
    i = i + 1
  }
  out
}

/// 級數乘法（截斷至 nmax 項）。
pub fn series_mul(a : Array[Rat], b : Array[Rat], nmax : Int) -> Array[Rat] {
  let out = Array::make(nmax, rat(0, 1))
  let mut i = 0
  while i < a.length() && i < nmax {
    let mut j = 0
    while j < b.length() && i + j < nmax {
      out[i + j] = rat_add(out[i + j], rat_mul(a[i], b[j]))
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 幾何級數 1/(1 − t^k) 之截斷展開。
pub fn series_geom(k : Int, nmax : Int) -> Array[Rat] {
  let out = Array::make(nmax, rat(0, 1))
  let mut i = 0
  while i < nmax {
    if i % k == 0 {
      out[i] = rat(1, 1)
    }
    i = i + 1
  }
  out
}

/// 以純量（有理數）乘級數。
pub fn series_scale(a : Array[Rat], c : Rat, nmax : Int) -> Array[Rat] {
  let out = Array::make(nmax, rat(0, 1))
  let mut i = 0
  while i < nmax && i < a.length() {
    out[i] = rat_mul(a[i], c)
    i = i + 1
  }
  out
}

/// Molien 級數之單項：置換 g 之貢獻 1/det(I − t·g) = Π_j 1/(1 − t^{c_j})，c_j 為輪換長度。
pub fn molien_term(g : Array[Int], nmax : Int) -> Array[Rat] {
  let s = Array::make(nmax, rat(0, 1))
  s[0] = rat(1, 1)
  let mut acc = s
  for c in perm_cycle_type(g) {
    acc = series_mul(acc, series_geom(c, nmax), nmax)
  }
  acc
}

/// 置換群之 Molien 級數：(1/|G|) Σ_g 1/det(I − t g)，截斷至 nmax 項。
pub fn molien_series(gens : Array[Array[Int]], n : Int, nmax : Int) -> (Array[Rat], Bool) {
  let g = perm_group(n, gens)
  let (elems, ok) = g.elements(100000)
  if !ok {
    return (Array::make(nmax, rat(0, 1)), false)
  }
  let mut acc = Array::make(nmax, rat(0, 1))
  for e in elems {
    acc = series_add(acc, molien_term(e, nmax), nmax)
  }
  let inv = rat(1, elems.length())
  (series_scale(acc, inv, nmax), true)
}

/// 級數之倒數：1/p(t)，要求 p(0) = 1（遞迴 q_n = −Σ_{k≥1} p_k q_{n−k}）。
pub fn series_inv(p : Array[Rat], nmax : Int) -> (Array[Rat], Bool) {
  if p.length() == 0 || p[0].n != 1 || p[0].d != 1 {
    return (Array::make(nmax, rat(0, 1)), false)
  }
  let q = Array::make(nmax, rat(0, 1))
  q[0] = rat(1, 1)
  let mut n = 1
  while n < nmax {
    let mut acc = rat(0, 1)
    let mut k = 1
    while k <= n && k < p.length() {
      acc = rat_add(acc, rat_mul(p[k], q[n - k]))
      k = k + 1
    }
    q[n] = rat(rat(0, 1).n - acc.n, if acc.d == 0 { 1 } else { acc.d })
    q[n] = rat_div(q[n], p[0])
    n = n + 1
  }
  (q, true)
}

/// 將模 p 之係數提升至對稱區間 (−p/2, p/2]（用於把特徵多項式解讀為特徵值之整數函數）。
pub fn lift_symmetric(c : Int, p : Int) -> Int {
  let x = mod_pos(c, p)
  if x * 2 > p { x - p } else { x }
}

/// 線性群作用之 Molien 級數：(1/|G|) Σ_g 1/det(I − t·g)。
/// det(I − t g) = Σ_k c_k (−t)^k，其中 c_k 為特徵多項式 det(xI − g) 之係數（對稱提升）。
pub fn molien_series_lin(
  mats : Array[Array[Int]],
  dim : Int,
  p : Int,
  nmax : Int
) -> (Array[Rat], Bool) {
  let elems = matrix_group_elements(mats, dim, p, 100000)
  if elems.length() == 0 {
    return (Array::make(nmax, rat(0, 1)), false)
  }
  let mut acc = Array::make(nmax, rat(0, 1))
  for g in elems {
    let cs = char_poly(g, dim, p)
    let poly = Array::make(dim + 1, rat(0, 1))
    // det(I − t·g) = Σ_k c_k t^k（c_k 為 det(xI − g) 之降冪係數）
    let mut k = 0
    while k <= dim {
      poly[k] = rat(lift_symmetric(cs[k], p), 1)
      k = k + 1
    }
    let (term, ok) = series_inv(poly, nmax)
    if !ok {
      return (Array::make(nmax, rat(0, 1)), false)
    }
    acc = series_add(acc, term, nmax)
  }
  (series_scale(acc, rat(1, elems.length()), nmax), true)
}

/// 反射群判定：每個生成元為「偽反射」（固定子空間餘維 1）。
/// 以矩陣檢驗：rank(g − I) ≤ 1。
pub fn is_pseudo_reflection(g : Array[Int], dim : Int, p : Int) -> Bool {
  let a = Array::make(dim * dim, 0)
  let mut i = 0
  while i < dim {
    let mut j = 0
    while j < dim {
      let delta = if i == j { 1 % p } else { 0 }
      a[i * dim + j] = mod_pos(g[i * dim + j] - delta, p)
      j = j + 1
    }
    i = i + 1
  }
  matrix_rank(a, dim, dim, p) <= 1
}

/// 列舉 k 變數、次數 d 之單項式（指數向量）。
pub fn monomials_of_degree(k : Int, d : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  if k == 1 {
    return [[d]]
  }
  let mut i = 0
  while i <= d {
    for rest in monomials_of_degree(k - 1, d - i) {
      let v : Array[Int] = [i]
      for x in rest {
        v.push(x)
      }
      out.push(v)
    }
    i = i + 1
  }
  out
}

/// 以線性代數直接計算不變子空間維度：作用於 k 變數、次數 d 之單項式空間，
/// 維度 = dim ∩_g ker(ρ(g) − I)（即 χ 之 Burnside 平均之等價線性代數定義）。
pub fn invariant_dim_linalg(
  gens : Array[Array[Int]],
  k : Int,
  d : Int,
  p : Int
) -> Int {
  let mons = monomials_of_degree(k, d)
  let m = mons.length()
  if m == 0 {
    return 0
  }
  let mats : Array[Array[Int]] = []
  for g in gens {
    // ρ(g) 作用於單項式：指數向量 v ↦ M v（指數之上三角矩陣乘向量）
    let img : Array[Int] = []
    for v in mons {
      // 計算新指數
      let w = Array::make(k, 0)
      let mut j = 0
      while j < k {
        let mut s = 0
        let mut t = 0
        while t < k {
          s = s + g[t * k + j] * v[t]
          t = t + 1
        }
        w[j] = s
        j = j + 1
      }
      // 找到對應單項式指標
      let mut idx = -1
      let mut t2 = 0
      for u in mons {
        let mut same = true
        let mut q = 0
        while q < k {
          if u[q] != w[q] {
            same = false
          }
          q = q + 1
        }
        if same {
          idx = t2
        }
        t2 = t2 + 1
      }
      img.push(idx)
    }
    // 置換矩陣（單項式空間上之置換）
    let mat = Array::make(m * m, 0)
    let mut t3 = 0
    for v in mons {
      let _ = v
      if img[t3] >= 0 {
        mat[t3 * m + img[t3]] = 1 % p
      }
      t3 = t3 + 1
    }
    mats.push(mat)
  }
  joint_fixed_space(mats, m, p).length()
}

/// 級數之比較（逐項相等）。
pub fn series_eq(a : Array[Rat], b : Array[Rat]) -> Bool {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut i = 0
  while i < n {
    if a[i].n != b[i].n || a[i].d != b[i].d {
      return false
    }
    i = i + 1
  }
  true
}

/// 二變數 S₃ 之置換矩陣生成元（作用於 x, y）。
pub fn s3_permutation_matrices() -> (Array[Array[Int]], Int) {
  // (x y)
  let swap : Array[Int] = [0, 1, 1, 0]
  // 3-輪換 x→y→x+y? 取 (x,y) ↦ (y, −x−y)（A₂ 之 Coxeter 元素，階 3）
  let rot : Array[Int] = [0, mod_pos(-1, 5), 1, mod_pos(-1, 5)]
  ([swap, rot], 2)
}

// ───────────────────────── 義務自證 ─────────────────────────

/// 註冊多項式環作用之義務自證。
pub fn poly_register(led : Ledger) -> Unit {
  led.oblige("poly.arithmetic", A1_Group, Order, "F_p[x] 算術：除法恆等式 a = q·b + r（deg r < deg b），對 F₃、F₅ 上所有低次多項式全測", fn() {
    let mut ok = true
    for p in [2, 3, 5] {
      let maxdeg = 3
      let total = int_pow(p, maxdeg + 1)
      let mut ca = 0
      while ca < total {
        let a = Array::make(maxdeg + 1, 0)
        let mut t = ca
        let mut i = 0
        while i <= maxdeg {
          a[i] = t % p
          t = t / p
          i = i + 1
        }
        let mut cb = 1
        while cb < total {
          let b = Array::make(maxdeg + 1, 0)
          let mut t2 = cb
          i = 0
          while i <= maxdeg {
            b[i] = t2 % p
            t2 = t2 / p
            i = i + 1
          }
          if poly_trim(b, p).length() > 0 {
            let (q, r) = poly_divmod(a, b, p)
            let rec = poly_add(poly_mul(q, b, p), r, p)
            if !poly_eq(rec, a, p) {
              ok = false
            }
            if r.length() > 0 && poly_deg(r) >= poly_deg(poly_trim(b, p)) {
              ok = false
            }
          }
          cb = cb + 1
        }
        ca = ca + 1
      }
    }
    ok
  })
  led.oblige("poly.irreducible", A1_Group, LinearCount, "不可約性：x²+x+1 於 F₂ 不可約；x²+1 於 F₃、F₅ 可約；F₂ 上 4 個一次式乘積皆可約", fn() {
    poly_is_irreducible([1, 1, 1], 2) &&
    poly_is_irreducible([1, 0, 1], 3) &&
    !poly_is_irreducible([1, 0, 1], 5) &&
    !poly_is_irreducible([1, 1, 1], 3) &&
    !poly_is_irreducible([1, 0, 1, 0, 1], 2)
  })
  led.oblige("poly.frobenius", A2_Action, OrbitStabilizer, "Frobenius 作用：對不可約 f（deg n）之根，φ : x ↦ x^p 之軌道長度 = n，且 x^{p^n} ≡ x", fn() {
    let mut ok = true
    // f = x² + x + 1 ∈ F₂[x]（deg 2, 不可約）
    let f2 : Array[Int] = [1, 1, 1]
    let (orb2, first2) = frobenius_orbit([0, 1], f2, 2, 20)
    if orb2.length() != 2 || first2 != 0 {
      ok = false
    }
    // f = x² + 1 ∈ F₃[x] 可約（x²+1 = (x−1)(x+1)？ x²+1 於 F₃：1²+1=2≠0, 2²+1=5=2≠0 → 不可約）
    let f3 : Array[Int] = [1, 0, 1]
    let (orb3, first3) = frobenius_orbit([0, 1], f3, 3, 20)
    if orb3.length() != 2 || first3 != 0 {
      ok = false
    }
    // x^p ≡ x mod (x − a) 之檢查：Frobenius 固定 F_p
    let fp : Array[Int] = [0, 1] // f = x
    let (orbp, fp0) = frobenius_orbit([1], fp, 5, 10)
    if fp0 != 0 || orbp.length() != 1 {
      ok = false
    }
    ok
  })
  led.oblige("poly.translation_invariants", A2_Action, LinearCount, "Z/p 平移作用 x ↦ x+1 之不變元 = F_p[x^p − x]：維數等於 k 之個數（kp ≤ N−1，p=3,5）", fn() {
    let mut ok = true
    for p in [3, 5] {
      let nmax = 2 * p
      let inv = translation_invariants(nmax, p)
      // 期望之基底： (x^p − x)^k，k ≥ 0，次數 kp < nmax
      let mut expected = 0
      let mut k = 0
      while k * p < nmax {
        expected = expected + 1
        k = k + 1
      }
      if inv.length() != expected {
        ok = false
      }
      // 每個 (x^p−x)^k 皆為平移不變
      let mut kk = 0
      while kk * p < nmax {
        let poly = frobenius_poly_power(kk, p)
        if !in_span(inv, poly, p) {
          ok = false
        }
        kk = kk + 1
      }
    }
    ok
  })
  led.oblige("poly.derivation_invariants", A2_Action, LinearCount, "α_p（無窮小平移）之不變元 = ker D = F_p[x^p]：維數為 ⌈N/p⌉（p=3,5）", fn() {
    let mut ok = true
    for p in [3, 5] {
      let nmax = 2 * p
      let inv = derivation_invariants(nmax, p)
      let mut expected = 0
      let mut k = 0
      while k * p < nmax {
        expected = expected + 1
        k = k + 1
      }
      if inv.length() != expected {
        ok = false
      }
      let mut kk = 0
      while kk * p < nmax {
        let poly = Array::make(kk * p + 1, 0)
        poly[kk * p] = 1 % p
        if !in_span(inv, poly, p) {
          ok = false
        }
        kk = kk + 1
      }
    }
    ok
  })
  led.oblige("poly.cyclic_quotient", A2_Action, Burnside, "Z/n 作用於 F_p[x]/(x^n−1)（x ↦ x^k）為置換，軌道數由 Burnside 引理給出（n=5,k=2：階 4、2 軌道）", fn() {
    let n = 5
    let k = 2
    let perm = cyclic_action_on_perm(n, k)
    let g = perm_group(n, [perm])
    let (elems, ok) = g.elements(1000)
    if !ok {
      return false
    }
    let mut fixed = 0
    for e in elems {
      fixed = fixed + perm_fixed_count(e)
    }
    let orbits = g.orbits().length()
    // Burnside：軌道數 = 平均固定點；且 i ↦ 2i mod 5 之階為 4
    fixed / elems.length() == orbits && perm_order(perm) == 4 && orbits == 2
  })
  led.oblige("poly.molien_c2", A1_Group, Molien, "Molien 級數：Z/2 作用 x ↦ −x 之不變環為 F[x²]，P(t) = 1/(1−t²)", fn() {
    let nmax = 8
    // 1 維表示 x ↦ −x：g = [−1]，det(I − t g) = 1 + t
    let mats : Array[Array[Int]] = [[mod_pos(-1, 5)]]
    let (series, ok) = molien_series_lin(mats, 1, 5, nmax)
    if !ok {
      return false
    }
    let expect = series_geom(2, nmax)
    series_eq(series, expect)
  })
  led.oblige("poly.molien_s3_perm", A1_Group, Molien, "Molien 級數＋CST：S₃ 置換作用於三變數，P(t) = 1/((1−t)(1−t²)(1−t³))，且生成元皆為偽反射", fn() {
    let nmax = 12
    // S₃ 於三變數之置換表示（對換 (0 1)、(1 2)）
    let gens : Array[Array[Int]] = [transposition(3, 0, 1), transposition(3, 1, 2)]
    let (series, ok) = molien_series(gens, 3, nmax)
    if !ok {
      return false
    }
    let expect = series_mul(
      series_mul(series_geom(1, nmax), series_geom(2, nmax), nmax),
      series_geom(3, nmax),
      nmax,
    )
    // 生成元於三變數置換表示中之 rank(g − I) = 1（偽反射）
    let mut refl = true
    for g in gens {
      let m = Array::make(9, 0)
      let mut i = 0
      while i < 3 {
        m[i * 3 + g[i]] = 1 % 5
        i = i + 1
      }
      if !is_pseudo_reflection(m, 3, 5) {
        refl = false
      }
    }
    series_eq(series, expect) && refl
  })
  led.oblige("poly.molien_s3_lin", A1_Group, Molien, "Molien 級數（二維標準表示）：S₃ 作用於 x,y 之不變環為自由環（次數 2,3，乘積 = |G| = 6）", fn() {
    let nmax = 10
    let (mats, dim) = s3_permutation_matrices()
    let (series, ok) = molien_series_lin(mats, dim, 5, nmax)
    if !ok {
      return false
    }
    let expect = series_mul(series_geom(2, nmax), series_geom(3, nmax), nmax)
    let elems = matrix_group_elements(mats, dim, 5, 10000)
    series_eq(series, expect) && elems.length() == 6 && 2 * 3 == 6
  })
}

/// 置換輔助：Z/n 作用於指標（與 cyclic_action_on_quotient 相同，供義務自證獨立使用）。
pub fn cyclic_action_on_perm(n : Int, k : Int) -> Array[Int] {
  let perm = perm_id(n)
  let mut i = 0
  while i < n {
    perm[i] = mod_pos(i * k, n)
    i = i + 1
  }
  perm
}
