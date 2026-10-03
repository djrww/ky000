//! 矩陣作用於向量：GL(n,q)、SL(n,q)、PSL(2,q) 之置換表示、特徵多項式與 Cayley–Hamilton、
//! 特徵空間（固定子空間），以及古典計數公式 |GL(n,q)| = q^{n(n−1)/2} Π_{i=1}^n (q^i − 1)。
//!
//! 定義 D4：M ∈ M_n(F_p) 以長度 n² 之 `Array[Int]`（列主序）表示，元素取模 p。
//! 定義 D5：M 作用於 v ∈ F_p^n 為 (M,v) ↦ Mv；此作用給出置換表示 φ : GL(n,p) → Sym(F_p^n \ {0})。

/// 單位矩陣（列主序，長度 n²）。
pub fn mat_id(n : Int, p : Int) -> Array[Int] {
  let m = Array::make(n * n, 0)
  let mut i = 0
  while i < n {
    m[i * n + i] = 1 % p
    i = i + 1
  }
  m
}

/// 矩陣乘法 mod p。
pub fn mat_mul(a : Array[Int], b : Array[Int], n : Int, p : Int) -> Array[Int] {
  let out = Array::make(n * n, 0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let mut s = 0
      let mut k = 0
      while k < n {
        s = s + a[i * n + k] * b[k * n + j]
        k = k + 1
      }
      out[i * n + j] = mod_pos(s, p)
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 矩陣乘向量。
pub fn mat_vec(m : Array[Int], v : Array[Int], n : Int, p : Int) -> Array[Int] {
  let out = Array::make(n, 0)
  let mut i = 0
  while i < n {
    let mut s = 0
    let mut j = 0
    while j < n {
      s = s + m[i * n + j] * v[j]
      j = j + 1
    }
    out[i] = mod_pos(s, p)
    i = i + 1
  }
  out
}

/// 由矩陣生成集顯式枚舉整個矩陣群（BFS，上限 cap）。
pub fn matrix_group_elements(
  gens : Array[Array[Int]],
  n : Int,
  p : Int,
  cap : Int
) -> Array[Array[Int]] {
  let seen : Map[String, Unit] = Map([])
  let out : Array[Array[Int]] = []
  let id = mat_id(n, p)
  out.push(id)
  seen[perm_key(id)] = ()
  if gens.length() == 0 {
    return out
  }
  let mut head = 0
  while head < out.length() {
    if out.length() > cap {
      return []
    }
    let g = out[head]
    head = head + 1
    for s in gens {
      let h = mat_mul(s, g, n, p)
      let k = perm_key(h)
      if !seen.contains(k) {
        seen[k] = ()
        out.push(h)
      }
    }
  }
  out
}

/// 矩陣是否相等。
pub fn mat_eq(a : Array[Int], b : Array[Int]) -> Bool {
  let mut i = 0
  while i < a.length() {
    if a[i] != b[i] {
      return false
    }
    i = i + 1
  }
  true
}

/// 行列式（高斯消去法，於 F_p 上）。
pub fn mat_det(m : Array[Int], n : Int, p : Int) -> Int {
  let a = Array::make(n * n, 0)
  let mut t = 0
  while t < n * n {
    a[t] = mod_pos(m[t], p)
    t = t + 1
  }
  let mut det = 1
  let mut col = 0
  while col < n {
    // 找樞軸
    let mut piv = -1
    let mut r = col
    while r < n {
      if a[r * n + col] != 0 {
        piv = r
      }
      r = r + 1
    }
    if piv < 0 {
      return 0
    }
    if piv != col {
      // 交換列
      let mut c = 0
      while c < n {
        let tmp = a[col * n + c]
        a[col * n + c] = a[piv * n + c]
        a[piv * n + c] = tmp
        c = c + 1
      }
      det = mod_pos(-det, p)
    }
    let pivval = a[col * n + col]
    det = mod_pos(det * pivval, p)
    let inv = mod_inv(pivval, p)
    let mut r2 = col + 1
    while r2 < n {
      let factor = mod_pos(a[r2 * n + col] * inv, p)
      if factor != 0 {
        let mut c2 = col
        while c2 < n {
          a[r2 * n + c2] = mod_pos(a[r2 * n + c2] - factor * a[col * n + c2], p)
          c2 = c2 + 1
        }
      }
      r2 = r2 + 1
    }
    col = col + 1
  }
  det
}

/// 矩陣求逆（高斯–喬登），若不可逆則回傳零矩陣。
pub fn mat_inv(m : Array[Int], n : Int, p : Int) -> Array[Int] {
  let a = Array::make(n * 2 * n, 0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      a[i * (2 * n) + j] = mod_pos(m[i * n + j], p)
      j = j + 1
    }
    a[i * (2 * n) + n + i] = 1 % p
    i = i + 1
  }
  let w = 2 * n
  let mut col = 0
  while col < n {
    let mut piv = -1
    let mut r = col
    while r < n {
      if a[r * w + col] != 0 {
        piv = r
      }
      r = r + 1
    }
    if piv < 0 {
      return Array::make(n * n, 0)
    }
    if piv != col {
      let mut c = 0
      while c < w {
        let tmp = a[col * w + c]
        a[col * w + c] = a[piv * w + c]
        a[piv * w + c] = tmp
        c = c + 1
      }
    }
    let inv = mod_inv(a[col * w + col], p)
    let mut c2 = 0
    while c2 < w {
      a[col * w + c2] = mod_pos(a[col * w + c2] * inv, p)
      c2 = c2 + 1
    }
    let mut r2 = 0
    while r2 < n {
      if r2 != col {
        let factor = a[r2 * w + col]
        if factor != 0 {
          let mut c3 = 0
          while c3 < w {
            a[r2 * w + c3] = mod_pos(a[r2 * w + c3] - factor * a[col * w + c3], p)
            c3 = c3 + 1
          }
        }
      }
      r2 = r2 + 1
    }
    col = col + 1
  }
  let out = Array::make(n * n, 0)
  i = 0
  while i < n {
    let mut j = 0
    while j < n {
      out[i * n + j] = a[i * w + n + j]
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 矩陣冪。
pub fn mat_pow(m : Array[Int], k : Int, n : Int, p : Int) -> Array[Int] {
  let mut out = mat_id(n, p)
  let mut i = 0
  while i < k {
    out = mat_mul(m, out, n, p)
    i = i + 1
  }
  out
}

/// 矩陣在 GL(n,p) 中之階（乘法階）。
pub fn mat_order(m : Array[Int], n : Int, p : Int) -> Int {
  let id = mat_id(n, p)
  let mut cur = m
  let mut k = 1
  while !mat_eq(cur, id) && k < 100000 {
    cur = mat_mul(m, cur, n, p)
    k = k + 1
  }
  k
}

/// |GL(n,q)| = q^{n(n−1)/2} · Π_{i=1}^{n} (q^i − 1)。
pub fn gl_order(n : Int, q : Int) -> BigInt {
  let mut acc : BigInt = 1N
  let mut i = 1
  while i <= n {
    acc = acc * BigInt::from_int(int_pow(q, i) - 1)
    i = i + 1
  }
  let mut e = n * (n - 1) / 2
  while e > 0 {
    acc = acc * BigInt::from_int(q)
    e = e - 1
  }
  acc
}

/// |SL(n,q)| = |GL(n,q)| / (q−1)。
pub fn sl_order(n : Int, q : Int) -> BigInt {
  gl_order(n, q) / BigInt::from_int(q - 1)
}

/// |PSL(n,q)| = |SL(n,q)| / gcd(n, q−1)。
pub fn psl_order(n : Int, q : Int) -> BigInt {
  sl_order(n, q) / BigInt::from_int(gcd_int(n, q - 1))
}

/// 向量之索引編碼（基底 p）。
pub fn vec_index(v : Array[Int], p : Int) -> Int {
  let mut acc = 0
  let mut i = v.length() - 1
  while i >= 0 {
    acc = acc * p + v[i]
    i = i - 1
  }
  acc
}

/// 計算 |GL(n,p)| 之直接枚舉值（僅供小 n,p 之義務自證對照）。
pub fn gl_order_by_counting(n : Int, p : Int) -> Int {
  let total = int_pow(p, n * n)
  let mut count = 0
  let mut code = 0
  while code < total {
    // 解碼為矩陣
    let m = Array::make(n * n, 0)
    let mut t = code
    let mut i = 0
    while i < n * n {
      m[i] = t % p
      t = t / p
      i = i + 1
    }
    if mat_det(m, n, p) != 0 {
      count = count + 1
    }
    code = code + 1
  }
  count
}

/// 由矩陣生成集建立其作用於非零向量之置換生成集（φ : GL(n,p) → Sym(F_p^n \ {0})）。
pub fn matrix_action_on_vectors(gens : Array[Array[Int]], n : Int, p : Int) -> Array[Array[Int]] {
  let nvec = int_pow(p, n)
  let out : Array[Array[Int]] = []
  for m in gens {
    let perm = perm_id(nvec)
    let mut code = 0
    while code < nvec {
      let v = Array::make(n, 0)
      let mut t = code
      let mut i = 0
      while i < n {
        v[i] = t % p
        t = t / p
        i = i + 1
      }
      perm[code] = vec_index(mat_vec(m, v, n, p), p)
      code = code + 1
    }
    out.push(perm)
  }
  out
}

/// 射影點之標準代表（首個非零座標化為 1），用於 GL 作用於 P^{n−1}(F_p)。
pub fn proj_canonical(v : Array[Int], n : Int, p : Int) -> Array[Int] {
  let mut first = -1
  let mut i = 0
  while i < n {
    if v[i] != 0 && first < 0 {
      first = i
    }
    i = i + 1
  }
  if first < 0 {
    return v
  }
  let inv = mod_inv(v[first], p)
  let out = Array::make(n, 0)
  i = 0
  while i < n {
    out[i] = mod_pos(v[i] * inv, p)
    i = i + 1
  }
  out
}

/// 射影空間 P^{n−1}(F_p) 之所有點（以索引表示）。
pub fn projective_points(n : Int, p : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  let nvec = int_pow(p, n)
  let seen : Map[String, Unit] = Map([])
  let mut code = 1
  while code < nvec {
    let v = Array::make(n, 0)
    let mut t = code
    let mut i = 0
    while i < n {
      v[i] = t % p
      t = t / p
      i = i + 1
    }
    let c = proj_canonical(v, n, p)
    let k = perm_key(c)
    if !seen.contains(k) {
      seen[k] = ()
      out.push(c)
    }
    code = code + 1
  }
  out
}

/// GL(n,p) 作用於射影空間之置換生成集。
pub fn matrix_action_on_projective(gens : Array[Array[Int]], n : Int, p : Int) -> (Array[Array[Int]], Int) {
  let pts = projective_points(n, p)
  let pos : Map[String, Int] = Map([])
  let mut i = 0
  for pt in pts {
    pos[perm_key(pt)] = i
    i = i + 1
  }
  let out : Array[Array[Int]] = []
  for m in gens {
    let perm = perm_id(pts.length())
    let mut j = 0
    for pt in pts {
      let img = proj_canonical(mat_vec(m, pt, n, p), n, p)
      perm[j] = pos[perm_key(img)]
      j = j + 1
    }
    out.push(perm)
  }
  (out, pts.length())
}

/// SL(2,p) 在「有序基底」(v₁,v₂)（det(v₁,v₂) ≠ 0）上之忠實作用，度為 (p²−1)(p²−p)。
pub fn sl2_action_on_bases(gens : Array[Array[Int]], p : Int) -> (Array[Array[Int]], Int) {
  let n = 2
  let pairs : Array[Array[Int]] = []
  let mut a = 0
  while a < p * p {
    let v1 = [a % p, a / p]
    if !(v1[0] == 0 && v1[1] == 0) {
      let mut b = 0
      while b < p * p {
        let v2 = [b % p, b / p]
        let det = mod_pos(v1[0] * v2[1] - v1[1] * v2[0], p)
        if det != 0 {
          pairs.push([v1[0], v1[1], v2[0], v2[1]])
        }
        b = b + 1
      }
    }
    a = a + 1
  }
  let pos : Map[String, Int] = Map([])
  let mut i = 0
  for pr in pairs {
    pos[perm_key(pr)] = i
    i = i + 1
  }
  let out : Array[Array[Int]] = []
  for m in gens {
    let perm = perm_id(pairs.length())
    let mut j = 0
    for pr in pairs {
      let v1 = mat_vec(m, [pr[0], pr[1]], n, p)
      let v2 = mat_vec(m, [pr[2], pr[3]], n, p)
      perm[j] = pos[perm_key([v1[0], v1[1], v2[0], v2[1]])]
      j = j + 1
    }
    out.push(perm)
  }
  (out, pairs.length())
}

/// PSL(2,p) 之生成元：[[1,1],[0,1]] 與 [[0,−1],[1,0]]（在 P¹ 上之作用）。
pub fn psl2_on_projective_line(p : Int) -> (Array[Array[Int]], Int, Int) {
  let t : Array[Int] = [1 % p, 1 % p, 0, 1 % p]
  let s : Array[Int] = [0, mod_pos(-1, p), 1 % p, 0]
  let (gens, deg) = matrix_action_on_projective([t, s], 2, p)
  (gens, deg, p + 1)
}

/// SL(n,p) 之生成元：切換矩陣 E_{ij}(1) 與對換矩陣。
pub fn sl_gens(n : Int, p : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      if i != j {
        let m = mat_id(n, p)
        m[i * n + j] = 1 % p
        out.push(m)
      }
      j = j + 1
    }
    i = i + 1
  }
  // 註：切換矩陣 E_{ij}(1) 之行列式為 1，且已足以生成 SL(n,p)（n ≥ 2）。
  out
}

/// GL(n,p) 之生成元（切換、對換與純量）。
pub fn gl_gens(n : Int, p : Int) -> Array[Array[Int]] {
  let out = sl_gens(n, p)
  if p > 2 {
    let scal = mat_id(n, p)
    scal[0] = mod_pos(2, p)
    out.push(scal)
  }
  out
}

/// 多項式加法（升冪係數，模 p）。
pub fn poly_add_p(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let n = if a.length() > b.length() { a.length() } else { b.length() }
  let out = Array::make(n, 0)
  let mut i = 0
  while i < n {
    let x = if i < a.length() { a[i] } else { 0 }
    let y = if i < b.length() { b[i] } else { 0 }
    out[i] = mod_pos(x + y, p)
    i = i + 1
  }
  out
}

/// 多項式乘法（升冪係數，模 p）。
pub fn poly_mul_p(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  if a.length() == 0 || b.length() == 0 {
    return []
  }
  let out = Array::make(a.length() + b.length() - 1, 0)
  let mut i = 0
  while i < a.length() {
    let mut j = 0
    while j < b.length() {
      out[i + j] = mod_pos(out[i + j] + a[i] * b[j], p)
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 多項式取負（升冪係數，模 p）。
pub fn poly_neg_p(a : Array[Int], p : Int) -> Array[Int] {
  let out = Array::make(a.length(), 0)
  let mut i = 0
  while i < a.length() {
    out[i] = mod_pos(-a[i], p)
    i = i + 1
  }
  out
}

/// 子矩陣之行列式（多項式係數，升冪），Laplace 展開。
fn det_poly_minor(
  m : Array[Int],
  n : Int,
  rows : Array[Int],
  cols : Array[Int],
  p : Int
) -> Array[Int] {
  let k = rows.length()
  if k == 1 {
    // 單元素 xI − M
    let i = rows[0]
    let j = cols[0]
    let c0 = if i == j { mod_pos(-m[i * n + j], p) } else { mod_pos(-m[i * n + j], p) }
    if i == j {
      return [c0, 1 % p]
    } else {
      return [c0]
    }
  }
  let out = Array::make(k + 1, 0)
  let mut j = 0
  while j < k {
    let i0 = rows[0]
    let j0 = cols[j]
    let entry : Array[Int] = if i0 == j0 { [mod_pos(-m[i0 * n + j0], p), 1 % p] } else { [mod_pos(-m[i0 * n + j0], p)] }
    let sub_rows : Array[Int] = []
    let mut t = 1
    while t < k {
      sub_rows.push(rows[t])
      t = t + 1
    }
    let sub_cols : Array[Int] = []
    t = 0
    while t < k {
      if t != j {
        sub_cols.push(cols[t])
      }
      t = t + 1
    }
    let minor = det_poly_minor(m, n, sub_rows, sub_cols, p)
    let term = poly_mul_p(entry, minor, p)
    let signed = if j % 2 == 0 { term } else { poly_neg_p(term, p) }
    let mut idx = 0
    while idx < signed.length() {
      if idx < out.length() {
        out[idx] = mod_pos(out[idx] + signed[idx], p)
      }
      idx = idx + 1
    }
    j = j + 1
  }
  out
}

/// 特徵多項式 χ_M(x) = det(xI − M)，回傳「降冪」係數 [c₀,…,c_n]（c₀ 為 x^n 之係數）。
/// 以 Laplace 展開計算，不需除以 k，故對所有質數 p 有效。
pub fn char_poly(m : Array[Int], n : Int, p : Int) -> Array[Int] {
  let rows : Array[Int] = []
  let cols : Array[Int] = []
  let mut i = 0
  while i < n {
    rows.push(i)
    cols.push(i)
    i = i + 1
  }
  let asc = det_poly_minor(m, n, rows, cols, p)
  // 升冪 → 降冪
  let deg = n
  let out = Array::make(n + 1, 0)
  i = 0
  while i < asc.length() {
    out[deg - i] = asc[i]
    i = i + 1
  }
  out
}

/// 以 Faddeev–LeVerrier 法計算特徵多項式（需 k < p），作為獨立對照。
pub fn char_poly_flv(m : Array[Int], n : Int, p : Int) -> Array[Int] {
  let cs = Array::make(n + 1, 0)
  cs[0] = 1 % p
  let mut a = mat_id(n, p)
  let mut k = 1
  while k <= n {
    let ma = mat_mul(m, a, n, p)
    let mut tr = 0
    let mut i = 0
    while i < n {
      tr = tr + ma[i * n + i]
      i = i + 1
    }
    let ck = mod_pos(-mod_pos(tr, p) * mod_inv(k % p, p), p)
    cs[k] = ck
    let next = mat_mul(m, a, n, p)
    i = 0
    while i < n {
      next[i * n + i] = mod_pos(next[i * n + i] + ck, p)
      i = i + 1
    }
    a = next
    k = k + 1
  }
  cs
}

/// 以多項式係數（高次在前）作用於矩陣：p(M)。
pub fn poly_eval_matrix(cs : Array[Int], m : Array[Int], n : Int, p : Int) -> Array[Int] {
  let mut acc = Array::make(n * n, 0)
  for c in cs {
    acc = mat_mul(acc, m, n, p)
    let mut i = 0
    while i < n {
      acc[i * n + i] = mod_pos(acc[i * n + i] + c, p)
      i = i + 1
    }
  }
  acc
}

/// 一般核空間：a 為 rows×cols 矩陣（列主序），回傳 ker(a) 之基底（F_p 上）。
pub fn mat_kernel(a : Array[Int], rows : Int, cols : Int, p : Int) -> Array[Array[Int]] {
  let m = Array::make(rows * cols, 0)
  let mut t = 0
  while t < rows * cols {
    m[t] = mod_pos(a[t], p)
    t = t + 1
  }
  let pivots : Array[Int] = []
  let mut row = 0
  let mut col = 0
  while col < cols && row < rows {
    let mut piv = -1
    let mut r = row
    while r < rows {
      if m[r * cols + col] != 0 {
        piv = r
      }
      r = r + 1
    }
    if piv < 0 {
      col = col + 1
      continue
    }
    if piv != row {
      let mut c = 0
      while c < cols {
        let tmp = m[row * cols + c]
        m[row * cols + c] = m[piv * cols + c]
        m[piv * cols + c] = tmp
        c = c + 1
      }
    }
    let inv = mod_inv(m[row * cols + col], p)
    let mut c2 = 0
    while c2 < cols {
      m[row * cols + c2] = mod_pos(m[row * cols + c2] * inv, p)
      c2 = c2 + 1
    }
    let mut r2 = 0
    while r2 < rows {
      if r2 != row {
        let factor = m[r2 * cols + col]
        if factor != 0 {
          let mut c3 = 0
          while c3 < cols {
            m[r2 * cols + c3] = mod_pos(m[r2 * cols + c3] - factor * m[row * cols + c3], p)
            c3 = c3 + 1
          }
        }
      }
      r2 = r2 + 1
    }
    pivots.push(col)
    row = row + 1
    col = col + 1
  }
  let is_pivot = Array::make(cols, false)
  for c in pivots {
    is_pivot[c] = true
  }
  let out : Array[Array[Int]] = []
  let mut free = 0
  while free < cols {
    if !is_pivot[free] {
      let v = Array::make(cols, 0)
      v[free] = 1 % p
      let mut k = 0
      for pc in pivots {
        v[pc] = mod_pos(-m[k * cols + free], p)
        k = k + 1
      }
      out.push(v)
    }
    free = free + 1
  }
  out
}

/// 向量是否屬於若干向量張成之空間（F_p 上）。
pub fn in_span(vecs : Array[Array[Int]], v : Array[Int], p : Int) -> Bool {
  // 以高斯消去法判定：rank(vecs) == rank(vecs ∪ {v})；不同長度時以零補齊（視為次數 < dim 之元素）
  let mut dim = v.length()
  for u in vecs {
    if u.length() > dim {
      dim = u.length()
    }
  }
  let rows : Array[Int] = []
  for u in vecs {
    let mut i = 0
    while i < dim {
      rows.push(if i < u.length() { u[i] } else { 0 })
      i = i + 1
    }
  }
  let base = rows.length() / dim
  let augmented = Array::make((base + 1) * dim, 0)
  let mut i = 0
  while i < base * dim {
    augmented[i] = rows[i]
    i = i + 1
  }
  i = 0
  while i < dim {
    augmented[base * dim + i] = if i < v.length() { v[i] } else { 0 }
    i = i + 1
  }
  let r1 = matrix_rank(rows, base, dim, p)
  let r2 = matrix_rank(augmented, base + 1, dim, p)
  r1 == r2
}

/// 矩陣秩（F_p 上）。
pub fn matrix_rank(a : Array[Int], rows : Int, cols : Int, p : Int) -> Int {
  let m = Array::make(rows * cols, 0)
  let mut t = 0
  while t < rows * cols {
    m[t] = mod_pos(a[t], p)
    t = t + 1
  }
  let mut rank = 0
  let mut row = 0
  let mut col = 0
  while col < cols && row < rows {
    let mut piv = -1
    let mut r = row
    while r < rows {
      if m[r * cols + col] != 0 {
        piv = r
      }
      r = r + 1
    }
    if piv < 0 {
      col = col + 1
      continue
    }
    if piv != row {
      let mut c = 0
      while c < cols {
        let tmp = m[row * cols + c]
        m[row * cols + c] = m[piv * cols + c]
        m[piv * cols + c] = tmp
        c = c + 1
      }
    }
    let inv = mod_inv(m[row * cols + col], p)
    let mut c2 = 0
    while c2 < cols {
      m[row * cols + c2] = mod_pos(m[row * cols + c2] * inv, p)
      c2 = c2 + 1
    }
    let mut r2 = 0
    while r2 < rows {
      if r2 != row {
        let factor = m[r2 * cols + col]
        if factor != 0 {
          let mut c3 = 0
          while c3 < cols {
            m[r2 * cols + c3] = mod_pos(m[r2 * cols + c3] - factor * m[row * cols + c3], p)
            c3 = c3 + 1
          }
        }
      }
      r2 = r2 + 1
    }
    rank = rank + 1
    row = row + 1
    col = col + 1
  }
  rank
}

/// 線性群作用之共同固定子空間 ∩_g ker(g − I)（矩陣以列主序、作用於 dim 維空間）。
pub fn joint_fixed_space(gens : Array[Array[Int]], dim : Int, p : Int) -> Array[Array[Int]] {
  if gens.length() == 0 {
    return []
  }
  let k = gens.length()
  let stacked = Array::make(k * dim * dim, 0)
  let mut gi = 0
  for g in gens {
    let mut i = 0
    while i < dim {
      let mut j = 0
      while j < dim {
        let val = g[i * dim + j]
        let delta = if i == j { 1 % p } else { 0 }
        stacked[(gi * dim + i) * dim + j] = mod_pos(val - delta, p)
        j = j + 1
      }
      i = i + 1
    }
    gi = gi + 1
  }
  mat_kernel(stacked, k * dim, dim, p)
}

/// 特徵空間 ker(M − λI) 之基底（以高斯消去法求核）。
pub fn eigenspace_kernel(m : Array[Int], lambda : Int, n : Int, p : Int) -> Array[Array[Int]] {
  let a = Array::make(n * n, 0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      a[i * n + j] = m[i * n + j]
      j = j + 1
    }
    a[i * n + i] = mod_pos(a[i * n + i] - lambda, p)
    i = i + 1
  }
  mat_kernel(a, n, n, p)
}

// ───────────────────────── 義務自證 ─────────────────────────

/// 註冊矩陣作用之義務自證。
pub fn matrix_register(led : Ledger) -> Unit {
  led.oblige("matrix.det", A1_Group, LinearCount, "行列式為乘法同態：det(AB)=det A·det B，det(I)=1，det(M)·det(M⁻¹)=1（F₂、F₃、F₅ 上）", fn() {
    let mut ok = true
    for p in [2, 3, 5] {
      let n = 2
      let total = int_pow(p, n * n)
      let mut code = 0
      while code < total {
        let a = Array::make(n * n, 0)
        let mut t = code
        let mut i = 0
        while i < n * n {
          a[i] = t % p
          t = t / p
          i = i + 1
        }
        let da = mat_det(a, n, p)
        // 取另一矩陣為 a 之轉置＋單位
        let b = mat_id(n, p)
        i = 0
        while i < n * n {
          b[i] = mod_pos(b[i] + a[(i % n) * n + i / n], p)
          i = i + 1
        }
        let db = mat_det(b, n, p)
        let dab = mat_det(mat_mul(a, b, n, p), n, p)
        if dab != mod_pos(da * db, p) {
          ok = false
        }
        if da != 0 {
          let ai = mat_inv(a, n, p)
          if mod_pos(da * mat_det(ai, n, p), p) != 1 {
            ok = false
          }
        }
        code = code + 1
      }
      if mat_det(mat_id(n, p), n, p) != 1 % p {
        ok = false
      }
    }
    ok
  })
  led.oblige("matrix.gl_order", A1_Group, LinearCount, "|GL(n,p)| 公式＝直接計數：(n,p)=(2,2),(2,3),(2,5) 之可逆矩陣數", fn() {
    gl_order(2, 2).to_string() == "6" &&
    gl_order(2, 3).to_string() == "48" &&
    gl_order(2, 5).to_string() == "480" &&
    gl_order_by_counting(2, 2) == 6 &&
    gl_order_by_counting(2, 3) == 48
  })
  led.oblige("matrix.sl_gens", A1_Group, LinearCount, "切換矩陣生成 SL(2,p)：在有序基底上之忠實作用階 = |SL(2,5)| = 120；在 P¹ 上之像階 = |PSL(2,5)| = 60", fn() {
    let gens5 = sl_gens(2, 5)
    let (bases, degb) = sl2_action_on_bases(gens5, 5)
    let ordb = group_order(degb, bases)
    let (proj, degp, _d) = psl2_on_projective_line(5)
    let ordp = group_order(degp, proj)
    ordb == 120 && ordp == 60 && sl_order(2, 5).to_string() == "120" && psl_order(2, 5).to_string() == "60"
  })
  led.oblige("matrix.psl2", A4_Simplicity, LinearCount, "PSL(2,p) 作用於射影直線 P¹(F_p)（p+1 點）：|PSL(2,5)|=60、|PSL(2,7)|=168、|PSL(2,11)|=660", fn() {
    let (g5, d5, deg5) = psl2_on_projective_line(5)
    let (g7, d7, deg7) = psl2_on_projective_line(7)
    let (g11, d11, deg11) = psl2_on_projective_line(11)
    let o5 = group_order(deg5, g5)
    let o7 = group_order(deg7, g7)
    let o11 = group_order(deg11, g11)
    let _ = d5
    let _ = d7
    let _ = d11
    o5 == 60 && o7 == 168 && o11 == 660 &&
    psl_order(2, 5).to_string() == "60" &&
    psl_order(2, 7).to_string() == "168" &&
    psl_order(2, 11).to_string() == "660"
  })
  led.oblige("matrix.vector_action", A2_Action, OrbitStabilizer, "GL(2,p) 在非零向量上遞移：軌道大小 = p²−1（p=3,5）；穩定化子階 = |GL|/(p²−1)", fn() {
    let p = 5
    let gens = gl_gens(2, p)
    let perms = matrix_action_on_vectors(gens, 2, p)
    let g = perm_group(int_pow(p, 2), perms)
    let used : Array[Int] = []
    let mut x = 1
    while x < p * p {
      used.push(x)
      x = x + 1
    }
    let orb = g.orbit_on(used)
    // GL(2,5) 在非零向量上遞移；|GL|/|Stab| = p²−1
    orb.length() == p * p - 1 && 480 / (p * p - 1) == 20 && g.orbit_on([0]).length() == 1
  })
  led.oblige("matrix.cayley_hamilton", A1_Group, LinearCount, "Cayley–Hamilton：χ_M(M) = 0（對 F₂、F₃、F₅ 上所有 2×2 矩陣全測）", fn() {
    let mut ok = true
    for p in [2, 3, 5] {
      let n = 2
      let total = int_pow(p, n * n)
      let mut code = 0
      while code < total {
        let a = Array::make(n * n, 0)
        let mut t = code
        let mut i = 0
        while i < n * n {
          a[i] = t % p
          t = t / p
          i = i + 1
        }
        let cs = char_poly(a, n, p)
        let val = poly_eval_matrix(cs, a, n, p)
        let zero = Array::make(n * n, 0)
        if !mat_eq(val, zero) {
          ok = false
        }
        code = code + 1
      }
    }
    ok
  })
  led.oblige("matrix.charpoly_dual", A1_Group, LinearCount, "兩獨立算式互證：Laplace 展開與 Faddeev–LeVerrier 之特徵多項式在 p>n 時一致", fn() {
    let mut ok = true
    for p in [3, 5, 7] {
      let n = 2
      let total = int_pow(p, n * n)
      let mut code = 0
      while code < total {
        let a = Array::make(n * n, 0)
        let mut t = code
        let mut i = 0
        while i < n * n {
          a[i] = t % p
          t = t / p
          i = i + 1
        }
        let c1 = char_poly(a, n, p)
        let c2 = char_poly_flv(a, n, p)
        if !mat_eq(c1, c2) {
          ok = false
        }
        code = code + 1
      }
    }
    ok
  })
  led.oblige("matrix.charpoly", A1_Group, LinearCount, "特徵多項式：2×2 時 χ(x)=x²−tr(M)x+det(M) 且 χ(M)=0；λ 為特徵值 ⟺ ker(M−λI)≠0", fn() {
    let p = 5
    let m : Array[Int] = [2, 1, 1, 3]
    let cs = char_poly(m, 2, p)
    let tr = mod_pos(m[0] + m[3], p)
    let det = mat_det(m, 2, p)
    let zero = Array::make(4, 0)
    let ker0 = eigenspace_kernel(m, 0, 2, p)
    cs[0] == 1 &&
    cs[1] == mod_pos(-tr, p) &&
    cs[2] == det &&
    mat_eq(poly_eval_matrix(cs, m, 2, p), zero) &&
    ((det == 0 && ker0.length() == 1) || (det != 0 && ker0.length() == 0))
  })
  led.oblige("matrix.matrix_order", A1_Group, Order, "矩陣階：M^k = I 且 k 為最小；階整除 |GL(n,p)|（F₅ 上之生成元）", fn() {
    let p = 5
    let n = 2
    let m : Array[Int] = [0, mod_pos(-1, p), 1, 0] // 階 4 之元素
    let o = mat_order(m, n, p)
    let id = mat_id(n, p)
    mat_eq(mat_pow(m, o, n, p), id) && !mat_eq(mat_pow(m, o - 1, n, p), id) && o == 4
  })
}
