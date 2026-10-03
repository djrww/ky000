//! 層作用（sheaf action）：Čech 複形、餘調維數、等變層與 Lefschetz 不動點公式。
//!
//! 以有限單純（胞腔）複形之細胞鏈複形實現層之餘調：
//!   * C⁰ →(δ⁰) C¹ →(δ¹) C²：δ⁰ 為頂點差分（0-上鏈 → 1-上鏈）、δ¹ 為多邊形面之邊界
//!     （對每個面之循環頂點序取帶號邊之和）。
//!   * H^i = ker δ^i / im δ^{i−1}，維數由矩陣秩計算（mod p 線性代數）。
//!   * 層作用：群 G 置換細胞，ρ^i 為 C^i 上之置換矩陣，且 δ∘ρ = ρ∘δ（等變層）。
//!   * Lefschetz 不動點公式（Hopf 跡公式）：Σ(−1)^i tr(ρ^i) = Σ(−1)^i dim Fix(H^i)，
//!     當 p ∤ |G| 時由不動子空間與餘調之分解計算。
//!
//! 三類標準例：圓周（三角形神經／3 胞腔）、球面（四面體邊界）、環面（3×3 方格胞腔複形）。

/// 一般（長方形）矩陣乘法：a 為 rows×mid，b 為 mid×cols（列主序）。
pub fn mat_mul_rect(
  a : Array[Int],
  b : Array[Int],
  rows : Int,
  mid : Int,
  cols : Int,
  p : Int
) -> Array[Int] {
  let out = Array::make(rows * cols, 0)
  let mut i = 0
  while i < rows {
    let mut j = 0
    while j < cols {
      let mut s = 0
      let mut k = 0
      while k < mid {
        s = s + a[i * mid + k] * b[k * cols + j]
        k = k + 1
      }
      out[i * cols + j] = mod_pos(s, p)
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// 尋找無向邊 {u,v} 於邊表之索引（-1 表不存在）。
pub fn edge_index(edges : Array[Array[Int]], u : Int, v : Int) -> Int {
  let mut i = 0
  while i < edges.length() {
    let e = edges[i]
    if (e[0] == u && e[1] == v) || (e[0] == v && e[1] == u) {
      return i
    }
    i = i + 1
  }
  -1
}

/// δ⁰：C⁰ → C¹（頂點差分；邊之定向取 u < v）。
pub fn cech_d0(nverts : Int, edges : Array[Array[Int]], p : Int) -> Array[Int] {
  let rows = edges.length()
  let m = Array::make(rows * nverts, 0)
  let mut r = 0
  while r < rows {
    let e = edges[r]
    let a = if e[0] < e[1] { e[0] } else { e[1] }
    let b = if e[0] < e[1] { e[1] } else { e[0] }
    m[r * nverts + a] = mod_pos(-1, p)
    m[r * nverts + b] = mod_pos(m[r * nverts + b] + 1, p)
    r = r + 1
  }
  m
}

/// δ¹：C¹ → C²（各面取循環頂點序之帶號邊界）。
pub fn cech_d1(edges : Array[Array[Int]], faces : Array[Array[Int]], p : Int) -> Array[Int] {
  let ne = edges.length()
  let nf = faces.length()
  let m = Array::make(nf * ne, 0)
  let mut f = 0
  while f < nf {
    let v = faces[f]
    let k = v.length()
    let mut t = 0
    while t < k {
      let u = v[t]
      let w = v[mod_pos(t + 1, k)]
      let sgn = if u < w { 1 } else { -1 }
      let a = if u < w { u } else { w }
      let b = if u < w { w } else { u }
      let idx = edge_index(edges, a, b)
      if idx >= 0 {
        m[f * ne + idx] = mod_pos(m[f * ne + idx] + sgn, p)
      }
      t = t + 1
    }
    f = f + 1
  }
  m
}

/// δ¹∘δ⁰ = 0 之檢查（複形條件）。
pub fn cech_is_complex(nverts : Int, edges : Array[Array[Int]], faces : Array[Array[Int]], p : Int) -> Bool {
  if faces.length() == 0 {
    return true
  }
  let d0 = cech_d0(nverts, edges, p)
  let d1 = cech_d1(edges, faces, p)
  let prod = mat_mul_rect(d1, d0, faces.length(), edges.length(), nverts, p)
  let mut i = 0
  while i < prod.length() {
    if prod[i] != 0 {
      return false
    }
    i = i + 1
  }
  true
}

/// Čech 餘調維數 (h⁰, h¹, h²)。
pub fn cech_betti(
  nverts : Int,
  edges : Array[Array[Int]],
  faces : Array[Array[Int]],
  p : Int
) -> (Int, Int, Int) {
  let ne = edges.length()
  let nf = faces.length()
  let d0 = cech_d0(nverts, edges, p)
  let r0 = matrix_rank(d0, ne, nverts, p)
  let h0 = nverts - r0
  if nf == 0 {
    (h0, ne - r0, 0)
  } else {
    let d1 = cech_d1(edges, faces, p)
    let r1 = matrix_rank(d1, nf, ne, p)
    (h0, ne - r1 - r0, nf - r1)
  }
}

/// 置換是否為複形之自同構（每一細胞皆有像，即無 -1）。
pub fn perm_total(perm : Array[Int]) -> Bool {
  let mut ok = true
  for x in perm {
    if x < 0 {
      ok = false
    }
  }
  ok
}

/// 定向鏈上之帶號置換矩陣 S₁：邊 e_i = (a,b) 之像為邊 (g a, g b)，
/// 若 g a < g b 取 +1、否則取 −1（上鏈作用即取此矩陣）。
pub fn signed_edge_matrix(edges : Array[Array[Int]], vperm : Array[Int], p : Int) -> Array[Int] {
  let ne = edges.length()
  let m = Array::make(ne * ne, 0)
  let mut i = 0
  while i < ne {
    let e = edges[i]
    let a = vperm[e[0]]
    let b = vperm[e[1]]
    let idx = edge_index(edges, a, b)
    if idx >= 0 {
      let sgn = if a < b { 1 } else { -1 }
      m[idx * ne + i] = mod_pos(sgn, p)
    }
    i = i + 1
  }
  m
}

/// 置換之矩陣（列主序：e_i ↦ e_{perm[i]}）。
pub fn perm_matrix(n : Int, perm : Array[Int]) -> Array[Int] {
  let m = Array::make(n * n, 0)
  let mut i = 0
  while i < n {
    m[perm[i] * n + i] = 1
    i = i + 1
  }
  m
}

/// 頂點置換所誘導之邊置換（無向邊）。
pub fn induced_edge_perm(edges : Array[Array[Int]], vperm : Array[Int]) -> Array[Int] {
  let ne = edges.length()
  let out = Array::make(ne, -1)
  let mut i = 0
  while i < ne {
    let e = edges[i]
    let idx = edge_index(edges, vperm[e[0]], vperm[e[1]])
    out[i] = idx
    i = i + 1
  }
  out
}

/// 置換所誘導之面置換（各面之循環頂點序）。
pub fn induced_face_perm(faces : Array[Array[Int]], vperm : Array[Int]) -> Array[Int] {
  let nf = faces.length()
  let out = Array::make(nf, -1)
  let mut f = 0
  while f < nf {
    let v = faces[f]
    let img : Array[Int] = []
    for x in v {
      img.push(vperm[x])
    }
    // 與既有面比對（循環序）
    let mut g = 0
    while g < nf {
      let w = faces[g]
      if w.length() == img.length() {
        let k = w.length()
        let mut s = 0
        while s < k {
          let mut same = true
          let mut t = 0
          while t < k {
            if w[mod_pos(s + t, k)] != img[t] {
              same = false
            }
            t = t + 1
          }
          if same {
            out[f] = g
          }
          s = s + 1
        }
      }
      g = g + 1
    }
    f = f + 1
  }
  out
}

/// 不動子空間之基底 = ker(ρ − I)。
pub fn fixed_basis(rho : Array[Int], n : Int, p : Int) -> Array[Array[Int]] {
  let id = mat_id(n, p)
  let m = Array::make(n * n, 0)
  let mut i = 0
  while i < n * n {
    m[i] = mod_pos(rho[i] - id[i], p)
    i = i + 1
  }
  mat_kernel(m, n, n, p)
}

/// 以「基底向量堆疊成矩陣再取秩」計算 δ 限制於不動子空間之秩。
pub fn restricted_rank(
  d : Array[Int],
  rows : Int,
  cols : Int,
  basis : Array[Array[Int]],
  p : Int
) -> Int {
  if basis.length() == 0 {
    return 0
  }
  // 對每個基底向量取其像，堆疊為 rows×basis.length 之矩陣
  let m = Array::make(rows * basis.length(), 0)
  let mut c = 0
  for v in basis {
    let img = mat_vec_rect(d, v, rows, cols, p)
    let mut r = 0
    while r < rows {
      m[r * basis.length() + c] = img[r]
      r = r + 1
    }
    c = c + 1
  }
  matrix_rank(m, rows, basis.length(), p)
}

/// 長方形矩陣 × 向量。
pub fn mat_vec_rect(m : Array[Int], v : Array[Int], rows : Int, cols : Int, p : Int) -> Array[Int] {
  let out = Array::make(rows, 0)
  let mut i = 0
  while i < rows {
    let mut s = 0
    let mut j = 0
    while j < cols {
      s = s + m[i * cols + j] * v[j]
      j = j + 1
    }
    out[i] = mod_pos(s, p)
    i = i + 1
  }
  out
}

/// 固定細胞數（置換矩陣之跡）＝單純 Lefschetz 數之組成。
pub fn fixed_cells(perm : Array[Int]) -> Int {
  let mut c = 0
  let mut i = 0
  while i < perm.length() {
    if perm[i] == i {
      c = c + 1
    }
    i = i + 1
  }
  c
}

// ───────────────── 標準複形 ─────────────────

/// 圓周：三角形神經（3 頂點、3 邊、無 2-胞腔）。
pub fn circle_cover() -> (Int, Array[Array[Int]], Array[Array[Int]]) {
  (3, [[0, 1], [1, 2], [0, 2]], [])
}

/// 球面：四面體邊界（4 頂點、6 邊、4 面）。
pub fn sphere_cover() -> (Int, Array[Array[Int]], Array[Array[Int]]) {
  (4, [[0, 1], [1, 2], [0, 2], [0, 3], [1, 3], [2, 3]], [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
}

/// 環面：3×3 方格胞腔複形（9 頂點、18 邊、9 個四邊形面）。
pub fn torus_cover() -> (Int, Array[Array[Int]], Array[Array[Int]]) {
  let edges : Array[Array[Int]] = []
  let faces : Array[Array[Int]] = []
  let idx = fn(i : Int, j : Int) -> Int { mod_pos(i, 3) * 3 + mod_pos(j, 3) }
  let mut i = 0
  while i < 3 {
    let mut j = 0
    while j < 3 {
      edges.push([idx(i, j), idx(i, j + 1)])
      edges.push([idx(i, j), idx(i + 1, j)])
      faces.push([idx(i, j), idx(i, j + 1), idx(i + 1, j + 1), idx(i + 1, j)])
      j = j + 1
    }
    i = i + 1
  }
  (9, edges, faces)
}

/// 層作用之義務自證。
pub fn sheaf_register(led : Ledger) -> Unit {
  led.oblige("sheaf.cech_complex", A2_Action, SheafCech, "Čech 複形條件 δ¹∘δ⁰ = 0：圓周（三角形神經）、球面（四面體邊界）、環面（3×3 方格）三者於 F_p（p=5,7）皆成立", fn() {
    let (nv1, e1, f1) = circle_cover()
    let (nv2, e2, f2) = sphere_cover()
    let (nv3, e3, f3) = torus_cover()
    cech_is_complex(nv1, e1, f1, 5) &&
    cech_is_complex(nv2, e2, f2, 5) &&
    cech_is_complex(nv3, e3, f3, 5) &&
    cech_is_complex(nv2, e2, f2, 7) &&
    cech_is_complex(nv3, e3, f3, 7)
  })
  led.oblige("sheaf.cohomology", A2_Action, SheafCech, "餘調維數與 Euler 特徵數：圓周 (h⁰,h¹,h²) = (1,1,0)、χ = 0；球面 (1,0,1)、χ = 2；環面 (1,2,1)、χ = 0；且 Σ(−1)^i dim C^i = Σ(−1)^i h^i", fn() {
    let (nv1, e1, f1) = circle_cover()
    let (h10, h11, h12) = cech_betti(nv1, e1, f1, 5)
    let (nv2, e2, f2) = sphere_cover()
    let (h20, h21, h22) = cech_betti(nv2, e2, f2, 5)
    let (nv3, e3, f3) = torus_cover()
    let (h30, h31, h32) = cech_betti(nv3, e3, f3, 5)
    let chi1 = nv1 - e1.length() + f1.length()
    let chi2 = nv2 - e2.length() + f2.length()
    let chi3 = nv3 - e3.length() + f3.length()
    h10 == 1 && h11 == 1 && h12 == 0 && chi1 == 0 &&
    h20 == 1 && h21 == 0 && h22 == 1 && chi2 == 2 &&
    h30 == 1 && h31 == 2 && h32 == 1 && chi3 == 0 &&
    h10 - h11 + h12 == chi1 &&
    h20 - h21 + h22 == chi2 &&
    h30 - h31 + h32 == chi3
  })
  led.oblige("sheaf.equivariant", A2_Action, SheafCech, "等變層（層作用）：Z/3 循環置換圓周之三圖卡、以及球面之 3 階旋轉，皆使 δ⁰∘ρ⁰ = ρ¹∘δ⁰（ρ⁰ 為頂點置換矩陣、ρ¹ 為定向邊上之帶號置換矩陣）；且有面之複形中面亦被置換（無 -1）", fn() {
    let (nv, edges, _faces) = circle_cover()
    let rho0perm = [1, 2, 0]
    let d0 = cech_d0(nv, edges, 5)
    let left = mat_mul_rect(d0, perm_matrix(nv, rho0perm), edges.length(), nv, nv, 5)
    let right = mat_mul_rect(signed_edge_matrix(edges, rho0perm, 5), d0, edges.length(), edges.length(), nv, 5)
    let circle_ok = mat_eq(left, right)
    let (nv2, e2, f2) = sphere_cover()
    let s0 = [0, 2, 3, 1]
    let s1 = induced_edge_perm(e2, s0)
    let s2 = induced_face_perm(f2, s0)
    let d0s = cech_d0(nv2, e2, 5)
    let sphere_ok = mat_eq(
      mat_mul_rect(d0s, perm_matrix(nv2, s0), e2.length(), nv2, nv2, 5),
      mat_mul_rect(signed_edge_matrix(e2, s0, 5), d0s, e2.length(), e2.length(), nv2, 5)
    )
    circle_ok &&
    sphere_ok &&
    perm_total(s1) &&
    perm_total(s2) &&
    fixed_cells(rho0perm) == 0 &&
    fixed_cells(induced_edge_perm(edges, rho0perm)) == 0
  })
  led.oblige("sheaf.lefschetz", A2_Action, SheafCech, "Lefschetz–Hopf 跡公式：Σ(−1)^i tr(ρ^i) = Σ(−1)^i dim Fix(H^i)。圓周旋轉 L = 0；球面 3 階旋轉（固定 1 頂點 1 面）L = 2；環面平移 L = 0", fn() {
    // 圓周：Z/3 旋轉，跡由固定細胞數給出
    let (nv1, e1, _f1) = circle_cover()
    let r0 = [1, 2, 0]
    let r1 = induced_edge_perm(e1, r0)
    let l_circle = fixed_cells(r0) - fixed_cells(r1)
    // 球面：旋轉 (1 2 3) 固定頂點 0 與面 {1,2,3}
    let (nv2, e2, f2) = sphere_cover()
    let s = [0, 2, 3, 1]
    let s1 = induced_edge_perm(e2, s)
    let s2 = induced_face_perm(f2, s)
    let l_sphere = fixed_cells(s) - fixed_cells(s1) + fixed_cells(s2)
    // 環面：平移 (i,j) ↦ (i+1,j+1)
    let (nv3, e3, f3) = torus_cover()
    let t : Array[Int] = []
    let mut i = 0
    while i < 9 {
      t.push(mod_pos(i / 3 + 1, 3) * 3 + mod_pos(i % 3 + 1, 3))
      i = i + 1
    }
    let t1 = induced_edge_perm(e3, t)
    let t2 = induced_face_perm(f3, t)
    let l_torus = fixed_cells(t) - fixed_cells(t1) + fixed_cells(t2)
    let total_actions = perm_total(r1) && perm_total(s1) && perm_total(s2) && perm_total(t1) && perm_total(t2)
    let _ = nv1
    let _ = nv2
    let _ = nv3
    // 餘調側：p=5（∤ |G|）以不動子空間分解計算 dim Fix(H^i)
    let d0c = cech_d0(3, e1, 5)
    let fb0c = fixed_basis(perm_matrix(3, r0), 3, 5)
    let fix_h0c = fb0c.length() - restricted_rank(d0c, 3, 3, fb0c, 5)
    let fb1c = fixed_basis(perm_matrix(3, r1), 3, 5)
    let fix_h1c = fb1c.length() - restricted_rank(d0c, 3, 3, fb0c, 5)
    let l_h_circle = fix_h0c - fix_h1c
    // 球面
    let d0s = cech_d0(4, e2, 5)
    let fb0s = fixed_basis(perm_matrix(4, s), 4, 5)
    let fix_h0s = fb0s.length() - restricted_rank(d0s, 6, 4, fb0s, 5)
    let d1s = cech_d1(e2, f2, 5)
    let fb1s = fixed_basis(perm_matrix(6, s1), 6, 5)
    let fix_h1s = fb1s.length() - restricted_rank(d1s, 4, 6, fb1s, 5) - restricted_rank(d0s, 6, 4, fb0s, 5)
    let fb2s = fixed_basis(perm_matrix(4, s2), 4, 5)
    let fix_h2s = fb2s.length() - restricted_rank(d1s, 4, 6, fb1s, 5)
    let l_h_sphere = fix_h0s - fix_h1s + fix_h2s
    l_circle == 0 && l_h_circle == 0 && l_sphere == 2 && l_h_sphere == 2 && l_torus == 0 && total_actions
  })
  led.oblige("sheaf.glue", A2_Action, SheafCech, "黏合與下降（descent）：圓周上 1-上鏈 f 為上邊界 ⟺ f₀₁+f₁₂ = f₀₂；餘圈 (1,0,0) 之同調類生成 H¹（非上邊界）", fn() {
    let (nv, edges, faces) = circle_cover()
    let d0 = cech_d0(nv, edges, 11)
    // im δ⁰ 由 δ⁰ 之各行（各頂點基底之像）張成
    let ne = edges.length()
    let rows : Array[Array[Int]] = []
    let mut c = 0
    while c < nv {
      let v : Array[Int] = []
      let mut r = 0
      while r < ne {
        v.push(d0[r * nv + c])
        r = r + 1
      }
      rows.push(v)
      c = c + 1
    }
    let f_coboundary : Array[Int] = [1, 1, mod_pos(1 + 1, 11)]
    let f_cocycle_class : Array[Int] = [1, 0, 0]
    let _ = faces
    in_span(rows, f_coboundary, 11) &&
    !in_span(rows, f_cocycle_class, 11) &&
    matrix_rank(d0, ne, nv, 11) == 2
  })
}
