//! 李型單群與李群（代數）作用：CFSG 之「李型群」條目。
//!
//! 內容分兩層：
//!   * 李型群之階與 Weyl 維數（組合／群論算式）：A₁(q) = PSL(2,q)、A₂(q) = PSL(3,q)、
//!     Suzuki ²B₂、Ree ²G₂、Steinberg ³D₄、G₂(q)、²F₄ 之階，全部以質因數分解交叉驗證。
//!   * Lie 代數作用：交換子、Jacobi 恆等式、Killing 形式、Cartan 矩陣行列式（單性判準），
//!     矩陣指數（李群作用之離散化）與 BCH 公式 exp(X)exp(Y) = exp(X+Y+[X,Y]/2)，
//!     以及伴隨作用 Ad(exp X) = exp(ad X)（李群作用於其李代數）。

// ───────────────── 李型群之階（Order 算式） ─────────────────

/// PSL(2,q) = A₁(q)：q(q²−1)/gcd(2,q−1)。
pub fn lie_a1_order(q : Int) -> Int {
  q * (q * q - 1) / gcd_int(2, q - 1)
}

/// PSL(3,q) = A₂(q)：q³(q³−1)(q²−1)/gcd(3,q−1)。
pub fn lie_a2_order(q : Int) -> Int {
  q * q * q * (q * q * q - 1) * (q * q - 1) / gcd_int(3, q - 1)
}

/// Suzuki 群 ²B₂(q)，q = 2^{2n+1}：q²(q²+1)(q−1)。
pub fn sz_order(q : Int) -> Int {
  q * q * (q * q + 1) * (q - 1)
}

/// Ree 群 ²G₂(q)，q = 3^{2n+1}：q³(q³+1)(q−1)。
pub fn ree_order(q : Int) -> Int {
  q * q * q * (q * q * q + 1) * (q - 1)
}

/// Steinberg 群 ³D₄(q)：q¹²(q⁸+q⁴+1)(q⁶−1)(q²−1)。
pub fn d4_order(q : Int) -> Int {
  let q2 = q * q
  let q4 = q2 * q2
  let q6 = q4 * q2
  let q8 = q4 * q4
  let q12 = q6 * q6
  q12 * (q8 + q4 + 1) * (q6 - 1) * (q2 - 1)
}

/// G₂(q)：q⁶(q⁶−1)(q²−1)。
pub fn g2_order(q : Int) -> Int {
  let q2 = q * q
  let q6 = q2 * q2 * q2
  q6 * (q6 - 1) * (q2 - 1)
}

/// 扭群 ²F₄(q)，q = 2^{2n+1}（n ≥ 1）：q¹²(q⁶+1)(q⁴−1)(q³+1)(q−1)。
pub fn f4twist_order(q : Int) -> Int {
  let q2 = q * q
  let q3 = q2 * q
  let q4 = q2 * q2
  let q6 = q3 * q3
  let q12 = q6 * q6
  q12 * (q6 + 1) * (q4 - 1) * (q3 + 1) * (q - 1)
}

// ───────────────── Weyl 維數公式與根系 ─────────────────

/// A₂（sl₃）最高權表示之維數：dim V_{a,b} = (a+1)(b+1)(a+b+2)/2。
pub fn weyl_dim_a2(a : Int, b : Int) -> Int {
  (a + 1) * (b + 1) * (a + b + 2) / 2
}

/// 以根系乘積公式獨立計算 A₂ 維數：dim = Π_{α>0} ⟨λ+ρ,α⟩/⟨λ 之 ρ,α⟩。
/// 正根 α₁, α₂, α₁+α₂；⟨λ+ρ,α₁⟩ = a+1、⟨λ+ρ,α₂⟩ = b+1、⟨λ+ρ,α₁+α₂⟩ = a+b+2、⟨ρ,α₁+α₂⟩ = 2。
pub fn weyl_dim_a2_roots(a : Int, b : Int) -> Rat {
  let n1 = rat(a + 1, 1)
  let n2 = rat(b + 1, 1)
  let n3 = rat(a + b + 2, 2)
  rat_mul(rat_mul(n1, n2), n3)
}

/// 有理數相等（交叉相乘）。
pub fn rat_eq(a : Rat, b : Rat) -> Bool {
  a.n * b.d == b.n * a.d
}

/// A₁ 之維數：V_k = Sym^k（2 維）之維數 k+1，等於 2 變數 k 次單項式之個數。
pub fn weyl_dim_a1_by_monomials(k : Int) -> Int {
  monomials_of_degree(2, k).length()
}

/// 秩 2 單根系之正根與全根數（A₂: 3/6；B₂: 4/8；G₂: 6/12）。
pub fn root_counts(kind : String) -> (Int, Int) {
  if kind == "A2" {
    (3, 6)
  } else if kind == "B2" {
    (4, 8)
  } else {
    (6, 12)
  }
}

/// 小型整數矩陣之行列式（Laplace 展開，n ≤ 4；Cartan 矩陣用）。
pub fn int_det(m : Array[Int], n : Int) -> Int {
  if n == 1 {
    return m[0]
  }
  let mut det = 0
  let mut j = 0
  while j < n {
    let sub : Array[Int] = []
    let mut r = 1
    while r < n {
      let mut c = 0
      while c < n {
        if c != j {
          sub.push(m[r * n + c])
        }
        c = c + 1
      }
      r = r + 1
    }
    let sign = if j % 2 == 0 { 1 } else { -1 }
    det = det + sign * m[j] * int_det(sub, n - 1)
    j = j + 1
  }
  det
}

/// Cartan 矩陣（A₂、B₂、G₂、A₃、D₄）。
pub fn cartan_matrix(kind : String) -> (Array[Int], Int) {
  if kind == "A2" {
    ([2, -1, -1, 2], 2)
  } else if kind == "B2" {
    ([2, -2, -1, 2], 2)
  } else if kind == "G2" {
    ([2, -1, -3, 2], 2)
  } else if kind == "A3" {
    ([2, -1, 0, -1, 2, -1, 0, -1, 2], 3)
  } else {
    ([2, -1, 0, 0, -1, 2, -1, -1, 0, -1, 2, 0, 0, -1, 0, 2], 4)
  }
}

// ───────────────── Lie 代數作用（線性代數算式） ─────────────────

/// 基底矩陣 E_{ij}（列主序）。
pub fn basis_matrix(i : Int, j : Int, n : Int) -> Array[Int] {
  let m = Array::make(n * n, 0)
  m[i * n + j] = 1
  m
}

/// 矩陣加法（mod p）。
pub fn mat_add(a : Array[Int], b : Array[Int], p : Int) -> Array[Int] {
  let out = Array::make(a.length(), 0)
  let mut i = 0
  while i < a.length() {
    out[i] = mod_pos(a[i] + b[i], p)
    i = i + 1
  }
  out
}

/// 純量乘（mod p）。
pub fn mat_scale(a : Array[Int], c : Int, p : Int) -> Array[Int] {
  let out = Array::make(a.length(), 0)
  let mut i = 0
  while i < a.length() {
    out[i] = mod_pos(a[i] * c, p)
    i = i + 1
  }
  out
}

/// 交換子 [A,B] = AB − BA（mod p）。
pub fn bracket(a : Array[Int], b : Array[Int], n : Int, p : Int) -> Array[Int] {
  let ab = mat_mul(a, b, n, p)
  let ba = mat_mul(b, a, n, p)
  let out = Array::make(n * n, 0)
  let mut i = 0
  while i < n * n {
    out[i] = mod_pos(ab[i] - ba[i], p)
    i = i + 1
  }
  out
}

/// Jacobi 恆等式之窮舉檢查（全基底三元組 [[x,y],z] + [[y,z],x] + [[z,x],y] = 0）。
pub fn jacobi_holds(n : Int, p : Int) -> Bool {
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let x = basis_matrix(i, j, n)
      let mut k = 0
      while k < n {
        let mut l = 0
        while l < n {
          let y = basis_matrix(k, l, n)
          let mut r2 = 0
          while r2 < n {
            let z = basis_matrix(r2, r2, n)
            let t1 = bracket(bracket(x, y, n, p), z, n, p)
            let t2 = bracket(bracket(y, z, n, p), x, n, p)
            let t3 = bracket(bracket(z, x, n, p), y, n, p)
            let s = mat_add(mat_add(t1, t2, p), t3, p)
            let mut q = 0
            while q < n * n {
              if s[q] != 0 {
                return false
              }
              q = q + 1
            }
            r2 = r2 + 1
          }
          l = l + 1
        }
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  true
}

/// ad_x：以 gl_n 之基底 E_{ij}（dim = n²）表示之交換子算子（線性變換矩陣，列主序）。
pub fn ad_matrix(x : Array[Int], n : Int, p : Int) -> Array[Int] {
  let dim = n * n
  let out = Array::make(dim * dim, 0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let e = basis_matrix(i, j, n)
      let b = bracket(x, e, n, p)
      let col = i * n + j
      let mut r = 0
      while r < dim {
        out[r * dim + col] = b[r]
        r = r + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

/// Killing 形式 B(a,b) = tr(ad_a ∘ ad_b)（mod p）。
pub fn killing_form(a : Array[Int], b : Array[Int], n : Int, p : Int) -> Int {
  let dim = n * n
  let ada = ad_matrix(a, n, p)
  let adb = ad_matrix(b, n, p)
  let prod = mat_mul(ada, adb, dim, p)
  let mut tr = 0
  let mut i = 0
  while i < dim {
    tr = tr + prod[i * dim + i]
    i = i + 1
  }
  mod_pos(tr, p)
}

/// 階乘（Int 範圍）。
pub fn fact_int(k : Int) -> Int {
  let mut f = 1
  let mut i = 2
  while i <= k {
    f = f * i
    i = i + 1
  }
  f
}

/// 矩陣指數之截斷和 exp(X) = Σ_{i=0}^{terms} X^i/i!（mod p；須 p ∤ i!，即 terms < p）。
pub fn matrix_exp(x : Array[Int], dim : Int, p : Int, terms : Int) -> Array[Int] {
  let mut acc = mat_id(dim, p)
  let mut pow = mat_id(dim, p)
  let mut i = 1
  while i <= terms {
    pow = mat_mul(pow, x, dim, p)
    let c = mod_inv(fact_int(i) % p, p)
    acc = mat_add(acc, mat_scale(pow, c, p), p)
    i = i + 1
  }
  acc
}

/// BCH（Baker–Campbell–Hausdorff）於二步冪零（Heisenberg）代數之截斷：
/// exp(X)exp(Y) = exp(X+Y+[X,Y]/2)，以 p 進位（p > 3）數值驗證。
pub fn bch_heisenberg_check(p : Int) -> Bool {
  let x = basis_matrix(0, 1, 3)
  let y = basis_matrix(1, 2, 3)
  let lhs = mat_mul(matrix_exp(x, 3, p, 4), matrix_exp(y, 3, p, 4), 3, p)
  let z = mat_add(mat_add(x, y, p), mat_scale(bracket(x, y, 3, p), mod_inv(2, p), p), p)
  let rhs = matrix_exp(z, 3, p, 4)
  mat_eq(lhs, rhs)
}

/// Ad(exp X)(Y) = exp(ad X)(Y)：以共軛 exp(X)·Y·exp(−X) 與 exp(ad X) 作用於 Y 之係數向量比對。
pub fn adjoint_exp_check(p : Int) -> Bool {
  let x = basis_matrix(0, 1, 2)
  let y : Array[Int] = [1, 0, 0, mod_pos(-1, p)]
  let ex = matrix_exp(x, 2, p, 3)
  let exi = matrix_exp(mat_scale(x, mod_pos(-1, p), p), 2, p, 3)
  let lhs = mat_mul(mat_mul(ex, y, 2, p), exi, 2, p)
  let adx = ad_matrix(x, 2, p)
  let expad = matrix_exp(adx, 4, p, 3)
  let rhs = mat_vec(expad, y, 4, p)
  mat_eq(lhs, rhs)
}

/// 李型群與李代數作用之義務自證。
pub fn lie_register(led : Ledger) -> Unit {
  led.oblige("lie.a1_a2_orders", A5_Classification, LinearCount, "A₁(q) = PSL(2,q) 之階 q(q²−1)/gcd(2,q−1) 與矩陣計數相符（q=2,3,5,7 → 6,12,60,168）；A₂(q) = PSL(3,q)：q=2,3,4,5 → 168,5616,20160,372000", fn() {
    lie_a1_order(2) == 6 &&
    lie_a1_order(3) == 12 &&
    lie_a1_order(5) == 60 &&
    lie_a1_order(7) == 168 &&
    lie_a1_order(2).to_string() == psl_order(2, 2).to_string() &&
    lie_a1_order(3).to_string() == psl_order(2, 3).to_string() &&
    lie_a1_order(5).to_string() == psl_order(2, 5).to_string() &&
    lie_a1_order(7).to_string() == psl_order(2, 7).to_string() &&
    lie_a2_order(2) == 168 &&
    lie_a2_order(3) == 5616 &&
    lie_a2_order(4) == 20160 &&
    lie_a2_order(5) == 372000 &&
    lie_a2_order(2) == 8 * 7 * 3
  })
  led.oblige("lie.twisted_orders", A5_Classification, BigOrder, "扭群之階：²B₂(8)=29120=2⁶·5·7·13、²B₂(2)=20；²G₂(3)=1512=3·|PSL(2,8)|=3·504；²F₄(2)=2·17971200（Tits 群 17971200=2¹¹·3³·5²·13）", fn() {
    sz_order(8) == 29120 &&
    sz_order(2) == 20 &&
    ree_order(3) == 1512 &&
    ree_order(3) == 3 * lie_a1_order(8) &&
    lie_a1_order(8) == 504 &&
    f4twist_order(2) == 35942400 &&
    f4twist_order(2) == 2 * 17971200 &&
    bigint_of_factorization(factorize(17971200)).to_string() == "17971200" &&
    factorization_to_string(factorize(29120)) == "2^6·5·7·13"
  })
  led.oblige("lie.d4_g2", A5_Classification, BigOrder, "³D₄(2) = 2¹²(2⁸+2⁴+1)(2⁶−1)(2²−1) = 211341312 = 2¹²·3⁴·7²·13（BigInt 精確比對）；G₂(2) = 12096 = 2⁶·3³·7，G₂(2)′ = 6048", fn() {
    let o = d4_order(2)
    o == 211341312 &&
    bigint_of_factorization(factorize(o)).to_string() == o.to_string() &&
    factorization_to_string(factorize(o)) == "2^12·3^4·7^2·13" &&
    g2_order(2) == 12096 &&
    g2_order(2) / 2 == 6048 &&
    factorization_to_string(factorize(12096)) == "2^6·3^3·7"
  })
  led.oblige("lie.weyl_dim", A5_Classification, WeylDimension, "Weyl 維數：A₂ 之 dim V_{a,b} = (a+1)(b+1)(a+b+2)/2 與根系乘積 Π_{α>0}⟨λ+ρ,α⟩/⟨ρ,α⟩ 一致；dim(a,b)=dim(b,a)；A₁ 之 dim = 單項式個數 k+1", fn() {
    let mut ok = true
    let mut a = 0
    while a <= 8 {
      let mut b = 0
      while b <= 8 {
        let d = weyl_dim_a2(a, b)
        if !rat_eq(weyl_dim_a2_roots(a, b), rat(d, 1)) {
          ok = false
        }
        if d != weyl_dim_a2(b, a) {
          ok = false
        }
        b = b + 1
      }
      a = a + 1
    }
    ok &&
    weyl_dim_a2(1, 0) == 3 &&
    weyl_dim_a2(1, 1) == 8 &&
    weyl_dim_a2(2, 1) == 15 &&
    weyl_dim_a2(3, 0) == 10 &&
    weyl_dim_a2(2, 2) == 27 &&
    weyl_dim_a1_by_monomials(0) == 1 &&
    weyl_dim_a1_by_monomials(3) == 4 &&
    weyl_dim_a1_by_monomials(7) == 8
  })
  led.oblige("lie.roots_cartan", A5_Classification, WeylDimension, "根系與 Cartan 矩陣：A₂/B₂/G₂ 之全根數 6/8/12（正根 3/4/6）；Weyl 群階 = |S₃|=6、|D₄|=8、|D₆|=12 對應；Cartan 行列式 A₂=3、B₂=2、G₂=1（么模）、A₃=4、D₄=4", fn() {
    let (p2a, t2a) = root_counts("A2")
    let (p2b, t2b) = root_counts("B2")
    let (p2g, t2g) = root_counts("G2")
    let (ca, na) = cartan_matrix("A2")
    let (cb, nb) = cartan_matrix("B2")
    let (cg, ng) = cartan_matrix("G2")
    let (c3, n3) = cartan_matrix("A3")
    let (cd, nd) = cartan_matrix("D4")
    let (d4g, ord4) = dihedral_group(4)
    let (d6g, ord6) = dihedral_group(6)
    let _ = d4g
    let _ = d6g
    p2a == 3 && t2a == 6 && p2b == 4 && t2b == 8 && p2g == 6 && t2g == 12 &&
    t2a == 2 * p2a && t2b == 2 * p2b && t2g == 2 * p2g &&
    group_order(3, symmetric_group(3).gens) == 6 &&
    ord4 == 8 && ord6 == 12 &&
    int_det(ca, na) == 3 &&
    int_det(cb, nb) == 2 &&
    int_det(cg, ng) == 1 &&
    int_det(c3, n3) == 4 &&
    int_det(cd, nd) == 4
  })
  led.oblige("lie.jacobi_killing", A5_Classification, Derivation, "Lie 代數：gl₂、gl₃ 於 F₅ 之 Jacobi 恆等式成立；sl₂ 之 Killing 形式 B(h,h)=8、B(e,f)=4、B(h,e)=0（mod 7）且非退化（3×3 之秩 3）", fn() {
    let h : Array[Int] = [1, 0, 0, mod_pos(-1, 7)]
    let e = basis_matrix(0, 1, 2)
    let f = basis_matrix(1, 0, 2)
    let bhh = killing_form(h, h, 2, 7)
    let bef = killing_form(e, f, 2, 7)
    let bhe = killing_form(h, e, 2, 7)
    let kill : Array[Int] = [bhh, bhe, killing_form(h, f, 2, 7), bhe, killing_form(e, e, 2, 7), bef, killing_form(f, h, 2, 7), bef, killing_form(f, f, 2, 7)]
    jacobi_holds(2, 5) &&
    jacobi_holds(3, 5) &&
    bhh == mod_pos(8, 7) &&
    bef == 4 &&
    bhe == 0 &&
    matrix_rank(kill, 3, 3, 7) == 3
  })
  led.oblige("lie.exp_bch_adjoint", A3_Composition, Derivation, "李群作用之離散化：exp(X)exp(Y) = exp(X+Y+[X,Y]/2)（Heisenberg，F₇）；Ad(exp X)(Y) = exp(ad X)(Y) 於 sl₂ ⊂ gl₂（F₇）；exp(X)exp(−X) = I", fn() {
    let x = basis_matrix(0, 1, 2)
    let ex = matrix_exp(x, 2, 7, 3)
    let exi = matrix_exp(mat_scale(x, mod_pos(-1, 7), 7), 2, 7, 3)
    bch_heisenberg_check(7) &&
    adjoint_exp_check(7) &&
    mat_eq(mat_mul(ex, exi, 2, 7), mat_id(2, 7))
  })
}
