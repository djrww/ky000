//! 數論引擎：階、分解、模算術。
//!
//! 這些是「衍生演算法」的最底層：所有群階計算、Sylow 分析、散在單群大整數階，
//! 皆化約至本檔之算式。每一函式皆對應 `num_register` 中之一條義務自證。

/// 正規化餘數：回傳 [0,m) 內之 r ≡ a (mod m)。（MoonBit 之 `%` 向零截斷，故需此式。）
pub fn mod_pos(a : Int, m : Int) -> Int {
  let r = a % m
  if r < 0 { r + m } else { r }
}

/// 快速冪取模：b^e mod m（e ≥ 0）。
pub fn pow_mod(b : Int, e : Int, m : Int) -> Int {
  let mut base = mod_pos(b, m)
  let mut exp = e
  let mut acc = 1 % m
  while exp > 0 {
    if exp % 2 == 1 {
      acc = mod_pos(acc * base, m)
    }
    base = mod_pos(base * base, m)
    exp = exp / 2
  }
  acc
}

/// 整數平方根（向下取整）。
pub fn int_sqrt(n : Int) -> Int {
  if n <= 0 {
    return 0
  }
  let mut r = 1
  while (r + 1) * (r + 1) <= n {
    r = r + 1
  }
  r
}

/// 質數判定（試除法，O(√n)）。
pub fn is_prime(n : Int) -> Bool {
  if n < 2 {
    return false
  }
  if n % 2 == 0 {
    return n == 2
  }
  let mut d = 3
  while d * d <= n {
    if n % d == 0 {
      return false
    }
    d = d + 2
  }
  true
}

/// 質因數分解：回傳 (p, e) 之列表，p 遞增，Π p^e = n。
pub fn factorize(n : Int) -> Array[(Int, Int)] {
  let out : Array[(Int, Int)] = []
  let mut m = if n < 0 { -n } else { n }
  let mut d = 2
  while d * d <= m {
    if m % d == 0 {
      let mut e = 0
      while m % d == 0 {
        m = m / d
        e = e + 1
      }
      out.push((d, e))
    }
    d = d + 1
  }
  if m > 1 {
    out.push((m, 1))
  }
  out
}

/// 分解之乘積（義務自證：須等於原數）。
pub fn factor_product(f : Array[(Int, Int)]) -> Int {
  let mut p = 1
  for __pr in f {
    let (b, e) = __pr
    let mut k = 0
    while k < e {
      p = p * b
      k = k + 1
    }
  }
  p
}

/// p 在 n 中之指數 v_p(n)。
pub fn valuation(n : Int, p : Int) -> Int {
  let mut m = if n < 0 { -n } else { n }
  let mut v = 0
  while m > 0 && m % p == 0 {
    m = m / p
    v = v + 1
  }
  v
}

/// Legendre 公式：v_p(n!) = Σ_{k≥1} ⌊n/p^k⌋。
pub fn legendre_factorial(n : Int, p : Int) -> Int {
  let mut v = 0
  let mut q = n / p
  while q > 0 {
    v = v + q
    q = q / p
  }
  v
}

/// 擴展歐幾里得：回傳 (g, x, y) 使 a·x + b·y = g = gcd(a,b)。
pub fn ext_gcd(a : Int, b : Int) -> (Int, Int, Int) {
  if b == 0 {
    let s = if a < 0 { -1 } else { 1 }
    return (if a < 0 { -a } else { a }, s, 0)
  }
  let (g, x, y) = ext_gcd(b, a % b)
  (g, y, x - (a / b) * y)
}

/// 最大公約數。
pub fn gcd_int(a : Int, b : Int) -> Int {
  let (g, _, _) = ext_gcd(a, b)
  g
}

/// 模反元素：a⁻¹ mod m（要求 gcd(a,m)=1），以擴展歐幾里得計算。
pub fn mod_inv(a : Int, m : Int) -> Int {
  let (g, x, _) = ext_gcd(a, m)
  if g != 1 {
    abort("mod_inv: 非互質")
  }
  mod_pos(x, m)
}

/// Euler φ 函式（由分解計算）。
pub fn totient(n : Int) -> Int {
  let mut r = n
  for __pr in factorize(n) {
    let (p, _e) = __pr
    r = r / p * (p - 1)
  }
  r
}

/// 正因數個數 τ(n)。
pub fn num_divisors(n : Int) -> Int {
  let mut t = 1
  for __pr in factorize(n) {
    let (_p, e) = __pr
    t = t * (e + 1)
  }
  t
}

/// 因數和 σ(n)。
pub fn sigma_sum(n : Int) -> Int {
  let mut t = 1
  for __pr in factorize(n) {
    let (p, e) = __pr
    let mut s = 1
    let mut pk = 1
    let mut k = 0
    while k < e {
      pk = pk * p
      s = s + pk
      k = k + 1
    }
    t = t * s
  }
  t
}

/// 所有正因數（遞增）。
pub fn divisors(n : Int) -> Array[Int] {
  let out : Array[Int] = []
  let mut d = 1
  while d * d <= n {
    if n % d == 0 {
      out.push(d)
      if d != n / d {
        out.push(n / d)
      }
    }
    d = d + 1
  }
  // 插入排序（避免額外依賴）
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] > key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

/// Möbius 函式 μ(n)。
pub fn mobius(n : Int) -> Int {
  if n == 1 {
    return 1
  }
  let f = factorize(n)
  for __pr in f {
    let (_p, e) = __pr
    if e > 1 {
      return 0
    }
  }
  if f.length() % 2 == 0 { 1 } else { -1 }
}

/// 二項係數 C(n,k)。
pub fn binomial(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    return 0
  }
  let mut r = 1
  let mut i = 0
  while i < k {
    r = r * (n - i) / (i + 1)
    i = i + 1
  }
  r
}

/// 前 n 個質數之列表。
pub fn primes_upto(n : Int) -> Array[Int] {
  let out : Array[Int] = []
  let mut k = 2
  while k <= n {
    if is_prime(k) {
      out.push(k)
    }
    k = k + 1
  }
  out
}

/// n 是否為質數冪（含 1 視為 p^0，回傳 (true, p, e)）。
pub fn prime_power(n : Int) -> (Bool, Int, Int) {
  if n <= 1 {
    return (true, n, 1)
  }
  let f = factorize(n)
  if f.length() == 1 {
    let (p, e) = f[0]
    (true, p, e)
  } else {
    (false, 0, 0)
  }
}

/// 由質因數分解建立大整數（散在單群階用）。
pub fn bigint_of_factorization(f : Array[(Int, Int)]) -> BigInt {
  let mut acc : BigInt = 1N
  for __pr in f {
    let (p, e) = __pr
    let mut k = 0
    while k < e {
      acc = acc * BigInt::from_int(p)
      k = k + 1
    }
  }
  acc
}

/// 指數比較：a 之分解是否整除 M 之分解（逐質數比較指數）。
pub fn factorization_divides(a : Array[(Int, Int)], m : Array[(Int, Int)]) -> Bool {
  let table : Map[Int, Int] = Map([])
  for __pr in m {
    let (p, e) = __pr
    table[p] = e
  }
  for __pr in a {
    let (p, e) = __pr
    match table.get(p) {
      Some(f) => if e > f { return false }
      None => return false
    }
  }
  true
}

/// 兩分解之乘積。
pub fn factorization_mul(a : Array[(Int, Int)], b : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let table : Map[Int, Int] = Map([])
  for __pr in a {
    let (p, e) = __pr
    table[p] = table.get_or_default(p, 0) + e
  }
  for __pr in b {
    let (p, e) = __pr
    table[p] = table.get_or_default(p, 0) + e
  }
  let keys : Array[Int] = []
  table.each(fn(k, _v) { keys.push(k) })
  keys.sort()
  let out : Array[(Int, Int)] = []
  for k in keys {
    out.push((k, table.get_or_default(k, 0)))
  }
  out
}

/// 分解之標準字串，如 "2^4·3^2·5·11"。
pub fn factorization_to_string(f : Array[(Int, Int)]) -> String {
  let mut s = ""
  let mut first = true
  for __pr in f {
    let (p, e) = __pr
    if !first {
      s = s + "·"
    }
    first = false
    s = if e == 1 { s + "\{p}" } else { s + "\{p}^\{e}" }
  }
  s
}

/// 從分解取 p-部分 p^e。
pub fn p_part(f : Array[(Int, Int)], p : Int) -> Int {
  for __pr in f {
    let (q, e) = __pr
    if q == p {
      let mut acc = 1
      let mut k = 0
      while k < e {
        acc = acc * p
        k = k + 1
      }
      return acc
    }
  }
  1
}

/// 階的「是否為偶數」＝ 分解含質數 2。
pub fn factorization_is_even(f : Array[(Int, Int)]) -> Bool {
  for __pr in f {
    let (p, _e) = __pr
    if p == 2 {
      return true
    }
  }
  false
}

/// 分解中相異質數個數 ω(n)。
pub fn omega(n : Int) -> Int {
  factorize(n).length()
}

// ───────────────────────── 義務自證 ─────────────────────────

/// 註冊數論引擎之義務自證。
pub fn num_register(led : Ledger) -> Unit {
  led.oblige("num.is_prime", A1_Group, Order, "質數判定與試除法一致：π(100)=25，且 561 非質數（Carmichael 數）", fn() {
    primes_upto(100).length() == 25 && !is_prime(561) && is_prime(7919) && !is_prime(1)
  })
  led.oblige("num.factorize", A1_Group, Order, "算術基本定理：Π p^e = n，且因子皆質數（對 n≤2000 全測）", fn() {
    let mut ok = true
    let mut n = 2
    while n <= 2000 {
      let f = factorize(n)
      if factor_product(f) != n {
        ok = false
      }
      for __pr in f {
        let (p, _e) = __pr
        if !is_prime(p) {
          ok = false
        }
      }
      n = n + 1
    }
    ok
  })
  led.oblige("num.pow_mod", A1_Group, Order, "快速冪正確性：b^e mod m 與費馬小定理 2^(p−1)≡1 (mod p)", fn() {
    let mut ok = true
    let mut e = 0
    while e <= 12 {
      if pow_mod(3, e, 1000000) != int_pow_small(3, e) {
        ok = false
      }
      e = e + 1
    }
    for p in primes_upto(60) {
      if p > 2 && pow_mod(2, p - 1, p) != 1 {
        ok = false
      }
    }
    ok
  })
  led.oblige("num.ext_gcd", A1_Group, Order, "Bézout 等式：a·x+b·y=gcd(a,b)（對 3≤a,b≤40 全測）", fn() {
    let mut ok = true
    let mut a = 3
    while a <= 40 {
      let mut b = 3
      while b <= 40 {
        let (g, x, y) = ext_gcd(a, b)
        if a * x + b * y != g || g != gcd_int(a, b) {
          ok = false
        }
        b = b + 1
      }
      a = a + 1
    }
    ok
  })
  led.oblige("num.mod_inv", A1_Group, Order, "模反元素：a·a⁻¹ ≡ 1 (mod m)（對所有與 m 互質之 a）", fn() {
    let m = 1009
    let mut ok = true
    let mut a = 1
    while a < m {
      if gcd_int(a, m) == 1 && mod_pos(a * mod_inv(a, m), m) != 1 {
        ok = false
      }
      a = a + 1
    }
    ok
  })
  led.oblige("num.totient", A1_Group, Order, "φ(n)=|(Z/n)^×|：直接計數與分解式一致（n≤300）", fn() {
    let mut ok = true
    let mut n = 1
    while n <= 300 {
      let mut c = 0
      let mut a = 1
      while a <= n {
        if gcd_int(a, n) == 1 {
          c = c + 1
        }
        a = a + 1
      }
      if c != totient(n) {
        ok = false
      }
      n = n + 1
    }
    ok
  })
  led.oblige("num.divisors", A1_Group, Order, "因數計數：τ(n)=|divisors(n)|，且 σ(n)=Σ 因數（n≤400）", fn() {
    let mut ok = true
    let mut n = 1
    while n <= 400 {
      let d = divisors(n)
      let mut sum = 0
      for x in d {
        sum = sum + x
      }
      if d.length() != num_divisors(n) || sum != sigma_sum(n) {
        ok = false
      }
      n = n + 1
    }
    ok
  })
  led.oblige("num.legendre", A1_Group, Order, "Legendre 公式：v_p(n!)=Σ⌊n/p^k⌋（n≤200, p≤13）", fn() {
    let mut ok = true
    for p in primes_upto(13) {
      let mut n = 1
      while n <= 300 {
        // 另一種算法：逐項累加 v_p(k)（避免階乘溢位）
        let mut acc = 0
        let mut k = 1
        while k <= n {
          acc = acc + valuation(k, p)
          k = k + 1
        }
        if acc != legendre_factorial(n, p) {
          ok = false
        }
        n = n + 1
      }
    }
    ok
  })
  led.oblige("num.mobius", A1_Group, Order, "Möbius 反演：Σ_{d|n} μ(d) = [n=1]（n≤500）", fn() {
    let mut ok = true
    let mut n = 1
    while n <= 500 {
      let mut s = 0
      for d in divisors(n) {
        s = s + mobius(d)
      }
      let want = if n == 1 { 1 } else { 0 }
      if s != want {
        ok = false
      }
      n = n + 1
    }
    ok
  })
  led.oblige("num.factorization_algebra", A1_Group, BigOrder, "分解代數：factorization_mul 之積等於大整數乘積，整除判定與 BigInt 餘數一致", fn() {
    let a : Array[(Int, Int)] = [(2, 4), (3, 2), (5, 1), (11, 1)]
    let b : Array[(Int, Int)] = [(2, 2), (7, 1), (11, 1)]
    let ab = factorization_mul(a, b)
    let ba = bigint_of_factorization(a) * bigint_of_factorization(b)
    if bigint_of_factorization(ab).to_string() != ba.to_string() {
      return false
    }
    factorization_divides(a, ab) && factorization_divides(b, ab) && !factorization_divides(ab, a)
  })
}

/// 小指數冪（僅供 pow_mod 之義務自證作對照）。
fn int_pow_small(b : Int, e : Int) -> Int {
  let mut acc = 1
  let mut k = 0
  while k < e {
    acc = acc * b
    k = k + 1
  }
  acc
}
