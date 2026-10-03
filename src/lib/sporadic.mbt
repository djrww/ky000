//! 26 個散在單群（sporadic simple groups）：CFSG 之「26 個例外」條目。
//!
//! 散在單群之階過大（|M| ≈ 8.08×10⁵³），故以「十進位字串 × 質因數指數表」雙表示法儲存，
//! 並以 BigInt 精確乘法互相驗證（不使用浮點、不使用第 3 方資料庫）：
//!   * 指數表乘積之十進位展開必等於所載之階；
//!   * 階 < 2³¹ 之群另以 Int 乘法獨立複核；
//!   * 與 Monster 階之整除結構：26 群中 24 群之階整除 |M|，例外為 J₄ 與 Ly；
//!   * |M| : |B| = 2⁵·3⁷·5³·7⁴·11·13²·29·41·59·71（以指數加法驗證）。

/// 散在單群之資料：名稱、階之十進位表示、質因數分解（指數表）。
pub struct Sporadic {
  name : String
  decimal : String
  fac : Array[(Int, Int)]
}

/// 26 個散在單群之名稱與階（依 CFSG 之分類列表）。
pub fn sporadics() -> Array[Sporadic] {
  [
    { name: "M11", decimal: "7920", fac: [(2, 4), (3, 2), (5, 1), (11, 1)] },
    { name: "M12", decimal: "95040", fac: [(2, 6), (3, 3), (5, 1), (11, 1)] },
    { name: "J1", decimal: "175560", fac: [(2, 3), (3, 1), (5, 1), (7, 1), (11, 1), (19, 1)] },
    { name: "M22", decimal: "443520", fac: [(2, 7), (3, 2), (5, 1), (7, 1), (11, 1)] },
    { name: "J2", decimal: "604800", fac: [(2, 7), (3, 3), (5, 2), (7, 1)] },
    { name: "M23", decimal: "10200960", fac: [(2, 7), (3, 2), (5, 1), (7, 1), (11, 1), (23, 1)] },
    { name: "HS", decimal: "44352000", fac: [(2, 9), (3, 2), (5, 3), (7, 1), (11, 1)] },
    { name: "J3", decimal: "50232960", fac: [(2, 7), (3, 5), (5, 1), (17, 1), (19, 1)] },
    { name: "M24", decimal: "244823040", fac: [(2, 10), (3, 3), (5, 1), (7, 1), (11, 1), (23, 1)] },
    { name: "McL", decimal: "898128000", fac: [(2, 7), (3, 6), (5, 3), (7, 1), (11, 1)] },
    { name: "He", decimal: "4030387200", fac: [(2, 10), (3, 3), (5, 2), (7, 3), (17, 1)] },
    { name: "Ru", decimal: "145926144000", fac: [(2, 14), (3, 3), (5, 3), (7, 1), (13, 1), (29, 1)] },
    { name: "Suz", decimal: "448345497600", fac: [(2, 13), (3, 7), (5, 2), (7, 1), (11, 1), (13, 1)] },
    { name: "ON", decimal: "460815505920", fac: [(2, 9), (3, 4), (5, 1), (7, 3), (11, 1), (19, 1), (31, 1)] },
    { name: "Co3", decimal: "495766656000", fac: [(2, 10), (3, 7), (5, 3), (7, 1), (11, 1), (23, 1)] },
    { name: "Co2", decimal: "42305421312000", fac: [(2, 18), (3, 6), (5, 3), (7, 1), (11, 1), (23, 1)] },
    { name: "Fi22", decimal: "64561751654400", fac: [(2, 17), (3, 9), (5, 2), (7, 1), (11, 1), (13, 1)] },
    { name: "HN", decimal: "273030912000000", fac: [(2, 14), (3, 6), (5, 6), (7, 1), (11, 1), (19, 1)] },
    { name: "Ly", decimal: "51765179004000000", fac: [(2, 8), (3, 7), (5, 6), (7, 1), (11, 1), (31, 1), (37, 1), (67, 1)] },
    { name: "Th", decimal: "90745943887872000", fac: [(2, 15), (3, 10), (5, 3), (7, 2), (13, 1), (19, 1), (31, 1)] },
    { name: "Fi23", decimal: "4089470473293004800", fac: [(2, 18), (3, 13), (5, 2), (7, 1), (11, 1), (13, 1), (17, 1), (23, 1)] },
    { name: "Co1", decimal: "4157776806543360000", fac: [(2, 21), (3, 9), (5, 4), (7, 2), (11, 1), (13, 1), (23, 1)] },
    { name: "J4", decimal: "86775571046077562880", fac: [(2, 21), (3, 3), (5, 1), (7, 1), (11, 3), (23, 1), (29, 1), (31, 1), (37, 1), (43, 1)] },
    { name: "Fi24'", decimal: "1255205709190661721292800", fac: [(2, 21), (3, 16), (5, 2), (7, 3), (11, 1), (13, 1), (17, 1), (23, 1), (29, 1)] },
    { name: "B", decimal: "4154781481226426191177580544000000", fac: [(2, 41), (3, 13), (5, 6), (7, 2), (11, 1), (13, 1), (17, 1), (19, 1), (23, 1), (31, 1), (47, 1)] },
    { name: "M", decimal: "808017424794512875886459904961710757005754368000000000", fac: [(2, 46), (3, 20), (5, 9), (7, 6), (11, 2), (13, 3), (17, 1), (19, 1), (23, 1), (29, 1), (31, 1), (41, 1), (47, 1), (59, 1), (71, 1)] },
  ]
}

/// 依名尋索引（-1 表示不存在）。
pub fn sporadic_index(name : String) -> Int {
  let arr = sporadics()
  let mut i = 0
  while i < arr.length() {
    if arr[i].name == name {
      return i
    }
    i = i + 1
  }
  -1
}

/// 最大質因數。
pub fn max_prime_factor(f : Array[(Int, Int)]) -> Int {
  let mut m = 1
  for pr in f {
    let (p, _e) = pr
    if p > m {
      m = p
    }
  }
  m
}

/// 指數表之整數乘積（僅供 < 2³¹ 之階使用）。
pub fn factor_product_small(f : Array[(Int, Int)]) -> Int {
  let mut acc = 1
  for pr in f {
    let (p, e) = pr
    acc = acc * int_pow(p, e)
  }
  acc
}

/// 散在單群之義務自證。
pub fn sporadic_register(led : Ledger) -> Unit {
  led.oblige("sporadic.count26", A5_Classification, BigOrder, "散在單群恰 26 個（CFSG 之 26 個例外）：名稱互異，且 J₁、J₃、J₄、Ly、O′N、Ru 六個 pariah 皆在列", fn() {
    let arr = sporadics()
    if arr.length() != 26 {
      return false
    }
    let seen : Map[String, Unit] = Map([])
    let mut distinct = true
    for s in arr {
      if seen.contains(s.name) {
        distinct = false
      }
      seen[s.name] = ()
    }
    distinct &&
    sporadic_index("J1") >= 0 &&
    sporadic_index("J3") >= 0 &&
    sporadic_index("J4") >= 0 &&
    sporadic_index("Ly") >= 0 &&
    sporadic_index("ON") >= 0 &&
    sporadic_index("Ru") >= 0
  })
  led.oblige("sporadic.factor_matches_decimal", A5_Classification, BigOrder, "26 個散在單群階之雙表示一致：質因數指數表之 BigInt 乘積 = 所載十進位階（逐群精確比對）", fn() {
    let arr = sporadics()
    let mut ok = true
    for s in arr {
      if bigint_of_factorization(s.fac).to_string() != s.decimal {
        ok = false
      }
    }
    ok
  })
  led.oblige("sporadic.small_int_products", A5_Classification, BigOrder, "階 < 2³¹ 之 10 群以 Int 乘法獨立複核：M₁₁=7920、M₁₂=95040、J₁=175560、M₂₂=443520、J₂=604800、M₂₃=10200960、HS=44352000、J₃=50232960、M₂₄=244823040、McL=898128000", fn() {
    let names = ["M11", "M12", "J1", "M22", "J2", "M23", "HS", "J3", "M24", "McL"]
    let arr = sporadics()
    let mut ok = true
    for nm in names {
      let i = sporadic_index(nm)
      if i < 0 {
        ok = false
      } else {
        let s = arr[i]
        if factor_product_small(s.fac).to_string() != s.decimal {
          ok = false
        }
      }
    }
    let hi = sporadic_index("He")
    hi >= 0 && factorization_to_string(arr[hi].fac) == "2^10·3^3·5^2·7^3·17" && ok
  })
  led.oblige("sporadic.monster_divisibility", A5_Classification, BigOrder, "與 Monster 階之整除結構：26 群中 24 群之階整除 |M|，例外恰為 J₄（11³、37、43）與 Ly（37、67）", fn() {
    let arr = sporadics()
    let mi = sporadic_index("M")
    let mfac = arr[mi].fac
    let mut divides = 0
    let exceptions : Array[String] = []
    for s in arr {
      if factorization_divides(s.fac, mfac) {
        divides = divides + 1
      } else {
        exceptions.push(s.name)
      }
    }
    divides == 24 &&
    exceptions.length() == 2 &&
    exceptions[0] == "Ly" &&
    exceptions[1] == "J4"
  })
  led.oblige("sporadic.monster_baby_index", A5_Classification, BigOrder, "|M| : |B| = 2⁵·3⁷·5³·7⁴·11·13²·29·41·59·71：以指數相加驗證 |B|·[M:B] = |M|（三者皆為 BigInt 精確值）", fn() {
    let arr = sporadics()
    let bfac = arr[sporadic_index("B")].fac
    let mfac = arr[sporadic_index("M")].fac
    let index : Array[(Int, Int)] = [(2, 5), (3, 7), (5, 3), (7, 4), (11, 1), (13, 2), (29, 1), (41, 1), (59, 1), (71, 1)]
    let prod = factorization_mul(bfac, index)
    let mut same = prod.length() == mfac.length()
    if same {
      let mut i = 0
      while i < mfac.length() {
        let (p1, e1) = prod[i]
        let (p2, e2) = mfac[i]
        if p1 != p2 || e1 != e2 {
          same = false
        }
        i = i + 1
      }
    }
    same &&
    max_prime_factor(mfac) == 71 &&
    max_prime_factor(bfac) == 47 &&
    mfac.length() == 15 &&
    bfac.length() == 11
  })
  led.oblige("sporadic.large_primes", A5_Classification, BigOrder, "大質數結構：26 群中 25 群之最大質因數 ≥ 11（例外恰為 J₂，最大質因數 7）；各群階皆為偶數", fn() {
    let arr = sporadics()
    let mut ge11 = 0
    let mut exception = ""
    let mut all_even = true
    for s in arr {
      let mp = max_prime_factor(s.fac)
      if mp >= 11 {
        ge11 = ge11 + 1
      } else {
        exception = s.name
      }
      let mut two = 0
      for pr in s.fac {
        let (p, _e) = pr
        if p == 2 {
          two = 1
        }
      }
      if two == 0 {
        all_even = false
      }
    }
    ge11 == 25 && exception == "J2" && all_even
  })
}
