//! 置換、有限群引擎（Schreier–Sims 基底與強生成集）。
//!
//! 定義 D1（置換）：p : {0,…,n−1} → {0,…,n−1} 為雙射，以 `Array[Int]` 表示，p[i] = p(i)。
//! 定義 D2（群作用）：G ≤ Sym(X) 以生成集 gens 給定；g·x = g[x]。
//! 定義 D3（基底與強生成集 BSGS）：base = (b₀,…,b_{k−1})，SGS 為一組對每層 i 生成
//!   G⁽ⁱ⁾ = Stab_G(b₀,…,b_{i−1}) 之元素；|G| = Π_i |Orb_{G⁽ⁱ⁾}(bᵢ)|  （軌道—穩定化連乘）。
//! 演算法採用 Schreier 引理（Schreier generators）與 sifting 反覆閉包。

/// 單位置換。
pub fn perm_id(n : Int) -> Array[Int] {
  let out : Array[Int] = []
  let mut i = 0
  while i < n {
    out.push(i)
    i = i + 1
  }
  out
}

/// 複合 (a∘b)(x) = a(b(x))。
pub fn perm_mul(a : Array[Int], b : Array[Int]) -> Array[Int] {
  let n = a.length()
  let out : Array[Int] = []
  let mut i = 0
  while i < n {
    out.push(a[b[i]])
    i = i + 1
  }
  out
}

/// 反置換。
pub fn perm_inv(a : Array[Int]) -> Array[Int] {
  let n = a.length()
  let out = perm_id(n)
  let mut i = 0
  while i < n {
    out[a[i]] = i
    i = i + 1
  }
  out
}

/// 共軛 g a g⁻¹。
pub fn perm_conj(a : Array[Int], g : Array[Int]) -> Array[Int] {
  perm_mul(g, perm_mul(a, perm_inv(g)))
}

/// 冪 a^k（k ≥ 0）。
pub fn perm_pow(a : Array[Int], k : Int) -> Array[Int] {
  let n = a.length()
  let mut out = perm_id(n)
  let mut i = 0
  while i < k {
    out = perm_mul(a, out)
    i = i + 1
  }
  out
}

/// 兩置換是否相等。
pub fn perm_eq(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let mut i = 0
  while i < a.length() {
    if a[i] != b[i] {
      return false
    }
    i = i + 1
  }
  true
}

/// 是否為單位置換。
pub fn perm_is_id(a : Array[Int]) -> Bool {
  let mut i = 0
  while i < a.length() {
    if a[i] != i {
      return false
    }
    i = i + 1
  }
  true
}

/// 是否固定 base[0..k−1]。
pub fn perm_fixes_prefix(a : Array[Int], base : Array[Int], k : Int) -> Bool {
  let mut i = 0
  while i < k {
    if a[base[i]] != base[i] {
      return false
    }
    i = i + 1
  }
  true
}

/// 由輪換表示法建立置換（輪換內元素為 0-based）。
pub fn perm_from_cycles(n : Int, cycles : Array[Array[Int]]) -> Array[Int] {
  let out = perm_id(n)
  for c in cycles {
    let m = c.length()
    if m > 0 {
      let mut i = 0
      while i < m {
        out[c[i]] = c[mod_pos(i + 1, m)]
        i = i + 1
      }
    }
  }
  out
}

/// 對換 (i j)。
pub fn transposition(n : Int, i : Int, j : Int) -> Array[Int] {
  let out = perm_id(n)
  out[i] = j
  out[j] = i
  out
}

/// 輪換 (c₀ c₁ … c_{m−1})。
pub fn cycle_perm(n : Int, c : Array[Int]) -> Array[Int] {
  perm_from_cycles(n, [c])
}

/// 輪換分解（長度 ≥ 2 之輪換，依最小元素排序）。
pub fn perm_cycles(a : Array[Int]) -> Array[Array[Int]] {
  let n = a.length()
  let seen : Array[Bool] = Array::make(n, false)
  let out : Array[Array[Int]] = []
  let mut i = 0
  while i < n {
    if !seen[i] {
      let cyc : Array[Int] = []
      let mut x = i
      while !seen[x] {
        seen[x] = true
        cyc.push(x)
        x = a[x]
      }
      if cyc.length() > 1 {
        out.push(cyc)
      }
    }
    i = i + 1
  }
  out
}

/// 置換階：輪換長度之最小公倍數。
pub fn perm_order(a : Array[Int]) -> Int {
  let mut ord = 1
  for c in perm_cycles(a) {
    ord = ord * c.length() / gcd_int(ord, c.length())
  }
  ord
}

/// 符號 sign(a) = (−1)^{n−#cycles}，回傳 ±1。
pub fn perm_sign(a : Array[Int]) -> Int {
  let n = a.length()
  let mut cycles = 0
  let seen : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if !seen[i] {
      cycles = cycles + 1
      let mut x = i
      while !seen[x] {
        seen[x] = true
        x = a[x]
      }
    }
    i = i + 1
  }
  if (n - cycles) % 2 == 0 { 1 } else { -1 }
}

/// 型（輪換型）：回傳各輪換長度之遞減列表（含固定點長度 1）。
pub fn perm_cycle_type(a : Array[Int]) -> Array[Int] {
  let n = a.length()
  let seen : Array[Bool] = Array::make(n, false)
  let lens : Array[Int] = []
  let mut i = 0
  while i < n {
    if !seen[i] {
      let mut len = 0
      let mut x = i
      while !seen[x] {
        seen[x] = true
        len = len + 1
        x = a[x]
      }
      lens.push(len)
    }
    i = i + 1
  }
  // 遞減插入排序
  let mut i2 = 1
  while i2 < lens.length() {
    let key = lens[i2]
    let mut j = i2 - 1
    while j >= 0 && lens[j] < key {
      lens[j + 1] = lens[j]
      j = j - 1
    }
    lens[j + 1] = key
    i2 = i2 + 1
  }
  lens
}

/// 固定點個數（作用於全集合 {0,…,n−1}）。
pub fn perm_fixed_count(a : Array[Int]) -> Int {
  let mut c = 0
  let mut i = 0
  while i < a.length() {
    if a[i] == i {
      c = c + 1
    }
    i = i + 1
  }
  c
}

/// 標準鍵（供集合去重）。
pub fn perm_key(a : Array[Int]) -> String {
  let mut s = ""
  let mut i = 0
  while i < a.length() {
    s = s + a[i].to_string()
    s = s + ","
    i = i + 1
  }
  s
}

// ───────────── 生成群：軌道、Schreier 引理、BSGS ─────────────

/// 以生成集表示之有限置換群。
pub struct PermGroup {
  n : Int
  gens : Array[Array[Int]]
}

/// 建立群（濾除單位置換與重複生成元）。
pub fn perm_group(n : Int, gens : Array[Array[Int]]) -> PermGroup {
  let out : Array[Array[Int]] = []
  let seen : Map[String, Unit] = Map([])
  for g in gens {
    if !perm_is_id(g) {
      let k = perm_key(g)
      if !seen.contains(k) {
        seen[k] = ()
        out.push(g)
      }
    }
  }
  { n, gens: out }
}

/// 軌道 Orb(x)（BFS）。
pub fn PermGroup::orbit(self : PermGroup, x : Int) -> Array[Int] {
  let seen : Map[Int, Unit] = Map([])
  let out : Array[Int] = []
  out.push(x)
  seen[x] = ()
  let mut head = 0
  while head < out.length() {
    let y = out[head]
    head = head + 1
    for s in self.gens {
      let z = s[y]
      if !seen.contains(z) {
        seen[z] = ()
        out.push(z)
      }
    }
  }
  out
}

/// 軌道分解（作用於 {0,…,n−1}）。
pub fn PermGroup::orbits(self : PermGroup) -> Array[Array[Int]] {
  let seen : Map[Int, Unit] = Map([])
  let out : Array[Array[Int]] = []
  let mut x = 0
  while x < self.n {
    if !seen.contains(x) {
      let o = self.orbit(x)
      for y in o {
        seen[y] = ()
      }
      out.push(o)
    }
    x = x + 1
  }
  out
}

/// 是否遞移。
pub fn PermGroup::is_transitive(self : PermGroup) -> Bool {
  self.orbit(0).length() == self.n
}

/// 傳遞給定集合上之軌道（作用限制）。
pub fn PermGroup::orbit_on(self : PermGroup, pts : Array[Int]) -> Array[Int] {
  if pts.length() == 0 {
    return []
  }
  let in_set : Map[Int, Unit] = Map([])
  for p in pts {
    in_set[p] = ()
  }
  let seen : Map[Int, Unit] = Map([])
  let out : Array[Int] = []
  out.push(pts[0])
  seen[pts[0]] = ()
  let mut head = 0
  while head < out.length() {
    let y = out[head]
    head = head + 1
    for s in self.gens {
      let z = s[y]
      if in_set.contains(z) && !seen.contains(z) {
        seen[z] = ()
        out.push(z)
      }
    }
  }
  out
}

/// Schreier 引理：穩定化子 Stab(x) 之生成元。
/// 設 T 為軌道 Orb(x) 之橫截面（u_y 滿足 u_y(x)=y），
/// 則 Stab(x) = ⟨ u_{s·y}⁻¹ · s · u_y : y ∈ Orb(x), s ∈ gens ⟩。
pub fn PermGroup::stabilizer_gens(self : PermGroup, x : Int) -> Array[Array[Int]] {
  let orb = self.orbit(x)
  let pos : Map[Int, Int] = Map([])
  let mut idx = 0
  for y in orb {
    pos[y] = idx
    idx = idx + 1
  }
  // 橫截面：trans[i] 把 x 送到 orb[i]
  let trans : Array[Array[Int]] = []
  let id = perm_id(self.n)
  let mut i = 0
  while i < orb.length() {
    trans.push(id)
    i = i + 1
  }
  // BFS 建橫截面
  let done : Map[Int, Unit] = Map([])
  done[x] = ()
  let mut head = 0
  let queue : Array[Int] = [x]
  while head < queue.length() {
    let y = queue[head]
    head = head + 1
    let uy = trans[pos[y]]
    for s in self.gens {
      let z = s[y]
      if !done.contains(z) {
        done[z] = ()
        queue.push(z)
        trans[pos[z]] = perm_mul(s, uy)
      }
    }
  }
  let out : Array[Array[Int]] = []
  for y in orb {
    let uy = trans[pos[y]]
    for s in self.gens {
      let sy = s[y]
      let usy = trans[pos[sy]]
      let w = perm_mul(perm_inv(usy), perm_mul(s, uy))
      if !perm_is_id(w) {
        out.push(w)
      }
    }
  }
  out
}

/// 由生成集顯式枚舉全部元素（BFS）。回傳 (元素, 是否未溢出)。
/// `cap` 為元素上限，超過則回傳空陣列與 false。
pub fn PermGroup::elements(self : PermGroup, cap : Int) -> (Array[Array[Int]], Bool) {
  let seen : Map[String, Unit] = Map([])
  let out : Array[Array[Int]] = []
  let id = perm_id(self.n)
  out.push(id)
  seen[perm_key(id)] = ()
  let mut head = 0
  while head < out.length() {
    if out.length() > cap {
      return ([], false)
    }
    let g = out[head]
    head = head + 1
    for s in self.gens {
      let h = perm_mul(s, g)
      let k = perm_key(h)
      if !seen.contains(k) {
        seen[k] = ()
        out.push(h)
      }
    }
  }
  (out, true)
}

/// BSGS 鏈（基底、強生成集、各層軌道與橫截面）。
pub struct Chain {
  n : Int
  base : Array[Int]
  sgs : Array[Array[Int]]
  orbits : Array[Array[Int]]
  trans : Array[Array[Array[Int]]]
  pos : Array[Map[Int, Int]]
}

/// 給定基底之前綴後，屬於該層穩定化子之強生成元。
fn chain_level_gens(ch : Chain, i : Int) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  for s in ch.sgs {
    if perm_fixes_prefix(s, ch.base, i) {
      out.push(s)
    }
  }
  out
}

/// 重算第 i 層之軌道與橫截面（第 i 層生成元＝固定 base[0..i−1] 之強生成元）。
fn chain_rebuild_level(ch : Chain, i : Int) -> Unit {
  let gens = chain_level_gens(ch, i)
  let b = ch.base[i]
  let orbit : Array[Int] = [b]
  let idx : Map[Int, Int] = Map([])
  idx[b] = 0
  let trans : Array[Array[Int]] = [perm_id(ch.n)]
  let mut head = 0
  while head < orbit.length() {
    let y = orbit[head]
    let uy = trans[head]
    head = head + 1
    for s in gens {
      let z = s[y]
      if !idx.contains(z) {
        idx[z] = orbit.length()
        orbit.push(z)
        trans.push(perm_mul(s, uy))
      }
    }
  }
  ch.orbits[i] = orbit
  ch.trans[i] = trans
  ch.pos[i] = idx
}

/// 重算 from 之後各層之軌道資料。
fn chain_rebuild(ch : Chain, from : Int) -> Unit {
  let mut i = from
  while i < ch.base.length() {
    chain_rebuild_level(ch, i)
    i = i + 1
  }
}

/// 把基底延長到足以容納新元素（取第一個被移動之點）。
fn chain_extend_base(ch : Chain, g : Array[Int]) -> Unit {
  let mut i = 0
  while i < ch.n {
    if g[i] != i {
      // 若該點已在基底中則不重複加入
      let mut found = false
      for b in ch.base {
        if b == i {
          found = true
        }
      }
      if !found {
        ch.base.push(i)
        ch.orbits.push([])
        ch.trans.push([])
        ch.pos.push(Map([]))
        return
      }
    }
    i = i + 1
  }
}

/// Sifting：把 g 逐層約化。回傳 (是否屬於群, 殘餘元素)。
pub fn chain_sift(ch : Chain, g : Array[Int]) -> (Bool, Array[Int]) {
  let mut cur = g
  let mut i = 0
  while i < ch.base.length() {
    let y = cur[ch.base[i]]
    match ch.pos[i].get(y) {
      Some(j) => {
        let u = ch.trans[i][j]
        cur = perm_mul(perm_inv(u), cur)
      }
      None => return (false, cur)
    }
    i = i + 1
  }
  (perm_is_id(cur), cur)
}

/// Schreier–Sims：建立基底與強生成集，滿足 |G| = Π 軌道長度。
pub fn schreier_sims(n : Int, gens : Array[Array[Int]]) -> Chain {
  let ch : Chain = {
    n,
    base: [],
    sgs: [],
    orbits: [],
    trans: [],
    pos: [],
  }
  let g0 = perm_group(n, gens)
  for g in g0.gens {
    ch.sgs.push(g)
    if ch.base.length() == 0 {
      chain_extend_base(ch, g)
    }
  }
  if ch.sgs.length() == 0 {
    return ch
  }
  let mut round = 0
  let mut added = true
  while added {
    added = false
    round = round + 1
    chain_rebuild(ch, 0)
    let mut i = 0
    while i < ch.base.length() {
      let lvl = chain_level_gens(ch, i)
      let orb = ch.orbits[i]
      let mut j = 0
      while j < orb.length() {
        let x = orb[j]
        let ux = ch.trans[i][j]
        for s in lvl {
          let sx = s[x]
          match ch.pos[i].get(sx) {
            Some(j2) => {
              let usx = ch.trans[i][j2]
              let w = perm_mul(perm_inv(usx), perm_mul(s, ux))
              if !perm_is_id(w) {
                let (mem, r) = chain_sift(ch, w)
                if !mem {
                  ch.sgs.push(r)
                  chain_extend_base(ch, r)
                  chain_rebuild(ch, i)
                  added = true
                }
              }
            }
            None => ()
          }
        }
        j = j + 1
      }
      i = i + 1
    }
    // 安全閥：避免病態輸入下無限迴圈
    if round > 4 * n + 16 {
      break
    }
  }
  ch
}

/// 群階 |G| = Π_i |Orb_{G⁽ⁱ⁾}(bᵢ)|。
pub fn chain_order(ch : Chain) -> Int {
  let mut ord = 1
  for o in ch.orbits {
    ord = ord * o.length()
  }
  ord
}

/// 計算生成群之階。
pub fn group_order(n : Int, gens : Array[Array[Int]]) -> Int {
  chain_order(schreier_sims(n, gens))
}

/// 成員判定。
pub fn group_member(n : Int, gens : Array[Array[Int]], g : Array[Int]) -> Bool {
  let ch = schreier_sims(n, gens)
  let (mem, _r) = chain_sift(ch, g)
  mem
}

/// 指數：群的元素階之最小公倍數（需要顯式元素表）。
pub fn exponent_of_elements(elems : Array[Array[Int]]) -> Int {
  let mut e = 1
  for g in elems {
    let o = perm_order(g)
    e = e * o / gcd_int(e, o)
  }
  e
}

// ───────────────────────── 義務自證 ─────────────────────────

/// 註冊置換與 BSGS 引擎之義務自證。
pub fn perm_register(led : Ledger) -> Unit {
  led.oblige("perm.basic", A1_Group, Order, "置換代數：結合律、反元素、階為輪換長度之 lcm（對 S₅ 全元素抽樣）", fn() {
    let n = 5
    let gens : Array[Array[Int]] = [transposition(n, 0, 1), cycle_perm(n, [0, 1, 2, 3, 4])]
    let g = perm_group(n, gens)
    let (elems, ok) = g.elements(200)
    if !ok {
      return false
    }
    let mut all = true
    for a in elems {
      if !perm_eq(perm_mul(a, perm_inv(a)), perm_id(n)) {
        all = false
      }
      let o = perm_order(a)
      // a^ord = e 且若 ord>1 則 a^(ord−1) ≠ e
      if !perm_is_id(perm_pow(a, o)) {
        all = false
      }
      if o > 1 && perm_is_id(perm_pow(a, o - 1)) {
        all = false
      }
      for b in elems {
        for c in elems {
          let l = perm_mul(perm_mul(a, b), c)
          let r = perm_mul(a, perm_mul(b, c))
          if !perm_eq(l, r) {
            all = false
          }
        }
      }
    }
    all
  })
  led.oblige("perm.schreier_sims", A2_Action, Order, "BSGS 與顯式枚舉互證：|S₄| = 24、|A₅| = 60、|GL(2,3) 之置換像| 等（兩演算法一致）", fn() {
    let s4 : Array[Array[Int]] = [transposition(4, 0, 1), transposition(4, 1, 2), transposition(4, 2, 3)]
    let a5 : Array[Array[Int]] = [
      cycle_perm(5, [0, 1, 2]),
      cycle_perm(5, [0, 1, 2, 3, 4]),
    ]
    let g4 = perm_group(4, s4)
    let g5 = perm_group(5, a5)
    let (e4, ok4) = g4.elements(100000)
    let (e5, ok5) = g5.elements(100000)
    if !ok4 || !ok5 {
      return false
    }
    e4.length() == 24 &&
    e5.length() == 60 &&
    group_order(4, s4) == 24 &&
    group_order(5, a5) == 60 &&
    chain_order(schreier_sims(4, s4)) == e4.length() &&
    chain_order(schreier_sims(5, a5)) == e5.length()
  })
  led.oblige("perm.stabilizer", A2_Action, OrbitStabilizer, "軌道—穩定化：|G| = |Orb(x)|·|Stab(x)|（S₅ 作用於 5 點；A₅ 作用於 5 點）", fn() {
    let s5 : Array[Array[Int]] = [transposition(5, 0, 1), cycle_perm(5, [0, 1, 2, 3, 4])]
    let g = perm_group(5, s5)
    if group_order(5, s5) != 120 {
      return false
    }
    let mut ok = true
    let mut x = 0
    while x < 5 {
      let orb = g.orbit(x)
      let st = g.stabilizer_gens(x)
      let nst = group_order(5, st)
      if orb.length() * nst != 120 {
        ok = false
      }
      x = x + 1
    }
    ok
  })
  led.oblige("perm.member", A2_Action, Order, "成員判定：偶置換屬於 A₅、對換不屬於 A₅（sifting 正確性）", fn() {
    let a5 : Array[Array[Int]] = [cycle_perm(5, [0, 1, 2]), cycle_perm(5, [0, 1, 2, 3, 4])]
    group_member(5, a5, cycle_perm(5, [0, 1, 3])) &&
    group_member(5, a5, perm_id(5)) &&
    !group_member(5, a5, transposition(5, 0, 1)) &&
    !group_member(5, a5, transposition(5, 2, 4))
  })
  led.oblige("perm.sign", A1_Group, Order, "符號為同態：sign(ab)=sign(a)sign(b)，sign(對換)=−1，sign(3-輪換)=+1", fn() {
    let n = 5
    let g = perm_group(n, [transposition(n, 0, 1), cycle_perm(n, [0, 1, 2, 3, 4])])
    let (elems, ok) = g.elements(200)
    if !ok {
      return false
    }
    let mut all = true
    for a in elems {
      for b in elems {
        if perm_sign(perm_mul(a, b)) != perm_sign(a) * perm_sign(b) {
          all = false
        }
      }
    }
    all && perm_sign(transposition(n, 0, 1)) == -1 && perm_sign(cycle_perm(n, [0, 1, 2])) == 1
  })
}
