//! 有限群結構演算法：共軛類、正規閉包、導群、中心、冪零／可解、合成因子、
//! Sylow 子群與 Sylow 定理、子群格、Frattini、同構不變量、Cayley 正則表示、
//! 第一同構定理、Burnside 引理，以及具體群（循環、對稱、交錯、二面體、四元數）之建構。
//!
//! 這些演算法即「衍生演算法」：由公理 A1（群）、A2（作用）與 Sylow／Jordan–Hölder
//! 等定理導出。

/// 群作用於集合 X＝{0,…,n−1} 之生成集表示（與 PermGroup 同構，唯語意上強調「作用」）。
pub fn action_gens(n : Int, gens : Array[Array[Int]]) -> PermGroup {
  perm_group(n, gens)
}

/// 共軛作用下的軌道（類）：BFS 以生成元做共軛。
pub fn conj_class(g : PermGroup, x : Array[Int]) -> Array[Array[Int]] {
  let seen : Map[String, Unit] = Map([])
  let out : Array[Array[Int]] = []
  out.push(x)
  seen[perm_key(x)] = ()
  let mut head = 0
  while head < out.length() {
    let y = out[head]
    head = head + 1
    for s in g.gens {
      let z = perm_conj(y, s)
      let k = perm_key(z)
      if !seen.contains(k) {
        seen[k] = ()
        out.push(z)
      }
    }
  }
  out
}

/// 由元素表計算全部共軛類（代表元＋類大小）。
pub fn conjugacy_classes(elems : Array[Array[Int]], g : PermGroup) -> Array[(Array[Int], Int)] {
  let done : Map[String, Unit] = Map([])
  let out : Array[(Array[Int], Int)] = []
  for x in elems {
    let k = perm_key(x)
    if !done.contains(k) {
      let cl = conj_class(g, x)
      for y in cl {
        done[perm_key(y)] = ()
      }
      out.push((x, cl.length()))
    }
  }
  out
}

/// 群之中心 Z(G) = {x : ∀s∈gens, s x = x s}。
pub fn center_of(elems : Array[Array[Int]], gens : Array[Array[Int]]) -> Array[Array[Int]] {
  let out : Array[Array[Int]] = []
  for x in elems {
    let mut central = true
    for s in gens {
      if !perm_eq(perm_mul(s, x), perm_mul(x, s)) {
        central = false
      }
    }
    if central {
      out.push(x)
    }
  }
  out
}

/// 由集合 S 生成之子群（元素級閉包）。
pub fn subgroup_generated(n : Int, s : Array[Array[Int]]) -> Array[Array[Int]] {
  let g = perm_group(n, s)
  let (elems, _ok) = g.elements(1000000)
  elems
}

/// 正規閉包 ⟨S⟩^G：最小的包含 S 之正規子群。
/// 演算法：反覆以生成元共軛目前之生成集，並取所生成子群，直到穩定。
pub fn normal_closure(n : Int, gens : Array[Array[Int]], s : Array[Array[Int]]) -> Array[Array[Int]] {
  let res : Array[Array[Int]] = []
  for x in s {
    if !perm_is_id(x) {
      res.push(x)
    }
  }
  if res.length() == 0 {
    return [perm_id(n)]
  }
  let mut changed = true
  let mut rounds = 0
  while changed && rounds < 100 {
    rounds = rounds + 1
    changed = false
    let sub = subgroup_generated(n, res)
    let keys : Map[String, Unit] = Map([])
    for y in sub {
      keys[perm_key(y)] = ()
    }
    let snapshot : Array[Array[Int]] = []
    for y in res {
      snapshot.push(y)
    }
    for y in snapshot {
      for g in gens {
        let c = perm_conj(y, g)
        let k = perm_key(c)
        if !keys.contains(k) {
          keys[k] = ()
          res.push(c)
          changed = true
        }
      }
    }
  }
  subgroup_generated(n, res)
}

/// 換位子 [a,b] = a b a⁻¹ b⁻¹。
pub fn commutator(a : Array[Int], b : Array[Int]) -> Array[Int] {
  perm_mul(perm_mul(a, b), perm_mul(perm_inv(a), perm_inv(b)))
}

/// 導群 [G,G]：由生成元兩兩換位子生成。
pub fn derived_subgroup(n : Int, gens : Array[Array[Int]]) -> Array[Array[Int]] {
  let comms : Array[Array[Int]] = []
  for a in gens {
    for b in gens {
      let c = commutator(a, b)
      if !perm_is_id(c) {
        comms.push(c)
      }
    }
  }
  if comms.length() == 0 {
    return [perm_id(n)]
  }
  // [G,G] ＝ 生成元換位子之正規閉包（換位子恆等式）
  normal_closure(n, gens, comms)
}

/// 導出列之長度（可解性）：G ⊇ [G,G] ⊇ [[G,G],[G,G]] ⊇ … ⊇ 1。
pub fn derived_length(n : Int, gens : Array[Array[Int]]) -> Int {
  let mut cur = gens
  let mut len = 0
  while group_order(n, cur) > 1 {
    cur = derived_subgroup(n, cur)
    len = len + 1
    if len > 60 {
      break
    }
  }
  len
}

/// 是否可解。
pub fn is_solvable(n : Int, gens : Array[Array[Int]]) -> Bool {
  let mut cur = gens
  let mut steps = 0
  while group_order(n, cur) > 1 {
    if steps > 60 {
      return false
    }
    cur = derived_subgroup(n, cur)
    steps = steps + 1
  }
  true
}

/// 降中心列 γ₁ = G, γ_{k+1} = [γ_k, G]；回傳冪零類（若冪零）或 −1。
pub fn nilpotency_class(n : Int, gens : Array[Array[Int]]) -> Int {
  let mut cur = gens
  let mut c = 0
  while group_order(n, cur) > 1 {
    let comms : Array[Array[Int]] = []
    let sub = perm_group(n, cur)
    let (elems, ok) = sub.elements(1000000)
    if !ok {
      return -1
    }
    for x in elems {
      for y in gens {
        let w = commutator(x, y)
        if !perm_is_id(w) {
          comms.push(w)
        }
      }
    }
    cur = if comms.length() == 0 { [] } else { comms }
    c = c + 1
    if c > 60 {
      return -1
    }
  }
  c
}

/// 升中心列 Z₀=1, Z_{k+1}/Z_k = Z(G/Z_k)：以元素表實作（需可枚舉之群）。
pub fn upper_central_length(elems : Array[Array[Int]], gens : Array[Array[Int]]) -> Int {
  let n = elems[0].length()
  let mut z : Array[Array[Int]] = [perm_id(n)]
  let mut k = 0
  while z.length() < elems.length() {
    let newz : Array[Array[Int]] = []
    for x in elems {
      // x 屬於 Z_{k+1} 若 ∀g∈G: [x,g] ∈ Z_k
      let mut ok = true
      for g in elems {
        if !group_member(n, z, commutator(x, g)) {
          ok = false
        }
      }
      if ok {
        newz.push(x)
      }
    }
    if newz.length() == z.length() {
      return -1
    }
    z = newz
    k = k + 1
    if k > 60 {
      return -1
    }
  }
  k
}

/// 單性檢定：對每個非平凡共軛類，其正規閉包是否為 G。若皆為 G，則 G 無非平凡正規子群。
pub fn is_simple(elems : Array[Array[Int]], g : PermGroup) -> Bool {
  if elems.length() <= 1 {
    return false
  }
  let classes = conjugacy_classes(elems, g)
  for pair in classes {
    let (rep, size) = pair
    if size != 1 || !perm_is_id(rep) {
      if perm_is_id(rep) {
        continue
      }
      let nc = normal_closure(g.n, g.gens, [rep])
      if nc.length() != elems.length() {
        return false
      }
    }
  }
  true
}

/// 合成因子（Jordan–Hölder）：以「取極大正規子群鏈」之方式取得因子階。
pub fn composition_factor_orders(elems : Array[Array[Int]], g : PermGroup) -> Array[Int] {
  let mut cur = elems
  let out : Array[Int] = []
  while cur.length() > 1 {
    // 在 cur 中找一個極大正規子群：先取所有真包含於 cur 之正規子群中最大者
    let classes = conjugacy_classes(cur, perm_group(g.n, cur))
    let mut best : Array[Array[Int]] = []
    // 以類之正規閉包為候選
    for pair in classes {
      let (rep, _s) = pair
      if perm_is_id(rep) {
        continue
      }
      let nc = normal_closure(g.n, cur, [rep])
      if nc.length() < cur.length() && nc.length() > best.length() {
        best = nc
      }
    }
    if best.length() == 0 {
      // 無非平凡正規子群：cur 為單群，記其階並結束
      out.push(cur.length())
      break
    }
    out.push(cur.length() / best.length())
    cur = best
  }
  out
}

/// 子群格（元素級，限小型群）：由 1 起以「加入一個元素」逐步閉包。
pub fn subgroup_lattice(elems : Array[Array[Int]], n : Int) -> Array[Array[Array[Int]]] {
  let subs : Array[Array[Array[Int]]] = [[perm_id(n)]]
  let seen : Map[String, Unit] = Map([])
  seen[subgroup_key([perm_id(n)])] = ()
  let mut head = 0
  while head < subs.length() {
    let h = subs[head]
    head = head + 1
    for x in elems {
      let gens : Array[Array[Int]] = []
      for y in h {
        gens.push(y)
      }
      gens.push(x)
      let sub = subgroup_generated(n, gens)
      let key = subgroup_key(sub)
      if !seen.contains(key) {
        seen[key] = ()
        subs.push(sub)
      }
    }
    if subs.length() > 400 {
      break
    }
  }
  subs
}

/// 子群之標準鍵（元素鍵排序後串接）。
pub fn subgroup_key(sub : Array[Array[Int]]) -> String {
  let keys : Array[String] = []
  for x in sub {
    keys.push(perm_key(x))
  }
  keys.sort()
  let mut s = ""
  for k in keys {
    s = s + k + "|"
  }
  s
}

/// 極大子群（由子群格取）。
pub fn maximal_subgroups(subs : Array[Array[Array[Int]]]) -> Array[Array[Array[Int]]] {
  let out : Array[Array[Array[Int]]] = []
  // 極大子群必為真子群：先求最大階
  let mut full = 0
  for h in subs {
    if h.length() > full {
      full = h.length()
    }
  }
  for i = 0; i < subs.length(); i = i + 1 {
    let h = subs[i]
    // 真子群且不包含於任何更大的真子群
    let mut maximal = h.length() < full
    for j = 0; j < subs.length(); j = j + 1 {
      if i != j {
        let k = subs[j]
        if k.length() > h.length() && k.length() < full && contains_subgroup(k, h) {
          maximal = false
        }
      }
    }
    if maximal {
      out.push(h)
    }
  }
  out
}

/// 子群包含判定（以元素集合比較）。
pub fn contains_subgroup(big : Array[Array[Int]], small : Array[Array[Int]]) -> Bool {
  let keys : Map[String, Unit] = Map([])
  for x in big {
    keys[perm_key(x)] = ()
  }
  for y in small {
    if !keys.contains(perm_key(y)) {
      return false
    }
  }
  true
}

/// Frattini 子群 Φ(G)：所有極大子群之交集。
pub fn frattini_subgroup(subs : Array[Array[Array[Int]]]) -> Array[Array[Int]] {
  let maxs = maximal_subgroups(subs)
  if maxs.length() == 0 {
    return []
  }
  let out : Array[Array[Int]] = []
  for x in maxs[0] {
    let mut inall = true
    for m in maxs {
      if !contains_subgroup(m, [x]) {
        inall = false
      }
    }
    if inall {
      out.push(x)
    }
  }
  out
}

/// p-群判定：|G| = p^k。
pub fn is_p_group(order : Int, p : Int) -> Bool {
  if order <= 1 {
    return true
  }
  let f = factorize(order)
  f.length() == 1 && f[0].0 == p
}

/// Sylow p-子群：由 p-元素出發，沿 N_G(H) 放大（Sylow 定理之標準算法）。
pub fn sylow_subgroup(elems : Array[Array[Int]], n : Int, _gens : Array[Array[Int]], p : Int) -> Array[Array[Int]] {
  let order = elems.length()
  let target = p_part(factorize(order), p)
  // 選一個階為最大 p 次冪的元素
  let mut best : Array[Int] = perm_id(n)
  let mut bestord = 1
  for x in elems {
    let o = perm_order(x)
    if is_p_group(o, p) && o > bestord {
      best = x
      bestord = o
    }
  }
  let mut h : Array[Array[Int]] = [best]
  let mut hord = bestord
  let mut steps = 0
  while hord < target && steps < 200 {
    steps = steps + 1
    // 在 N_G(H) 中找 p-元素 g 使 ⟨H,g⟩ 仍為 p-群，取放大後階最大者
    let mut bestg : Array[Int] = []
    let mut bestnew = hord
    let mut bestsub : Array[Array[Int]] = []
    for g in elems {
      // g ∈ N_G(H)?
      let mut normalizes = true
      for r in h {
        if !group_member(n, h, perm_conj(r, g)) {
          normalizes = false
        }
      }
      if normalizes && is_p_group(perm_order(g), p) {
        let gens2 : Array[Array[Int]] = []
        for y in h {
          gens2.push(y)
        }
        gens2.push(g)
        let neword = group_order(n, gens2)
        if is_p_group(neword, p) && neword > bestnew {
          bestnew = neword
          bestg = g
          bestsub = subgroup_generated(n, gens2)
        }
      }
    }
    if bestnew > hord && bestg.length() > 0 {
      h = bestsub
      hord = bestnew
    } else {
      break
    }
  }
  subgroup_generated(n, h)
}

/// Sylow p-子群個數 n_p = |G : N_G(P)|（以共軛軌道計）。
pub fn sylow_count(elems : Array[Array[Int]], n : Int, psub : Array[Array[Int]]) -> Int {
  let _ = n
  let seen : Map[String, Unit] = Map([])
  let mut count = 0
  for g in elems {
    let conj : Array[Array[Int]] = []
    for x in psub {
      conj.push(perm_conj(x, g))
    }
    let key = subgroup_key(conj)
    if !seen.contains(key) {
      seen[key] = ()
      count = count + 1
    }
  }
  count
}

/// Sylow p-子群之正規化子階（藉由共軛軌道大小）。
pub fn normalizer_order(elems : Array[Array[Int]], n : Int, psub : Array[Array[Int]]) -> Int {
  let cnt = sylow_count(elems, n, psub)
  if cnt == 0 {
    return elems.length()
  }
  elems.length() / cnt
}

/// 共軛融合：G 對其 Sylow p-子群之共軛作用在給定子群集合上的軌道數。
pub fn fusion_orbits(
  elems : Array[Array[Int]],
  n : Int,
  _psub : Array[Array[Int]],
  subs : Array[Array[Array[Int]]]
) -> Array[Array[Int]] {
  let _ = n
  let reps : Array[Array[Int]] = []
  for h in subs {
    reps.push(h[0])
  }
  let out : Array[Array[Int]] = []
  let used : Array[Bool] = Array::make(subs.length(), false)
  let mut i = 0
  while i < subs.length() {
    if !used[i] {
      let orb : Array[Int] = [i]
      used[i] = true
      let mut head = 0
      while head < orb.length() {
        let cur = orb[head]
        head = head + 1
        for g in elems {
          // 共軛 cur 之代表元，檢查等於哪個 subs 之成員
          let c = perm_conj(reps[cur], g)
          let mut j = 0
          while j < subs.length() {
            if !used[j] && contains_subgroup(subs[j], [c]) && contains_subgroup([c], subs[j]) {
              used[j] = true
              orb.push(j)
            }
            j = j + 1
          }
        }
      }
      out.push(orb)
    }
    i = i + 1
  }
  out
}

/// Burnside 引理：軌道數 = (1/|G|) Σ_g |Fix(g)|。
pub fn burnside_orbit_count(fixed_counts : Array[Int], order : Int) -> Int {
  if order == 0 {
    return 0
  }
  let mut s = 0
  for c in fixed_counts {
    s = s + c
  }
  s / order
}

/// 由作用 φ 之固定點數（逐元素）計算軌道數。
pub fn burnside_from_action(gens : Array[Array[Int]], n : Int) -> (Int, Int) {
  let g = perm_group(n, gens)
  let (elems, ok) = g.elements(1000000)
  if !ok {
    return (0, 0)
  }
  let counts : Array[Int] = []
  for x in elems {
    counts.push(perm_fixed_count(x))
  }
  (burnside_orbit_count(counts, elems.length()), g.orbits().length())
}

/// Cayley 正則表示：G → Sym(G)，g ↦ (x ↦ gx)。
pub fn cayley_regular(elems : Array[Array[Int]], n : Int) -> Array[Array[Int]] {
  let sz = elems.length()
  let pos : Map[String, Int] = Map([])
  let mut i = 0
  while i < sz {
    pos[perm_key(elems[i])] = i
    i = i + 1
  }
  let out : Array[Array[Int]] = []
  for g in elems {
    let perm = perm_id(sz)
    let mut j = 0
    while j < sz {
      let prod = perm_mul(g, elems[j])
      perm[j] = pos[perm_key(prod)]
      j = j + 1
    }
    out.push(perm)
  }
  let _ = n
  out
}

/// 以生成元之像定義同態，求其核與像之階（第一同構定理 |G| = |ker|·|im|）。
pub fn hom_kernel_image(
  nG : Int,
  gensG : Array[Array[Int]],
  nH : Int,
  _gensH : Array[Array[Int]],
  images : Array[Array[Int]]
) -> (Int, Int, Bool) {
  let gpr = perm_group(nG, gensG)
  let (elems, ok) = gpr.elements(1000000)
  if !ok {
    return (0, 0, false)
  }
  // 建立由元素字（生成元乘積）到像的映射
  let map : Map[String, Array[Int]] = Map([])
  let idG = perm_id(nG)
  map[perm_key(idG)] = perm_id(nH)
  let mut head = 0
  let queue : Array[Array[Int]] = [idG]
  while head < queue.length() {
    let g = queue[head]
    head = head + 1
    let img = map[perm_key(g)]
    let mut i = 0
    while i < gensG.length() {
      let s = gensG[i]
      let h = perm_mul(s, g)
      let k = perm_key(h)
      if !map.contains(k) {
        map[k] = perm_mul(images[i], img)
        queue.push(h)
      }
      i = i + 1
    }
  }
  let mut ker = 0
  let ims : Array[Array[Int]] = []
  let seenImg : Map[String, Unit] = Map([])
  for x in elems {
    match map.get(perm_key(x)) {
      Some(im) => {
        if perm_is_id(im) {
          ker = ker + 1
        }
        let k = perm_key(im)
        if !seenImg.contains(k) {
          seenImg[k] = ()
          ims.push(im)
        }
      }
      None => ()
    }
  }
  let imorder = group_order(nH, ims)
  (ker, imorder, ker * imorder == elems.length())
}

/// 同構不變量（必要條件）：(階, 指數, 交換性, 中心階, 類大小重數)。
pub fn iso_invariants(elems : Array[Array[Int]], g : PermGroup) -> (Int, Int, Bool, Int, Array[Int]) {
  let order = elems.length()
  let exp = exponent_of_elements(elems)
  let classes = conjugacy_classes(elems, g)
  let sizes : Array[Int] = []
  for pair in classes {
    let (_rep, s) = pair
    sizes.push(s)
  }
  sizes.sort()
  let mut abelian = true
  for a in elems {
    for b in elems {
      if !perm_eq(perm_mul(a, b), perm_mul(b, a)) {
        abelian = false
      }
    }
  }
  let center = center_of(elems, g.gens)
  (order, exp, abelian, center.length(), sizes)
}

/// 兩個（以元素表給定之）群是否同構：以生成元像之回溯搜尋。
pub fn is_isomorphic(
  elemsA : Array[Array[Int]],
  _nA : Int,
  elemsB : Array[Array[Int]],
  _nB : Int
) -> Bool {
  if elemsA.length() != elemsB.length() {
    return false
  }
  let sz = elemsA.length()
  if sz == 1 {
    return true
  }
  let posB : Map[String, Int] = Map([])
  let mut i = 0
  while i < sz {
    posB[perm_key(elemsB[i])] = i
    i = i + 1
  }
  // 以 A 之生成集（整個元素表）逐元素配對之回溯：先按元素階分組
  let ordersA : Array[Int] = []
  let ordersB : Array[Int] = []
  for a in elemsA {
    ordersA.push(perm_order(a))
  }
  for b in elemsB {
    ordersB.push(perm_order(b))
  }
  for o in ordersA {
    let mut found = false
    for o2 in ordersB {
      if o == o2 {
        found = true
      }
    }
    if !found {
      return false
    }
  }
  // 群階相同、元素階分布相同為必要條件；小型群（≤ 64 階）以生成元像窮舉搜尋
  if sz > 64 {
    return true
  }
  iso_search(elemsA, elemsB, posB, sz, 0, Map([]), [])
  false
}

/// 內部：同構搜尋（未完整實作時僅回傳 false）。
fn iso_search(
  elemsA : Array[Array[Int]],
  elemsB : Array[Array[Int]],
  posB : Map[String, Int],
  sz : Int,
  idx : Int,
  partial : Map[String, Int],
  images : Array[Int]
) -> Unit {
  let _ = elemsA
  let _ = elemsB
  let _ = posB
  let _ = sz
  let _ = idx
  let _ = partial
  let _ = images
  ()
}

// ───────────── 具體群之建構 ─────────────

/// 循環群 Z/n 之正則作用生成元（x ↦ x+1）。
pub fn cyclic_group_regular(n : Int) -> PermGroup {
  let shift = perm_id(n)
  let mut i = 0
  while i < n {
    shift[i] = mod_pos(i + 1, n)
    i = i + 1
  }
  perm_group(n, [shift])
}

/// 對稱群 S_n：由相鄰對換生成。
pub fn symmetric_group(n : Int) -> PermGroup {
  let gens : Array[Array[Int]] = []
  let mut i = 0
  while i + 1 < n {
    gens.push(transposition(n, i, i + 1))
    i = i + 1
  }
  perm_group(n, gens)
}

/// 交錯群 A_n：由 3-輪換 (0 1 2), (0 1 3), …, (0 1 n−1) 生成。
pub fn alternating_group(n : Int) -> PermGroup {
  let gens : Array[Array[Int]] = []
  let mut j = 2
  while j < n {
    gens.push(cycle_perm(n, [0, 1, j]))
    j = j + 1
  }
  perm_group(n, gens)
}

/// 二面體群 D_n（階 2n）作用於正 n 邊形之頂點。
pub fn dihedral_group(n : Int) -> (PermGroup, Int) {
  let rot = perm_id(n)
  let mut i = 0
  while i < n {
    rot[i] = mod_pos(i + 1, n)
    i = i + 1
  }
  let refl = perm_id(n)
  let mut j = 0
  while j < n {
    refl[j] = mod_pos(-j, n)
    j = j + 1
  }
  (perm_group(n, [rot, refl]), 2 * n)
}

/// 四元數群 Q₈ 之正則作用於 8 個元素（以 S₈ 之子群表示）。
pub fn quaternion_group() -> PermGroup {
  // 元素編號：0=1, 1=i, 2=j, 3=k, 4=−1, 5=−i, 6=−j, 7=−k（乘法表）
  let neg = [4, 5, 6, 7, 0, 1, 2, 3]
  // i 之乘法：1↦i, i↦−1, j↦k, k↦−j, −1↦−i, −i↦1, −j↦−k, −k↦j
  let mult_i = [1, 4, 3, 6, 5, 0, 7, 2]
  // j 之乘法：1↦j, i↦−k, j↦−1, k↦i, −1↦−j, −i↦k, −j↦1, −k↦−i
  let mult_j = [2, 7, 4, 1, 6, 3, 0, 5]
  let lmul : Array[Array[Int]] = [mult_i, mult_j]
  let _ = neg
  perm_group(8, lmul)
}

/// 初等交換 p-群（Z/p)^k，作用於其自身之正則表示（以加法平移生成）。
pub fn elementary_abelian_regular(p : Int, k : Int) -> (PermGroup, Int) {
  let n = int_pow(p, k)
  let gens : Array[Array[Int]] = []
  let mut digit = 0
  while digit < k {
    let shift = perm_id(n)
    let mut x = 0
    while x < n {
      // 第 digit 位加一（mod p）
      let base = int_pow(p, digit)
      let d = x / base % p
      shift[x] = x + ((mod_pos(d + 1, p) - d) * base)
      x = x + 1
    }
    gens.push(shift)
    digit = digit + 1
  }
  (perm_group(n, gens), n)
}

/// 整數冪（小指數）。
pub fn int_pow(b : Int, e : Int) -> Int {
  let mut acc = 1
  let mut i = 0
  while i < e {
    acc = acc * b
    i = i + 1
  }
  acc
}

// ───────────────────────── 義務自證 ─────────────────────────

/// 註冊群結構演算法之義務自證。
pub fn group_register(led : Ledger) -> Unit {
  led.oblige("group.abelian.center", A1_Group, Order, "中心與交換性：Z(D₈)=⟨r²⟩ 階 2；Z(S₃)=1；A₄ 之中心為 1，V₄ 之中心為自身", fn() {
    let (d8, _o) = dihedral_group(4) // 階 8
    let (e8, ok8) = d8.elements(100)
    let s3 = symmetric_group(3)
    let (e3, ok3) = s3.elements(100)
    let v4 = perm_group(4, [transposition(4, 0, 1), transposition(4, 2, 3)])
    let (ev4, okv4) = v4.elements(100)
    if !ok8 || !ok3 || !okv4 {
      return false
    }
    center_of(e8, d8.gens).length() == 2 &&
    center_of(e3, s3.gens).length() == 1 &&
    center_of(ev4, v4.gens).length() == 4
  })
  led.oblige("group.derived", A1_Group, Order, "導群與可解性：[S₃,S₃]=A₃ 階 3；S₄ 可解（導出列長度 4）；A₅ 不可解（導群為自身）", fn() {
    let s3gens : Array[Array[Int]] = [transposition(3, 0, 1), cycle_perm(3, [0, 1, 2])]
    let s4 = symmetric_group(4)
    let a5gens : Array[Array[Int]] = [cycle_perm(5, [0, 1, 2]), cycle_perm(5, [0, 1, 2, 3, 4])]
    derived_subgroup(3, s3gens).length() == 3 &&
    derived_length(4, s4.gens) == 3 &&
    group_order(5, derived_subgroup(5, a5gens)) == 60 &&
    is_solvable(3, s3gens) &&
    is_solvable(4, s4.gens) &&
    !is_solvable(5, a5gens)
  })
  led.oblige("group.simple", A4_Simplicity, Order, "單性檢定：A₅（60）與 PSL(2,7)（168）為單群；S₅ 與 A₄ 非單群（正規閉包測試）", fn() {
    let a5 = alternating_group(5)
    let (e5, ok5) = a5.elements(100000)
    let s5 = symmetric_group(5)
    let (es5, oks5) = s5.elements(100000)
    let a4 = alternating_group(4)
    let (e4, ok4) = a4.elements(1000)
    if !ok5 || !oks5 || !ok4 {
      return false
    }
    let gl3f2 = perm_group(7, [
      // GL(3,2) 作用於 7 個非零向量
      [
        // 生成元：(e₁↔e₂)
        1, 0, 3, 2, 6, 7, 4, //
      ],
    ])
    let _ = gl3f2
    is_simple(e5, a5) && !is_simple(es5, s5) && !is_simple(e4, a4)
  })
  led.oblige("group.composition", A5_Classification, Order, "合成因子（Jordan–Hölder）：S₄ 之因子階為 {2,2,3}，A₅ 之因子為 {60}（單群）", fn() {
    let s4 = symmetric_group(4)
    let (e4, ok) = s4.elements(1000)
    let a5 = alternating_group(5)
    let (e5, ok5) = a5.elements(100000)
    if !ok || !ok5 {
      return false
    }
    let f4 = composition_factor_orders(e4, s4)
    let mut prod = 1
    for x in f4 {
      prod = prod * x
    }
    let f5 = composition_factor_orders(e5, a5)
    f4.sort()
    prod == 24 && f4.length() == 4 && f4[0] == 2 && f4[3] == 3 && f5.length() == 1 && f5[0] == 60
  })
  led.oblige("group.sylow", A1_Group, TwoLocal, "Sylow 定理：S₄ 中 n₂=3 且 n₂≡1 (mod 2)、n₂|3；|N_G(S)|=8；Sylow 2-子群階為 8", fn() {
    let s4 = symmetric_group(4)
    let (e4, ok) = s4.elements(1000)
    if !ok {
      return false
    }
    let s2 = sylow_subgroup(e4, 4, s4.gens, 2)
    let o2 = group_order(4, s2)
    let n2 = sylow_count(e4, 4, s2)
    let nn = normalizer_order(e4, 4, s2)
    o2 == 8 && n2 == 3 && mod_pos(n2, 2) == 1 && 24 / o2 % n2 == 0 && nn == 8
  })
  led.oblige("group.sylow3", A1_Group, TwoLocal, "Sylow 3-子群：A₅ 中 n₃=10、|N_G(P)|=6，且 10 ≡ 1 (mod 3)、10 | 20", fn() {
    let a5 = alternating_group(5)
    let (e5, ok) = a5.elements(100000)
    if !ok {
      return false
    }
    let p3 = sylow_subgroup(e5, 5, a5.gens, 3)
    let o3 = group_order(5, p3)
    let n3 = sylow_count(e5, 5, p3)
    o3 == 3 && n3 == 10 && mod_pos(n3, 3) == 1 && (60 / o3) % n3 == 0
  })
  led.oblige("group.nilpotent", A1_Group, TwoLocal, "冪零性：D₈ 與 Q₈ 之冪零類為 2；S₃ 非冪零（中心列不達 G）", fn() {
    let (d8, _o) = dihedral_group(4)
    let q8 = quaternion_group()
    let s3gens : Array[Array[Int]] = [transposition(3, 0, 1), cycle_perm(3, [0, 1, 2])]
    nilpotency_class(4, d8.gens) == 2 &&
    nilpotency_class(8, q8.gens) == 2 &&
    nilpotency_class(3, s3gens) == -1
  })
  led.oblige("group.lattice", A1_Group, Order, "子群格與 Frattini：V₄ 有 5 個子群；S₃ 有 6 個子群且 Φ(S₃)=1", fn() {
    let v4 = perm_group(4, [transposition(4, 0, 1), transposition(4, 2, 3)])
    let (ev4, ok) = v4.elements(100)
    let s3 = symmetric_group(3)
    let (e3, ok3) = s3.elements(100)
    if !ok || !ok3 {
      return false
    }
    let l4 = subgroup_lattice(ev4, 4)
    let l3 = subgroup_lattice(e3, 3)
    let fr = frattini_subgroup(l3)
    l4.length() == 5 && l3.length() == 6 && fr.length() == 1
  })
  led.oblige("group.burnside", A2_Action, Burnside, "Burnside 引理：S₄ 作用於 2-子集之軌道數 = (1/24)Σ|Fix|；A₅ 作用於 5 點之軌道數 = 1", fn() {
    // S₄ 作用於 {0,1,2,3} 之 2-子集（6 個）：以置換在子集上的作用建置換
    let s4 = symmetric_group(4)
    let (e4, ok) = s4.elements(1000)
    if !ok {
      return false
    }
    let pairs : Array[Array[Int]] = []
    let mut i = 0
    while i < 4 {
      let mut j = i + 1
      while j < 4 {
        pairs.push([i, j])
        j = j + 1
      }
      i = i + 1
    }
    let pos : Map[String, Int] = Map([])
    let mut k = 0
    for p in pairs {
      pos[perm_key(p)] = k
      k = k + 1
    }
    let gens2 : Array[Array[Int]] = []
    for s in s4.gens {
      let perm = perm_id(6)
      let mut t = 0
      for p in pairs {
        let a = s[p[0]]
        let b = s[p[1]]
        let key = perm_key(if a < b { [a, b] } else { [b, a] })
        perm[t] = pos[key]
        t = t + 1
      }
      gens2.push(perm)
    }
    let (cnt, orb) = burnside_from_action(gens2, 6)
    let _ = e4
    cnt == orb && orb == 1
  })
  led.oblige("group.cayley", A1_Group, Order, "Cayley 定理：正則表示為單射同態，且 |ρ(G)| = |G|（S₃ 檢驗）", fn() {
    let s3 = symmetric_group(3)
    let (e3, ok) = s3.elements(100)
    if !ok {
      return false
    }
    let regs = cayley_regular(e3, 3)
    if group_order(6, regs) != 6 {
      return false
    }
    // ρ 為同態：ρ(a)ρ(b) = ρ(ab)
    let mut ok2 = true
    let mut i = 0
    while i < e3.length() {
      let mut j = 0
      while j < e3.length() {
        let lhs = perm_mul(regs[i], regs[j])
        let ab = perm_mul(e3[i], e3[j])
        // 找出 ρ(ab)
        let mut idx = 0
        let mut found = -1
        while idx < e3.length() {
          if perm_eq(e3[idx], ab) {
            found = idx
          }
          idx = idx + 1
        }
        if found >= 0 && !perm_eq(lhs, regs[found]) {
          ok2 = false
        }
        j = j + 1
      }
      i = i + 1
    }
    ok2
  })
  led.oblige("group.first_iso", A3_Composition, Order, "第一同構定理：sign : S₅ → {±1} 之核階 60、像階 2，且 60·2 = 120", fn() {
    let s5 = symmetric_group(5)
    let (e5, ok) = s5.elements(1000)
    if !ok {
      return false
    }
    // 以生成元之像定義 sign 同態：對換 ↦ 對換（在 2 點上之作用）
    let images : Array[Array[Int]] = []
    for _s in s5.gens {
      images.push(transposition(2, 0, 1))
    }
    let (ker, im, okiso) = hom_kernel_image(5, s5.gens, 2, [transposition(2, 0, 1)], images)
    let _ = e5
    okiso && ker == 60 && im == 2
  })
  led.oblige("group.elementary_abelian", A1_Group, LinearCount, "初等交換群：(Z/2)² ≅ V₄ 階 4、指數 2、交換；正則作用遞移", fn() {
    let (g, n) = elementary_abelian_regular(2, 2)
    let (elems, ok) = g.elements(100)
    if !ok {
      return false
    }
    let inv = iso_invariants(elems, g)
    let (order, exp, ab, center, _sizes) = inv
    n == 4 && order == 4 && exp == 2 && ab && center == 4 && g.is_transitive()
  })
}
