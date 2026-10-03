//! 2-群作用（2-group action）與 2-局部結構：Sylow 定理、融合、中心列、Frattini、
//! 以及 2-群對其子群／陪集之共軛作用。
//!
//! 由 Sylow 定理（A4 公理之推論）導出之算式：
//!   * n_p ≡ 1 (mod p)、n_p | [G : P]，且全部 Sylow 子群共軛（融合為單一軌道）。
//!   * 冪零 2-群之中心列：γ₂(G) = [G,G]、Z(G)、Φ(G) 之階；超特殊 2-群（D₈、Q₈）滿足
//!     Z(G) = [G,G] = Φ(G)，且 G/Φ ≅ (Z/2)²（指數 2 子群恰 3 個）。
//!   * 2-群作用：正則作用遞移、對指數 2 子群之共軛作用平凡（因指數 2 子群必正規）。

/// 子群格中階為 ord 之子群。
pub fn subgroups_of_order(lat : Array[Array[Array[Int]]], ord : Int) -> Array[Array[Array[Int]]] {
  let out : Array[Array[Array[Int]]] = []
  for h in lat {
    if h.length() == ord {
      out.push(h)
    }
  }
  out
}

/// 子群是否正規：|N_G(H)| = |G|。
pub fn is_normal_subgroup(elems : Array[Array[Int]], n : Int, h : Array[Array[Int]]) -> Bool {
  normalizer_order(elems, n, h) == elems.length()
}

/// Frattini 子群之階（由子群格之極大子群交得）。
pub fn frattini_order(lat : Array[Array[Array[Int]]]) -> Int {
  frattini_subgroup(lat).length()
}

/// 兩子群（以元素表給定）是否相同。
pub fn same_subgroup(a : Array[Array[Int]], b : Array[Array[Int]]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let keys : Map[String, Unit] = Map([])
  for x in a {
    keys[perm_key(x)] = ()
  }
  for y in b {
    if !keys.contains(perm_key(y)) {
      return false
    }
  }
  true
}

/// 單一元素之共軛（輔助）。
pub fn conj_of(x : Array[Int], g : Array[Int]) -> Array[Int] {
  perm_conj(x, g)
}

/// 共軛融合軌道：直接以 G 之元素逐一共軛子群（不依賴任何未驗證假設）。
pub fn conjugation_orbits(
  elems : Array[Array[Int]],
  n : Int,
  subs : Array[Array[Array[Int]]]
) -> Array[Array[Int]] {
  let _ = n
  let out : Array[Array[Int]] = []
  let used = Array::make(subs.length(), false)
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
          let conj : Array[Array[Int]] = []
          for x in subs[cur] {
            conj.push(perm_conj(x, g))
          }
          let mut j = 0
          while j < subs.length() {
            if !used[j] && same_subgroup(subs[j], conj) {
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

/// n_p ≡ 1 (mod p) 且 n_p | [G:P] 之檢查（Sylow 定理之兩個同餘條件）。
pub fn sylow_congruence_ok(order : Int, psub_order : Int, p : Int, np : Int) -> Bool {
  np % p == mod_pos(1, p) &&
  (order / psub_order) % np == 0
}

/// 2-群作用之義務自證。
pub fn twolocal_register(led : Ledger) -> Unit {
  led.oblige("2lc.sylow_theorem", A1_Group, TwoLocal, "Sylow 定理：S₃ 之 Sylow 2-子群階 2、個數 3；S₄ 階 8、個數 3；S₅ 階 8、個數 15；A₅ 階 4、個數 5——皆滿足 n_p ≡ 1 (mod p) 且 n_p | [G:P]", fn() {
    let s3 = symmetric_group(3)
    let (e3, _o1) = s3.elements(1000)
    let s4 = symmetric_group(4)
    let (e4, _o2) = s4.elements(1000)
    let s5 = symmetric_group(5)
    let (e5, _o3) = s5.elements(1000)
    let a5 = alternating_group(5)
    let (ea5, _o4) = a5.elements(1000)
    let p3 = sylow_subgroup(e3, 3, s3.gens, 2)
    let p4 = sylow_subgroup(e4, 4, s4.gens, 2)
    let p5 = sylow_subgroup(e5, 5, s5.gens, 2)
    let p5a = sylow_subgroup(ea5, 5, a5.gens, 2)
    let n3 = sylow_count(e3, 3, p3)
    let n4 = sylow_count(e4, 4, p4)
    let n5 = sylow_count(e5, 5, p5)
    let n5a = sylow_count(ea5, 5, p5a)
    p3.length() == 2 && n3 == 3 &&
    p4.length() == 8 && n4 == 3 &&
    p5.length() == 8 && n5 == 15 &&
    p5a.length() == 4 && n5a == 5 &&
    sylow_congruence_ok(6, 2, 2, 3) &&
    sylow_congruence_ok(24, 8, 2, 3) &&
    sylow_congruence_ok(120, 8, 2, 15) &&
    sylow_congruence_ok(60, 4, 2, 5)
  })
  led.oblige("2lc.fusion_conjugate", A2_Action, TwoLocal, "融合：S₄ 之三個 Sylow 2-子群互相共軛（單一融合軌道，長度 3），|N_G(P)| = |G|/n₂ = 8 = |P|（P 自正規化）", fn() {
    let s4 = symmetric_group(4)
    let (e4, _o1) = s4.elements(1000)
    let lat = subgroup_lattice(e4, 4)
    let syl = subgroups_of_order(lat, 8)
    let psub = sylow_subgroup(e4, 4, s4.gens, 2)
    let orbits = conjugation_orbits(e4, 4, syl)
    syl.length() == 3 &&
    orbits.length() == 1 &&
    orbits[0].length() == 3 &&
    normalizer_order(e4, 4, psub) == 8 &&
    normalizer_order(e4, 4, psub) * sylow_count(e4, 4, psub) == 24
  })
  led.oblige("2lc.extraspecial", A1_Group, TwoLocal, "超特殊 2-群：D₈ 與 Q₈ 皆滿足 |Z(G)| = |[G,G]| = |Φ(G)| = 2，且指數 2 子群恰 3 個（G/Φ ≅ (Z/2)²）", fn() {
    let (d8g, _ord) = dihedral_group(4)
    let (ed8, _o1) = d8g.elements(100)
    let q8 = quaternion_group()
    let (eq8, _o2) = q8.elements(100)
    let latd = subgroup_lattice(ed8, 4)
    let latq = subgroup_lattice(eq8, 8)
    let idxd = subgroups_of_order(latd, 4)
    let idxq = subgroups_of_order(latq, 4)
    ed8.length() == 8 &&
    eq8.length() == 8 &&
    center_of(ed8, d8g.gens).length() == 2 &&
    center_of(eq8, q8.gens).length() == 2 &&
    derived_subgroup(4, d8g.gens).length() == 2 &&
    derived_subgroup(8, q8.gens).length() == 2 &&
    frattini_order(latd) == 2 &&
    frattini_order(latq) == 2 &&
    idxd.length() == 3 &&
    idxq.length() == 3
  })
  led.oblige("2lc.central_series", A1_Group, TwoLocal, "冪零性與中心列：D₈、Q₈ 之冪零類 = 2、上中心列長 = 2（皆為 2-群：|G| = 8 = 2³）；Q₈ 為 Hamiltonian（全部子群正規），D₈ 之反射子群不正規", fn() {
    let (d8g, _ord) = dihedral_group(4)
    let (ed8, _o1) = d8g.elements(100)
    let q8 = quaternion_group()
    let (eq8, _o2) = q8.elements(100)
    let latd = subgroup_lattice(ed8, 4)
    let latq = subgroup_lattice(eq8, 8)
    let mut all_normal_q = true
    for h in latq {
      for g in eq8 {
        let conj : Array[Array[Int]] = []
        for x in h {
          conj.push(perm_conj(x, g))
        }
        if !same_subgroup(h, conj) {
          all_normal_q = false
        }
      }
    }
    let mut non_normal_d = 0
    for h in latd {
      if !is_normal_subgroup(ed8, 4, h) {
        non_normal_d = non_normal_d + 1
      }
    }
    nilpotency_class(4, d8g.gens) == 2 &&
    nilpotency_class(8, q8.gens) == 2 &&
    upper_central_length(ed8, d8g.gens) == 2 &&
    upper_central_length(eq8, q8.gens) == 2 &&
    is_p_group(8, 2) &&
    all_normal_q &&
    latq.length() == 6 &&
    non_normal_d == 4
  })
  led.oblige("2lc.two_group_action", A2_Action, TwoLocal, "2-群作用：D₈ 作用於正方形 4 頂點遞移（軌道 4、穩定子階 2，8 = 4·2）；V₄ 之正則作用遞移於 4 點；指數 2 子群必正規故共軛作用平凡（3 個固定點）", fn() {
    let (d8g, _ord) = dihedral_group(4)
    let (ed8, _o1) = d8g.elements(100)
    let (v4, v4ord) = elementary_abelian_regular(2, 2)
    let _ = v4ord
    let (ev4, _o2) = v4.elements(100)
    let latd = subgroups_of_order(subgroup_lattice(ed8, 4), 4)
    let mut all_fixed = true
    for h in latd {
      let x = h[h.length() - 1]
      for g in ed8 {
        let conj : Array[Array[Int]] = []
        for y in h {
          conj.push(perm_conj(y, g))
        }
        if !same_subgroup(h, conj) || !group_member(4, h, conj_of(x, g)) {
          all_fixed = false
        }
      }
    }
    let (orbits, order) = burnside_from_action(d8g.gens, 4)
    orbits == 1 &&
    order == 1 &&
    d8g.orbits().length() == 1 &&
    v4.orbits().length() == 1 &&
    ev4.length() == 4 &&
    latd.length() == 3 &&
    all_fixed
  })
}
