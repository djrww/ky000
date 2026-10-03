//! 公理系統、定義、命題／定理、算式與義務自證（proof obligation）之登記與稽核。
//!
//! 本檔案建立整份形式化的「骨架」：
//!   * `Axiom`  —— 五條公理 A1–A5（群、作用、作用之複合、單性、有限單群分類）。
//!   * `Calculus` —— 算式種類（階、軌道—穩定化、Burnside、Molien、Weyl 維數、特徵標正交…）。
//!   * `Obligation` —— 義務自證：每條演算法必須附上可執行之布林驗證。
//!   * `Ledger`  —— 稽核帳本：執行全部義務並輸出總表。
//!
//! 設計原則：所有命題皆以「可執行算式」證明；任何一項失敗即記為未證，並在總表中列出。

/// 五條公理（本形式化之出發點）。
pub enum Axiom {
  /// A1 群公理：結合律、單位元、反元素。
  A1_Group
  /// A2 作用公理：e·x = x，(gh)·x = g·(h·x)。
  A2_Action
  /// A3 複合公理：作用沿同態與等變映射複合，仍為作用（作用之組合）。
  A3_Composition
  /// A4 單性公理：|G|>1 且僅有平凡正規子群。
  A4_Simplicity
  /// A5 分類公理（CFSG）：有限單群必同構於 Z/p、交錯群、李型群或 26 個散在單群之一。
  A5_Classification
}

/// 公理代號。
pub fn Axiom::tag(self : Axiom) -> String {
  match self {
    A1_Group => "A1"
    A2_Action => "A2"
    A3_Composition => "A3"
    A4_Simplicity => "A4"
    A5_Classification => "A5"
  }
}

/// 公理之數學陳述（算式化）。
pub fn Axiom::statement(self : Axiom) -> String {
  match self {
    A1_Group => "∀a,b,c∈G: (ab)c=a(bc)；∃e: ea=ae=a；∀a∃a⁻¹: aa⁻¹=e ｜ 群之簽章 (G,·,e,(-)⁻¹) 滿足三式"
    A2_Action => "G×X→X 滿足 e·x=x 與 (gh)·x=g·(h·x) ｜ 等價於同態 φ:G→Sym(X)，kerφ 為作用之核"
    A3_Composition => "α:G→Sym(X)、β:H→Sym(Y)、同態 f:G→H、等變 φ:X→Y ⟹ g↦β(f(g)) 為 G 在 Y 上之作用；乘積作用 g·(x,y)=(αg x, βg y) 亦為作用"
    A4_Simplicity => "G 有限且 |G|>1；若 N ⊴ G 則 N=1 或 N=G ｜ 等價：每一非單位元素之正規閉包皆為 G"
    A5_Classification => "G 有限單群 ⟹ G ≅ Z/p (p 質數) ∨ A_n (n≥5) ∨ 李型群 ∨ 26 個散在單群之一 ｜ 推論（階反查）：|G| 之質因數分解決定候選"
  }
}

/// 列出全部公理（供報表與文件遍歷）。
pub fn all_axioms() -> Array[Axiom] {
  [A1_Group, A2_Action, A3_Composition, A4_Simplicity, A5_Classification]
}

/// 算式種類（每一步證明都必須歸結到可計算之算式）。
pub enum Calculus {
  /// 結構常數與階：|G|、指數、階之質因數分解。
  Order
  /// 軌道—穩定化：|G| = |Orb(x)|·|Stab(x)|。
  OrbitStabilizer
  /// Burnside 引理：|X/G| = (1/|G|)Σ_g |Fix(g)|。
  Burnside
  /// 型（模式）與計數：矩陣／置換群之線性代數計數（如 |GL(n,q)|）。
  LinearCount
  /// Molien 級數與不變環：P(t) = (1/|G|)Σ_g 1/det(I−t·g)。
  Molien
  /// 導子與 Leibniz 律：[D,x]=1、D(fg)=Df·g+f·Dg。
  Derivation
  /// Weyl 維數公式：dim V_λ = Π_{α>0} ⟨λ+ρ,α⟩/⟨ρ,α⟩。
  WeylDimension
  /// 群概型：Hopf 公理 (Δ,ε,S) 與不變子空間之線性代數。
  HopfScheme
  /// 2-群局部結構：Sylow 定理、融合引理、中心列。
  TwoLocal
  /// 散在單群：大整數階之分解與整除關係。
  BigOrder
  /// 層與餘調：Čech 複形之核與像，Euler 特徵數。
  SheafCech
  /// 富化作用：Maschke／Reynolds 平均算子與特徵標內積。
  Enriched
}

/// 算式代號。
pub fn Calculus::tag(self : Calculus) -> String {
  match self {
    Order => "ORD"
    OrbitStabilizer => "ORB"
    Burnside => "BUR"
    LinearCount => "LIN"
    Molien => "MOL"
    Derivation => "DER"
    WeylDimension => "WEY"
    HopfScheme => "HOP"
    TwoLocal => "2LC"
    BigOrder => "BIG"
    SheafCech => "SHE"
    Enriched => "ENR"
  }
}

/// 義務自證：一條演算法／定理搭配一個可執行之驗證。
pub struct Obligation {
  id : Int
  algo : String
  axiom : String
  calculus : String
  statement : String
  run : () -> Bool
}

/// 稽核帳本：收集義務、執行、統計。
pub struct Ledger {
  items : Array[Obligation]
  mut counter : Int
}

/// 建立空帳本。
pub fn Ledger::new() -> Ledger {
  { items: [], counter: 0 }
}

/// 登記一條義務自證。`run` 回傳 true 表示該命題於本次執行中被驗證通過。
pub fn Ledger::oblige(
  self : Ledger,
  algo : String,
  ax : Axiom,
  calc : Calculus,
  statement : String,
  run : () -> Bool
) -> Unit {
  self.counter = self.counter + 1
  self.items.push({
    id: self.counter,
    algo,
    axiom: ax.tag(),
    calculus: calc.tag(),
    statement,
    run,
  })
}

/// 執行全部義務，輸出逐條結果與總表。回傳 (通過數, 總數)。
pub fn Ledger::audit(self : Ledger, verbose : Bool) -> (Int, Int) {
  let mut passed = 0
  let failed : Array[String] = []
  for it in self.items {
    let ok = (it.run)()
    if ok {
      passed = passed + 1
    } else {
      failed.push("[\{it.id}] \{it.algo} — \{it.statement}")
    }
    if verbose {
      let mark = if ok { "✓" } else { "✗" }
      println("\{mark} [\{it.id}] \{it.algo} (\{it.axiom}/\{it.calculus})")
    }
  }
  println("")
  println("──────── 義務自證總表 ────────")
  println("通過：\{passed} / \{self.items.length()}")
  if failed.length() > 0 {
    println("未通過：")
    for f in failed {
      println("  ✗ \{f}")
    }
  } else {
    println("全部義務皆由可執行算式驗證通過。")
  }
  (passed, self.items.length())
}

/// 依演算法名稱聚合，回傳不重複之演算法數（用於「演算法必須多於 48 條」之計數稽核）。
pub fn Ledger::algorithm_count(self : Ledger) -> Int {
  let seen : Map[String, Unit] = Map([])
  for it in self.items {
    seen[it.algo] = ()
  }
  let mut c = 0
  seen.each(fn(_k, _v) { c = c + 1 })
  c
}

/// 依算式種類統計義務數。
pub fn Ledger::calculus_histogram(self : Ledger) -> Array[(String, Int)] {
  let h : Map[String, Int] = Map([])
  for it in self.items {
    h[it.calculus] = h.get_or_default(it.calculus, 0) + 1
  }
  let out : Array[(String, Int)] = []
  h.each(fn(k, v) { out.push((k, v)) })
  out
}
