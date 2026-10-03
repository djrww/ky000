# ky000 ｜ 有限單群 × 作用：以 MoonBit 寫成之公理與義務自證

> **命題即算式，算式即證明。**
> 本專案以五條公理 (A1–A5) 為骨架，將「有限單群」(Z/p、交錯群、李型群、26 個散在單群)
> 與「作用」(矩陣、多項式環、微分算子、層、李群、群概型、2-群、富化……) 組合起來，
> 全部以 MoonBit 實作，並用**可執行之義務自證 (proof obligation)** 逐條稽核。
> 無第 3 方依賴形式化：每一條命題都由本專案自己的算式算出來。

```
[建置成功 9872 ms] djrww/ky000/lib -> djrww/ky000/main
══════════════════════════════════════════════════════════════
 有限單群 × 作用：公理 A1–A5 與義務自證稽核（MoonBit）
══════════════════════════════════════════════════════════════

──────── 義務自證總表 ────────
通過：82 / 82
全部義務皆由可執行算式驗證通過。

演算法條數（去重）=82
算式分布：
  ORD : 24    LIN : 13    BIG : 9     2LC : 8     SHE : 5
  ORB : 4     MOL : 4     HOP : 4     ENR : 4     BUR : 3
  WEY : 2     DER : 2
```

* **82 / 82 條義務**全部通過（最後一次實跑輸出，見下節完整引錄）。
* **82 條去重演算法 > 48**（需求門檻），由 `Ledger::algorithm_count()` 以雜湊去重統計。
* 12 種算式 (calculus) 家族，對應 5 條公理。

---

## 1. 五條公理（Axiom）

| 代號 | 公理 | 陳述（`Axiom::statement()` 原文） |
|---|---|---|
| **A1** | 群公理 | ∀a,b,c∈G: (ab)c=a(bc)；∃e: ea=ae=a；∀a∃a⁻¹: aa⁻¹=e ⟹ 群之簽章 (G,·,e,(-)⁻¹) 滿足三式 |
| **A2** | 作用公理 | G×X→X 滿足 e·x=x 與 (gh)·x=g·(h·x) ⟺ 同態 φ:G→Sym(X)，kerφ 為作用之核 |
| **A3** | 複合公理 | α:G→Sym(X)、β:H→Sym(Y)、同態 f:G→H、等變 φ:X→Y ⟹ g↦β(f(g)) 為 G 在 Y 上之作用；乘積作用 g·(x,y)=(αg x, βg y) 亦為作用 |
| **A4** | 單性公理 | G 有限且 \|G\|>1；若 N ⊴ G 則 N=1 或 N=G ⟺ 每一非單位元素之正規閉包皆為 G |
| **A5** | 分類公理 (CFSG) | G 有限單群 ⟹ G ≅ Z/p ∨ A_n (n≥5) ∨ 李型群 ∨ 26 個散在單群之一 ⟹ 推論：\|G\| 之質因數分解決定候選 |

12 種算式家族（`Calculus`）：`ORD` 階 / `ORB` 軌道—穩定化 / `BUR` Burnside /
`LIN` 線性代數計數 / `MOL` Molien 級數 / `DER` 導子與 Leibniz / `WEY` Weyl 維數 /
`HOP` 群概型 (Hopf) / `2LC` 2-局部 / `BIG` 大整數階 / `SHE` 層與餘調 / `ENR` 富化 (Reynolds、特徵標)。

---

## 2. 作用域對照表（13 個領域 × 模組 × 代表義務）

| # | 作用域 | 模組 | 條數 | 代表義務（實跑通過之具體結論） |
|---|---|---|---|---|
| 1 | 數字論引擎（群階之基礎算式） | `num.mbt` | 10 | 質數表、質因數分解與大整數、`pow_mod`、Euler φ、Möbius、二項係數 |
| 2 | 置換作用（Schreier–Sims） | `perm.mbt` | 5 | BSGS 鏈、置換群階、軌道—穩定化、置換之階與符號 |
| 3 | 群結構作用（共軛、Sylow、子群格） | `group.mbt` | 12 | 正規閉包、導群／導長、合成因子、Sylow 子群與個數、Frattini、Cayley 正則 |
| 4 | 矩陣作用於向量 | `matrix.mbt` | 9 | \|GL(2,2)\|=6、\|GL(2,3)\|=48、\|GL(2,5)\|=480；GL(2,5) 在 F₅²∖{0} 之軌道 = 24 = p²−1；\|PSL(2,5)\|=60 作用於 P¹(F₅) |
| 5 | 多項式環作用 | `poly.mbt` | 9 | 多項式算術環、不可約性、Frobenius 軌道、平移不變元 = F_p[x^p−x]、導子不變元 = ker D = F_p[x^p]、Z/n 對 F_p[x]/(xⁿ−1) 之 Burnside 軌道數、Molien 級數（S₃ 反射表示 = 1/((1−t²)(1−t³))） |
| 6 | 微分算子作用 | `poly.mbt` `lie.mbt` | 2 | d/dx 之 Leibniz 律與核；Jacobi 恆等式；Killing 形式 B(h,h)=8、B(e,f)=4、非退化 |
| 7 | 循環單群 Z/p | `num.mbt` `group.mbt` `hopf.mbt` | — | 質數篩、Z/p 正則作用遞移、μ_p 群概型之不變量 = ⌈N/p⌉ |
| 8 | 交錯群 A_n | `alt.mbt` | 5 | \|A₅\|=60=\|PSL(2,5)\|=\|PSL(2,4)\|；A₅ 之 5 個共軛類 {1,15,20,12,12}；A₅ 單（每個非單位元素之正規閉包 = 60）；A₄ 不單（V₄、導長 2）；3-輪換生成 |
| 9 | 李型單群 | `lie.mbt` | 7 | PSL(2,q) 階 q(q²−1)/d；PSL(3,q) 階（q=2,3,4,5 → 168, 5616, 20160, 372000）；²B₂(8)=29120=2⁶·5·7·13；²G₂(3)=1512=3·\|PSL(2,8)\|；³D₄(2)=211341312=2¹²·3⁴·7²·13；G₂(2)=12096；Weyl 維數 dim V_{a,b}=(a+1)(b+1)(a+b+2)/2 與根系乘積一致；Cartan 行列式 A₂=3、B₂=2、G₂=1 |
| 10 | 26 個散在單群（例外） | `sporadic.mbt` | 6 | 恰 26 群、名稱互異；26 個階之「十進位 × 質因數指數表」以 BigInt 精確互相驗證（含 \|M\| = 808017424794512875886459904961710757005754368000000000）；階 < 2³¹ 之 10 群另以 Int 乘法複核；**24/26 之階整除 \|M\|，例外恰為 J₄ 與 Ly**；\|M\|:\|B\| = 2⁵·3⁷·5³·7⁴·11·13²·29·41·59·71 |
| 11 | 層作用 | `sheaf.mbt` | 5 | Čech 複形 δ¹∘δ⁰ = 0（圓周／球面／環面，F₅ 與 F₇）；餘調 (h⁰,h¹,h²)：圓周 (1,1,0)、球面 (1,0,1)、環面 (1,2,1) 且 χ 相符；等變層 δ⁰∘ρ⁰ = ρ¹∘δ⁰；**Lefschetz–Hopf 跡公式**（圓周旋轉 L=0、球面三階旋轉 L=2）；黏合／下降：f₀₁+f₁₂=f₀₂ ⟺ 上邊界 |
| 12 | 李群作用（離散化） | `lie.mbt` | 7 | exp(X)exp(Y)=exp(X+Y+[X,Y]/2)（Heisenberg，F₇）；**Ad(exp X) = exp(ad X)**；exp(X)exp(−X)=I |
| 13 | 群概型作用 | `hopf.mbt` | 4 | 群代數 Hopf 公理（Δ(g)=g⊗g、ε、S、Δ 乘法性、S 反乘法性）；共作用之不變量 = 餘不變量 = 軌道數；μ_p 不變量 = ⌈N/p⌉；F₃[x]/(x³−1) 非半單（(x−1)³=0、冪等元僅 0 與 1） |
| 14 | 2-群作用 | `twolocal.mbt` | 5 | Sylow：S₄ 階 8／個數 3、S₅ 階 8／個數 15、A₅ 階 4／個數 5，皆滿足 n_p ≡ 1 (mod p)；S₄ 之三個 Sylow 共軛融合為單一軌道、\|N(P)\|=8；D₈、Q₈ 超特殊：\|Z\|=\|[G,G]\|=\|Φ\|=2、指數 2 子群恰 3 個；Q₈ 全部子群正規 |
| 15 | 富化作用 | `enrich.mbt` | 5 | Reynolds 算子 π²=π、rank π = dim V^G；Maschke 互補構造（F₅）；S₃ 特徵標兩兩正交、Σd²=6；A₅ 次數 1,3,3,4,5 之 Σd²=60；Frobenius 互反 χ_perm=χ_triv+χ_std；不變維數 = Molien 係數 = 分拆數 (1,1,2,3,4,5,7) |

> 條數合計：10+5+12+9+9+5+7+6+5+4+5+5 = **82** 條義務，**82** 條去重演算法名稱。

---

## 3. 實跑稽核（原文引錄）

```
公理：
  A1 — ∀a,b,c∈G: (ab)c=a(bc)；∃e: ea=ae=a；∀a∃a⁻¹: aa⁻¹=e ｜ 群之簽章 (G,·,e,(-)⁻¹) 滿足三式
  A2 — G×X→X 滿足 e·x=x 與 (gh)·x=g·(h·x) ｜ 等價於同態 φ:G→Sym(X)，kerφ 為作用之核
  A3 — α:G→Sym(X)、β:H→Sym(Y)、同態 f:G→H、等變 φ:X→Y ⟹ g↦β(f(g)) 為 G 在 Y 上之作用
  A4 — G 有限且 |G|>1；若 N ⊴ G 則 N=1 或 N=G ｜ 等價：每一非單位元素之正規閉包皆為 G
  A5 — G 有限單群 ⟹ G ≅ Z/p ∨ A_n (n≥5) ∨ 李型群 ∨ 26 個散在單群之一

義務自證逐條：
  ✓ [1] num.is_prime (A1/ORD)          ✓ [2] num.factorize (A1/ORD)
  ✓ [3] num.pow_mod (A1/ORD)           ✓ [4] num.ext_gcd (A1/ORD)
  …
  ✓ [80] enrich.char_orthogonality (A1/ENR)
  ✓ [81] enrich.frobenius (A3/ENR)
  ✓ [82] enrich.invariants_molien (A2/MOL)

──────── 義務自證總表 ────────
通過：82 / 82
全部義務皆由可執行算式驗證通過。
通過 82/82
演算法條數（去重）=82
```

完整 82 條清單（含每條之命題陳述）見 **[`docs/obligations.md`](docs/obligations.md)**。

---

## 4. 建置與執行

```bash
# 一次性：安裝固定版本之 moonc-worker 與 playground payload（需網路）
node dev/setup.mjs

# 建置 + 執行稽核（預設 JS 後端）
node dev/build.mjs

# 只建置、不執行
node dev/build.mjs --build

# 換後端（wasm-gc / wasm / native 皆可，需 payload 支援）
MBT_TARGET=wasm-gc node dev/build.mjs
```

`dev/build.mjs` 直接讀取 `moon.pkg.json` 之 `import` 關係建立套件圖，逐套件呼叫
`buildPackage`，最後 `linkCore` 產出 `dev/.out/ky000.mjs`；核心套件為
`@moonbit/moonc-worker@0.1.202608195`（與 payload `@marianoguerra/tutuca-playground-payload@0.23.0`
嚴格配對——版本不符會出現誤導性的 `E4018` 錯誤）。

官方工具鏈路徑亦已備齊：`moon.mod.json`（模組 `djrww/ky000`，源碼目錄 `src`）與
`src/lib/moon.pkg.json`、`src/main/moon.pkg.json`（`is-main` + 依賴 lib）。

---

## 5. 檔案地圖

```
src/lib/axiom.mbt      公理 A1–A5、算式家族、Obligation、Ledger（oblige / audit / algorithm_count / histogram）
src/lib/num.mbt        數字論引擎（階、質因數分解、BigInt 階乘、模運算）       10 條
src/lib/perm.mbt       置換群與 Schreier–Sims BSGS                               5 條
src/lib/group.mbt      群結構：共軛、Sylow、子群格、Frattini、合成因子           12 條
src/lib/matrix.mbt     GL/SL/PSL、向量與射影作用、特徵多項式、核與秩             9 條
src/lib/poly.mbt       多項式環、Frobenius、平移／導子不變元、Molien、有理數級數  9 條
src/lib/alt.mbt        交錯群 A_n：階、共軛類、單性、A₅ ≅ PSL(2,5)                5 條
src/lib/lie.mbt        李型群階公式、Weyl 維數、根系／Cartan、Jacobi/Killing、exp/BCH/Ad  7 條
src/lib/sporadic.mbt   26 個散在單群：BigInt 雙表示、Monster 整除結構             6 條
src/lib/sheaf.mbt      層作用：Čech 複形、餘調、等變、Lefschetz、下降             5 條
src/lib/hopf.mbt       群概型作用：Hopf 公理、共作用不變量、μ_p 非半單            4 條
src/lib/twolocal.mbt   2-群作用：Sylow、融合、超特殊 2-群、中心列                5 條
src/lib/enrich.mbt     富化作用：Reynolds、Maschke、特徵標正交、Frobenius         5 條
src/main/main.mbt      入口：登記全部義務 → 稽核 → 總表
dev/build.mjs          建置驅動（套件圖、buildPackage、linkCore）
dev/setup.mjs          固定版本依賴安裝與 payload 佈署
docs/obligations.md    82 條義務自證總表（自動彙整）
docs/pali.md           巴利文（拉丁轉寫）對照與偈頌
```

---

## 6. 巴利文（拉丁轉寫）對照

完整對照與偈頌見 [`docs/pali.md`](docs/pali.md)；此處僅列核心語彙：

| 中文 | Pali (拉丁轉寫) | 說明 |
|---|---|---|
| 公理 | *mūla-ñāya* | 根本法則 |
| 定義 | *lakkhaṇa* | 相、特徵 |
| 定理 | *siddhanta* | 已成立之義 |
| 算式 | *gaṇana* | 計算 |
| 義務自證 | *sādhana-kicca* | 應成辦之事 |
| 演算法 | *vidhi* | 方法、行相 |
| 作用 | *kamma* | 業、作 |
| 有限單群 | *paricchinna-kevala-gaṇa* | 有量而無雜之眾 |
| 散在單群 | *ākiṇṇa-kevala-gaṇa* | 散布而無雜之眾 |

偈（gāthā）：

> *Sabbe sādhana-kiccā sādhitā honti;*
> *aṭṭhāsīti vidhī, dvāsīti kiccāni,*
> *sabbaṃ gaṇanāya siddhaṃ.*
>
> （一切義務自證皆已成辦；八十二演算法、八十二義務，一切由算式而成立。）
