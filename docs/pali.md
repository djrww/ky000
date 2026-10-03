# 巴利文（拉丁轉寫）對照 ｜ *Pāli-pāṭha*（Pali Text）

> **說明（誠實聲明）**：巴利三藏中並無「矩陣作用」「群概型」「Čech 餘調」等現代數學名相；
> 本檔之巴利語為**為本專案而造之新複合詞（nava-samāsa）**，用作文件與偈頌之莊嚴語，
> 數學內容一律以每一條之 MoonBit 算式 (`*algo*`) 與中文陳述為準。
> 巴利語在此扮演「同義莊嚴層」：可誦、可對照、不充作證明。

---

## 一、五公理 ｜ *pañca mūla-ñāyā*

| 代號 | Pali（拉丁轉寫） | 中文 |
|---|---|---|
| A1 | *Gaṇa-ñāya* — sabbe a, b, c gaṇe: (ab)c = a(bc); atthi e: ea = ae = a; sabbassa a atthi a⁻¹: a a⁻¹ = e. | 群公理：結合律、單位元、反元素 |
| A2 | *Kamma-ñāya* — G × X → X: e · x = x; (gh) · x = g · (h · x). Samo ca: homaṇḍalaṃ φ : G → Sym(X), tassa khandhaṃ (ker φ) kammassa mūlaṃ. | 作用公理：單位與結合，等價於同態 φ:G→Sym(X) |
| A3 | *Saṃyoga-ñāya* — α : G → Sym(X), β : H → Sym(Y), homaṃ f : G → H, anuvattī φ : X → Y ⟹ g ↦ β(f(g)) kammaṃ G yā Y; guṇakammaṃ g · (x, y) = (αg x, βg y) pi kammaṃ. | 複合公理：作用沿同態與等變映射複合 |
| A4 | *Kevala-ñāya* — G paricchinnaṃ, \|G\| > 1; sace N ⊴ G, atha N = 1 vā N = G. Samo: sabbassa x ≠ e, tassa pakati-parikkhepo (normal closure) = G. | 單性公理：僅平凡正規子群 |
| A5 | *Vibhajja-ñāya* (CFSG) — Paricchinna-kevala-gaṇo G ⟹ G ≅ Z/p (p paviṭṭha-pañha) vā A_n (n ≥ 5) vā Lī-sadisa-gaṇo vā chabbīsatiyā ākiṇṇa-gaṇānaṃ aññataro. | 分類公理：Z/p、交錯群、李型、26 散在群 |

---

## 二、十二算式家族 ｜ *dvādasa gaṇana-kulāni*

| 代號 | Pali | 中文 |
|---|---|---|
| ORD | *Kama-vidhi*（階與次第） | 結構常數與階 |
| ORB | *Sayaṃ-kamma-vidhi*（自行與止住） | 軌道—穩定化 |
| BUR | *Vīthī-gaṇana*（軌道數） | Burnside 引理 |
| LIN | *Paṭṭa-gaṇana*（矩陣計算） | 線性代數計數 |
| MOL | *Moliena-saraṇī*（Molien 級數） | Molien 級數與不變環 |
| DER | *Vibhajana-kāraka*（微分作者） | 導子與 Leibniz 律 |
| WEY | *Weyla-vitthāra*（Weyl 展開） | Weyl 維數公式 |
| HOP | *Gaṇa-vidhāna*（群概型） | Hopf 公理 (Δ, ε, S) |
| 2LC | *Dvi-gaṇa-vidhi*（2-群法） | Sylow 與 2-局部融合 |
| BIG | *Mahā-saṅkhyā*（大數） | 散在單群之大整數階 |
| SHE | *Paṭala-vidhi*（層法） | 層與 Čech 餘調 |
| ENR | *Upavaddhita-kamma*（富化作用） | Reynolds 平均與特徵標 |

---

## 三、作用域（十五處） ｜ *kamma-ṭṭhānāni*

| 中文 | Pali | 算式錨點（實跑通過） |
|---|---|---|
| 矩陣作用於向量 | *paṭṭaṃ sarassa upari kammaṃ* | GL(2,5) 在 F₅²∖{0} 之軌道 = 24 = p²−1 |
| 置換作用 | *vipallāsa-kammaṃ* | Schreier–Sims BSGS，置換群階 |
| 群結構作用 | *gaṇa-saṇṭhāna-kammaṃ* | 正規閉包、Sylow 個數、Frattini |
| 多項式環作用 | *nānāpada-cakka-kammaṃ* | 平移不變元 = F_p[x^p−x]；Molien 級數 |
| 微分算子作用 | *vibhajana-kāraka-kammaṃ* | ker D = F_p[x^p]；Jacobi 恆等式 |
| 循環單群 | *vatta-kevala-gaṇo* | Z/p 正則作用遞移；μ_p 不變量 = ⌈N/p⌉ |
| 交錯群 | *parivatta-gaṇo* | \|A₅\| = 60；類 {1,15,20,12,12}；單性 |
| 李型單群 | *Lī-sadisa-kevala-gaṇo* | PSL(2,q)、PSL(3,q)、²B₂(8)、²G₂(3)、³D₄(2)、G₂(2) |
| 二十六散在單群 | *chabbīsati ākiṇṇa-gaṇā* | BigInt 雙表示一致；24/26 整除 \|M\|（例外 J₄、Ly） |
| 層作用 | *paṭala-kammaṃ* | Čech 餘調 (1,1,0)、(1,0,1)、(1,2,1)；Lefschetz L=0,2,0 |
| 李群作用 | *Lī-gaṇa-kammaṃ* | exp(X)exp(Y)=exp(X+Y+[X,Y]/2)；Ad(exp X)=exp(ad X) |
| 群概型作用 | *gaṇa-vidhāna-kammaṃ* | Δ(g)=g⊗g 之 Hopf 公理；不變量 = 餘不變量 = 軌道數 |
| 2-群作用 | *dvi-gaṇa-kammaṃ* | n₂(S₄)=3、n₂(S₅)=15、n₂(A₅)=5；融合為單一軌道 |
| 富化作用 | *upavaddhita-kammaṃ* | π² = π；⟨χᵢ,χⱼ⟩ = δᵢⱼ；不變維數 = 分拆數 |
| 數字論引擎 | *saṅkhyā-vidhi* | 質因數分解、模反元素、Legendre、Möbius |

---

## 四、偈頌 ｜ *Gāthā*

> **Ācariya-gāthā**
>
> *Pañca ñāyā mūlaṃ katvā,*
> *gaṇaṃ kammañ ca saṃyujja;*
> *asīti-dve vidhī katvā,*
> *saccaṃ gaṇanāya sādhiya.*
>
> *Sīha-Lī-sadisā gaṇā,*
> *ākiṇṇa-gaṇā chabbīsati;*
> *paṭalaṃ gaṇa-vidhānañ ca,*
> *dvi-gaṇaṃ upavaddhitaṃ.*
>
> *Sabbe sādhana-kiccā sādhitā honti:*
> *dvāsīti kiccāni, dvāsīti vidhiyo,*
> *ekūnāsīti na atthi —*
> *sabbaṃ gaṇanāya siddhaṃ.*
>
> （以五公理為根，結合群與作用；造八十二法，以算式而立於真實。
> 李型與獅子之相似群、二十六散布之群；層與群概型、二群與富化。
> 一切義務自證皆已成辦：八十二義務、八十二演算法，無一未成——
> 一切由算式而成立。）

### 逐句對照

| Pali | 中文 |
|---|---|
| *pañca ñāyā mūlaṃ katvā* | 以五條公理為根本 |
| *gaṇaṃ kammañ ca saṃyujja* | 把群與作用組合起來 |
| *asīti-dve vidhī katvā* | 造作八十二條演算法 |
| *saccaṃ gaṇanāya sādhiya* | 以可執行之算式證成真實 |
| *sādhana-kicca* | 義務自證（應成辦之事） |
| *sabbe sādhana-kiccā sādhitā honti* | 一切義務自證皆已成辦（= 「全部義務皆由可執行算式驗證通過」） |

---

## 五、程式對照 ｜ *Vidhi-paṭisandhi*

| Pali | MoonBit 對應 |
|---|---|
| *mūla-ñāya* | `pub enum Axiom { A1_Group, A2_Action, A3_Composition, A4_Simplicity, A5_Classification }` |
| *lakkhaṇa* | `Obligation.statement`（每條之命題／定理／定義／算式陳述） |
| *sādhana-kicca* | `Ledger::oblige(algo, ax, calc, statement, run)` |
| *siddhanta-ñāṇa* | `Ledger::audit(verbose) -> (通過, 總數)` |
| *vidhi-gaṇanā* | `Ledger::algorithm_count()`（去重；本專案 = 82 > 48） |
| *gaṇana-kulāni* | `Ledger::calculus_histogram()` |

> *Sādhu! Sādhu! Sādhu!* — 善哉，善哉，善哉。
