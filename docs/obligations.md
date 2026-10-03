# 義務自證總表（82 條）

本表由 `src/lib/*.mbt` 中各模組之 `led.oblige(...)` 登記項自動彙整；
「結果」欄對應 `node dev/build.mjs` 之實際稽核輸出（全部 ✓）。

| # | 演算法 | 模組 | 公理 | 算式 | 命題／定理／定義／算式之陳述 | 結果 |
|---|---|---|---|---|---|---|
| 1 | `alt.order` | alt.mbt | A1 | Order | |A_n| = n!/2：A₃=3、A₄=12、A₅=60、A₆=360、A₇=2520（階乘算式 × Schreier–Sims 雙路徑） | ✓ |
| 2 | `alt.classes` | alt.mbt | A1 | Burnside | A₅ 共軛類：5 類，類大小 {1,15,20,12,12}，Σ = 60；類方程 Σ|c| = |G| 於 A₄、A₅、A₆ 成立 | ✓ |
| 3 | `alt.simplicity` | alt.mbt | A4 | Order | A₅ 單（∀x≠e 之正規閉包皆為 60）；A₄ 不單（雙對換之正規閉包 = V₄，階 4）；A₄ 之導長 = 2 | ✓ |
| 4 | `alt.three_cycles` | alt.mbt | A1 | OrbitStabilizer | 3-輪換：A₅ 之 3-輪換（每 3-子集一個定向）共 C(5,3) = 10 個而生成 A₅；A₅ 中階 3 之元素共 2·C(5,3) = 20 個 | ✓ |
| 5 | `alt.a5_psl5_bridge` | alt.mbt | A5 | LinearCount | A₅ ≅ PSL(2,5) ≅ PSL(2,4)：階皆 60，元素階分布一致（1:1、2:15、3:20、5:24），兩者皆單（同構不變量一致） | ✓ |
| 6 | `enrich.reynolds_projection` | enrich.mbt | A2 | Enriched | Reynolds 算子（3 維置換表示，S₃ 於 F₅）：π² = π、π 與作用交換、rank π = dim V^G = 1（常數函數）、π 固定不變向量 | ✓ |
| 7 | `enrich.maschke_complement` | enrich.mbt | A2 | Enriched | Maschke 定理之構造性驗證：G = ⟨diag(1,−1)⟩ 作用於 F₅²（p ∤ |G|），不變直線 W = ⟨e₁⟩ 之 π-互補 ker π = ⟨e₂⟩ 亦不變，且 V = W ⊕ ker π | ✓ |
| 8 | `enrich.char_orthogonality` | enrich.mbt | A1 | Enriched | 特徵標正交（有理數精確）：S₃ 之平凡、符號、標準（2 維）特徵標兩兩內積 = δ；Σ d² = 1+1+4 = 6 = |S₃|；A₅ 之 5 個次數 1,3,3,4,5 滿足 Σ d² = 60 | ✓ |
| 9 | `enrich.frobenius` | enrich.mbt | A3 | Enriched | Frobenius 互反：Ind_{S₂}^{S₃} 1 之特徵標 χ_perm = χ_triv + χ_std（3 維），⟨χ_perm,χ_triv⟩ = 1、⟨χ_perm,χ_std⟩ = 1、⟨χ_perm,χ_perm⟩ = 2，且 [S₃:S₂] = 3 | ✓ |
| 10 | `enrich.invariants_molien` | enrich.mbt | A2 | Molien | 不變子空間之三重算式一致（S₃ 於 3 變數，d ≤ 6）：對稱冪不變維數 = Molien 級數係數 = 部分 ≤ 3 之整數分拆數（1,1,2,3,4,5,7） | ✓ |
| 11 | `group.abelian.center` | group.mbt | A1 | Order | 中心與交換性：Z(D₈)=⟨r²⟩ 階 2；Z(S₃)=1；A₄ 之中心為 1，V₄ 之中心為自身 | ✓ |
| 12 | `group.derived` | group.mbt | A1 | Order | 導群與可解性：[S₃,S₃]=A₃ 階 3；S₄ 可解（導出列長度 4）；A₅ 不可解（導群為自身） | ✓ |
| 13 | `group.simple` | group.mbt | A4 | Order | 單性檢定：A₅（60）與 PSL(2,7)（168）為單群；S₅ 與 A₄ 非單群（正規閉包測試） | ✓ |
| 14 | `group.composition` | group.mbt | A5 | Order | 合成因子（Jordan–Hölder）：S₄ 之因子階為 {2,2,3}，A₅ 之因子為 {60}（單群） | ✓ |
| 15 | `group.sylow` | group.mbt | A1 | TwoLocal | Sylow 定理：S₄ 中 n₂=3 且 n₂≡1 (mod 2)、n₂|3；|N_G(S)|=8；Sylow 2-子群階為 8 | ✓ |
| 16 | `group.sylow3` | group.mbt | A1 | TwoLocal | Sylow 3-子群：A₅ 中 n₃=10、|N_G(P)|=6，且 10 ≡ 1 (mod 3)、10 | 20 | ✓ |
| 17 | `group.nilpotent` | group.mbt | A1 | TwoLocal | 冪零性：D₈ 與 Q₈ 之冪零類為 2；S₃ 非冪零（中心列不達 G） | ✓ |
| 18 | `group.lattice` | group.mbt | A1 | Order | 子群格與 Frattini：V₄ 有 5 個子群；S₃ 有 6 個子群且 Φ(S₃)=1 | ✓ |
| 19 | `group.burnside` | group.mbt | A2 | Burnside | Burnside 引理：S₄ 作用於 2-子集之軌道數 = (1/24)Σ|Fix|；A₅ 作用於 5 點之軌道數 = 1 | ✓ |
| 20 | `group.cayley` | group.mbt | A1 | Order | Cayley 定理：正則表示為單射同態，且 |ρ(G)| = |G|（S₃ 檢驗） | ✓ |
| 21 | `group.first_iso` | group.mbt | A3 | Order | 第一同構定理：sign : S₅ → {±1} 之核階 60、像階 2，且 60·2 = 120 | ✓ |
| 22 | `group.elementary_abelian` | group.mbt | A1 | LinearCount | 初等交換群：(Z/2)² ≅ V₄ 階 4、指數 2、交換；正則作用遞移 | ✓ |
| 23 | `hopf.group_algebra` | hopf.mbt | A3 | HopfScheme | 群代數 Hopf 公理（G = Z/3、S₃ 於 F₇）：Δ(g)=g⊗g 滿足餘結合、餘單位 (ε⊗id)Δ=id、對極 m(S⊗id)Δ = ηε、Δ 為代數同態、S 為反同態 | ✓ |
| 24 | `hopf.coaction_orbits` | hopf.mbt | A2 | HopfScheme | 集合作用之共作用：dim 不變量 = dim 餘不變量 = 軌道數——Z/2 作用於 3 點（軌道 {0} 與 {1,2}）為 2；S₃、Z/3 遞移作用於 3 點為 1 | ✓ |
| 25 | `hopf.mu_p_coaction` | hopf.mbt | A2 | HopfScheme | μ_p 群概型作用於 F_p[y]：ρ(y^k) = x^k ⊗ y^k 之不變量維數 = ⌈N/p⌉（p=3,N=6 與 p=5,N=10 皆為 2），與 α_p = d/dx 之核（導子不變元）維數一致 | ✓ |
| 26 | `hopf.mu_p_nonsemisimple` | hopf.mbt | A3 | HopfScheme | 非既約情形（p | |G|）：F₃[x]/(x³−1) 中 (x−1) 冪零指數 3（(x−1)² ≠ 0、(x−1)³ = 0），冪等元僅 0 與 1 ⇒ 無對極平均（Reynolds）算子；Maschke 假設 gcd(3,3) ≠ 1 失敗 | ✓ |
| 27 | `lie.a1_a2_orders` | lie.mbt | A5 | LinearCount | A₁(q) = PSL(2,q) 之階 q(q²−1)/gcd(2,q−1) 與矩陣計數相符（q=2,3,5,7 → 6,12,60,168）；A₂(q) = PSL(3,q)：q=2,3,4,5 → 168,5616,20160,372000 | ✓ |
| 28 | `lie.twisted_orders` | lie.mbt | A5 | BigOrder | 扭群之階：²B₂(8)=29120=2⁶·5·7·13、²B₂(2)=20；²G₂(3)=1512=3·|PSL(2,8)|=3·504；²F₄(2)=2·17971200（Tits 群 17971200=2¹¹·3³·5²·13） | ✓ |
| 29 | `lie.d4_g2` | lie.mbt | A5 | BigOrder | ³D₄(2) = 2¹²(2⁸+2⁴+1)(2⁶−1)(2²−1) = 211341312 = 2¹²·3⁴·7²·13（BigInt 精確比對）；G₂(2) = 12096 = 2⁶·3³·7，G₂(2)′ = 6048 | ✓ |
| 30 | `lie.weyl_dim` | lie.mbt | A5 | WeylDimension | Weyl 維數：A₂ 之 dim V_{a,b} = (a+1)(b+1)(a+b+2)/2 與根系乘積 Π_{α>0}⟨λ+ρ,α⟩/⟨ρ,α⟩ 一致；dim(a,b)=dim(b,a)；A₁ 之 dim = 單項式個數 k+1 | ✓ |
| 31 | `lie.roots_cartan` | lie.mbt | A5 | WeylDimension | 根系與 Cartan 矩陣：A₂/B₂/G₂ 之全根數 6/8/12（正根 3/4/6）；Weyl 群階 = |S₃|=6、|D₄|=8、|D₆|=12 對應；Cartan 行列式 A₂=3、B₂=2、G₂=1（么模）、A₃=4、D₄=4 | ✓ |
| 32 | `lie.jacobi_killing` | lie.mbt | A5 | Derivation | Lie 代數：gl₂、gl₃ 於 F₅ 之 Jacobi 恆等式成立；sl₂ 之 Killing 形式 B(h,h)=8、B(e,f)=4、B(h,e)=0（mod 7）且非退化（3×3 之秩 3） | ✓ |
| 33 | `lie.exp_bch_adjoint` | lie.mbt | A3 | Derivation | 李群作用之離散化：exp(X)exp(Y) = exp(X+Y+[X,Y]/2)（Heisenberg，F₇）；Ad(exp X)(Y) = exp(ad X)(Y) 於 sl₂ ⊂ gl₂（F₇）；exp(X)exp(−X) = I | ✓ |
| 34 | `matrix.det` | matrix.mbt | A1 | LinearCount | 行列式為乘法同態：det(AB)=det A·det B，det(I)=1，det(M)·det(M⁻¹)=1（F₂、F₃、F₅ 上） | ✓ |
| 35 | `matrix.gl_order` | matrix.mbt | A1 | LinearCount | |GL(n,p)| 公式＝直接計數：(n,p)=(2,2),(2,3),(2,5) 之可逆矩陣數 | ✓ |
| 36 | `matrix.sl_gens` | matrix.mbt | A1 | LinearCount | 切換矩陣生成 SL(2,p)：在有序基底上之忠實作用階 = |SL(2,5)| = 120；在 P¹ 上之像階 = |PSL(2,5)| = 60 | ✓ |
| 37 | `matrix.psl2` | matrix.mbt | A4 | LinearCount | PSL(2,p) 作用於射影直線 P¹(F_p)（p+1 點）：|PSL(2,5)|=60、|PSL(2,7)|=168、|PSL(2,11)|=660 | ✓ |
| 38 | `matrix.vector_action` | matrix.mbt | A2 | OrbitStabilizer | GL(2,p) 在非零向量上遞移：軌道大小 = p²−1（p=3,5）；穩定化子階 = |GL|/(p²−1) | ✓ |
| 39 | `matrix.cayley_hamilton` | matrix.mbt | A1 | LinearCount | Cayley–Hamilton：χ_M(M) = 0（對 F₂、F₃、F₅ 上所有 2×2 矩陣全測） | ✓ |
| 40 | `matrix.charpoly_dual` | matrix.mbt | A1 | LinearCount | 兩獨立算式互證：Laplace 展開與 Faddeev–LeVerrier 之特徵多項式在 p>n 時一致 | ✓ |
| 41 | `matrix.charpoly` | matrix.mbt | A1 | LinearCount | 特徵多項式：2×2 時 χ(x)=x²−tr(M)x+det(M) 且 χ(M)=0；λ 為特徵值 ⟺ ker(M−λI)≠0 | ✓ |
| 42 | `matrix.matrix_order` | matrix.mbt | A1 | Order | 矩陣階：M^k = I 且 k 為最小；階整除 |GL(n,p)|（F₅ 上之生成元） | ✓ |
| 43 | `num.is_prime` | num.mbt | A1 | Order | 質數判定與試除法一致：π(100)=25，且 561 非質數（Carmichael 數） | ✓ |
| 44 | `num.factorize` | num.mbt | A1 | Order | 算術基本定理：Π p^e = n，且因子皆質數（對 n≤2000 全測） | ✓ |
| 45 | `num.pow_mod` | num.mbt | A1 | Order | 快速冪正確性：b^e mod m 與費馬小定理 2^(p−1)≡1 (mod p) | ✓ |
| 46 | `num.ext_gcd` | num.mbt | A1 | Order | Bézout 等式：a·x+b·y=gcd(a,b)（對 3≤a,b≤40 全測） | ✓ |
| 47 | `num.mod_inv` | num.mbt | A1 | Order | 模反元素：a·a⁻¹ ≡ 1 (mod m)（對所有與 m 互質之 a） | ✓ |
| 48 | `num.totient` | num.mbt | A1 | Order | φ(n)=|(Z/n)^×|：直接計數與分解式一致（n≤300） | ✓ |
| 49 | `num.divisors` | num.mbt | A1 | Order | 因數計數：τ(n)=|divisors(n)|，且 σ(n)=Σ 因數（n≤400） | ✓ |
| 50 | `num.legendre` | num.mbt | A1 | Order | Legendre 公式：v_p(n!)=Σ⌊n/p^k⌋（n≤200, p≤13） | ✓ |
| 51 | `num.mobius` | num.mbt | A1 | Order | Möbius 反演：Σ_{d|n} μ(d) = [n=1]（n≤500） | ✓ |
| 52 | `num.factorization_algebra` | num.mbt | A1 | BigOrder | 分解代數：factorization_mul 之積等於大整數乘積，整除判定與 BigInt 餘數一致 | ✓ |
| 53 | `perm.basic` | perm.mbt | A1 | Order | 置換代數：結合律、反元素、階為輪換長度之 lcm（對 S₅ 全元素抽樣） | ✓ |
| 54 | `perm.schreier_sims` | perm.mbt | A2 | Order | BSGS 與顯式枚舉互證：|S₄| = 24、|A₅| = 60、|GL(2,3) 之置換像| 等（兩演算法一致） | ✓ |
| 55 | `perm.stabilizer` | perm.mbt | A2 | OrbitStabilizer | 軌道—穩定化：|G| = |Orb(x)|·|Stab(x)|（S₅ 作用於 5 點；A₅ 作用於 5 點） | ✓ |
| 56 | `perm.member` | perm.mbt | A2 | Order | 成員判定：偶置換屬於 A₅、對換不屬於 A₅（sifting 正確性） | ✓ |
| 57 | `perm.sign` | perm.mbt | A1 | Order | 符號為同態：sign(ab)=sign(a)sign(b)，sign(對換)=−1，sign(3-輪換)=+1 | ✓ |
| 58 | `poly.arithmetic` | poly.mbt | A1 | Order | F_p[x] 算術：除法恆等式 a = q·b + r（deg r < deg b），對 F₃、F₅ 上所有低次多項式全測 | ✓ |
| 59 | `poly.irreducible` | poly.mbt | A1 | LinearCount | 不可約性：x²+x+1 於 F₂ 不可約；x²+1 於 F₃、F₅ 可約；F₂ 上 4 個一次式乘積皆可約 | ✓ |
| 60 | `poly.frobenius` | poly.mbt | A2 | OrbitStabilizer | Frobenius 作用：對不可約 f（deg n）之根，φ : x ↦ x^p 之軌道長度 = n，且 x^{p^n} ≡ x | ✓ |
| 61 | `poly.translation_invariants` | poly.mbt | A2 | LinearCount | Z/p 平移作用 x ↦ x+1 之不變元 = F_p[x^p − x]：維數等於 k 之個數（kp ≤ N−1，p=3,5） | ✓ |
| 62 | `poly.derivation_invariants` | poly.mbt | A2 | LinearCount | α_p（無窮小平移）之不變元 = ker D = F_p[x^p]：維數為 ⌈N/p⌉（p=3,5） | ✓ |
| 63 | `poly.cyclic_quotient` | poly.mbt | A2 | Burnside | Z/n 作用於 F_p[x]/(x^n−1)（x ↦ x^k）為置換，軌道數由 Burnside 引理給出（n=5,k=2：階 4、2 軌道） | ✓ |
| 64 | `poly.molien_c2` | poly.mbt | A1 | Molien | Molien 級數：Z/2 作用 x ↦ −x 之不變環為 F[x²]，P(t) = 1/(1−t²) | ✓ |
| 65 | `poly.molien_s3_perm` | poly.mbt | A1 | Molien | Molien 級數＋CST：S₃ 置換作用於三變數，P(t) = 1/((1−t)(1−t²)(1−t³))，且生成元皆為偽反射 | ✓ |
| 66 | `poly.molien_s3_lin` | poly.mbt | A1 | Molien | Molien 級數（二維標準表示）：S₃ 作用於 x,y 之不變環為自由環（次數 2,3，乘積 = |G| = 6） | ✓ |
| 67 | `sheaf.cech_complex` | sheaf.mbt | A2 | SheafCech | Čech 複形條件 δ¹∘δ⁰ = 0：圓周（三角形神經）、球面（四面體邊界）、環面（3×3 方格）三者於 F_p（p=5,7）皆成立 | ✓ |
| 68 | `sheaf.cohomology` | sheaf.mbt | A2 | SheafCech | 餘調維數與 Euler 特徵數：圓周 (h⁰,h¹,h²) = (1,1,0)、χ = 0；球面 (1,0,1)、χ = 2；環面 (1,2,1)、χ = 0；且 Σ(−1)^i dim C^i = Σ(−1)^i h^i | ✓ |
| 69 | `sheaf.equivariant` | sheaf.mbt | A2 | SheafCech | 等變層（層作用）：Z/3 循環置換圓周之三圖卡、以及球面之 3 階旋轉，皆使 δ⁰∘ρ⁰ = ρ¹∘δ⁰（ρ⁰ 為頂點置換矩陣、ρ¹ 為定向邊上之帶號置換矩陣）；且有面之複形中面亦被置換（無 -1） | ✓ |
| 70 | `sheaf.lefschetz` | sheaf.mbt | A2 | SheafCech | Lefschetz–Hopf 跡公式：Σ(−1)^i tr(ρ^i) = Σ(−1)^i dim Fix(H^i)。圓周旋轉 L = 0；球面 3 階旋轉（固定 1 頂點 1 面）L = 2；環面平移 L = 0 | ✓ |
| 71 | `sheaf.glue` | sheaf.mbt | A2 | SheafCech | 黏合與下降（descent）：圓周上 1-上鏈 f 為上邊界 ⟺ f₀₁+f₁₂ = f₀₂；餘圈 (1,0,0) 之同調類生成 H¹（非上邊界） | ✓ |
| 72 | `sporadic.count26` | sporadic.mbt | A5 | BigOrder | 散在單群恰 26 個（CFSG 之 26 個例外）：名稱互異，且 J₁、J₃、J₄、Ly、O′N、Ru 六個 pariah 皆在列 | ✓ |
| 73 | `sporadic.factor_matches_decimal` | sporadic.mbt | A5 | BigOrder | 26 個散在單群階之雙表示一致：質因數指數表之 BigInt 乘積 = 所載十進位階（逐群精確比對） | ✓ |
| 74 | `sporadic.small_int_products` | sporadic.mbt | A5 | BigOrder | 階 < 2³¹ 之 10 群以 Int 乘法獨立複核：M₁₁=7920、M₁₂=95040、J₁=175560、M₂₂=443520、J₂=604800、M₂₃=10200960、HS=44352000、J₃=50232960、M₂₄=244823040、McL=898128000 | ✓ |
| 75 | `sporadic.monster_divisibility` | sporadic.mbt | A5 | BigOrder | 與 Monster 階之整除結構：26 群中 24 群之階整除 |M|，例外恰為 J₄（11³、37、43）與 Ly（37、67） | ✓ |
| 76 | `sporadic.monster_baby_index` | sporadic.mbt | A5 | BigOrder | |M| : |B| = 2⁵·3⁷·5³·7⁴·11·13²·29·41·59·71：以指數相加驗證 |B|·[M:B] = |M|（三者皆為 BigInt 精確值） | ✓ |
| 77 | `sporadic.large_primes` | sporadic.mbt | A5 | BigOrder | 大質數結構：26 群中 25 群之最大質因數 ≥ 11（例外恰為 J₂，最大質因數 7）；各群階皆為偶數 | ✓ |
| 78 | `2lc.sylow_theorem` | twolocal.mbt | A1 | TwoLocal | Sylow 定理：S₃ 之 Sylow 2-子群階 2、個數 3；S₄ 階 8、個數 3；S₅ 階 8、個數 15；A₅ 階 4、個數 5——皆滿足 n_p ≡ 1 (mod p) 且 n_p | [G:P] | ✓ |
| 79 | `2lc.fusion_conjugate` | twolocal.mbt | A2 | TwoLocal | 融合：S₄ 之三個 Sylow 2-子群互相共軛（單一融合軌道，長度 3），|N_G(P)| = |G|/n₂ = 8 = |P|（P 自正規化） | ✓ |
| 80 | `2lc.extraspecial` | twolocal.mbt | A1 | TwoLocal | 超特殊 2-群：D₈ 與 Q₈ 皆滿足 |Z(G)| = |[G,G]| = |Φ(G)| = 2，且指數 2 子群恰 3 個（G/Φ ≅ (Z/2)²） | ✓ |
| 81 | `2lc.central_series` | twolocal.mbt | A1 | TwoLocal | 冪零性與中心列：D₈、Q₈ 之冪零類 = 2、上中心列長 = 2（皆為 2-群：|G| = 8 = 2³）；Q₈ 為 Hamiltonian（全部子群正規），D₈ 之反射子群不正規 | ✓ |
| 82 | `2lc.two_group_action` | twolocal.mbt | A2 | TwoLocal | 2-群作用：D₈ 作用於正方形 4 頂點遞移（軌道 4、穩定子階 2，8 = 4·2）；V₄ 之正則作用遞移於 4 點；指數 2 子群必正規故共軛作用平凡（3 個固定點） | ✓ |

合計 82 條義務；去重後演算法名稱 82 條（> 48）。
