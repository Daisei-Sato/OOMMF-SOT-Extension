# OOMMF-SOT-Extension
Custom SOT (Spin-Orbit Torque) extension module for OOMMF and simulation MIF files.
# OOMMF-SOTEvolve

スピンホールナノオシレータ（SHNO）のマイクロマグネティクスシミュレーションにおいて、スピン軌道トルク（SOT）の影響を正確に評価するため、OOMMFの標準エボルバーを拡張したカスタムC++モジュールです。

⚠️ **【注意書き】**
本リポジトリで公開している `Oxs_SOTEvolve` モジュール（C++）およびLLG方程式の4項独立出力ロジックは、**現在、物理モデルの妥当性検証およびデバッグを行っている最中（進行中）**のものです。理論解や先行研究の挙動との完全な一致を確認するためのテストパターンを実行中であり、検証が完了した段階で順次アップデートを行います。

## 概要 (Overview)
本プロジェクトは、マイクロマグネティクスシミュレータ **OOMMF** の標準エボルバー（`Oxs_SpinXferEvolve`）をベースに、物理モデルをC++で拡張し、**SOT（ダンピングライク項およびフィールドライク項）**を計算可能にした `Oxs_SOTEvolve` モジュールと、それを用いたSHNOシミュレーション用のMIFファイル一式を格納しています。

## 主な実装・改良点 (Key Contributions)
1. **SOT物理モデルのC++実装 (`sotevolve.cc`)**
   - 従来のSTT（スピン注入トルク）モデルを拡張し、LLG方程式にダンピングライクSOT項、フィールドライクSOT項を組み込みました。
   - 空間的な電流密度分布やスピン偏極方向の異方性を動的に反映できるロジックを実装しています。

2. **LLG方程式の各項の独立出力機能（新規追加）**
   - 物理現象の深い詳細解析やデバッグを容易にするため、**LLG方程式の4つの項（プリセッション項、ダンピング項、ダンピングライクSOT項、フィールドライクSOT項）をそれぞれ独立したトルクとして個別に出力・保存できるロジック**を独自に実装しました。これにより、各トルクが磁気反転や発振ダイナミクスに与える寄与度を定量的に評価可能です。

3. **OOMMFアーキテクチャへの適合**
   - `OXS_EXT_REGISTER` を用いてOOMMFの拡張プラグイン（Oxs_Ext）として正しく登録。
   - `field_like_ratio` や `sigma_rule` などのカスタムパラメータをMIFファイル側から動的に制御できる柔軟な設計にしました。

4. **パラメータスイープを最適化したマルチフィジックスデータ連携 (`simulation_sot_run.mif`)**
   - 外部ツール（COMSOL等）で事前計算された不均一な電流密度分布（J）およびエルステッド磁場分布（H）のデータ（`.ovf`形式）をインポートし、実験環境に近い高度な結合シミュレーション環境を構築しました。
   - **高度なシミュレーション制御の自動化:** 効率的なパラメータスイープ（実験自動化）を行うため、COMSOL由来の電流密度分布とエルステッド磁場分布に共通の倍率パラメータ（`J_multiplier`）を連動させつつ、OOMMF側で印加する均一外部磁場には独立した別の倍率（`H_multiplier`）を適用できる変数を設計しました。これにより、物理的な整合性を保ったまま効率的な自動大量計算（スイープ検証）を可能にしています。

## 🔬 物理モデル (Physics Model)

### LLG方程式

`Oxs_SpinXferEvolve` をベースとした標準的なLLG方程式は以下の通りです。

$$
\frac{d\mathbf{m}}{dt} = -|\gamma|\\mathbf{m}\times\mathbf{H}_{\rm eff} + \alpha\left(\mathbf{m}\times\frac{d\mathbf{m}}{dt}\right) + |\gamma|\beta\theta_{SH}\(\mathbf{m}\times\mathbf{m}_p\times\mathbf{m}) - |\gamma|\beta\epsilon'\(\mathbf{m}\times\mathbf{m}_p)
$$

`Oxs_SOTEvolve` では、これをSOT (spin-orbit torque) 用に書き換え、以下の形式で実装しています。

$$
\frac{d\mathbf{m}}{dt} = -|\gamma|\\mathbf{m}\times\mathbf{H}_{\rm eff} + \alpha\left(\mathbf{m}\times\frac{d\mathbf{m}}{dt}\right) + |\gamma|\beta_{DL}\(\mathbf{m}\times\boldsymbol{\sigma}\times\mathbf{m}) - |\gamma|\beta_{FL}\(\mathbf{m}\times\boldsymbol{\sigma})
$$

第一項はプリセッション項、第二項はダンピング項で、これらはオリジナルの式と同一です。第三項・第四項がSOTトルクで、それぞれダンピングライク項 (DL)、フィールドライク項 (FL) に対応します。

### パラメータの対応関係

オリジナルの式では $\theta_{SH}$ と $\epsilon'$ という二つの独立した物理パラメータでDL項・FL項を記述していますが、本モジュールでは以下の通り再定義しています。

**$\boldsymbol{\sigma}$ （スピン偏極方向）**

$$
\boldsymbol{\sigma} = s\cdot(\hat{z}\times\hat{\mathbf{J}})
$$

$s=\pm1$ は `sigma_rule`（`z_cross_j` / `j_cross_z`）によって選択されます。

**$\beta_{DL}$ （ダンピングライク係数）**

$$
\beta_{DL} = \frac{\hbar\\theta_{SH}}{2\mu_0 et_F}\|\mathbf{J}(\mathbf{r})| \cdot J_{mult}\cdot M_s^{-1}
$$

**$\beta_{FL}$ （フィールドライク係数）**

$$
\beta_{FL} = r_{FL}\times\beta_{DL}
$$

ここで $r_{FL}$ は MIFファイルで指定する `field_like_ratio` パラメータに対応します。入力パラメータとしては `theta_SH`（SOT全体の強さ）と `field_like_ratio`（FL/DL比）の2つを与えることで $\beta_{DL}, \beta_{FL}$ が決定されます。
$\mathbf{J}(\mathbf{r})$ は空間分布を持つ電流密度場です（COMSOL等から `J_field` としてインポート）。

$t_F$ は強磁性層の膜厚です（`fm_thickness`）。

入力パラメータとしては `theta_SH`（SOT全体の強さ）と `field_like_ratio`（FL/DL比）の2つを与えることで、$\beta_{DL}$, $\beta_{FL}$ が決定されます。
### プリセッション補正項との結合について

LLG方程式を陽関数形に整理する過程で、$\beta_{DL}$ と $\beta_{FL}$ は最終的に以下のように交差して現れます（`do_precess = 1` の場合）。

$$
\mathrm{coef}_{(\mathbf{m}\times\boldsymbol{\sigma})} = \alpha\beta_{DL} - \beta_{FL}, \qquad
\mathrm{coef}_{(\mathbf{m}\times\boldsymbol{\sigma}\times\mathbf{m})} = \alpha\beta_{FL} + \beta_{DL}
$$

これは標準的なSlonczewski型STTの式変形に由来するものです。

$\beta_{DL}$ はdamping-likeトルク項だけでなくfield-likeトルク項にも寄与します（ダンピング定数 $\alpha$ を介して）。

逆に $\beta_{FL}$ もdamping-likeトルク項に寄与します。

$\alpha \ll 1$ の典型的な強磁性体では、この交差項は小さく、 $\theta_{SH}$ と `field_like_ratio` によるDL/FL項の制御が支配的です。

### 出力される各トルク項

デバッグ・解析用に、LLG方程式の各項を個別に出力できます。

- `Precession term` : $-|\gamma|\,\mathbf{m}\times\mathbf{H}_{\rm eff}$
- `Damping term` : $\alpha\left(\mathbf{m}\times d\mathbf{m}/dt\right)$
- `DL SOT term` : $|\gamma|\beta_{DL}\(\mathbf{m}\times\boldsymbol{\sigma}\times\mathbf{m})$ を中心とする項
- `FL SOT term` : $-|\gamma|\beta_{FL}\(\mathbf{m}\times\boldsymbol{\sigma})$ を中心とする項

## 今後の展望 (Future Work / Migration)
**GPU駆動シミュレータ（mumax3）への移行**
本モジュールにより、CPU環境（OOMMF）での物理ロジックの正当性検証は完了しました。現在は、より大規模なパラメータスイープおよび計算の高速化を目指し、このSOTモデルのロジックを、Go言語およびCUDA（GPU）ベースで動作する **mumax3** 環境へ移植・移行する作業を進めています。

## ディレクトリ構成 (Repository Structure)
* `sotevolve.cc` / `sotevolve.h`：自作したSOT拡張エボルバーのソースコード
* `sotEvolve_beta.mif`：SHNOシミュレーション用の設定ファイル（Tclスクリプト）

## 📖 謝辞・引用元 (References & Acknowledgments)
本プロジェクトは、米国国立標準技術研究所（NIST）が開発・提供しているオープンソースのマイクロマグネティクスシミュレータ「OOMMF」を基に拡張を行いました。素晴らしいソフトウェアを開発・維持されているコミュニティに深く感謝いたします。

* **OOMMF 公式サイト:** [NIST OOMMF Project](https://math.nist.gov/oommf/)
* **OOMMF 引用情報:**
  > M. J. Donahue and D. G. Porter, "OOMMF User's Guide, Version 1.0", NISTIR 6376, National Institute of Standards and Technology, Gaithersburg, MD (1999).

*ベースとした標準モジュール: `Oxs_SpinXferEvolve`*

Author:[佐藤　大生]
所属/専攻:[近畿大学大学院 産業理工学学研究科電子情報コース / 2028年卒予定]
