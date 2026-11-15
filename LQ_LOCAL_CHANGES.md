# LQ 项目对这份 sirocco clone 的本地改动

> 这是一份**上游 sirocco v1.3 的 clone**（`e3a8c4d`）。下面列出我们所有偏离上游的
> 改动、理由、验收和回滚方法。**在这棵树里改任何东西，都要回来补一行。**
>
> 完整背景在 `../Workspace/sirocco/README.md`，项目权威记忆在
> `../Workspace/LQ_RT_PROJECT_MEMORY.md`。这里只做索引，不重复论证。

用 `git status` 随时看当前偏离；每项都有 `*.lq_backup` 或 `*_backup` 原件。

---

## 1. `xdata/atomic_macro/h20_hetop_standard80_elem_ions.dat` — 金属丰度（2026-08-19）

**改了什么**：14 行 active `Element` 从 Verner+94 换成 AGSS09（Pasta 的 FreeEOS/OPTAB
表所用的同一套），12 个共有金属再整体 +0.0037 dex 补回 sirocco 没有的 6 个元素
（P, Cl, Ti, Cr, Mn, Ni，仅占金属质量 0.87%）。

```
              X        Y         Z
改之前     0.7069   0.2743   0.01874
改之后     0.7374   0.2492   0.01336   = Pasta，4 位一致
```

**为什么**：两个码要逐格对照，成分不一致会污染每一个结论。动 sirocco 而不动 Pasta
是因为这边是一个文本文件，那边要重跑整个 OPTAB（单色谱已删）。

**⚠️ 被注释掉的 `#Element` 行原样保留，仍是 Verner+94。** 取消注释而不换算会把
不一致带回来。

**验收**：`sirocco --dry-run` 退出 0，无 `Element line incorrectly formatted`。

**复现 / 回滚**：
```
python ../Workspace/analysis/lq_abundance_alignment.py --apply
cp xdata/atomic_macro/h20_hetop_standard80_elem_ions.dat.verner94_backup \
   xdata/atomic_macro/h20_hetop_standard80_elem_ions.dat
```

---

## 2. `source/agn.c` — 归一化窗口不再硬编码 2–10 keV（2026-08-19）

**改了什么**：`agn_init` 里 `emit_2_10 <= 0` 时不再 `exit(0)`，改用**模型自身的
频率范围**当归一化锚点，并打印说明。另加一行 `#include "models.h"`。

**为什么**：上游把 `models` 路径的锚点写死为 2–10 keV，那是 **AGN 惯例**（AGN 有
X 射线冕，2–10 keV 是标准归一化）。我们的源是 TDE 风光球，该带**按构造是空的** ——
普朗克在 2 keV 之上的份额在最热处是 `1.1e-16`，交接半径处是 `4.8e-48`。
而我们的频率网格顶面是 920 eV，所以任何当前 SED 都会让上游版本直接退出。

**只有本来会 `exit(0)` 的分支改变行为** ⟹ 有 2–10 keV 功率的模型（即所有 AGN
用法）**逐位不变**，`.pf` 格式也不动（无新增输入键）。

**验收**（4 rank / 2e4 photons）：注入光度 `2.416295883857e+45` vs 目标
`2.4163e+45`，比值 **0.99999825**。改前需要一个 1.5% 的经验标定因子，改后不需要。

**回滚**：`cp source/agn.c.lq_backup source/agn.c && cd source && make sirocco`

**值得向上游提**：硬编码一个 X 射线带来归一化任意 `models` 输入是真实的设计缺陷，
我们的修法对 AGN 零影响。真要提的话，把那段 TDE 专用的注释压成三行讲通用原则。

---

## 3. `source/sirocco.h` — `MAXSCAT` 从 2000 提到 20000/200000（2026-08-23）

⚠️ **本条补记于 2026-08-25。** `sirocco.h` 一直在 `git status` 的改动列表里，
但这份文件此前**没有记录它** —— 十条生产臂的 provenance 因此有一个缺口，
现在补上。（发现方式：为重跑规格核对"九条臂是否来自同一构建"时，
逐条对照 `sirocco_status.txt` 的四个改动文件与本文件的条目，`sirocco.h` 对不上。）

**改了什么**：

```c
#ifndef MAXSCAT
#define MAXSCAT    20000
#endif
```
（原为 `#define MAXSCAT 2000`）。`#ifndef` 包裹是为了能用 `-DMAXSCAT=200000`
在编译期覆盖，于是有两个二进制：`bin/sirocco_maxscat20k` 与 `bin/sirocco_maxscat200k`。
**十条生产臂全部用 `sirocco_maxscat200k`**（`checksums.sha256` 里二进制哈希
`d49dac9b`，十条臂一致）。

**为什么**：这些臂的边界光深到 τ_es = 362，散射次数分布是**幂律**（α≈1.1–2.1，
无特征标度），所以 MAXSCAT 是一个会截断真实分布的硬上限，不是安全网。

**验收（2026-08-25 实测）**：最深臂 R2c（τ=362）的日志
`luminosity lost by too many scatters = 0.000000000000e+00` ⟹ **零截断，200k 够用。**

**回滚**：`source/sirocco.h.lq_backup` 是原件。

---

## 构建注意

- `make prototypes` **跑不了**，本机没装 `cproto`。⟹ **不要新增 `.c` 文件**，
  除非自带头文件绕开 `templates.h`（那是生成文件，不要手改）。
- 构建结尾报 `could not find a gnu version of indent` 是**无害的**，二进制已经
  拷进 `bin/`。

---

## 4. `windsave2table`: the master table now carries the radial radiative flux

**File:** `source/windsave2table_sub.c` (backup `*.lq_backup`)
**Date:** 2026-08-23

**Why.** The LQ project's first deliverable is `L(r)`, and nothing
`windsave2table` wrote was a flux. `cell_spec_flux` is named flux but
`estimators_simple.c:160,726` accumulates `w_ave*ds/(4 pi V dnu)` — the mean
intensity `J_nu`. `dmo_dt` is momentum from resonant scattering only
(`resonate.c:1249`). `cool_adiab` is the gas pdV term (`cooling.c`,
`n k T V div v`), not the radiation's. So `L(r)` was simply unavailable, and the
only fallback — `4 pi r^2 * 4 pi J_direct <= L <= 4 pi r^2 * 4 pi J` — is
useless inside `tau ~ 1` because the direct component is `exp(-tau)` and the
lower end collapses.

**What changed.** Sirocco already computes the flux and never wrote it out:
`F_vis` / `F_UV` / `F_Xray` accumulate `w_ave*ds*lmn` and are divided by the
cell volume (`estimators_simple.c:330-365,762`) — the Monte Carlo flux
estimator — and element `[0]` is the radial component in spherical geometry.
The three bands **partition** the spectrum (`< UV_low = 7.4e14`,
`> UV_hi = 3e16`, else), so their sum is the total with nothing left over.
Two columns were added to `create_master_table` (23 -> 25 columns; the array
was already sized 51) plus two cases in `get_one`:

    F_rad_r     F_vis[0] + F_UV[0] + F_Xray[0]                    last cycle
    F_rad_per   the same from the *_persistent (cycle-smoothed) arrays

Then `L(r) = 4 pi r^2 * F_rad_r`. Both are written so they can be checked
against each other: on a converged run they agree, and disagreement is the
diagnostic. Measured on the R1 tau=10 pilot they agree to 1-6%.

**Acceptance.** No behaviour change to `sirocco` itself — only the table writer
was touched, and only by addition; every pre-existing column keeps its name,
order and value. **The flux is already inside every `.wind_save`** (the whole
plasma struct is written as a binary blob, `windsave.c:98`), so runs made before
this change need only `windsave2table` re-run, not `sirocco` re-run.

**Rollback.** `git checkout source/windsave2table_sub.c && (cd source && make windsave2table)`.
