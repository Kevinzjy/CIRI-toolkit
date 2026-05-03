# 基于 Scan1/Scan2 的 circRNA 全长结构识别设计

本文档定义 CIRI-toolkit 在 **不引入 remap** 的前提下，如何在现有 `Scan1 -> Scan2 -> Summary` baseline 输出基础上，进一步识别同一 BSJ 下的 circRNA 全长结构与候选 isoform 序列。

当前目标不是 genome-wide 地重新发现一批新的 circRNA，也不是用新的启发式去改变 baseline `.bsj` / `.out` 结果；目标是：

- 先用 Rust CIRI3 baseline 找到可信 BSJ。
- 再围绕这些已确认 BSJ 收集局部 read evidence。
- 在局部范围内构建 circRNA splice graph。
- 输出同一 BSJ 下的 anchored isoform 与 candidate side path。
- 同时输出结构表 `.isoforms` 与序列文件 `.fa`。

从算法形态上看，这条路线整体上更接近 `docs/CIRI-AS.md` 描述的思路，而不是依赖 RO remap 的 `CIRI-full RO2` 路线：

- 都以 confirmed circRNA / BSJ 作为后处理锚点。
- 都需要在 `Summary` 之后重扫一次原始比对结果。
- 都需要抽取 circ span 内的 internal splice junction、coverage 与 exon-like 结构。
- 都需要把 read-level evidence 进一步汇总成 circ-local path / exon chain。

因此，后续 full-length 实现时，`CIRI-AS` 文档中的 circ cluster、internal splice 提取、splice signal 校验、cirexon 构建与 exon path 规则应作为第一参考；`CIRI-full` 的 RO1/RO2 逻辑则主要保留为未来 side evidence 增强方案。

当前实现进度（2026-04-30）：

- `--as` sidecar 已生成 `<as_prefix>_splice.list` 与 `<as_prefix>.list`。
- 已基于 validated cirexon graph 输出 `<as_prefix>.isoforms`、`<as_prefix>.isoform_summary` 和 `<as_prefix>.fa`。
- 当前不做 ES/A5SS/A3SS/IR 类型判断；本阶段目标是回答同一 BSJ 下有几条 full-length isoform path，以及每条 path 的 exon chain / junction chain / sequence。
- `<as_prefix>_AS.list` 仅作为历史兼容 header 保留，不再作为当前开发主线。

## 1. 设计定位

本阶段采用：

```text
baseline BSJ detection
  -> Summary confirmed BSJ
  -> second BAM/SAM sweep
  -> BSJ-local evidence bundle
  -> circular splice graph
  -> isoform path ranking
  -> .isoforms + .fa
```

设计边界：

- 不修改 baseline `Scan1 -> Scan2 -> Summary` 的判定逻辑。
- 不改变 baseline `.bsj1`、`.bsj`、`.out` 的内容与正确性。
- full-length 模块只作为 `Summary` 之后的独立后处理阶段运行。
- 第一版不依赖 `bwa mem`、RO remap、RO BAM 或 integrated Scan2。
- 第一版优先做 **BSJ-anchored full-length reconstruction**，不做无锚点的 novo circRNA 发现。

## 1.1 与 CIRI-AS 的对应关系

当前 full-length 路线与 `docs/CIRI-AS.md` 的对应关系可以概括为：

### 直接复用的思路

- `load_circ_records(result_path)`
  - 直接对应本文档的 `load_bsj_seeds`
- circRNA 区间聚类
  - 对应本路线 second BAM/SAM sweep 前的 seed 分桶与局部扫描范围建立
- `mapping_check1 / mapping_check1_add`
  - 对应从 seed read 与 circ span 内普通 read 中提取 internal splice / coverage / range evidence
- `splice_loci_check`
  - 对应 de novo internal junction 的 motif / annotation-assisted 校验
- `cluster_reads`
  - 对应 internal splice candidate 聚类，避免单 read 抖动直接生成大量假边界
- cirexon 构造
  - 对应 anchored path / candidate side path 使用的 candidate exon 生成
- exon path 构建
  - 对应本路线的 circular splice graph path 枚举

### 需要保留但不必第一版完整照搬的部分

- coverage validation
  - 第一版建议保留“candidate 内 coverage 是否连续”“junction 两侧是否有明显支撑”的轻量规则；Mann-Whitney U-test 风格统计可后置
- intron retention / AS taxonomy
  - `ES / A5SS / A3SS / IR` 属于 path 解释层；第一版全长结构输出不要求同步完整复刻
- corrected PSI
  - 第一版不是必须项，先把结构和 sequence 做对，再讨论 insert length normalization

### 当前 full-length 路线对 CIRI-AS 的收敛

与原始 CIRI-AS 相比，本路线做了两点收敛：

1. 主目标从 “AS event classification” 转成 “full-length circRNA path / sequence reconstruction”。
2. 输出优先级从 `_splice.list / _AS.list / PSI` 转成：
   - `<prefix>.isoforms`
   - `<prefix>.isoform_summary`
   - `<prefix>.fa`
   - `<prefix>.full.bundle.tsv`（后续调试型 evidence bundle）

因此，实现时不应机械地复刻 CIRI-AS 的历史输出命名，而应优先复用其：

- read-level evidence 抽取方式
- splice candidate 聚类方式
- exon/path 构建方式

再把结果落到当前 full-length 专用协议里。

## 2. 为什么第一版不做 remap

当前目标是“快速识别全长结构”，而不是“重新做一轮全基因组候选发现”。对这个目标来说，全局 `bwa mem` 有两个问题：

- 计算开销高，尤其是 genome-wide merged read remap。
- 大多数真正有价值的信息已经存在于 ordinary BAM/SAM 中，包括：
  - baseline BSJ supporting reads
  - circ span 内的 split alignments
  - supplementary / soft-clipped alignments
  - mate orientation 与 insert 关系
  - annotation exon boundary

因此第一版 full-length 模块采用：

- **alignment-based**
- **annotation-assisted**
- **BSJ-anchored**

后续如果某些 unresolved case 确实需要 remap，再把 remap 作为单独补救路径讨论，而不是第一版主流程依赖。

## 3. 输入与输出边界

### 3.1 输入

第一版 full-length 模块只依赖以下输入：

- `<prefix>.out`
  - Summary 最终 circRNA 表，提供 confirmed BSJ seed
- `<prefix>.bsj`
  - read-level BSJ support 记录，提供每个 circ 的 direct BSJ evidence
- input BAM/SAM
  - 第二遍扫描读取原始 alignment evidence
- reference FASTA
  - 从最终 exon chain 提取 circRNA sequence
- annotation GTF
  - 提供 exon boundary、gene exon 与 strand 辅助

### 3.2 输出

第一版新增以下 sidecar 输出：

```text
<prefix>.full.bundle.tsv
<prefix>.isoforms
<prefix>.isoform_summary
<prefix>.fa
```

约束：

- baseline `.bsj` / `.out` 保持不变。
- `.isoforms` / `.isoform_summary` / `.fa` 是新模块输出，不回写 baseline 结果。
- 所有 anchored / candidate path 必须显式标注证据层级，不能混成一个最终置信度口径。

## 4. 术语与路径分层

### 4.1 anchored path

anchored path 指：

- 以 baseline confirmed BSJ 为锚点；
- path 的边界严格使用该 circ 的 `BSJ start` 与 `BSJ end`；
- 所有 exon / splice junction 均落在 `BSJ start .. end` 内；
- 具备直接 BSJ support。

它是第一版主输出，也是默认优先级最高的 path。

### 4.2 candidate side path

candidate side path 指：

- 仍然挂靠在某个 confirmed BSJ 下；
- 没有直接 BSJ support；
- 但有 circ span 内的 side evidence，例如：
  - `5' RO-like` ordinary paired-end pattern
  - mate orientation 异常但仍位于 circ span 内的 read pair
  - internal split read
  - annotation / motif 辅助边界

这类 path 的角色是：

- 作为同一 BSJ 下的候补内部结构；
- 先保留、先输出；
- 后续再通过 simulator / truth / real sample 评估其正确率。

第一版要求所有 candidate side path 必须显式标注：

- `path_tier = candidate`
- `bsj_supported = no`

### 4.3 mate orientation 异常 read

本文档中的 “mate orientation 异常 read” 指：

- paired-end 两端都能在 ordinary BAM/SAM 中找到 alignment；
- 落在某个 confirmed circ 的 `BSJ start .. end` 内或近邻窗口；
- 配对方向 / 插入关系不符合该局部线性转录本的常规预期；
- 但本身没有进入 baseline `junction_reads_ID`。

这类 read 不直接当作新的 BSJ 证据，而是作为：

- `candidate side path` 的 side support；
- 或 annotation / split evidence 不足时的弱结构辅助。

## 5. 模块拆分建议

建议新增：

```text
src/full_length.rs
src/full_bundle.rs
src/full_graph.rs
src/full_fasta.rs
```

职责：

- `full_length.rs`
  - full-length runner
  - 负责阶段调度、I/O 和输出
- `full_bundle.rs`
  - second BAM/SAM sweep
  - 构建每个 circ 的局部 evidence bundle
- `full_graph.rs`
  - 从 bundle 构建 circular splice graph
  - 枚举 path、打分和排序
- `full_fasta.rs`
  - 按最终 exon chain 提取 circRNA sequence

第一版不要把这些逻辑塞回 `Scan1` / `Scan2` 热路径，避免影响 baseline parity 与性能边界。

## 6. 阶段 1：加载 BSJ seeds

从 `<prefix>.out` 读取 confirmed circRNA，构建 `CircSeed`：

```rust
struct CircSeed {
    circ_id: String,
    chr: String,
    start: i32,
    end: i32,
    strand: char,
    gene_id: String,
    circ_type: String,
    junction_reads: usize,
    junction_read_ids: Vec<String>,
    score: i32,
}
```

约束：

- `CircSeed` 只来自 Summary confirmed circ。
- 第一版不从 RO-only、candidate-only 或 sidecar-only 结果中额外新增 seed。
- `start/end` 直接视为 anchored path 的硬边界。

## 7. 阶段 2：second BAM/SAM sweep 与 evidence bundle

第二遍扫描不重新判 BSJ，而是围绕已知 `CircSeed` 收集局部结构证据。

### 7.1 evidence 类型

建议收集以下 evidence：

- `bsj_seed_read`
  - 来自 `.out` 的 `junction_reads_ID`
- `internal_split_read`
  - CIGAR 含 `N`，且 splice junction 落在 circ span 内
- `boundary_clip_read`
  - soft-clip / supplementary 断点靠近 exon boundary 或 BSJ boundary
- `mate_link_read`
  - seed read 的 mate，或 circ span 内普通覆盖 read
- `mate_orientation_anomalous`
  - 位于 circ span 内、配对方向异常、但无直接 BSJ support 的 paired-end read
- `coverage_read`
  - 支持 exon body 覆盖的普通 read

### 7.2 bundle 中间产物

建议先落地为：

```text
<prefix>.full.bundle.tsv
```

第一版 schema 建议：

```text
circ_id
read_id
mate
evidence_type
chr
aln_start
aln_end
cigar
flag
mapq
strand
junction_start
junction_end
is_seed_bsj_read
within_bsj_span
boundary_support
annotation_support
motif_support
```

字段说明：

- `circ_id`
  - 该 read 当前归属的 confirmed BSJ
- `evidence_type`
  - `bsj_seed_read / internal_split_read / boundary_clip_read / mate_link_read / mate_orientation_anomalous / coverage_read`
- `junction_start/junction_end`
  - 若该 read 提供内部 splice junction，则记录其边界；否则记 `NA`
- `within_bsj_span`
  - 当前 alignment 或推断 junction 是否落在 `BSJ start .. end` 内
- `boundary_support`
  - 当前 read 是否支持某个 exon boundary / BSJ boundary
- `annotation_support`
  - 该 boundary 或 junction 是否与 annotation 一致
- `motif_support`
  - de novo junction 是否具备 canonical motif 辅助

设计原则：

- bundle 是 full-length 阶段的正式输入，不再依赖从 `.out` / `.bsj` 逆向推断所有 read-level 结构细节。
- bundle 可以同时保存 anchored path 与 candidate side path 所需证据，但必须保留 `evidence_type` 区分。
- 这一步的实现应优先参考 `docs/CIRI-AS.md` 中的 `record_mapping_detail`、`mapping_check1`、`mapping_check1_add` 与 circ cluster 扫描思路，而不是另起一套全新的 read-level 协议。

## 8. 候选 boundary 生成

第一版 candidate boundary 来源分四类：

1. annotation exon boundary
2. split-read confirmed boundary
3. de novo boundary with motif / annotation-nearby support
4. side-evidence-only boundary

优先级建议：

```text
annotation exon boundary
  > split-read confirmed boundary
  > de novo boundary with motif/annotation support
  > side-evidence-only boundary
```

为了减少单条 read 带来的边界抖动，第一版建议直接参考 `CIRI-AS cluster_reads` 的两层聚类思路：

- 先按 `site2` 近邻聚类；
- 再按 `site1` 近邻聚类；
- cluster 坐标取中位 read 的 boundary；
- 只有进入 splice cluster 的 candidate 才能进入后续 exon/path 构造。

### 8.1 de novo boundary 纳入条件

用户要求第一版允许：

- `1 条高质量 split read`
- 再加 `annotation` 或 `motif` 至少一个辅助条件

因此第一版 de novo junction 的推荐 gate 为：

```text
split_read_count >= 1
AND high_mapq = true
AND (annotation_nearby = true OR canonical_motif = true)
```

其中：

- `annotation_nearby`
  - junction boundary 与 annotation exon boundary 的距离在允许窗口内
- `canonical_motif`
  - donor/acceptor 满足 canonical splice signal

### 8.2 side-evidence-only boundary

仅由 `mate_orientation_anomalous` 或类似弱 side evidence 支持的 boundary：

- 不得单独升级成 anchored path 的主 boundary；
- 只能作为 candidate side path 的弱辅助；
- 必须在输出中显式体现其较低证据层级。

## 9. circular splice graph

第一版建议使用 boundary graph，而不是直接用 read segment 作节点。

### 9.1 节点

节点类型：

- `BSJ_START`
- `BSJ_END`
- `donor`
- `acceptor`

所有 anchored path 节点必须满足：

```text
BSJ start <= node <= BSJ end
```

### 9.2 边

边类型：

- `exon_edge`
  - 表示 exon body
- `splice_edge`
  - 表示 donor -> acceptor 的内部 splice
- `closure_edge`
  - 表示 `BSJ_END -> BSJ_START`

### 9.3 anchored path 约束

anchored path 必须：

- 经过 `closure_edge`
- 边界严格使用 `CircSeed.start/end`
- 所有 exon / splice junction 落在 circ span 内
- 具备 direct BSJ support

### 9.4 candidate side path 约束

candidate side path：

- 仍然挂在某个 confirmed BSJ 下
- 不要求 direct BSJ support
- 可以使用 side evidence 补充内部结构
- 但必须在输出中标记：
  - `path_tier = candidate`
  - `bsj_supported = no`

## 10. path 打分与排序

第一版目标是“可解释、可调试”，而不是一开始就做复杂统计模型。

建议 path score 由以下部分组成：

- `+` baseline BSJ seed support
- `+` internal split read support
- `+` annotation-supported boundary
- `+` canonical motif support
- `+` mate-link consistency
- `+` circ span 内连续 coverage 支持
- `-` 仅靠弱 side evidence 的 boundary
- `-` 不被 annotation 支持且无 motif 的 de novo junction
- `-` 需要长距离、低支持的跳跃连接

排序规则建议：

1. anchored path 优先于 candidate side path
2. 在 anchored path 内按 score 排序
3. 在 candidate side path 内按 score 排序

每个 circ 第一版可保留：

```text
top-k = 3
```

后续再视真实数据复杂度调整。

## 11. 输出协议

### 11.1 `<prefix>.isoforms`

第一版建议字段：

```text
circ_id
isoform_id
chr
start
end
strand
path_tier
bsj_supported
gene_id
exon_chain
junction_chain
exon_count
isoform_len
bsj_seed_reads
internal_split_reads
boundary_clip_reads
anomalous_pair_reads
annotation_supported_junctions
de_novo_supported_junctions
score
rank
```

字段说明：

- `path_tier`
  - `anchored` 或 `candidate`
- `bsj_supported`
  - `yes` 或 `no`
- `exon_chain`
  - 建议复用 simulator 风格：
    `start:end!strand,start:end!strand,`
- `junction_chain`
  - 内部 splice junction 序列，便于调试与后续比较

### 11.2 `<prefix>.fa`

### 11.2 `<prefix>.isoform_summary`

每个 Summary confirmed circRNA 一行，用于直接查看同一 BSJ 下当前识别出几种 full-length isoform：

```text
circ_id
chr
start
end
strand
gene_id
junction_reads
cirexon_count
isoform_count
```

`isoform_count = 0` 表示当前已确认 BSJ 存在，但 validated cirexon graph 还没有形成从 circ start 到 circ end 的完整 anchored path。

### 11.3 `<prefix>.fa`

两类 path 都输出，但 header 必须显式标注证据层级：

```text
>circ_id|isoform_id|tier=anchored|bsj_supported=yes|score=...
>circ_id|isoform_id|tier=candidate|bsj_supported=no|score=...
```

序列生成规则：

- 按最终 `exon_chain` 从 reference FASTA 提取 sequence
- 若为负链，按 transcript 方向 reverse-complement
- 不依赖 remap 或额外 merged read FASTQ

## 12. 与现有 RO 设计的关系

第一版 full-length reconstruction 与 `docs/05-ro-feature-plan.md` 的关系是：

- baseline `.bsj` / `.out` 仍然是主锚点
- full-length 先走 **ordinary alignment-only** 路线
- RO remap / integrated Scan2 / RO-enhanced BSJ rescue 暂不作为第一版前置依赖
- 若 future 需要，可把 RO sidecar evidence 作为第二阶段增量接入

因此，当前 full-length 路线与 RO 规划是并行关系，不是严格串行依赖。

## 13. 第一版刻意不做的事情

为保证边界清晰，第一版明确不做：

- 修改 `.bsj` / `.out`
- 依赖 `bwa mem` remap
- 从 candidate-only / RO-only 结果新增 circ seed
- 无 BSJ 锚点的 genome-wide circRNA 发现
- PSI / usage correction
- 纯 de novo transcript reconstruction

第一版要先回答的问题是：

- 对于已确认 BSJ，能否快速、稳定、可解释地输出 anchored isoform 与 candidate side path？

## 14. 验证建议

第一版验证建议分三层：

### 14.1 bundle layer

- `<prefix>.full.bundle.tsv` 是否覆盖了所有 `junction_reads_ID`
- circ span 内 internal split read 是否被完整抓取
- anomalous pair 是否按定义落在 circ span / 近邻窗口内

### 14.2 path layer

- anchored path 是否始终严格限制在 `BSJ start .. end`
- candidate side path 是否都带有 `bsj_supported=no`
- annotation-only、annotation+split、de novo-assisted 的 path 是否可区分

### 14.3 sequence layer

- `.fa` sequence 是否与 `exon_chain` 一致
- 负链 reverse-complement 是否正确
- simulator truth 中已知 isoform 的 exon chain / length 是否可回归

## 15. 推荐实现顺序

1. 读取 `.out`，建立 `CircSeed`
2. 第二遍 BAM/SAM sweep，输出 `<prefix>.full.bundle.tsv`
3. 只实现 annotation-assisted anchored path
4. 加入 de novo internal junction
5. 加入 `mate_orientation_anomalous` candidate side path
6. 输出 `.isoforms` 与 `.fa`

第一版建议先把 collector 做对，再做 graph；因为如果没有稳定的 bundle，后续图和 path 排序都无法调试。

---
最后更新：2026-04-30
