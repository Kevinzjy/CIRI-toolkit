# 基于 Summary 后处理的 circRNA segments 识别设计

本文档定义 CIRI-toolkit 当前阶段在 `Scan1 -> Scan2 -> Summary` 之后的默认后处理目标：先稳定识别 circRNA 相关 reads 的 read-level segments，并输出 `<prefix>.segments`。  
直接从现有 evidence 一步到位做 genome-wide full-length isoform reconstruction，目前实践上不可控，容易把 read-level 归属错误放大成 path-level 假阳性，因此当前路线先退回到更可验证的 segments 层。

当前目标不是重新发现新的 circRNA，也不是引入新启发式改变 `priority=1` 主流程证据或 `.out`。当前阶段要先回答的问题是：

- 对于已识别到的 circRNA 相关 reads，能否稳定恢复 read-level segment chain；
- 能否先把 `circ_id`、`is_circular` 和 `r1/r2_segments` 判定做对；
- 能否让这些输出直接和 simulator 的 `.reads.tsv` 对照，作为后续 isoform reconstruction 的可信输入层。

## 1. 当前阶段定位

当前主线改为：

```text
Scan1   -> <prefix>.bsj1 + <prefix>.segments1
Scan2   -> <prefix>.bsj2 + <prefix>.segments2
Summary -> <prefix>.out
segments -> read <prefix>.out + <prefix>.segments1/2
         -> confirmed BSJ read segment chains
         -> supplemental backward read scan
         -> <prefix>.segments
         -> future full-length reconstruction
```

设计边界：

- 不修改 baseline `Scan1 -> Scan2 -> Summary` 的判定逻辑。
- 不改变 `priority=1` 主流程证据和 `.out` 的判定结果。
- `<prefix>.segments1` / `<prefix>.segments2` 是 sidecar evidence，不作为 Summary 输入。
- 不再使用 `--as` / `--as-out` 作为额外入口。
- 主流程在写完 `<prefix>.out` 后默认继续执行 segments 后处理。
- 第一阶段只输出 confirmed BSJ reads 和 backward/circular candidate reads 的 read-level segments，不直接输出 `.isoforms`、`.isoform_summary`、`.fa`。
- full-length reconstruction 暂时降级为后续阶段，必须建立在 read-level segments 已经稳定可回归的前提上。

## 2. 为什么先做 segments，不直接做 full-length

直接 full-length reconstruction 当前行不通，核心原因不是“缺少图模型”，而是 read-level 输入层还不够稳定：

- 同一 read pair 可能有 `primary / supplementary / secondary` 多套可选比对链；
- circRNA read 的最佳解释需要同时满足 topology、BSJ 结构和 splice 校正，而不是简单读取某一条主比对；
- 一旦 read-level `circ_id`、`is_circular` 或 segment chain 错了，后续 exon graph / isoform path 会系统性放大错误；
- simulator 目前最直接、最严格的 oracle 也是 `.reads.tsv`，而不是 path-level event taxonomy。

因此当前阶段必须先把“单条 read 应该被解释成什么样”做成正式输出，再谈 path 构建。

## 3. 与 CIRI-AS 的关系

这条新路线仍然是 CIRI-AS-style 的后处理，只是阶段目标从“直接产出 full-length path”收缩为“先产出可靠的 segments”。

当前仍然直接保留的 CIRI-AS 核心思路：

- 以 Summary confirmed circRNA / BSJ 为后处理锚点；
- 在 `Scan1` / `Scan2` 处理 read group 时同步保存 BSJ reads 的 alignment evidence；
- 使用 circ-local read grouping 和 alignment-chain 解释；
- 使用 annotation-guided + de novo splice signal 逻辑修正 splice site；
- 用局部 topology 解释 circRNA 相关 reads，而不是把每条 supplementary 直接当成独立证据。

当前明确暂缓的 CIRI-AS / full-length 输出：

- `_splice.list`
- `.list`
- `.isoforms`
- `.isoform_summary`
- `.fa`
- `_AS.list`
- PSI / usage correction
- ES / A5SS / A3SS / IR 等 event taxonomy

这些能力不是被否定，而是顺序后移。只有当 `<prefix>.segments` 已经能稳定对齐 simulator truth，才值得继续往 path 层推进。

## 4. 触发方式与输入边界

### 4.1 触发方式

segments 后处理是主流程默认阶段：

```text
ciri \
  -i <input.sam|input.bam> \
  -o <prefix> \
  -r <ref.fa> \
  -a <anno.gtf>
```

不再使用：

- `--as`
- `--as-out`

### 4.2 输入

当前阶段只依赖以下输入：

- `<prefix>.out`
  - Summary 最终 circRNA 表，提供 confirmed BSJ 和 `junction_reads_ID`
- `<prefix>.segments1`
  - Scan1 对 BSJ read groups 捕获的 mapper alignment evidence
- `<prefix>.segments2`
  - Scan2 对 rescued / mate-level BSJ read groups 捕获的 mapper alignment evidence
- input BAM/SAM
  - Summary 后补充扫描非 BSJ read groups，用于识别 `type=backward` circular/backward candidates
- annotation GTF
  - 提供 exon boundary、gene span、strand 辅助
- reference FASTA
  - 当前阶段不是必须依赖，但保留给后续 splice signal / path 扩展

### 4.3 输出

当前阶段只新增一个正式输出：

```text
<prefix>.segments
```

约束：

- `priority=1` 主流程证据和 `.out` 判定结果保持不变；
- `<prefix>.segments` 是后续 full-length reconstruction 的正式输入层；
- `<prefix>.segments1` / `<prefix>.segments2` 是临时 sidecar 协议，不直接作为用户最终解释结果；
- 第一阶段不再把 full-length 路径文件作为主输出协议。

### 4.4 sidecar evidence 协议

`<prefix>.segments1` / `<prefix>.segments2` 是内部临时文件，当前列顺序为：

```text
read_id
stage
mate
flag
chrom
pos
mapq
cigar
read_len
clips
```

语义：

- `stage` 取 `scan1` / `scan2`，只用于追踪 evidence 来源；
- `mate` 取 `R1` / `R2`，方便人工排查；
- `flag/chrom/pos/mapq/cigar/read_len` 是最终 chain reconstruction 需要的 mapper block；
- `clips` 只在 CIGAR 含 soft clip 时写入 `L:<seq>` / `R:<seq>` clipped subsequence，否则写 `*`，用于追踪 validator 接受时可定位的 clip 序列；
- validator 接受 BSJ 后，长度 `>=10bp` 且能在 circ 区间 exact match 的 soft clip 会在同一 sidecar 中追加为 `scan1_local` / `scan2_local` pseudo-alignment row；若整段 clip 无法 exact match，则允许记录最长 prefix/suffix partial exact match，但仍要求 retained match 长度 `>=10bp`；
- 最终 `<prefix>.segments` 会用 `<prefix>.out` 的 `junction_reads_ID` 过滤 confirmed BSJ reads；`type=backward` / `type=outward` rows 来自 Summary 后额外扫描的非 BSJ read groups；
- sidecar 文件不得反向影响 `.bsj1` / `.bsj2` / `.out`。

## 5. `<prefix>.segments` 协议

### 5.1 列顺序

```text
read_id
type
circ_id
chrom
start
end
strand
is_circular
is_r1_bsj
is_r2_bsj
r1_cigar
r1_segments
r2_cigar
r2_segments
```

最终 `<prefix>.segments` 是 read-level row，但写出顺序按 genomic position 排序，便于后续 circRNA-level graph construction 顺序读取同一局部区域的 evidence。排序键为：

```text
chrom, start, end, type_rank, circ_id, read_id
```

其中 `type_rank` 当前为 `bsj < backward < outward < forward`。行内 `r1_segments / r2_segments` 仍保持 read-chain order，不改成 genomic order。

### 5.2 `type` 枚举

`type` 的固定枚举值为：

- `bsj`
- `backward`
- `outward`
- `forward`

其中当前 `v1` 只实际输出：

- `bsj`
- `backward`
- `outward`

`forward` 保留给后续 circ span 内部普通线性 reads / linear-compatible reads；当前版本暂不输出。

### 5.3 `type` 与 mate-level BSJ 标记

pair-level `is_bsj` 不再单独输出，因为它和 `type` 完全重复：

- `type=bsj` 本身表示该 read pair 是 confirmed BSJ read；
- `type=backward` 表示非 confirmed BSJ、但 topology 呈 circular/backward 的 read；
- `type=outward` 表示非 confirmed BSJ、R1/R2 各自 linear-compatible，但 pair orientation 呈 outward circular-compatible 的 read；
- `type=forward` 为后续 linear-compatible reads 预留。

mate-level BSJ 只保留：

- `is_r1_bsj`
- `is_r2_bsj`

这两个字段表示对应 mate 的最佳 chain 是否跨过 BSJ / back-splice boundary，或是否被 post-Summary mate-level BSJ evidence 标记为 BSJ。

### 5.4 `circ_id` 约定

- `type=bsj`
  - 必须写唯一 `circ_id`
- `type=backward`
  - 当前写 `NA`
- `type=outward`
  - 当前写 `NA`，因为同一 read 可以在 circRNA-level graph 阶段兼容多个 circ span
- `type=forward`
  - 未来定义，当前不输出

`chrom / start / end` 的含义按 `type` 区分：

- `type=bsj`
  - `chrom / start / end` 来自 Summary-confirmed circRNA BSJ locus
- `type=backward`
  - `chrom` 来自 selected R1/R2 chains 的唯一 chromosome
  - `start / end` 是该 read pair 所有 retained alignment segments 的最小 / 最大 genomic position
- `type=outward`
  - `chrom` 来自 selected R1/R2 primary chains 的唯一 chromosome
  - `start / end` 是该 read pair 所有 retained alignment segments 的最小 / 最大 genomic position，不表示 candidate BSJ boundary，也不做 read-level circRNA 唯一归属
  - 这两个坐标只用于后续 genomic-position 排序和局部建图，不表示已确定 candidate BSJ boundary

### 5.5 `is_circular` 约定

`is_circular` 是基于当前 alignment-chain topology 的判定结果，而不是简单的来源标签：

- `1`
  - 该 read pair 的最佳解释体现出 circRNA topology
- `0`
  - 该 read pair 当前仍可被更合理地解释为线性 topology

当前阶段尤其关注：

- `type=bsj` 的 read 是否能稳定识别为 `is_circular=1`
- `type=backward` 的 read 是否能在不强行分配 `circ_id` 的情况下，被稳定识别为 `is_circular=1`

## 6. segment token 与 extended CIGAR 协议

`r1_segments` / `r2_segments` 直接对齐 simulator `.reads.tsv` 的协议，便于逐条比较。

单个 genomic segment token 格式为：

```text
start-end:strand
```

多个 segment 使用 `|` 分隔：

```text
10500-10550:+|10600-10698:+
```

如果跨过 BSJ，则必须在两个 BSJ 分区之间插入 `<bsj>`：

```text
10500-10550:+|<bsj>|10000-10098:+
```

固定规则：

- segment 顺序必须按 read-chain order 输出，不能为了坐标排序破坏 circRNA 拓扑；
- 对正链 BSJ read，典型链为 `circ-end-side|<bsj>|circ-start-side`；
- 对负链 read，先按 reverse-complement 后的 RNA/read-chain 方向解释，再输出 `strand=-`，但每个 token 内坐标仍保持 `start <= end`；
- supplementary 若属于同一最佳 chain，必须拼接进最终 segment chain；
- 对于 confirmed BSJ reads，Scan1/Scan2 sidecar 中的 `scan*_local` pseudo rows 可作为 local clip block 纳入最终 chain；
- 当前要输出的是类似 STAR 的 splice-aware chain，而不是把 BWA chimeric 记录原样逐条外泄。

### 6.1 `r1_cigar` / `r2_cigar`

`r1_cigar` / `r2_cigar` 是与 `r1_segments` / `r2_segments` 同步的 chain-level CIGAR。

它不是原始 SAM/BAM CIGAR，也不要求满足 SAM 规范；它的目标是把 CIRI-AS 已经解释过的 mate-level alignment chain 表达成一条 read-chain、splice-aware 的内部结构记录。

当前约定：

- `M` 表示一个输出 segment 覆盖的 genomic match block；
- `S` 表示 mate read 在当前输出 chain 两端没有纳入任何 segment 的 read bases，用于判断 terminal soft clip；
- `N` 表示 read-chain 中相邻两个 segment 之间的普通 skipped interval；
- `B` 是 CIRI-specific operation，表示 read-chain 中相邻两个分区之间跨过 BSJ/back-splice boundary；
- `N` / `B` 的长度都按相邻两个输出 segment 的 genomic interval gap 计算；若 read-chain 是高坐标到低坐标，则使用高坐标 block 的 `start` 与低坐标 block 的 `end` 之间的距离；
- 如果 splice/BSJ 校正后两个分区在 genomic 坐标上紧邻，允许输出 `0N` 或 `0B`，以保留显式边界；
- `S` 只在 chain 开头或结尾输出，不表示内部 splice / BSJ 结构，也不生成对应 segment token；
- CIGAR 和 segments 都必须按相同的 read-chain order 输出；
- `<bsj>` / `B` 的位置必须是 read-chain 中真实发生 circular wrap 的位置，不能把 chain 重排成 genomic order 后再插入；
- `type=backward` rows 也必须在触发 backward topology 的 read-chain wrap 处输出 `<bsj>` / `B`，但这只是 read-level topology marker，不代表 confirmed BSJ evidence；
- 对于已经完整 materialize 的 BSJ chain，`is_r*_bsj=1` 应该等价于对应 `r*_cigar` 中存在 `B`。

示例：

```text
r1_segments = 240-260:+|<bsj>|100-149:+|180-199:+
r1_cigar    = 5S21M90B50M30N20M4S
```

当前过渡期仍保留 `is_r1_bsj` / `is_r2_bsj`，原因是：

- 方便验证脚本在不解析 CIGAR 的情况下做快速统计；
- 兼容已有 simulator truth 和当前评估口径；
- 后续确认所有 BSJ-local alignment 都能 materialize 到 `B` 后，可以再讨论是否移除。

### 6.2 internal splice boundary 校正

最终 `<prefix>.segments` 写出前，会对 read-chain 相邻 segment 之间的普通 `N` junction 做边界校正：

- 只校正 internal splice boundary，不校正 `B` 对应的 BSJ 分区边界；
- annotation first：优先在小窗口内选择同一 gene、同一 RNA strand 的 exon end / exon start；
- read-specific validated junction hints 次之：复用 CIRI-AS candidate validation 已经校正过的 `(site2, site1)`；
- de novo splice signal fallback：没有 annotation / read-specific hint 时，根据 RNA strand 检查 canonical / semi-canonical splice dinucleotide；
- 普通 mapper block 保持保守窗口；validator 保存的 local clip block 或很短的 retained block 可对 annotation / read-specific hint 使用更宽窗口；
- `B` 对应的 circ boundary gap 不参与 internal correction，避免把 read-level segments 修正反向污染 BSJ 拓扑解释；
- 校正后的坐标同时写入 `r*_segments` 和 `r*_cigar`。

当前实现把这些信号合并进同一个候选排序，而不是命中 annotation 后立即停止：

- read-specific validated junction hint 优先级最高；
- transcript-consistent splice pair 高于普通 exon-boundary pair；
- confirmed BSJ reads 的 preliminary segment junction support 可作为第二轮 support-aware tie breaker；
- 第二轮不是全量重建：只有 preliminary chain 中某个 internal `N` junction 附近存在不同的 supported splice pair 时，才重新 materialize 该 read；其他 confirmed BSJ reads 直接复用第一轮结果；
- splice motif score 参与排序，但 motif-only 校正仍限制在保守窗口内；
- offset movement 只作为上述信号同级时的最后 tie breaker。

### 6.3 approximate clip rescue 的阶段边界

approximate local junction rescue 暂不进入当前 `<prefix>.segments` strong chain。

当前 read-level segments 只表达该 read 自身已经能高置信 materialize 的 alignment chain：

- mapper block、supplementary / secondary chain、以及 validator 已经定位到 circ 区间的 local clip block 可以进入 `r*_segments`；
- `<10bp` 的短片段或短 clip 不作为 strong exon / junction evidence；
- terminal soft clip 推断出的 approximate junction 不直接写成 `r*_segments` 中的新 segment；
- approximate clip 也不参与当前 strong `junction_chain` 评估口径。

原因是 approximate clip 本质上是 path-level 证据，而不是单 read level 的强模板。它需要在同一个 circRNA 内结合更多信息共同解释：

- 其他 BSJ reads 是否支持相同 junction；
- circRNA 内部 linear junction 是否支持相同 splice pair；
- annotation junction 是否支持该 exon boundary；
- splice motif / strand 是否一致；
- short clip 是否只提供 junction path hint，而不足以定义 exon boundary。

因此 approximate clip rescue 应后移到 circRNA-level full-length assembly 阶段。后续可以把这类 clip 序列作为 `boundary_clip_read` / weak junction path evidence 使用，但必须经过 circRNA-level graph 或 path voting 后，才允许影响 isoform/path 输出。

## 7. 第一阶段覆盖范围

当前阶段只输出两类明确属于 circRNA 解释域的 reads：

### 7.1 `type=bsj`

来源：

- Summary `junction_reads_ID` 中的 reads

要求：

- 必须挂到唯一 `circ_id`
- 必须输出最终 `r1_cigar / r1_segments / r2_cigar / r2_segments`
- 至少一个 mate 通常应有 `is_r*_bsj=1`

### 7.2 `type=backward`

来源：

- 当前 CIRI-AS 扫描里进入识别逻辑、并表现出 backward / circular topology 的 reads

要求：

- 当前先不强行分配 `circ_id / strand`，统一写 `NA`
- `chrom` 必须能由 selected R1/R2 chain 唯一确定；跨染色体或无法确定单一 chromosome 的 read 不进入最终 `<prefix>.segments`
- `start / end` 表示该 read pair 所有 retained alignment segments 覆盖到的最小 / 最大 genomic position，不表示 candidate BSJ boundary
- `is_circular=1`
- 触发 backward topology 的 mate chain 必须在对应 `r*_segments` 中写出 `<bsj>` marker，并在 `r*_cigar` 中用 `B` 标记该 read-order wrap；`is_r*_bsj` 仍保持 0，因为它不是 Summary-confirmed BSJ mate
- 必须输出最终 `r1_cigar / r1_segments / r2_cigar / r2_segments`
- `is_r1_bsj=0` 且 `is_r2_bsj=0`，因为该 read pair 不是 Summary-confirmed BSJ read

### 7.3 `type=outward`

来源：

- Summary 后额外扫描到的非 BSJ read pair；
- R1/R2 各自都能形成 linear-compatible primary chain；
- 两个 mate 的 primary alignment 呈 5' RO-like outward geometry：一个 mate 为 reverse、另一个为 forward，reverse mate 的 5' 端与 forward mate 的 5' 端存在 overlap，同时两个 3' 端向外展开；普通 proper pair 的 contained overlap 不属于 `outward`。

要求：

- 不强行分配 `circ_id`，统一写 `NA`；`strand` 默认写 `NA`，但如果 mate 内部已有普通 `N` junction 且 annotation / splice signal 给出唯一不冲突的 RNA strand，则写入该 strand，并用同一 strand 重新 materialize `r1_segments / r2_segments`；
- `chrom` 必须能由 selected R1/R2 primary chains 唯一确定；跨染色体或无法确定单一 chromosome 的 read 不进入最终 `<prefix>.segments`；
- `start / end` 表示该 read pair 所有 retained alignment segments 覆盖到的最小 / 最大 genomic position；
- `is_circular=1`，因为它是 pair-level circular-compatible topology support；
- 不写 `<bsj>` marker，也不在 `r*_cigar` 中写 `B`，因为该 read 没有明确的 mate-chain backward junction；
- `is_r1_bsj=0` 且 `is_r2_bsj=0`；
- graph completeness 评估中，`type=outward` 只能贡献自身 CIGAR 中已经存在的普通 `N` junction；没有 `N` 的 outward read 只作为 circRNA-level pair/path support，不能凭空生成 exon-exon junction。

### 7.4 当前暂缓

以下范围明确不进入 `v1`：

- circ span 内部普通线性 reads
- future `type=forward` / FSJ reads
- 因阈值或证据不足未被当前识别逻辑稳定接纳的弱支持 reads

## 8. alignment-chain 选择规则

每个 `read_id` 最多输出一行。

当前规则：

- 默认优先 `primary + supplementary`
- 如果这套链不能形成合法 circRNA 解释，则允许 `secondary` 进入候选
- supplementary 应该被拼接成同一条 splice-aware chain，而不是拆成独立输出行
- secondary 不是默认优先，只在 primary 路径不成立时才作为回退候选

“合法 chain”的判断必须同时满足：

1. 当前 topology 判定成立
2. 与当前 circ / BSJ 结构兼容
3. splice site 经过 annotation-guided + de novo splice signal 逻辑校正后仍成立

在所有合法 chain 中按以下顺序择优：

1. topology 更匹配
2. 更符合当前 circ / BSJ 结构
3. 更符合 annotation-guided + de novo splice signal 校正
4. splice offset 总量更小
5. insert size 更符合整体分布
6. 若仍并列，使用稳定 tie-break

推荐稳定 tie-break：

- `secondary` 使用更少
- `supplementary` 使用更少
- 最终 segment token 字典序更小

这个顺序的目标不是“最像 mapper 原始输出”，而是“最稳定、最符合当前 circRNA 解释”。

### 8.1 XA / SA tag 的阶段边界

`sim:433108` 暴露了一个重要问题：BWA `XA` tag 里可能包含比已经展开成 supplementary record 的位置更合理的替代比对。例如短的低 MAPQ supplementary block 被放到远端基因组位置时，会在 `type=backward` chain 中形成很长的假 `N`，而同一 read 的 `XA` 候选可能落在 circRNA 局部区域内。

源码确认后的当前边界：

- Java CIRI3 Scan1 / Scan2 不解析 SAM optional tags；它只从每条 SAM record 中取 `flag / chr / pos / MAPQ / CIGAR / SEQ`，因此没有显式处理 `XA`；
- Java CIRI3 也不解析 `SA:Z` 字符串；`SA` 的作用只来自 aligner 已经 materialize 成独立 SAM/BAM record 的 supplementary alignment；
- Rust CIRI3 parity 主流程保持同一边界，不在 Scan1 / Scan2 / Summary 中读取 aux tag；
- 当前 CIRI-AS / segments chain builder 会把 primary + supplementary record 作为优先链，并只在这条链不成立时考虑 secondary fallback，但它还看不到未展开的 `XA` alternative。

当前已在 post-Summary segments / CIRI-AS sidecar 中实现 XA-aware 修正与 negative filter：

- 将 `XA` 解析成低优先级 alternative alignment block，只用于 read-level sidecar chain 修正或拒绝弱 `type=backward`，不是 BSJ strong evidence；
- 对 selected chain 中 `MAPQ=0` 的原始 alignment record，包括 primary 和 supplementary，寻找 `XA` 中 read-coordinate 几乎一致、strand 一致、anchor 长度达到最小可信阈值的无缝替代；
- 多个 selected records 可以联合替换：单独替换 R1 或 R2 可能无法降低 pair-level span，但联合替换后能暴露更合理的 local linear / circular chain；
- 对所有无缝 `XA` 用同一套 topology-neutral ranking 比较 linear / circular 替代：query coverage 基本不降低、single chromosome、selected span 更短且优先回到 CIRI3 `max_span=200000` 范围内；
- 如果 best XA-derived chain 仍保持 backward topology，则用它修正 `type=backward` read-level chain；如果 best chain 变成 non-circular / linear-compatible，则拒绝该 `type=backward` read；
- 当前过滤暂不要求 candidate BSJ position / circ boundary 匹配，目标是先消除 `sim:433108`、`sim:432652`、`sim:432594`、`sim:4599148` 这类由重复 alignment hit 造成的假 backward，同时保留 `sim:3966179` 这类可以由本地 `XA` 修成合理 backward/BSJ-like chain 的 read；
- `XA` 的 `NM` 不作为 hard filter，也不参与排序。RNA editing / mutation 场景下 mismatch 不应被当成负证据；sidecar 只用 same-strand / same-read-slice / coverage / span 这些结构条件决定 read-level chain；
- `XA-only` 不得增加 `.out` junction read count，也不得改变 `.bsj1 / .bsj2 / Summary` 行为；
- 后续若要让 `XA` 进入 Scan1 / Scan2 主 BSJ 识别，必须先作为“超越 CIRI3”的独立阶段验证；在 CIRI3 parity 模式下，主流程继续不读取 aux tag。

已知 CIRI3 主流程缺陷需要留档：BWA-MEM 本身不做 splice-aware/circ-aware chain 选择，可能把短 supplementary block materialize 到远端重复位点，同时在 `XA` 中保留更合理的本地替代。例如 `sim:3966179` 的 R2 有 `chr1:90982 23M127H` 远端 supplementary，但同一 record 的 `XA:Z` 包含 `chr1,-203745224,22M128S,0`，该替代与 read slice 基本一致且更符合同一 read pair 的局部 spanning。Java CIRI3 不解析 `XA`，因此 Scan1 / Scan2 无法利用这个更合理解释；Rust 为保持 parity 暂不修改主流程，只在 `<prefix>.segments` 后处理阶段使用该信息修正 read-level chain。若后续该策略在 simulator 与真实数据上稳定，可考虑作为非 parity 模式的 Scan1 / Scan2 增强。

## 9. splice site 解释规则

当前阶段不应引入新的 splice 识别逻辑，而应继续复用现有 CIRI-AS 实现的解释规则：

- annotation-guided motif / strand / offset tie-break
- de novo canonical splice signal 检查
- circ-local candidate 聚类后再决定边界

换言之，`<prefix>.segments` 不是“原始 alignment dump”，而是：

- 已经过 circRNA-specific topology 过滤；
- 已经过 splice-aware chain 重组；
- 已经过 annotation + motif 校正后的 read-level 解释结果。

## 10. 验证口径

当前阶段验证必须优先对齐 simulator 的 `.reads.tsv`。

### 10.1 `type=bsj`

优先比较：

- `circ_id`
- `is_circular`
- `r1_cigar`
- `r1_segments`
- `r2_cigar`
- `r2_segments`
- `is_r1_bsj`
- `is_r2_bsj`

### 10.2 `type=backward`

优先比较：

- `is_circular`
- `r1_cigar`
- `r1_segments`
- `r2_cigar`
- `r2_segments`

当前不强制比较：

- `circ_id`

因为 `backward` 当前明确允许 `circ_id=NA`。

### 10.3 为什么先不比较 isoform

如果 `circ_id` 或 `is_circular` 已经错了，那么后续 isoform reconstruction 的错误只会被放大。  
因此当前阶段必须先把 read-level truth 做对，再谈 path-level truth。

## 11. 后续 full-length 方向中值得保留的内容

虽然“直接 full-length reconstruction”当前暂停，但旧方案里以下内容仍值得保留，作为下一阶段的设计补充。

### 11.1 保留的总体边界

- 仍然以 Summary confirmed BSJ 为硬锚点
- 仍然不修改 `priority=1` 主流程证据和 `.out` 判定结果
- 仍然不把 remap 作为第一依赖
- full-length 仍然应作为 `Summary` 之后的独立后处理阶段

### 11.2 保留的 evidence 思路

后续 path 层仍可继续使用以下 evidence 分类：

- `bsj_seed_read`
- `internal_split_read`
- `boundary_clip_read`
- `mate_link_read`
- `mate_orientation_anomalous`
- `coverage_read`

但这些 evidence 必须建立在 `<prefix>.segments` 已经稳定的前提上，再进一步汇总。

其中 `boundary_clip_read` 只作为后续 circRNA-level path assembly 的候选证据来源。short / approximate clip 不能在 read-level `<prefix>.segments` 阶段升级成 strong segment template；只有当同一 circRNA 内的其他 BSJ reads、内部 linear junction 或 annotation junction 共同支持同一 path 时，才可在 full-length assembly 中作为弱证据参与路径选择。

### 11.3 保留的 path 层约束

后续若重新启动 full-length reconstruction，仍建议保留：

- anchored path 与 candidate path 的证据分层
- circ-local graph，而不是 genome-wide novo transcript graph
- 不把 side-evidence-only boundary 直接升级成强 path 证据
- path 输出必须显式区分 `anchored` 与 `candidate`

### 11.4 保留的“暂不做”原则

在 segments 阶段和 future path 阶段都继续成立：

- 不修改 `.bsj` / `.out`
- 不依赖 `bwa mem` remap 作为第一前提
- 不从 sidecar 结果反向新增 circ seed
- 不做无 BSJ 锚点的 genome-wide circRNA 发现
- 不把 PSI / usage correction 提前到 read-level 基础仍不稳定的阶段

## 12. 推荐实现顺序

1. 去掉 `--as` / `--as-out`，让后处理默认运行
2. 在 Scan1 写 `<prefix>.bsj1` 的同时写 `<prefix>.segments1`
3. 在 Scan2 写 `<prefix>.bsj2` 的同时写 `<prefix>.segments2`
4. Summary 只读取 `.bsj1/.bsj2` 并输出 `<prefix>.out`
5. segments 阶段读取 `<prefix>.out + <prefix>.segments1/2`，生成 confirmed `type=bsj` rows
6. 用 simulator `.reads.tsv` 做 read-level 对照
7. Summary 后额外扫描非 BSJ read groups，补充 mate-chain wrap 型 `type=backward` rows
8. 同一扫描中补充 pair-orientation 型 `type=outward` rows；这些 rows 不写 `<bsj>` / `B`，只作为 circRNA-level graph support
9. 在 sidecar chain selection 中评估 XA-aware alternative alignment，只修复 read-level segments，不改变 CIRI3 parity 主流程
10. 后续再通过 circRNA-level region extraction 补 `forward` / internal linear reads，并重新评估 full-length path 层

---
最后更新：2026-05-07
