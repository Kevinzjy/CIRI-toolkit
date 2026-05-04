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
Scan1 -> Scan2 -> Summary
  -> second BAM/SAM sweep
  -> circ-related read collection
  -> best alignment-chain selection
  -> <prefix>.segments
  -> future full-length reconstruction
```

设计边界：

- 不修改 baseline `Scan1 -> Scan2 -> Summary` 的判定逻辑。
- 不改变 `priority=1` 主流程证据和 `.out` 的判定结果。
- 不再使用 `--as` / `--as-out` 作为额外入口。
- 主流程在写完 `<prefix>.out` 后默认继续执行 segments 后处理。
- 第一阶段只输出 read-level segments，不直接输出 `.isoforms`、`.isoform_summary`、`.fa`。
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
- 在 `Summary` 之后重扫原始 BAM/SAM；
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
- input BAM/SAM
  - 第二遍扫描读取原始 alignment evidence
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
- 第一阶段不再把 full-length 路径文件作为主输出协议。

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
is_bsj
r1_segments
r1_is_bsj
r2_segments
r2_is_bsj
```

### 5.2 `type` 枚举

`type` 的固定枚举值为：

- `bsj`
- `backward`
- `forward`

其中当前 `v1` 只实际输出：

- `bsj`
- `backward`

`forward` 仅作为未来兼容 linear / linear-compatible reads 的预留枚举值，当前版本不输出。

### 5.3 `type` 与 `is_bsj` 的固定关系

当前约定写死为：

- `type=bsj` 时，`is_bsj=1`
- `type=backward` 时，`is_bsj=0`
- `type=forward` 时，`is_bsj=0`

### 5.4 `circ_id` 约定

- `type=bsj`
  - 必须写唯一 `circ_id`
- `type=backward`
  - 当前写 `NA`
- `type=forward`
  - 未来定义，当前不输出

### 5.5 `is_circular` 约定

`is_circular` 是基于当前 alignment-chain topology 的判定结果，而不是简单的来源标签：

- `1`
  - 该 read pair 的最佳解释体现出 circRNA topology
- `0`
  - 该 read pair 当前仍可被更合理地解释为线性 topology

当前阶段尤其关注：

- `type=bsj` 的 read 是否能稳定识别为 `is_circular=1`
- `type=backward` 的 read 是否能在不强行分配 `circ_id` 的情况下，被稳定识别为 `is_circular=1`

## 6. segment token 协议

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
10000-10098:+|<bsj>|10500-10550:+
```

固定规则：

- segment 顺序必须按 genomic order 输出；
- 对负链 read，`strand=-`，但坐标仍保持 `start <= end`；
- supplementary 若属于同一最佳 chain，必须拼接进最终 segment chain；
- 当前要输出的是类似 STAR 的 splice-aware chain，而不是把 BWA chimeric 记录原样逐条外泄。

## 7. 第一阶段覆盖范围

当前阶段只输出两类明确属于 circRNA 解释域的 reads：

### 7.1 `type=bsj`

来源：

- Summary `junction_reads_ID` 中的 reads

要求：

- 必须挂到唯一 `circ_id`
- 必须输出最终 `r1_segments / r2_segments`
- `is_bsj=1`

### 7.2 `type=backward`

来源：

- 当前 CIRI-AS 扫描里进入识别逻辑、并表现出 backward / circular topology 的 reads

要求：

- 当前先不强行分配 `circ_id`，统一写 `NA`
- 必须输出最终 `r1_segments / r2_segments`
- `is_bsj=0`

### 7.3 当前暂缓

以下范围明确不进入 `v1`：

- circ span 内部普通线性 reads
- future `type=forward`
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
- `r1_segments`
- `r2_segments`
- `r1_is_bsj`
- `r2_is_bsj`

### 10.2 `type=backward`

优先比较：

- `is_circular`
- `r1_segments`
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
2. 复用当前 CIRI-AS second BAM/SAM sweep 与识别逻辑
3. 只收敛输出到 `<prefix>.segments`
4. 先实现 `type=bsj`
5. 再实现 `type=backward`
6. 用 simulator `.reads.tsv` 做 read-level 对照
7. 等 `circ_id / is_circular / segments` 稳定后，再重新评估 full-length path 层

---
最后更新：2026-05-03
