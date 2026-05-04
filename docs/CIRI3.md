# CIRI3 主流程拆解与 Rust 对齐笔记

本文档用于整理 `vendor/CIRI3` 的主流程实现思路，作为 CIRI-toolkit Rust 版本的核心留档。和 `docs/CIRI-AS.md`、`docs/CIRI-full.md` 一样，这份文档的目标不是复述用户手册，而是记录上游 Java 的模块边界、数据契约、关键启发式和当前 Rust 必须保留的 parity 约束。

## 0. 当前结论（2026-05）

当前 CIRI-toolkit 主流程已经完成 CIRI3 parity，对齐结论如下：

- circRNA-level：100%
- read ID-level：100%
- read-assignment-level：100%
- FSJ count-level：100%

这意味着本项目当前阶段对 `vendor/CIRI3` 的理解已经不再停留在“能跑出相近结果”，而是已经验证到：

- Scan1 候选发现逻辑已对齐；
- Scan2 rescue 与 FSJ 统计逻辑已对齐；
- Summary 合并与 stringency 过滤逻辑已对齐；
- 最终注释输出与列契约已对齐。

因此，本文档除了记录 Java 规范本身，也要记录这次实现过程中已经验证过的经验结论，避免后续重构、扩展或优化重新引入曾经排掉的问题。

## 1. 定位

在当前项目里，`vendor/CIRI3` 是唯一的主流程行为规范。

- `Scan1 -> Scan2 -> Summary` 的判定结果必须以 `vendor/CIRI3` 为准。
- `CIRI-AS` 与 `CIRI-full` 都是下游结构分析或 side evidence，不得反向改写主 CIRI3 的 circRNA 判定。
- Rust 可以重构代码结构、I/O 路径和性能实现，但不能改变 Java 已知的 read-level 决策顺序、过滤条件和输出契约。

推荐把几个上游文档的关系理解为：

```text
CIRI3      = 主 BSJ/FSJ 检测规范
CIRI-AS    = circ 内部结构与 splice evidence 参考
CIRI-full  = RO/remap/full-length 重构参考
```

本文档的覆盖范围限定为：

- BWA-MEM 口径下的主 CIRI3 检测流程；
- Rust 当前已经完成 parity 的主链行为；
- 后续维护时必须继续保留的实现约束。

本文档不打算完整覆盖：

- DE 子命令的统计流程；
- STAR 三文件路径的全部细节；
- 多样本/RNase-R wrapper 的工程性分发逻辑；
- CIRI-AS / CIRI-full / RO feature 的下游结构分析细节。

这些内容要么不是当前 Rust 功能主线，要么已经有独立文档承接。

## 2. 上游源码入口

`vendor/CIRI3` 中与主 circRNA 检测最相关的类主要是：

- `src/com/zx/test/TestParameters.java`
  - CLI 参数解析与运行模式分发入口。
- `src/com/zx/test/SingleTest.java`
  - 单样本 BWA-MEM SAM/BAM 主流程。
- `src/com/zx/findcircrna/FindCircRNAScan1.java`
  - 第一遍扫描，收集 BSJ1 候选。
- `src/com/zx/findcircrna/IsBSJScan1.java`
  - Scan1 核心 read/candidate 判定逻辑。
- `src/com/zx/findcircrna/FindCircRNAScan2.java`
  - 第二遍扫描，基于 Scan1 唯一位点索引补救 BSJ 并统计 FSJ。
- `src/com/zx/findcircrna/IsBSJScan2.java`
  - Scan2 候选遍历、FSJ 统计和二次判定逻辑。
- `src/com/zx/findcircrna/Summary.java`
  - 合并 BSJ1/BSJ2、按 stringency 过滤、生成 Summary circ 列表。
- `src/com/zx/hg38/Annotation.java`
  - 把 Summary 结果写成最终 `.out` 表。
- `src/com/zx/findcircrna/GetAnnotationInformation.java`
  - GTF/GFF 注释读取。
- `src/com/zx/findcircrna/ReadFaFile.java`
  - 全参考 FASTA 读取。
- `src/com/zx/findcircrna/Misd.java`
  - CIGAR 分类器。
- `src/com/zx/findcircrna/IsStand.java`
  - 从 flag 里提取当前比对链方向和 pair 内 segment 编号。

此外还存在与主流程平行的变体：

- `Bam*`：BAM 读取版本。
- `*STAR*`：STAR 三文件输入版本。
- `Mut*` / `File*` / `Tsv*`：多线程、多样本、RNase-R 和 user circ collection 版本。

当前 Rust 主流程主要对齐的是 BWA-MEM 的 `SingleTest` / `Bam*` 这条路径。

## 3. CLI 与运行模式

`TestParameters.java` 把 CIRI3 分成两大类入口：

1. circRNA 检测主流程
2. DE 子命令
   - `DE_BSJ`
   - `DE_Ratio`
   - `DE_Relative`

对 Rust 主流程最重要的是 circRNA 检测参数默认值：

- `-Max / --max_span = 200000`
- `-Min / --min_span = 140`
- `-S / --strigency = 2`
- `-U / --mapq_uni = 10`
- `-E / --rel_exp = 0`
- `-Mc / --mitochondria = 0`
- `-M / --chrM = chrM`
- `-T / --thread_num = 1`
- `-It / --intron = 0`
- `-Sp / --splicing_signals = 0`
- `-Ma / --mapper = 0`
- `-W / --way = 0`

其中需要特别注意：

- Java 使用拼写 `strigency`，不是 `stringency`。Rust 对外可以用更正常的拼写，但行为必须对齐 Java 的 `0/1/2` 三档规则。
- `-Ma 0` 表示 BWA-MEM，`-Ma 1` 表示 STAR。
- `-W 0/1/2` 分别表示单样本、样本列表、带 RNase-R 信息的样本列表。
- `-UC` 会绕过常规 Scan1 候选发现，改为“只对用户给定 circ 位点做 BSJ/FSJ 计数”。

## 4. 主流程总览

`SingleTest.CIRI3(...)` 的主线可以概括为：

```text
1. 读取注释
2. 读取参考 FASTA
3. Scan1: 从 SAM/BAM 中发现 BSJ1 候选
4. 用 Scan1 唯一 circ 位点建立 Scan2 索引
5. Scan2: 二次补救 BSJ，并统计 FSJ
6. Summary: 合并、去重、过滤
7. Annotation: 输出最终 circRNA 结果表
```

对应的中间数据流是：

```text
SAM/BAM
  -> <input>BSJ1
  -> chrCircSiteMap
  -> Scan2 index + circFSJMap
  -> SummaryCircList
  -> final .out
```

当前 Rust 对齐阶段反复验证后，可以把主流程的真正“不可随意重构”点再压缩成一句话：

```text
不是“先找 circ 位点，再补点计数”这么简单，而是“先由 read 形态决定候选，再由候选集合反过来约束 Scan2 和 Summary 的解释空间”。
```

很多看似局部的小改动，例如候选去重时机、桶内遍历方向、辅助序列的取向规则，都会顺着这条链影响最终 `.out`。

## 5. 注释与参考序列加载

### 5.1 注释

`GetAnnotationInformation.hand(...)` 会根据 `--intron` 分成两种模式：

- 非 intron 模式：
  - 记录 `chr + exon_start -> gene_id + strand`
  - 记录 `chr + exon_end -> gene_id + strand`
- intron 模式：
  - 记录 `chr + exon_start/end -> [gene_id, strand, transcript_id]`

同时两种模式都会维护：

- `geneExonMap`
  - `chr -> 按 gene_start 排序的 gene interval`
- `exonListMap`
  - `chr + gene_id -> exon 区间列表`

需要保留的 Java 特征：

- GTF/GFF 解析用正则直接抓 `gene_id` / `transcript_id`，不是通用 parser。
- `geneExonMap` 和 `exonListMap` 的构建依赖输入注释里同一 `gene_id` 的 exon 记录近似连续出现。
- 最后按 `gene_start` 排序，但 exon 本身不做复杂归并。

### 5.2 FASTA

`ReadFaFile.readFa(...)` 会把整条染色体序列读入内存并转成大写，保存为：

- `chrTCGAMap`
- `chrLenMap`

这是 Java 后续所有 motif、边界和 Smith-Waterman 合并的基础。Rust 当前也应保留“按染色体整段访问”的语义，而不是把 reference 查询变成带缓存抖动的随机 I/O。

## 6. Scan1：BSJ1 候选发现

### 6.1 read group 组织方式

`FindCircRNAScan1.findCircRNAScan1(...)` 按 query name 顺序读取 SAM/BAM，并为每个 read group 维护两套结构：

- `readsMap`
  - `segment_id -> 该 segment 的所有 alignment 片段`
- `standMap`
  - `segment_id -> alignment_strand + 原始 read sequence`

这里的 `segment_id` 来自 `IsStand.stand7(flag)`，本质上是从 flag 中取 pair 内 segment 位。

`standMap` 的前缀方向来自 `IsStand.stand5(flag)`，即当前 alignment 的 strand bit。Java 不是按 read1/read2 固定方向解释序列，而是按“当前比对链方向”逐条决定是否需要 reverse-complement。这个细节对 parity 很重要。

### 6.2 CIGAR 分类

`Misd.misd(cigar, seq_len)` 会把 CIGAR 映射成一个四元组，Rust 需要保留这套语义：

- `flag = -1`
  - `SM` 风格，左软截断。
- `flag = 1`
  - `MS` 风格，右软截断。
- `flag = 10`
  - `SMS` 风格，中间匹配、两端软截断。
- `flag = 0`
  - 连续匹配或普通 `M/D/I/M` 组合，不直接构成 BSJ 候选。
- `ref_len = -2`
  - 视为无效/不支持的 CIGAR。

Java 会先把 `H` 全部替换成 `S`，因此硬截断在主逻辑里被软截断语义吸收。Rust 如果直接保留 `H`，很容易造成 Scan1/Scan2 候选分叉。

### 6.3 BSJ1 判定骨架

`IsBSJScan1.isBSJScan1(...)` 会在每个 segment 的 alignment 之间两两配对，寻找 potential PCC signal。它主要处理两类模式：

1. two-segment
   - 典型 `SM` 与 `MS` 配对。
2. multiple-segment
   - 典型 `SMS` 与 `SM` / `MS` 配对。

核心步骤是：

1. 要求同染色体、同比对方向、CIGAR 有效，且至少一个片段 `MAPQ >= minMapqUni`。
2. 根据 CIGAR 类型计算 circ span，并要求：
   - `minCircle <= span <= maxCircle`
   - 两段关键 clipped 长度差不超过 6。
3. 计算 `end_adjustment1/2`，并要求绝对值不超过 4。
4. 从当前 read 的原始序列中抽取：
   - `str1`
   - `str2`
   - 必要时 `str3`
5. 若是 paired-end，再检查另一端是否完整落在候选 circ 区间内部，并生成：
   - `str4`
   - `str4_ok`
6. 把这些信息送给 `IsBSJHg1` / `IsBSJIntronHg1` 做 splice-signal 与 repeat 级别判断。

Java 返回值会直接写入 `<input>BSJ1`，格式是：

```text
read_id \t <CIGAR payload> \t <判定结果 payload>
```

Rust 当前持久化的 `<prefix>.bsj1` 在 Java payload 前增加 `mate_label` 和 `priority` 两列，用于后续 mate-level evidence 展示；`Summary` 会只消费 `priority=1` 行，并在内部转换回 Java payload 口径。

### 6.4 Scan1 的几个硬约束

- 候选发现按 queryname 分组顺序运行，默认读到新 read id 才结算上一组。
- Java 对一个 read group 只要找到第一个通过 `IsBSJHg1` 的候选就直接返回，不会继续选“更优候选”。
- `readLen` 取所有输入 read 中出现过的最大序列长度。
- `Mapped_Reads` 统计是“该 read group 至少有一个非 `*` CIGAR”的 reads 数，不是 alignment 数。

## 7. Scan2 索引：基于唯一 circ 位点而不是原始 BSJ1 行

`SingleTest` 在 Scan1 结束后会先读取 `<input>BSJ1`，构造：

- `scan1IdMap`
  - 所有已在 Scan1 命中的 read id。
- `chrCircSiteMap`
  - `chr -> 唯一 circ 位点集合`

这里的关键不是原始 BSJ1 行，而是按 `chr + start + end` 去重后的唯一 circ 位点。Rust 当前已经把这条规则写成硬约束，原因就在于 Java 的 Scan2 也是这样做的。

之后 Java 用 `seqLen = readLen - 12` 建索引，并生成：

- `circFSJMap`
  - 每个 circ 位点的 FSJ 计数初始值。
- `siteArrayMap1/siteArrayMap2`
  - 起点/终点桶门控数组。
- `chrSiteMap1/chrSiteMap2`
  - 桶内的候选位点列表，元素类型是 `SiteSort`。

`SiteSort.length` 实际承载的是候选 payload 数组。当前 Rust 需要保持与 Java 一致的字段语义：

```text
[site1, site2, strand, signal1, signal2, sum_q]
```

这也是 Scan2 不能直接拿原始 BSJ1 文本做索引的原因。

## 8. Scan2：候选补救与 FSJ 计数

### 8.1 输入与输出

`FindCircRNAScan2.findCircRNAScan2(...)` 会再次完整扫描输入文件：

- Java 如果 read id 已在 `scan1IdMap` 中，直接跳过，不再重复补救 BSJ。
- 如果 read id 已在 `scan1IdMap` 中，直接跳过，不再重复补救 BSJ。
- 否则调用 `IsBSJScan2.isCandidate(...)`。

Scan2 若成功补救出 BSJ，会把结果追加写回同一个 `<input>BSJ1` 文件。也就是说，Java 的 `BSJ1` 文件实际上同时承载了：

- Scan1 原始候选
- Scan2 rescued 候选

### 8.1.1 Rust mate-level BSJ 中间文件协议

Rust 的用户版 `.bsj` 需要显示具体是哪一端 mate 支持 BSJ，同时不能让这类展示信息反向改变 `Summary`。当前实现将中间文件收敛为两个文件：

- `<prefix>.bsj1`
  - Scan1 输出。
- `<prefix>.bsj2`
  - Scan2 输出。

两者都使用同一字段布局：

```text
read_id \t mate_label \t priority \t <CIGAR payload> \t <判定结果 payload>
```

字段语义：

- `mate_label`
  - `R1` 或 `R2`，表示这条 BSJ 证据来自哪一端 mate。
- `priority`
  - `1`：主流程证据。该行等价于旧 CIRI3/Rust parity 路径会交给 `Summary` 的 BSJ 行。
  - `0`：mate-level 附加证据。该行保留到最终用户版 `.bsj`，不参与 `Summary`、主 Scan2 candidate index 或主 read-count 统计。

`priority` 不是质量分，也不表示该 mate 一定比另一端更好。它只是内部路由字段，用来保证一个文件可以同时承载主流程输入和 mate-level 展示输入。

Scan1 写入规则：

- 对每个 read group，先按当前 CIRI3 first-hit 语义得到唯一主流程 BSJ1；这条写为 `priority=1`。
- 若同一 read group 的另一端 mate 也能独立通过 Scan1 BSJ 判定，额外写为 `priority=0`。
- 当 R1/R2 都是 BSJ1 时，哪一条是 `priority=1` 由 CIRI3-compatible first-hit 顺序决定，不做质量重选。

Scan2 写入规则：

- 主流程仍按 `read_id` 判断 Scan1 是否已 claim；只要 `.bsj1` 中存在该 read id 的 `priority=1` 行，Scan2 不再为 `Summary` rescue 该 read group。
- 对未被 Scan1 主流程 claim 的 read group，Scan2 first-hit rescue 行写为 `priority=1`。
- mate-level 展示路径可以继续检查未被同 mate claim 的另一端；额外命中的 mate 写为 `priority=0`。

消费者规则：

- `Scan2` 主 candidate index 只读取 `.bsj1` 中 `priority=1` 的唯一 circ 位点，保持 Java `chrCircSiteMap` 语义。
- `Summary` 按 `.bsj1 -> .bsj2` 的原始写入顺序读取，并且只消费 `priority=1` 行；解析时跳过 `mate_label` 与 `priority` 两列。
- 最终用户版 `<prefix>.bsj` 在 `Summary` 结束后由 `.bsj1 + .bsj2` 的全部行生成，保留 `priority` 字段，按 `read_id, mate_label, priority, scan_stage` 排序，并在尾列保留 `scan1/scan2` 来源。

最终 `.bsj` 保留 `priority=0` 是有意设计：后续内部 splice site 识别需要看到同一 read pair 的另一端 mate-level BSJ evidence。`priority=0` 行不能增加 `Summary` junction read count，但可用于 mate consistency、内部 splice 辅助证据、边界解释和 R1/R2 discordant case 标记。任何后续结构识别模块若用 `.bsj` 计数，都必须显式区分 read-pair-level 计数和 mate-level evidence 计数。

该设计已经替代早期的额外展示临时文件：当前不再单独维护 `<prefix>.scan1.tmp`、`<prefix>.scan2.tmp` 或 `<prefix>.bsj.raw.tmp`，但仍保持主流程和 mate-level 展示输出的逻辑隔离。

### 8.2 FSJ 统计

`IsBSJScan2` 在尝试 rescue 之前，还会通过 `GetFSJClass` 在 circ 区间内统计 forward splice junction 证据，并累加到 `circFSJNewMap`。

这部分逻辑与 BSJ rescue 共享同一套候选索引，因此 bucket gating 和位点遍历顺序都不能随便改。

### 8.3 候选遍历顺序

对单端 clipped 候选，Java 会按两种 anchor 分开查：

- `SM` 风格：
  - 用对齐起点附近去查 `site1`。
- `MS` 风格：
  - 用对齐终点附近去查 `site2`。

并且每次都分成两个桶：

- `num1`
- `num2`

需要保留的顺序规则是：

- `num1` 桶逆序遍历。
- `num2` 桶正序遍历。

这直接影响“先命中哪个候选就返回哪个候选”的行为，是当前 Rust parity 的关键点之一。

### 8.4 Scan2 tag 语义

`IsBSJHg2.isBSJHg2(...)` 返回的 tag 在 Scan2 里有三类语义：

- `"0"`
  - 不是 BSJ，但可计为某个 circ 的 FSJ 证据。
- `"2"`
  - 无效候选，继续看下一个候选。
- 其他非 `2` tag
  - 视为有效 BSJ，立即返回。

因此 Java 的正常语义不是“穷举所有候选后择优”，而是：

```text
遇到首个有效非 2 tag -> 立即返回
```

Rust 必须保留这一点。

## 9. Summary：合并、去重、过滤

`Summary.summary(...)` 会读取 BSJ 候选行，把候选按 `chr/start/end` 聚合，并执行三层处理。Java 读取的是 Scan1 加 Scan2 追加后的 `BSJ` 文件；Rust 当前按 `.bsj1 -> .bsj2` 顺序读取，并且只消费 `priority=1` 行，在进入 Summary merge/stringency 前转换回 Java payload 口径。

### 9.1 以 `chr/start/end` 为主键聚合 read

`circMap`：

- key: `chr \t start \t end`
- value: 支持该 circ 的原始候选行集合

同时建立：

- `circStartMap`
- `circEndMap`

分别按共享 start 或共享 end 对 circ 分组。

### 9.2 start-shared / end-shared 的 Smith-Waterman 合并

如果同一 start 或同一 end 下有多个 circ，Java 会进一步比较支持性更强和支持性更弱的候选：

- 抽取边界附近序列
- 用 `SmithWaterman` 比较
- 达到阈值后把弱候选合并进强候选

这是 Java Summary 消除近邻重复 circ 的重要步骤。当前 Rust 中凡是看起来“多余”的局部对齐和容器遍历，通常都是为了保留这一行为。

### 9.3 TP / non / FP 计数与 stringency

每个 circ 最终会统计：

- `TPReads`
- `nonReads`
- `FPReads`
- `TPReads3`
- `nonReads3`
- `FPReads3`
- `CIGARCount3 = [SM, MS, SMS]`
- `tag`

然后按 `strigency` 过滤：

- `strigency == 2`
  - 要求更严格的 distinct PCC signal 约束。
- `strigency == 1`
  - 放宽到主要看 junction reads。
- `strigency == 0`
  - 输出所有满足基本 read evidence 的 circ。

Java 代码里的精确规则是：

```text
S=2:
((TPReads > 19 * FPReads || FPReads <= 1)
 && TPReads > nonReads + FPReads
 && CIGARSet.size() >= 3
 && TPReads >= 2)
||
(tag > 0 && falseCIGARSet3.size() == 0 && CIGARSet3.size() >= 3 && TPReads3 >= 2)

S=1:
((TPReads > 19 * FPReads || falseCIGARSet.size() <= 2)
 && TPReads > nonReads + FPReads
 && TPReads >= 2)
||
(tag > 0 && falseCIGARSet3.size() == 0 && TPReads3 >= 2)

S=0:
((TPReads > 19 * FPReads || falseCIGARSet.size() <= 2)
 && TPReads > nonReads + FPReads)
||
(tag > 0 && falseCIGARSet3.size() == 0 && TPReads3 >= 2)
```

这些条件必须按 Java 原式对齐，不能改写成“看起来等价”的更现代逻辑。

## 10. Annotation：最终结果表

`Summary` 产出的 `SummaryCircList` 只是中间表，最终输出由 `Annotation.java` 或 `AnnotationIntron.java` 完成。

非 intron 模式的标准输出 header 是：

```text
circRNA_ID
chr
circRNA_start
circRNA_end
#junction_reads
SM_MS_SMS
#non_junction_reads
junction_reads_ratio
circRNA_type
gene_id
strand
junction_reads_ID
Score
```

`Annotation.java` 会根据：

- 起点/终点是否正好落在 exon boundary
- 是否落在同一基因区间
- 是否落在 exon 内

把 circ 标成：

- `exon`
- `intron`
- `intergenic_region`
- `NA`

同时保留 `junction_reads_ID` 和 Summary 里的分数/标签信息。当前 Rust 主流程输出 13 列，也是沿着这里对齐的。

## 11. 当前 Rust 必须保留的 Java 细节

下面这些是已经反复证明会影响 parity 的点：

- Scan2 索引必须来自唯一 circ 位点集合，也就是 Java 的 `chrCircSiteMap` 语义。
- Scan2 candidate payload 字段布局必须保持 `[site1, site2, strand, signal1, signal2, sum_q]`。
- Scan2 桶门控遍历顺序必须是：
  - `num1` 逆序
  - `num2` 正序
- Scan2 正常模式保持“首个有效非 `2` tag 即返回”。
- read 的序列方向必须按当前 alignment strand 计算，不能按 pair 固定方向复用。
- `H` 必须先按 `S` 处理，再进入 CIGAR 分类。
- `seqLen` 的 Scan2 索引尺度是 `readLen - 12`，不是 readLen 本身。
- `Summary` 的 stringency 条件必须按 Java 原式保留。
- FASTA 必须统一大写，否则 motif/比对结果会漂移。
- 注释层不只是加 gene name，而是会改变 `circRNA_type` 和最终行是否保留。

## 12. 已验证的实现经验总结

下面这些经验不是抽象建议，而是本次完成 100% parity 过程中已经验证过、确实会影响结果或维护成本的结论。

### 12.1 先定 spec，再定实现

- `vendor/CIRI3` 的 Java 代码本身才是规范，论文描述和 README 只能作为背景。
- 如果 Java 代码、文档说明和直觉理解不一致，以 Java 实际行为为准。
- 遇到“这段 Java 看起来不够优雅”的代码时，默认先假设它承载了行为约束，而不是先重写成更 Rust 风格。

### 12.2 不要把“无序容器”当作自由优化空间

- Java 里的 `HashMap` / `HashSet` / 分组后首命中行为，很多时候会间接影响 Summary 合并和 Scan2 early-return。
- Rust 即使不必逐字复刻 Java 容器，也必须把“最终观测到的遍历语义”恢复出来。
- 一旦某段逻辑依赖“先遇到哪个候选就停”，遍历方向本身就是行为，不再只是实现细节。

### 12.3 输入上下文本身就是算法的一部分

- Scan2 不是只看当前 read，还看由 Scan1 唯一位点集合构成的候选空间。
- Summary 不是只汇总 surviving circ，还依赖同 family 的竞争候选与局部边界序列。
- 因此 subset 复现失败时，不能只怀疑单条 read 的判定函数，往往是中间上下文被删掉了。

### 12.4 注释版本和参考版本必须锁定

- whole-genome parity 不是“只要算法一样就一定一致”。
- FASTA 版本不同、GTF 版本不同、甚至某个位点 exon 覆盖是否存在，都会制造表面上的 parity gap。
- 当前 whole-genome 推荐复核口径是：
  - FASTA：`/data/public/database/gencode/hg38/_BWAindex/hg38.fa`
  - GTF：`/data/public/database/gencode/hg38/gencode.v44.annotation.gtf`

换句话说，环境锁定本身就是 parity 工作的一部分，不是外围杂务。

### 12.5 方向语义必须分层

- 服务主 CIRI3 parity 的 `standMap` 语义是 alignment-oriented，不是原始测序方向。
- 任何需要“原始 read1/read2 测序方向”的扩展功能，都不能复用这套语义偷懒。
- 这条经验对后续 RO、full-length、sidecar 模块尤其重要，否则很容易一边修扩展一边破坏主流程。

### 12.6 parity 调试必须分层做

- 先看 circRNA-level。
- 再看 read ID-level。
- 最后看 read-assignment-level 和 FSJ。

如果一开始就对着最终 `.out` 逐行肉眼 diff，效率很低，也容易把 Summary 问题误判成 Scan1/Scan2 问题。

### 12.7 release 才是有效性能与真实流程口径

- 真实 SAM/BAM fixture、whole-genome parity 和性能观察都应使用 `--release`。
- debug build 适合做局部单元调试，不适合拿来判断主流程耗时、候选密度或整体资源行为。
- 性能优化必须建立在 parity 已锁定的前提下进行，否则很难区分“变快了”还是“少做了事”。

### 12.8 热路径优化要优先做低风险、可验证的改动

已验证有效的优化方式通常是：

- 减少热循环分配；
- 复用 buffer；
- 延迟或避免中间字符串物化；
- 让 SAM/BAM 两条读取路径收敛到同一判定语义；
- 在不改变输出的前提下释放大工作集，控制 RSS。

而下面这类改动风险更高：

- 改候选集合语义；
- 改遍历顺序；
- 改 read 分组边界；
- 改 Summary 聚合条件；
- 把启发式“抽象重写”为自以为等价的新逻辑。

### 12.9 扩展能力必须 sidecar 隔离

- CIRI-AS、CIRI-full、RO、full-length reconstruction 都只能建立在 confirmed BSJ 之上。
- 它们可以消费 `.out`、最终 `.bsj` 的 mate-level evidence、原始 BAM/SAM 和额外 evidence，但默认不得改动 `priority=1` 主流程证据和 `.out` 判定结果。
- 这是本项目能在保持主流程稳定的同时继续扩展功能的根本前提。

## 13. Rust 实现时的边界建议

为了后续维护清晰，建议继续把 Rust 主流程边界固定为：

```text
annotation/fasta load
-> Scan1
-> Scan2
-> Summary
-> final annotation output
```

并把以下能力明确隔离在主流程之外：

- CIRI-AS sidecar
- CIRI-full / RO remap
- full-length isoform reconstruction
- 任何新的 re-scoring / ranking 试验

理由很直接：

- `vendor/CIRI3` 的行为规范已经足够复杂。
- 主流程 parity 已完成，不应被下游扩展反向污染。
- 一旦 `priority=1` 主证据或 `.out` 结果变化，就必须先按 CIRI3 口径解释，而不是用扩展模块兜底。

## 14. 文档完备性检查

从当前维护需求看，这份文档已经覆盖了主流程留档最需要的内容：

- 规范来源：哪些 Java 类才是主流程 spec。
- 参数口径：哪些默认值必须与 CIRI3 保持一致。
- 执行顺序：`Scan1 -> Scan2 -> Summary -> Annotation` 如何串起来。
- 数据契约：哪些中间结构与字段布局不能乱改。
- 判定细节：哪些局部启发式真正影响 read-level 结果。
- parity 经验：哪些坑已经踩过、以后不该再踩。
- 维护边界：哪些扩展必须隔离在主流程之外。

仍然刻意没有展开的部分有：

- STAR 路径的全部实现细节；
- 多样本与 RNase-R 模式的 wrapper 差异；
- DE 子命令；
- CIRI-AS / CIRI-full 下游流程。

这不是缺漏，而是刻意控制作用域。对当前 Rust 主线维护来说，继续把这份文档聚焦在“主 CIRI3 parity 规范”上，比把所有上游模块都塞进来更有用。

## 15. 结论

`vendor/CIRI3` 的核心不是一个“先找 BSJ 再数数”的简单流程，而是一套强依赖顺序、容器语义、字符串协议和局部启发式的两遍扫描系统：

- Scan1 负责发现高置信 BSJ read 形态；
- Scan2 负责在唯一 circ 位点集合上做补救和 FSJ 计数；
- Summary 负责近邻合并、read 级 evidence 汇总和 stringency 过滤；
- Annotation 负责输出最终 circ 类型与 gene/strand 注释。

对 Rust 而言，真正的规范不是某个论文描述，而是 `vendor/CIRI3` 这套 Java 实现本身。现在既然已经完成 100% parity，这份文档的职责就不仅是解释 Java 做了什么，更是明确说明：

- 哪些行为已经被验证为必须保留；
- 哪些“优化”其实会破坏 parity；
- 哪些扩展必须严格隔离在主流程之外。

后续任何重构、并行化、I/O 优化或 sidecar 扩展，都应继续以这份主流程留档为参照。
