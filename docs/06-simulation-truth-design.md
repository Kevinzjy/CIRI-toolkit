# `ciri-simulator` answer 输出设计

本文档定义 CIRI-toolkit 中独立模拟器 `ciri-simulator` 的输入、随机采样模型、FASTQ 输出和结构化 answer TSV 协议。模拟器直接读取 FASTA/GTF，生成压缩 paired FASTQ，并输出 isoform 层与 read-pair 层的标准答案，服务后续 CIRI-AS-style internal structure / full-length isoform reconstruction 验证。

## 1. 设计目标

`ciri-simulator` 的 answer 协议按 full-length reconstruction 的验证需求设计：

- 两张答案文件承担不同职责：circRNA/isoform 汇总一张表，read pair 逐条来源一张表。
- circRNA 层一行一个 circRNA，直接包含该 circ 下所有模拟 isoform 结构与 reads/BSJ reads 计数。
- read 层一行一个模拟 read pair，记录对应 circRNA、isoform、是否具有 circRNA 证据、是否跨 BSJ，以及 R1/R2 的真实基因组来源片段。
- 所有随机选择由 `--seed` 控制，便于复现。
- circRNA locus、exon window、coverage 和 fragment offset 都默认随机抽样，尽量避免只覆盖固定构造。
- answer 文件聚焦 circ locus、isoform exon chain、isoform read count 和 read pair 的真实来源片段；sequence-overlap / RO 评估使用独立 debug 输出协议。

## 2. 工具入口

模拟器的正式开发入口是 `ciri-simulator`，薄 wrapper 位于
`src/bin/ciri-simulator.rs`，核心实现位于 `src/simulator/`：

```bash
cargo run --bin ciri-simulator -- \
  -r <reference.fa> \
  -a <annotation.gtf> \
  -o tmp/sim \
  --circ-count 100 \
  --circ-coverage 10 \
  --linear-coverage 0.1 \
  --scale 0.5 \
  --read-len 150 \
  --insert-len 260 \
  --insert-sd 40 \
  --seed 5
```

默认不设置 `--chrom`，因此会在 GTF/FASTA 共有的所有染色体上选择候选 transcript。需要复现小型 chr1 fixture 时再显式传 `--chrom chr1`。

核心参数：

```text
--chrom                   可选染色体过滤；默认不限制
--circ-count              模拟 circRNA 数量，默认 100
--circ-coverage           circRNA coverage 均值，默认 10
--linear-coverage         线性 transcript 背景 coverage 均值，默认 0.1
--scale                   coverage Gaussian 相对标准差，默认 0.5；设为 0 时固定使用均值
--read-len                普通 paired-end read 长度，默认 150
--insert-len              主要 fragment span 均值，默认 260
--insert-sd               主要 fragment span 标准差，默认 40
--insert-len-minor        minor 长片段分布均值，默认 420
--insert-sd-minor         minor 长片段分布标准差，默认 60
--minor-insert-fraction   minor 长片段分布比例，默认 0.10
--error-rate              普通 read substitution rate，默认 0.002
--exon-exclusive-rate     每个去冗余 exon 坐标被标记为 circRNA-exclusive 的概率，范围 [0, 1]，默认 0.25
```

默认 PE150 模型：

```text
read_len = 150
insert_len ~ 90% N(260, 40) + 10% N(420, 60)
ordinary read substitution rate = 0.002
sampled_circ_coverage = max(0, N(circ_coverage, (circ_coverage * scale)^2))
sampled_linear_coverage = max(0, N(linear_coverage, (linear_coverage * scale)^2))
circ read pairs = max(1, round(circ_len * sampled_circ_coverage / (2 * read_len)))
linear read pairs = round(summed transcript exon length * sampled_linear_coverage / (2 * read_len))
sequence-overlap side metadata = optional check over sampled circ read pairs
```

所有 truth TSV 中的 genomic coordinates 使用 GTF 风格的 `1-based, closed interval`。Rust 内部取序列、切片和长度计算可以转换为 0-based half-open，但写出 `.isoforms.tsv`、`.reads.tsv` 和 `.annotation.gtf` 时必须回到 1-based closed coordinates。

read pair 数量整数化规则固定为：

- circRNA read pairs 使用 `max(1, round(circ_len * sampled_circ_coverage / (2 * read_len)))`，确保每个模拟 circRNA 至少有一个 circular read pair。
- linear background read pairs 使用 `round(summed transcript exon length * sampled_linear_coverage / (2 * read_len))`。

### 2.1 Linear annotation sampling model

`ciri-simulator` 使用输入 GTF 生成一个面向 CIRI 识别流程的 linear annotation：

```text
<prefix>.annotation.gtf
```

该文件表示模拟数据中的线性 RNA 注释空间；它由输入 GTF 过滤得到，并用于 linear background read simulation。circRNA simulation 仍使用完整输入 GTF，因此 circRNA 可以来自输入注释中的任意 exon。

`--exon-exclusive-rate` 控制 circRNA-exclusive exon sampling，取值范围为 `[0, 1]`。sampling 单位是去冗余后的 exon 坐标，而不是原始 GTF exon records。exon identity 定义为：

```text
(chrom, start, end, strand)
```

处理流程：

1. 从输入 GTF 收集所有 exon records。
2. 按 `(chrom, start, end, strand)` 去重，得到唯一 exon 坐标集合。
3. 对每个唯一 exon 坐标生成一个由 `--seed` 控制的 `[0, 1]` 随机数；随机数 `<= exon_exclusive_rate` 时，该坐标标记为 circRNA-exclusive；随机数 `> exon_exclusive_rate` 时，该坐标保留在线性注释空间。
4. 写出 `<prefix>.annotation.gtf` 时，删除所有坐标命中 exclusive set 的 exon records；如果多个 transcript 或 gene 共享同一 exon 坐标，这些 records 会一起从 linear annotation 中移除。
5. `<prefix>.annotation.gtf` 只保留 `gene`、`transcript` 和 `exon` records。过滤后没有任何非 exclusive exon 的 transcript 不保留；所有 transcript 都被删除的 gene 也不保留。输出 record 顺序必须保持 raw input GTF 中 surviving records 的原始顺序。
6. internal exon 被移除后的 transcript 可以继续作为 linear simulation candidate，但只能使用剩余的非 exclusive exons。linear background reads 不得使用任何 circRNA-exclusive exon。
7. circRNA reads 继续从完整输入 GTF 采样；被抽为 circRNA-exclusive 的 exon 不保证一定出现在最终模拟出的 circRNA isoforms 中。

如果 transcript 在过滤后没有可用 exon，或剩余 exon chain 已不足以模拟有效 linear fragment，该 transcript 不进入 linear background read candidate pool。若过滤后没有任何可用 linear candidate，模拟器应给出明确错误，而不是退回完整 GTF。

circRNA simulation 产生的一部分 reads 可能不包含 BSJ 或 outward-facing pair 特征，因此即使它们覆盖 circRNA-exclusive exon，也可能无法仅从比对结果确定其 circRNA 来源。这种覆盖是模拟设计的一部分，不需要由 linear background reads 补足。

这个设计让 `.isoforms.tsv` 和 `.reads.tsv` 保持自明：truth 文件直接记录实际模拟出的 circRNA isoforms、read 来源和 BSJ 标签；评估时优先使用这两张 truth 表判断 circRNA BSJ reads 与 isoform 识别是否准确。

在 **circular reads 的采样阶段**，短环使用 full-length 友好的 insert bias：

- 当 `circ_len <= 300` 时，约 70% 的片段长度会围绕 `circ_len` 本身做更窄的 Gaussian 采样；
- 该 bias 只作用于 circRNA reads，linear background 仍使用普通 insert 分布；
- 该 bias 用于提高短环 full-length / closed path 验证样本中覆盖完整环长的正例比例。

上述默认值用于生成规模适中的 PE150 风格 fixture；大规模 benchmark 可提高 `--circ-count`、`--circ-coverage` 和 `--linear-coverage`。

## 3. 输出文件

当 `-o <prefix>` 时，模拟器输出 5 个正式文件：

```text
<prefix>_1.fq.gz
<prefix>_2.fq.gz
<prefix>.annotation.gtf
<prefix>.isoforms.tsv
<prefix>.reads.tsv
```

其中：

- `_1.fq.gz` / `_2.fq.gz` 是面向 `ciri` / aligner 的真实 paired FASTQ 输入。
- `.annotation.gtf` 是由输入 GTF 过滤得到的 linear RNA annotation，作为后续 CIRI 识别使用的注释输入。
- `.isoforms.tsv` 是 circRNA / isoform 层标准答案。
- `.reads.tsv` 是 read-pair 层标准答案。

主 paired FASTQ 默认压缩为 `.fq.gz`：运行开始时先检查压缩工具，优先使用 `pigz -n -c`；如果 `pigz` 不存在但 `gzip` 可用，则输出提示并 fallback 到 `gzip -n -c`。如果二者都不存在，模拟器在生成 reads 前直接报错。

模拟结束时标准错误输出摘要：

```text
Simulated <circ_count> circRNAs, <isoform_count> circular isoforms
Total <read_pair_count> read pairs, <circ_read_pair_count> circRNA read pairs, <linear_read_pair_count> linear read pairs
BSJ feature: <junction_spanning_mate_count> reads / <bsj_read_pair_count> read pairs
```

这里 `Total/circRNA/linear` 的计数口径是 paired-end read pair；`BSJ feature` 的第一个数字是 individual mate 级别，只统计真正跨过 circular boundary 的 R1/R2，第二个数字是至少一个 mate 跨过 BSJ 的 read pair 数。来自 circRNA 的普通 read pair 不会自动计为 BSJ feature。

## 4. `<prefix>.isoforms.tsv`

`.isoforms.tsv` 一行一个模拟 circRNA。目标是直接回答：

- 模拟了哪些 circRNA；
- 每个 circRNA 有几种 isoform；
- 每种 isoform 的 exon chain 是什么；
- 每种 isoform 模拟出多少 read pairs；
- 每种 isoform 模拟出多少 BSJ read pairs。

字段顺序：

```text
circ_id
chrom
start
end
strand
gene_id
transcript_id
coverage
read_cnt
bsj_read_cnt
isoform_cnt
isoform_exons
isoform_len
isoform_read_cnt
isoform_bsj_read_cnt
```

字段说明：

- `circ_id` 使用 `chr:start|end`，其中 `start/end` 是该模拟 circRNA exon window 的 genomic span。
- `coverage` 是该 circRNA 实际抽样到的 coverage 值。
- `read_cnt` 是该 circRNA 下所有 circular read pairs 数。
- `bsj_read_cnt` 是该 circRNA 下至少一个 mate 跨 BSJ 的 read pairs 数。
- `isoform_cnt` 是该 circRNA 下模拟出的 isoform 数。

### 4.1 isoform compact fields

`isoform_exons`、`isoform_len`、`isoform_read_cnt`、`isoform_bsj_read_cnt` 都采用分号分隔的 `key=value` 列表，key 为 isoform ID。

示例：

```text
isoform_exons =
isoform1=100-200:+,300-400:+;isoform2=100-200:+,500-600:+

isoform_len =
isoform1=202;isoform2=202

isoform_read_cnt =
isoform1=120;isoform2=80

isoform_bsj_read_cnt =
isoform1=35;isoform2=20
```

`isoform_exons` 中单个 exon token 的格式为：

```text
start-end:strand
```

同一 isoform 内的 exon token 使用 `,` 分隔，不同 isoform 的 `key=value` 使用 `;` 分隔。这里不使用 `|`，因为 `|` 保留给 `circ_id` 的 `start|end` 和 `.reads.tsv` 的 segment chain。

如果某个 exon 坐标被标记为 circRNA-exclusive，并且出现在实际模拟出的 circRNA isoform 中，`isoform_exons` 可以在该 exon token 后追加 `*`：

```text
isoform1=100-200:+*,300-400:+
```

`*` 不是新增 truth 字段，只是 `isoform_exons` 内的结构标记，用于快速查看该 exon 不存在于 `<prefix>.annotation.gtf`。下游解析 exon 坐标时应先去掉可选的 `*` 后缀。

### 4.2 isoform generation model

初始实现使用固定、可测试的简单 isoform 模型：

- 每个 circRNA 至少生成一个 isoform，`isoform1` 使用该 circRNA 抽样到的连续 exon window。
- 如果 exon window 至少包含 3 个 exons，且该 circRNA 的 `read_cnt >= 2`，则以固定概率 `0.5` 生成 `isoform2`；`isoform2` 从同一 circRNA exon window 中均匀随机跳过一个内部 exon。
- `isoform2` 只能跳过内部 exon，不能改变 circRNA 的 `start/end` back-splice boundary。
- 多 isoform 的随机选择、跳过哪个内部 exon、每个 isoform 的 read allocation 都由 `--seed` 控制，必须可复现。
- circRNA read pairs 在该 circRNA 的 isoforms 之间均匀随机分配；生成出的每个 isoform 至少分配 1 个 read pair，因此不会输出没有 read 支持的 isoform。

注意：

- 同一个 circRNA 下的 isoform 共享 `circ_id`，但 exon chain 可以不同。
- `isoform1` 使用一个连续 exon window。
- `isoform2` 在满足生成条件时跳过一个内部 exon，用于多 isoform reconstruction 测试。
- `isoform_read_cnt` 与 `isoform_bsj_read_cnt` 必须能按 isoform ID 与 `.reads.tsv` 汇总结果对齐。

## 5. `<prefix>.reads.tsv`

`.reads.tsv` 一行一个模拟 read pair。它记录 read pair 的真实来源，以及该 read pair 在判定层面是否应视为 circRNA read、是否包含 BSJ read。

这个文件同时覆盖 circRNA reads 与 linear background reads。linear background reads 的 `circ_id`、`start`、`end` 和 `isoform_id` 使用 `NA`，`chrom` 与 `strand` 记录该 read pair 的线性来源坐标系，`is_circular=0`，`is_bsj=0`。

字段顺序：

```text
read_id
circ_id
chrom
start
end
strand
isoform_id
is_circular
is_bsj
r1_segments
r1_is_bsj
r2_segments
r2_is_bsj
type
```

字段说明：

- `read_id` 不带 `/1` 或 `/2`，对应 FASTQ 中同一 read pair 的 shared ID。
- 对 circRNA reads，`circ_id`、`chrom`、`start`、`end`、`strand` 和 `isoform_id` 指明该 read pair 的真实 circRNA / isoform 来源，并与 `.isoforms.tsv` 中对应 circRNA 的字段一致。
- 对 linear background reads，`circ_id`、`start`、`end` 和 `isoform_id` 使用 `NA`；`chrom` 保留线性来源 chromosome，`strand` 保留线性来源 strand，使不带 chromosome 的 segment token 仍可独立还原 genomic interval。
- `is_circular=1` 是 observable topology rule，而不是来源标签；只有该 read pair 的 mate 自身跨 BSJ，或两个 mate 构成 outward-facing circular-compatible pair 时才标记为 1。
- `is_circular=0` 表示该 read pair 的来源片段可以被解释为完美线性 reads；即使该 read pair 实际来自 circRNA，只要 topology 上可线性解释，也标记为 0。
- `is_bsj=1` 表示 R1 或 R2 至少一个 mate 跨过 BSJ 位点。
- `is_bsj=0` 表示两个 mate 都不跨 BSJ 位点。
- `r1_is_bsj` / `r2_is_bsj` 分别表示 R1 / R2 是否为 BSJ read。
- `r1_segments` / `r2_segments` 是该 mate 的真实模拟来源 read-chain 坐标片段。
- `type` 是 read-pair 的主 truth 类型，当前取值为 `bsj` / `outward` / `forward`：
  - `bsj` 表示至少一个 mate 真实跨 BSJ；
  - `outward` 表示两个 mate 都不跨 BSJ，但 primary pair 呈严格 outward-facing geometry，即 reverse-oriented mate 位于 forward-oriented mate 左侧，且首个 read-chain block 的 start/end 两个边界相对位移都至少为 19 bp；
  - `forward` 表示 read pair 可按线性片段解释。

Simulator truth 不再区分 `bsj` / `backward`。在 truth 侧，`backward` 必须意味着有 read 跨 BSJ，因此与 `bsj` 是同一个可观测事件；`backward` 只保留在 CIRI 预测结果中，用来表示有 BSJ-like / circular chain 结构但无法定位具体 BSJ 位点的 read。

### 5.1 segment protocol

`r1_segments` / `r2_segments` 使用 read-chain order 输出 genomic segments。每个 genomic segment 的格式为：

```text
start-end:strand
```

该格式与 `.isoforms.tsv` 的 `start-end:strand` exon token 保持一致；read 行已经显式包含 `chrom`，所以 segment token 不再重复写 chromosome。

多个 segment 使用 `|` 分隔：

```text
10500-10550:+|10600-10698:+
```

对于 BSJ read，segments 必须按 read-chain order 输出，并在 read chain 跨过 circular boundary 的位置插入 `<bsj>` token：

```text
10000-10098:+|10500-10550:+
10500-10550:+|<bsj>|10000-10098:+
```

BSJ read 必须使用第二种显式格式，即在跨 BSJ 位置加入 `<bsj>`，这样可以直接从 answer 文件中看出该 read 的来源跨过 back-splice junction。

约定：

- segment 顺序必须按 read-chain order 输出，而不是按 genomic coordinate 重排。
- 对负链 isoform，`strand=-`，但每个 token 内部的 `start <= end` 仍使用 genomic coordinate 的自然顺序。
- 如果 read 横跨多个 exon，需要拆成多个 genomic segments。
- 如果 read 在 circular boundary 回绕，必须保留回绕前后的 read-chain 顺序，并在 BSJ 分区之间插入 `<bsj>`。
- 对不跨 BSJ 的 circular-origin read，如果它来源于一个线性上连续可解释的 exon 片段，则 `is_circular=0` 且 `type=forward`；如果 mate 自身不跨 BSJ 但 read pair 呈严格 outward-facing geometry，且 reverse/forward 首个 block 的两个边界位移都至少为 19 bp，则 `is_circular=1` 且 `type=outward`。其他 circular-origin fragments 不写成 truth-side `backward`。

### 5.2 设计动机

`.reads.tsv` 记录 simulator 的真实 read-pair 来源：

- CIRI-AS/full-length reconstruction 可以按 `read_id` 定位具体不一致 read。
- 对每个 read pair，可以直接查看它来自哪个 circRNA 和哪个 isoform。
- 可以判断某个 mate 是否真的跨 BSJ，以及跨 BSJ 的两段分别来自哪些 genomic intervals。
- 可以把 CIRI 输出的 `junction_reads_ID` 与 simulator answer 做 read-level 对照。

## 6. RO / sequence-overlap debug 输出

主 answer 协议聚焦 full-length isoform reconstruction 所需的两类标准答案：

- circRNA / isoform 结构与 read count；
- read pair 的真实来源片段和 circRNA/BSJ 标签。

RO detector 评估使用独立 debug 输出，例如：

```text
<prefix>.ro_debug.tsv
```

该 debug 输出与主 answer 协议解耦，适合保存 sequence-overlap、merged sequence、candidate-entry 和 final strict evidence 等调试字段。

RO debug 输出如果恢复，应继续遵守三层口径：

- raw overlap layer：sequence overlap 本身；
- candidate / evidence layer：是否进入 downstream candidate；
- strict final layer：是否成为最终 full-length evidence。

这些层与 circ/BSJ answer 分别评估。

## 7. 回归检查口径

`src/simulator/tests.rs` 中的 simulator contract test 应固定检查以下 invariant：

- 输出目录中只出现第 3 节列出的 5 个正式文件。
- 所有 TSV header 与本文档列出的字段顺序一致。
- `<prefix>.annotation.gtf` 存在，且不包含被抽为 circRNA-exclusive 的 exon coordinates。
- `--exon-exclusive-rate` 的抽样基于去冗余 exon 坐标集合；共享同一 `(chrom, start, end, strand)` 的 exon records 应作为同一个 sampling unit。
- `--exon-exclusive-rate` 取值范围为 `[0, 1]`；每个去冗余 exon 独立生成随机数，随机数 `<= exon_exclusive_rate` 时标记为 circRNA-exclusive。
- 在固定 `--seed` 下，exclusive exon set、`.annotation.gtf`、`.isoforms.tsv` 和 `.reads.tsv` 必须稳定复现。
- `<prefix>.annotation.gtf` 只包含 `gene`、`transcript` 和 `exon` records；没有剩余 exon 的 transcript 以及没有剩余 transcript 的 gene 必须删除。
- `<prefix>.annotation.gtf` 中 surviving `gene`、`transcript` 和 `exon` records 的输出顺序必须与 raw input GTF 保持一致。
- 所有 truth TSV 坐标必须使用 1-based closed interval。
- circRNA read pair 数必须按 `max(1, round(circ_len * coverage / (2 * read_len)))` 整数化；linear read pair 数必须按 `round(total_linear_len * coverage / (2 * read_len))` 整数化。
- 初始 isoform 模型固定为 `isoform1` 连续 exon window；当 exon window 至少 3 个 exons 且 `read_cnt >= 2` 时，以概率 `0.5` 生成一个 skip-one-internal-exon 的 `isoform2`。
- circRNA read pairs 在 isoforms 间均匀随机分配，并且每个输出 isoform 至少有 1 个 read pair。
- linear background reads 只能从 `<prefix>.annotation.gtf` 的可用 transcript / exon chain 采样。
- circular reads 可以来自完整输入 GTF，因此可以包含 `<prefix>.annotation.gtf` 中不存在的 circRNA-exclusive exons。
- `<prefix>.isoforms.tsv` 的 `isoform_exons` 可用 `*` 标记实际使用的 circRNA-exclusive exon，标记不得改变 exon 坐标解析结果。
- `<prefix>.isoforms.tsv` 的数据行数等于模拟 circRNA 数。
- `<prefix>.reads.tsv` 的数据行数等于解压后 `<prefix>_1.fq.gz` 与 `<prefix>_2.fq.gz` 的 FASTQ record 数。
- `<prefix>.isoforms.tsv` 中 `isoform_cnt` 之和等于所有 circRNA 的 isoform 数。
- `<prefix>.isoforms.tsv` 中 `read_cnt` 之和等于 `<prefix>.reads.tsv` 中 `circ_id != NA` 的行数。
- `<prefix>.isoforms.tsv` 中 `bsj_read_cnt` 之和等于 `<prefix>.reads.tsv` 中 `is_bsj=1` 的行数，也等于控制台 `BSJ feature: <read pairs>`。
- `<prefix>.isoforms.tsv` 中每个 `isoform_read_cnt` 的和等于该 circRNA 的 `read_cnt`。
- `<prefix>.isoforms.tsv` 中每个 `isoform_bsj_read_cnt` 的和等于该 circRNA 的 `bsj_read_cnt`。
- `<prefix>.reads.tsv` 中每行的 `is_bsj` 必须等于 `r1_is_bsj OR r2_is_bsj`。
- `<prefix>.reads.tsv` 中 `is_circular` 必须按 topology rule 标注；来自 circRNA 但可完美线性解释的 read pair 应为 `is_circular=0`。
- `<prefix>.reads.tsv` 中 `type=outward` 必须满足 `is_circular=1` 且 `is_bsj=0`，并具有严格 outward-facing pair geometry：reverse-oriented 首个 block 在 forward-oriented 首个 block 左侧，且 start/end 两个边界相对位移都至少为 19 bp；`type=backward` 不再由 simulator truth 输出。
- `<prefix>.reads.tsv` 中 `r1_is_bsj=1` 的行必须在 `r1_segments` 中包含 `<bsj>` token；`r2_is_bsj=1` 同理。
- 对每个 circRNA，按 `.reads.tsv` 汇总 `isoform_id` 的 `circ_id != NA` 行数必须等于 `.isoforms.tsv` 中对应 `isoform_read_cnt`。
- 对每个 circRNA，按 `.reads.tsv` 汇总 `isoform_id` 的 `is_bsj=1` 行数必须等于 `.isoforms.tsv` 中对应 `isoform_bsj_read_cnt`。

## 8. 实现边界

`ciri-simulator` 的基础实现边界：

- base quality 使用固定 `I`。
- 普通 reads 只模拟 substitution error。
- indel、复杂测序错误和 read-position-specific quality profile 属于后续误差模型扩展。
- SAM/BAM 由外部 aligner 生成。
- `<prefix>.annotation.gtf` 是 linear RNA annotation；后续运行 `ciri` 时应使用该文件作为 `-a` 输入，以评估 circRNA-specific exon / isoform 的识别效果。
- `.reads.tsv` 记录所有模拟 read pairs。
- linear background read pairs 使用 `circ_id=NA`、`start=NA`、`end=NA`、`isoform_id=NA`、`is_circular=0`、`is_bsj=0`，并保留 `chrom` 与 `strand` 以解释不带 chromosome 的 segment tokens。
- circular origin read pairs 根据可观测 topology rule 标注 `is_circular`，并根据 R1/R2 是否真实跨 BSJ 标注 `is_bsj`；没有 BSJ 或 outward-facing pair 证据的 circular-origin reads 仍写作 `type=forward`。

## 9. 当前实现状态

已实现：

- 正式输出文件：`<prefix>_1.fq.gz`、`<prefix>_2.fq.gz`、`<prefix>.annotation.gtf`、`<prefix>.isoforms.tsv`、`<prefix>.reads.tsv`。
- 去冗余 exon-coordinate Bernoulli sampling 与 `--exon-exclusive-rate`。
- 只包含 `gene` / `transcript` / `exon` 且保持 raw GTF surviving record 顺序的 linear annotation。
- circRNA 从完整输入 GTF 采样，linear background 从 `<prefix>.annotation.gtf` 对应的 filtered transcript/exon chain 采样。
- `start-end:strand` segment protocol、`<bsj>` token、`is_circular` topology rule 和 read-level truth。
- simulator contract test，锁定正式文件集合、TSV header、read count 汇总、BSJ 标签和 segment chromosome 省略规则。

后续扩展：

1. 增加更小的 synthetic reference/GTF fixture，用于专门覆盖负链、共享 exon、全 transcript 被过滤、跨多个 exon 和跨 BSJ 的边界 case。
2. 按第 6 节补充 `<prefix>.ro_debug.tsv`，用于 RO detector 评估；该文件仍保持 debug side output，不进入正式 answer 文件集合。

---
最后更新：2026-05-03
