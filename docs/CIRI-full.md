# CIRI-full 功能拆解与历史参考

本文档记录 `vendor/CIRI-full` Java 模块的功能边界、数据契约和算法阶段。当前 CIRI-toolkit 已经实现基于 `<prefix>.segments` 的 major isoform 重构，后续活跃路线是 multi-isoform usage 和 multi-sample integration；不再沿 CIRI-full RO remap / Merge pipeline 继续开发，也不计划完整复刻 Java 历史输出。

本文件仅作为历史参考和风险清单保留。若将来重新评估 RO side evidence，应先作为独立 sidecar 设计验证，不能反向改变 CIRI3 parity 主流程。

## 1. 定位

CIRI-full 的目标不是重新识别 BSJ，而是在已有 circRNA 候选基础上重构全长 circRNA isoform。它有两条信息来源：

- CIRI-AS：提供 circRNA 内部 cirexon、junction read mapping 和可变剪接候选。
- RO reads：从 paired-end RNA-seq 中识别 5' reverse overlap，将 read pair 合并成长单端 read，再通过 BWA-MEM 比对和 splice motif 检查推断全长或部分 circRNA 结构。

历史 CIRI-full pipeline 边界为：

```text
Scan1 -> Scan2 -> Summary -> CIRI-AS -> CIRI-full
                         \                  ^
                          \-> RO1 -> BWA -> RO2
```

其中：

- `Scan1 -> Scan2 -> Summary` 继续负责 BSJ 识别、FSJ 计数、stringency 过滤和主结果输出。
- `CIRI-AS` 作为下游结构分析，为 CIRI-full Merge 提供 `_jav.list`。
- `RO1/RO2` 是 CIRI-full 自己的 RO read 识别与验证路径。
- `Merge` 只整合结构信息，不应反向改变主 `.ciri/.out` 的 circRNA 判定。

当前不再把该 pipeline 作为 CIRI-toolkit 的实现路线。CIRI-toolkit 的当前边界是：

```text
Scan1 -> Scan2 -> Summary -> segments -> major isoform -> future usage/multi-sample
```

当前实现口径为：

- 不要求 CIRI-full 与 Java 模块 100% parity。
- 不复刻历史 `Merge` detail annotation、旧 header、未验证的 RO-only 聚类输出和辅助 FASTA 提取功能。
- 不再优先实现 RO1/RO2 或 CIRI-full Merge。
- 仍可保留 RO1、RO2、Merge 章节作为术语和风险参考：例如 sequence-level RO overlap 不能直接等价于 full-length structure evidence。
- 后续 multi-isoform usage / multi-sample integration 的活跃设计以 `docs/07-full-length-reconstruction.md` 为准。

## 2. 上游 Java 模块概览

源码入口为 `vendor/CIRI-full/src/CIRI_Full2.java`，支持 4 个模式：

- `Pipeline`：自动运行 CIRI、CIRI-AS、RO1、BWA、RO2、Merge。
- `RO1`：从 paired FASTQ 中识别 5' RO，并输出候选合并 read。
- `RO2`：分析 RO1 FASTQ 经 BWA-MEM 后的 SAM，判断 RO read 是否支持 Full/Part circRNA 重构。
- `Merge`：整合 CIRI-AS `_jav.list` 和 RO2 `_ro2_info.list`，输出最终 full-length circRNA detail annotation。

Java Pipeline 的执行顺序：

```text
1. 创建 workdir:
   CIRI_output/
   CIRI-AS_output/
   CIRI-full_output/
   sam/

2. 启动 tmp1.sh:
   bwa mem -T 19 -t <threads> <ref> <read1> <read2> > sam/<prefix>_ciri.sam
   perl CIRI2.pl ...
   perl CIRI_AS_v1.2.pl ... -D yes

3. 主线程同时运行 RO1:
   RO1_single(read1, read2, CIRI-full_output/<prefix>, 95, 13)

4. 运行 BWA 比对 RO1 输出:
   bwa mem -T 19 <ref> CIRI-full_output/<prefix>_ro1.fq > sam/<prefix>_ro.sam

5. 运行 RO2:
   RO2(ref, sam/<prefix>_ro.sam, read_length, 100000, CIRI-full_output/<prefix>)

6. 等待 CIRI-AS 结束后运行 Merge:
   Merge3(
     CIRI_output/<prefix>.ciri,
     CIRI-AS_output/<prefix>_jav.list,
     CIRI-full_output/<prefix>_ro2_info.list,
     annotation.gtf,
     CIRI-full_output/<prefix>
   )
```

需要注意两个 Java 行为：

- Pipeline 的 RO2 `range` 固定传 `100000`；RO2 独立模式 CLI 默认值代码中是 `200000`，help 文本写的是 `100000`。Rust 对齐阶段应优先明确选择并用基线验证。
- Merge CLI 解析了 `-r reference.fa`，但 `Merge3` 构造函数没有使用 reference；真正输入只有 CIRI、CIRI-AS、RO2、GTF 和输出前缀。

## 3. 输入输出契约

### 3.1 Pipeline 输入

Java Pipeline 主要参数：

- `-1`：read1 FASTQ。
- `-2`：read2 FASTQ。
- `-r`：reference FASTA。
- `-a`：可选 GTF annotation。
- `-o`：输出 prefix，默认 `out`。
- `-d`：workdir，必填，且 Java 要求该目录运行前不存在。
- `-t`：线程数，默认 `1`。
- `-0`：透传给 CIRI2 的选项。

### 3.2 RO1 输出

RO1 生成：

- `<prefix>_ro1_align.txt`
- `<prefix>_ro1.fq`

`_ro1_align.txt` 每条候选一行：

```text
read_name    identity%    r1_overlap_start    r1_overlap_end    r2_overlap_start    r2_overlap_end    read_length
```

坐标为 Java 输出的 1-based 闭区间。

`_ro1.fq` 为合并后的候选 RO read FASTQ：

```text
@read_id
merged_sequence
+
merged_quality
```

其中 read name 直接使用 R1 的 FASTQ header。

### 3.3 RO2 输出

RO2 生成：

- `<prefix>_ro2.sam`
- `<prefix>_ro2_info.list`

`_ro2_info.list` header：

```text
Read_ID    Chr    BSJ_position    Strand    Reconstructed_state    Cirexon    Mapping_order    Splice_site_state+    Splice_site_state-
```

字段语义：

- `Read_ID`：RO1 FASTQ read id。
- `Chr`：参考染色体。
- `BSJ_position`：`chr:start|end`，或 `No BSJ`。
- `Strand`：`+`、`-` 或 `0`。
- `Reconstructed_state`：`Full` 或 `Part`。
- `Cirexon`：RO2 推断出的 circRNA exon 区间，格式如 `start-end,start-end,`。
- `Mapping_order`：RO read 上各段映射到基因组的顺序，格式同样为 `start-end,`。
- `Splice_site_state+`：正链 splice motif 搜索记录，字符串以 `readposBound+:` 开头。
- `Splice_site_state-`：负链 splice motif 搜索记录，字符串以 `readposBound-:` 开头。

`_ro2.sam` 保存被 RO2 接受并用于重构的 SAM 记录，包括原始 BWA 记录和局部 remap 生成的 `Index_find_1/2` 伪记录。

### 3.4 Merge 输出

Merge 生成：

- `<prefix>_merge_circRNA_detail.anno`

Java header：

```text
BSJ    Chr    Start    End    GTF-annotated_exon    Cirexon    Coveage    BSJ_reads_information RO_reads_information Original_gene    strain
```

注意：

- `Coveage` 是 Java 原始拼写。
- header 中 `BSJ_reads_information RO_reads_information Original_gene` 被写成一个带空格的字段名，但数据行实际按 tab 写出 `Coverage`、`BSJ_reads_information`、`RO_reads_information`、`Original_gene`、`strain` 等字段。Rust 若要完全复刻旧输出，需要保留该 header；若要提供新格式，应另设版本化输出。
- `strain` 是 Java 原始拼写，语义上对应 strand。

## 4. RO1：5' RO 识别与 read pair 合并

实现类：

- `RO1_single.java`
- `index_head_finder_more.java`
- `reversc.java`

### 4.1 FASTQ 读取

RO1 逐对读取 R1/R2 FASTQ：

- 判断是否 gzip 仅依据 R1 文件名最后一段是否为 `gz`。
- gzip 输入使用 `GBK` 编码读取。
- 非 gzip 输入使用默认 `FileReader`。
- read length 使用最后读取到的 R1 长度返回给 RO2。

Rust 实现建议：

- 支持 plain FASTQ 和 gzip FASTQ。
- 明确使用字节级 FASTQ parser，避免编码行为影响序列和质量值。
- 对齐阶段仍需保留 Java 对 read name、序列大写化和质量反转的输出行为。

### 4.2 R2 反向互补

Java `reversc` 的逻辑：

1. `A -> t`
2. `T -> a`
3. `G -> c`
4. `C -> g`
5. 转大写
6. 字符串 reverse

因此只显式处理 `A/T/G/C`，其他字符在互补阶段保持原样，然后整体反转。R2 quality 只 reverse，不做其他处理。

### 4.3 重叠识别

`index_head_finder_more(seq1, seq2_rc, q1, q2_rc, identity, min_len)`：

- `seq1` 为 R1 sequence。
- `seq2` 为 R2 reverse-complement sequence。
- 默认 identity 为 `95`。
- 默认最小 overlap 长度为 `13`。

搜索策略：

```text
for indexn in [0, 1, 2]:
  seed = seq1[indexn..indexn+8]
  for startpoint in 0..seq2.len()-8:
    if seed 与 seq2[startpoint..startpoint+8] 至少 7bp 相同且非 N:
      overlap_r1 = seq1[indexn..seq2.len()-startpoint]
      overlap_r2 = seq2[startpoint..]
      identity = matching_non_N_bases * 100 / overlap_len
      if identity >= threshold and overlap_len >= min_len:
        accept first hit
```

重要 parity 点：

- 只尝试 R1 起点 `0/1/2` 的 8bp seed。
- seed match 要求 `>=7` 个非 N 相同碱基。
- overlap identity 使用整数除法语义：Java 表达式 `right * 100 / i` 先做整数除法，再赋给 double。
- 一旦找到第一个通过阈值的 hit，立即停止搜索，不选择最高 identity。
- 合并序列为 `seq2_rc[0..startpoint] + seq1[indexn..]`。
- 合并质量为 `q2_rc[0..startpoint] + q1[indexn..]`。

输出坐标：

```text
r1_overlap_start = indexn + 1
r1_overlap_end   = indexn + overlap_len
r2_overlap_start = startpoint + 1
r2_overlap_end   = startpoint + overlap_len
```

## 5. RO2：RO read 比对验证与重构

实现类：

- `RO2.java`
- `sam5.java`
- `index_head_finder2.java`

RO2 输入为 reference FASTA、RO1 FASTQ 经 BWA-MEM 后的 SAM、read length、range 和输出 prefix。其核心任务是：

1. 按染色体拆分 SAM。
2. 逐染色体加载 reference sequence。
3. 按 read id 聚合多段比对。
4. 过滤低质量或非候选 read。
5. 对未对齐首尾片段做局部 remap。
6. 在相邻 mapping segment 间搜索 GTAG/CTAC splice motif。
7. 输出 Full/Part RO 重构结果。

### 5.1 SAM 预处理

RO2 首先读取 SAM header：

- 只处理 `@SQ` 中长度 `>500000` 的 contig。
- 每个 contig 创建临时文件。
- 非 `*` 比对按 `RNAME` 写入对应临时文件。
- Java 使用当前工作目录下的临时文件 `ciri_tmp_<chr>.sam`。
- 写入时尽量保持同一个 read id 的记录不被拆散。

Rust 实现建议不落地临时文件，改为按 contig 分桶或外部排序；但 parity 阶段要保证同一 read group 的遍历顺序与 Java 一致。

### 5.2 CIGAR 解析

Java `sam5.getpos()` 返回 4 个整数：

```text
pos[0] = left soft/hard clip length
pos[1] = query mapped length，M 和 I 都计入
pos[2] = right soft/hard clip length
pos[3] = reference indel adjustment，D 增加，I 减少
```

行为细节：

- 先把 `S` 替换为 `H`，因此 soft/hard clip 在逻辑中等价。
- `M` 增加 `pos[1]`。
- `I` 增加 `pos[1]`，同时 `pos[3] -= I`。
- 其他非数字 CIGAR op 统一走 `else`，会把长度加到 `pos[3]`。这意味着 `D/N/=/X` 等并没有被区分。

Rust 版本如果用 `noodles-sam` 的结构化 CIGAR，需要显式复刻上述分类，否则 RO2 的 segment 坐标会偏离 Java。

### 5.3 read group 初筛

对同一 read id 的多条 SAM 记录，Java 逻辑大致如下：

- 跳过 `chrM`。
- 以第一条有效记录确定 chr 和 strand。
- 只保留同 chr、同 strand、与首条记录起点距离 `< range/2`、MAPQ `>5` 的记录。
- 按 read-coordinate 起点排序。
- 计算：
  - `total = first_left_clip + first_mapped_len + first_right_clip`
  - `mapped_bases = sum(segment_mapped_len)`
  - `per = mapped_bases / total`
- 如果 `per < 0.5`，丢弃。
- 如果最佳 MAPQ `<15`，丢弃。

这里的 strand 来自 SAM flag：

- flag `< 16` 视为正向。
- flag `>= 16` 视为反向。

这是 Java 的简化判断，不等价于完整 SAM bit flag 解析。Rust parity 阶段应保持这个行为。

### 5.4 首尾局部 remap

当首端或尾端有未对齐片段长度 `>=10` 时，RO2 尝试在参考序列局部窗口中重新定位：

- 首端 remap 使用 read 的前端未对齐序列。
- 尾端 remap 使用 read 的末端未对齐序列。
- 搜索窗口围绕现有 segment 的边界，窗口大小约 20kb。
- `index_head_finder2` 使用短序列中间 8bp seed：
  - `sp = (seq1.len() - 8) / 2`
  - 要求 seed 至少 7bp 非 N 相同。
  - 全长 identity 要求 `>=95`。
  - Java 变量 `idmax` 在循环内重置，实际行为是保存最后一个符合条件的位置，而不是全局最高 identity。
- remap 成功后会追加伪 SAM 记录，read name 为原 read id，字段中包含 `Index_find_1` 或 `Index_find_2`。

这些伪记录进入后续 splice motif 和输出 `_ro2.sam` 的流程，因此 Rust 版本需要显式建模，不能只作为临时判断。

### 5.5 Full/Part 分类

RO2 通过 segment 间 gap/overlap 和原始 read length 判断候选状态：

- `Full`：RO merged read 覆盖了完整 circRNA 结构，Java 内部 `yes=1`。
- `Part`：RO merged read 只支持部分结构，Java 内部 `yes=2`。
- 其他情况丢弃或只计入日志。

其中 `readlength == 0` 被视为 long-read 兼容模式，会跳过部分基于原始 read length 的边界过滤。

### 5.6 splice motif 搜索

对相邻 mapping segment，RO2 搜索 canonical splice motif：

- 正链 motif：`GTAG`。
- 负链 motif：`CTAC`。
- 先从 `mv = 5` 到 `1` 向左侧回扫，再从 `mv = 0` 到 `gap + 5` 向右侧扫描。
- 对候选 motif，还会检查 read/reference 两侧 10bp 相似性：
  - 某些边界组合要求 `>=9` 个匹配。
  - 其他组合要求 `>=10` 个匹配。
- 分别记录 plus/minus motif 候选到 `Splice_site_state+` 和 `Splice_site_state-`。

strand 选择：

- 如果 plus/minus 都不能覆盖所有 segment junction，则 strand 为 `0`。
- 否则优先选择 motif 数量更多的一侧。
- 数量相同则使用 Java 的 point sum tie-break，选择得分更小的一侧。

### 5.7 exon 区间与 BSJ 输出

RO2 将 splice motif 位置转换为 exon intervals：

- 输出区间 start 通常为 Java 0-based 内部坐标 `+1`。
- end 保持 Java substring 风格的右边界结果。
- overlapping exon 会被合并。
- 如果 junction 不完整或 segment overlap 不满足规则，会丢弃。

接受候选时：

- 如果 `junction == 0 && state == Part`，`BSJ_position = No BSJ`。
- 否则 `BSJ_position = chr:first_exon_start|last_exon_end`。
- `Cirexon` 输出合并后的 circ exon。
- `Mapping_order` 输出 RO read mapping segment 在基因组上的顺序。

## 6. Merge：整合 CIRI-AS 与 RO2

实现类：

- `Merge3.java`
- `Exon_gtf.java`

Merge 输入：

- CIRI circRNA 输出。
- CIRI-AS `_jav.list`，要求 CIRI-AS 运行时开启 `-D yes`。
- RO2 `_ro2_info.list`。
- 可选 GTF。
- 输出 prefix。

### 6.1 GTF 读取

Merge 只读取 GTF 中 `exon` 行：

- 保存 `chr, start, end, gene_id`。
- `gene_id` 从 attributes 中第一个 key 为 `gene_id` 的字段提取。
- 后续输出 `GTF-annotated_exon` 时，只选完全位于 circ 区间内的 exon：`exon.start >= circ.start && exon.end <= circ.end`。
- `Original_gene` 则取第一个与 circ span 有 overlap 的 GTF exon 的 gene_id。

Java 没有按 transcript 分组，也没有使用 GTF strand。

### 6.2 CIRI-AS 与 CIRI 结果载入

CIRI-AS `_jav.list` 中 Merge 只保留前 6 列：

```text
BSJ, Chr, Start, End, GTF-annotated_exon, predicted_cirexon
```

同时在第二遍读取 `_jav.list` 时，Merge 还使用：

- coverage 字段。
- BSJ reads information 字段。
- junction read mapping 字段。

具体列位依赖 CIRI-AS Java/Perl 输出格式，Rust 实现应先和 `docs/CIRI-AS.md` 的 `_jav.list` 契约对齐，禁止用脆弱的尾列推断。

CIRI `.ciri` 结果中 Merge 使用：

- `tmpc[0]`：BSJ id。
- `tmpc[4]`：junction read count，作为表达量 tie-break。
- `tmpc[10]`：strand。

### 6.3 RO read 关联到 CIRI-AS circRNA

Merge 首先尝试把每条 RO2 记录分配给 CIRI-AS circRNA：

1. 解析 RO2 `Mapping_order` 为多个 interval。
2. 如果 RO2 有 BSJ：
   - 找同 chr 的 CIRI-AS circ。
   - circ start/end 与 RO2 BSJ start/end 的差值都在 3bp 内则匹配。
3. 如果 RO2 是 `No BSJ`：
   - 找同 chr 且 circ span 覆盖所有 RO interval 的 CIRI-AS circ。
   - tie-break 顺序：
     1. 与 AS cirexon overlap 数更多。
     2. CIRI junction read count 更高。
     3. GTF exon overlap 数更多。
4. 成功分配后，RO read 以 `read_id##mapping_order&&` 追加到该 BSJ 的 RO read 信息，并从 RO-only 待聚类列表中移除。

### 6.4 用 RO junction 补全 CIRI-AS cirexon

对每条 CIRI-AS circ：

- 读取原始 AS predicted cirexon。
- 从 AS junction read mapping 中抽取 junction start/end。
- 从分配来的 RO reads 的 consecutive mapping intervals 中抽取 junction。
- 加入 circ BSJ end/start 作为闭环 junction。
- 对正向距离 `<=400` 的 junction 生成缺失 cirexon segment。
- 如果新 segment 不在已有 predicted cirexon 中，则追加。

最终输出该 circ 的：

- BSJ。
- Chr/Start/End。
- GTF-annotated exon。
- 合并后的 Cirexon。
- Coverage。
- BSJ reads information。
- RO reads information。
- Original_gene。
- strain。

### 6.5 RO-only 聚类输出

未被 CIRI-AS 吸收的 RO records 会被单独聚类成 RO-only circRNA：

Full RO 聚类：

- 只处理 `Reconstructed_state != Part`。
- 要求 strand 不为 `0`。
- 按同 chr、首 exon start 与 cluster start 差值 `<4`、末 exon end 与 cluster end 差值 `<4` 聚类。
- exon 也按 start/end 差值 `<4` 合并。

Part RO 聚类：

- 只处理 `BSJ_position != No BSJ`。
- 以 BSJ start/end 差值 `<4` 和同 chr 聚类。
- 使用 `exonstate` 标记部分 exon：
  - `0`：完整 exon。
  - `1`：左侧 partial。
  - `2`：右侧 partial。
- 对首末 partial exon 有多处特殊合并逻辑，包括距离 `<300` 和环形首尾修正。

RO-only 输出：

- `Coverage` 写 `n/a`。
- `BSJ_reads_information` 写 `n/a`。
- `RO_reads_information` 写聚类 reads。
- `Original_gene` 从 GTF overlap 推断。
- `strain` 优先使用 CIRI 中相同 BSJ 的 strand，否则为 `none`。

## 7. 辅助类与非主流程工具

### 7.1 Exon_gtf

`Exon_gtf` 只做一件事：从给定 GTF exon list 中返回完全落在 circ 区间内的 exon。注释掉的近似合并逻辑没有生效。

### 7.2 Exon_length

`Exon_length` 会把 exon list 裁剪到指定边界内，合并 overlap 后计算总长度。当前 CIRI-full 主入口未直接调用，可能是历史辅助代码。

### 7.3 Exon_cov2

`Exon_cov2` 根据 coverage 字符串用滑窗和阈值推断 exon-like 区间。当前 CIRI-full 主入口未直接调用。若后续需要复刻旧辅助功能，应单独做 parity，不应混入 Merge 第一版。

### 7.4 Seq

`Seq` 从 detail annotation 和 reference FASTA 中提取 `Full` circRNA 的重构序列，输出 `<newlist>_circle.fa`。当前 `CIRI_Full2` 没有暴露该模式，Java Pipeline 也没有调用它。Rust 第一版可以暂不实现，后续作为 `ciri-full seq` 或 `--emit-circle-fa` 可选功能讨论。

## 8. Rust 数据模型建议

为了和 Java 逻辑对齐，建议先按模块建立直接映射的数据结构，避免过早抽象。

```rust
struct Ro1Candidate {
    read_name: String,
    identity_percent: i32,
    r1_overlap_start: i32,
    r1_overlap_end: i32,
    r2_overlap_start: i32,
    r2_overlap_end: i32,
    original_read_len: i32,
    merged_seq: Vec<u8>,
    merged_qual: Vec<u8>,
}

struct Ro2Segment {
    read_id: String,
    chr: String,
    flag: i32,
    strand: i8,
    ref_start: i32,
    mapq: i32,
    cigar: String,
    left_clip: i32,
    mapped_query_len: i32,
    right_clip: i32,
    ref_adjust: i32,
    seq: Vec<u8>,
    source: Ro2SegmentSource,
}

enum Ro2SegmentSource {
    Bwa,
    IndexFind1,
    IndexFind2,
}

struct Ro2Record {
    read_id: String,
    chr: String,
    bsj_position: Option<(i32, i32)>,
    strand: char,
    state: Ro2State,
    cirexon: Vec<(i32, i32)>,
    mapping_order: Vec<(i32, i32)>,
    splice_plus_state: String,
    splice_minus_state: String,
}

enum Ro2State {
    Full,
    Part,
}

struct FullMergeRecord {
    bsj: String,
    chr: String,
    start: i32,
    end: i32,
    gtf_exons: Vec<(i32, i32)>,
    cirexons: Vec<(i32, i32)>,
    coverage: Option<String>,
    bsj_reads_info: Option<String>,
    ro_reads_info: Vec<String>,
    original_gene: Option<String>,
    strand: Option<String>,
}
```

关键点：

- CIRI-full 使用大量 1-based 输出区间，内部计算却混合 Java substring 的 0-based end-exclusive 习惯。实现时建议用专门类型区分 internal coordinate 和 output coordinate。
- RO2 中 `Mapping_order` 和 `Cirexon` 均保留尾随逗号。若要逐字节对齐 Java 输出，writer 必须保留。
- `read_id##mapping_order&&` 是 Merge 的字符串协议，应作为兼容 writer 的显式格式，不要临时拼接散落在代码中。

## 9. Rust CLI 接入建议

建议将 CIRI-full 暴露为独立子命令或主命令可选阶段。为了不干扰当前 CIRI3 parity，优先采用独立子命令：

```text
ciri full pipeline \
  -1 <read1.fq.gz> \
  -2 <read2.fq.gz> \
  -r <ref.fa> \
  -a <anno.gtf> \
  -o <prefix> \
  -d <workdir> \
  -t <threads>

ciri full ro1 \
  -1 <read1.fq.gz> \
  -2 <read2.fq.gz> \
  -o <prefix> \
  --min-identity 95 \
  --min-overlap 13

ciri full ro2 \
  -r <ref.fa> \
  -s <ro.sam|ro.bam> \
  -l <read_length> \
  -o <prefix> \
  --range 100000

ciri full merge \
  -c <ciri.out> \
  --as <prefix_jav.list> \
  --ro <prefix_ro2_info.list> \
  -a <anno.gtf> \
  -o <prefix>
```

后续可在主 pipeline 增加：

```text
--ciri-as
--ciri-full
--full-workdir <dir>
--full-range <int>
--full-keep-intermediate
```

第一版不建议自动调用外部 `bwa mem`，而是支持两种路径：

- pipeline convenience：检测并调用 `bwa`，行为接近 Java。
- reproducible mode：用户显式提供 RO1 FASTQ 的 BWA SAM/BAM，Rust 只做 RO2/Merge。

## 10. 分阶段实现计划

### 阶段 1：文档与 golden 数据

- 固化 Java CIRI-full 在小数据上的输入输出。
- 保存 RO1、RO2、Merge 每一步中间文件作为 golden。
- 明确 `_jav.list` 列契约，与 `docs/CIRI-AS.md` 对齐。
- 记录 Java 中 `range` 默认值差异，选择 Rust 默认值。

### 阶段 2：RO1 parity

- 实现 FASTQ pair reader。
- 实现 `reversc` 行为。
- 实现 `index_head_finder_more` 的 first-hit 搜索。
- 输出 `_ro1_align.txt` 和 `_ro1.fq`。
- 与 Java 对比候选 read id、merged sequence、merged quality、坐标和 identity。

RO1 是最适合先落地的子模块，因为不依赖 CIRI-AS 或 CIRI 主流程。

### 阶段 3：RO2 CIGAR 与 read group parity

- 实现 Java `sam5` 兼容 CIGAR 分类。
- 实现 read group 聚合和排序。
- 复刻 chr/strand/range/MAPQ/filter 逻辑。
- 先输出调试结构，与 Java `_ro2.sam` 候选记录对齐。

### 阶段 4：RO2 remap 与 splice motif

- 实现 `index_head_finder2`。
- 建模 `Index_find_1/2` 伪 segment。
- 实现 GTAG/CTAC 搜索和 strand tie-break。
- 对齐 `_ro2_info.list`。

这是 CIRI-full 中最高风险阶段，建议加入 read-level trace，能打印每个 read 的 segment、motif candidates、chosen strand、Full/Part 判定原因。

### 阶段 5：Merge parity

- 实现 GTF exon loader。
- 实现 CIRI-AS/CIRI/RO2 载入。
- 实现 RO-to-AS assignment。
- 实现 AS cirexon 补全。
- 实现 RO-only Full/Part 聚类。
- 对齐 `_merge_circRNA_detail.anno`。

### 阶段 6：集成主流程

- 将 CIRI-AS 和 CIRI-full 作为可选后处理阶段接入。
- 保证默认 CIRI3 输出完全不变。
- 增加 end-to-end smoke test。
- 如果需要自动调用 BWA，单独封装外部命令 runner，并记录命令和版本。

## 11. Parity 风险清单

后续实现时需要重点保护以下 Java 行为：

- RO1 overlap identity 使用整数除法语义。
- RO1 只取 first hit，不取 best hit。
- RO1 只搜索 R1 起点 `0/1/2` 的 seed。
- `reversc` 只显式互补 `A/T/G/C`。
- RO2 CIGAR 中 `S` 与 `H` 等价，`I` 同时影响 query length 和 ref adjustment，其他 op 都按 ref adjustment 增加。
- RO2 SAM flag strand 使用 `flag < 16` 简化判断。
- RO2 只处理 header 中长度 `>500000` 的 contig。
- RO2 默认跳过 `chrM`。
- RO2 局部 remap 的 `idmax` 实际不是全局最大值。
- RO2 splice motif 搜索顺序和 tie-break 会影响 strand 和 exon 边界。
- RO2 output interval 有大量 `+1`，必须逐项对照 Java。
- Merge 的 3bp/4bp 近似匹配阈值分别用于不同阶段，不可统一。
- Merge 的 `read_id##mapping_order&&` 字符串协议被后续输出直接消费。
- Merge header 存在拼写和 tab/space 历史问题，是否保留需要提前决定。

## 12. 测试建议

以下测试分层是历史 CIRI-full parity 路线的建议，当前不作为活跃实现任务：

- `full_ro1_overlap_unit`：覆盖 perfect overlap、1 mismatch、N、min overlap、first-hit。
- `full_cigar_sam5_unit`：覆盖 `S/H/M/I/D/N/=/X` 在 Java 兼容分类下的结果。
- `full_index_find_unit`：覆盖局部 remap seed 和最后符合位置行为。
- `full_ro2_read_trace`：固定一组 RO SAM，对齐 Full/Part/No BSJ 输出。
- `full_merge_assignment_unit`：覆盖 RO with BSJ、RO No BSJ、表达量 tie-break、GTF overlap tie-break。
- `full_merge_ro_only_unit`：覆盖 Full 聚类和 Part 聚类的 4bp 阈值。
- `full_end_to_end_golden`：对齐 Java Pipeline 的最终 `_merge_circRNA_detail.anno`。

对比策略：

- RO1：逐行比较 `_ro1_align.txt`，并比较 `_ro1.fq`。
- RO2：先比较 `_ro2_info.list`，再抽样比较 `_ro2.sam`。
- Merge：逐行比较 `_merge_circRNA_detail.anno`；如果后续决定修正 header，应同时保留 Java-compatible writer 用于 parity。

## 13. 与现有 CIRI-toolkit 的关系

当前 Rust CIRI3 主流程已经达到 parity，且 CIRI-toolkit 已经通过 `<prefix>.segments` 实现 major isoform 输出。CIRI-full 不再作为下一阶段接入路线。若未来重新评估本文件中的 RO side evidence，应遵守：

- 默认不开启 CIRI-full/RO remap。
- 不修改 Scan1/Scan2/Summary 判定逻辑。
- 不要求 CIRI-AS `_jav.list` 或 CIRI-full Merge 兼容输出。
- RO1/RO2 只能作为 sidecar evidence 重新设计，不能作为主流程前置依赖。
- 当前优先级仍是 multi-isoform usage 与 multi-sample integration，见 `docs/07-full-length-reconstruction.md`。

历史建议优先顺序如下，当前已归档：

```text
RO1 -> RO2 -> CIRI-AS _jav.list writer -> Merge -> full pipeline orchestration
```

这样可以把 CIRI-full 中最独立、最容易建立 golden 的 RO 路径先做出来，再整合 CIRI-AS 和主流程。

---
最后更新：2026-05-19
