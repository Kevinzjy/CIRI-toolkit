# CIRI-AS 功能拆解与历史参考

本文档记录 `vendor/CIRI-AS/CIRI_AS_v1.2.pl` 的功能边界、数据契约和算法阶段。当前 CIRI-toolkit 已经实现自己的 `<prefix>.segments -> major isoform` 后处理路线，不再按 CIRI-AS 输出协议继续开发，也不计划完整复刻 Perl 脚本的历史文件。

本文件作为历史参考保留，主要用于解释仍被当前 segments/isoform 逻辑复用的局部思想：circ-local read 重扫、alignment-chain 解释、annotation-aware splice boundary 校正，以及 read-level evidence 如何进入后续 isoform usage。

## 1. 定位

CIRI-AS 是 CIRI3 之后的扩展分析阶段，用于识别 circRNA 内部组成和可变剪接事件。它不替代 `Scan1 -> Scan2 -> Summary` 的 BSJ 检测逻辑，而是基于最终 circRNA 结果和同一份比对文件重新扫描 reads。

当前实现口径为：

- 不要求 CIRI-AS 与 Perl 脚本 100% parity。
- CIRI-AS v1.2 未按 CIRI3 主流程同等级别严格验证；Rust 实现允许在文档记录清楚的前提下修正不稳定或生物学上较弱的启发式。
- 不复刻未被后续使用的历史中间文件。
- 只保留对当前 segments、major isoform 和后续 usage 有价值的思路：
  - circ 区间内 read 重新扫描。
  - 内部 splice junction 识别。
  - exon boundary / candidate junction 解释。
  - GTF annotation 辅助。
  - read-level evidence 到 isoform usage 的映射。
- 当前活跃开发路线以 `docs/07-full-length-reconstruction.md` 为准：已实现 read-level segments 与单样本 major isoform，下一阶段是 multi-isoform usage 和 multi-sample integration。

当前 Rust pipeline 边界为：

```text
Scan1 -> Scan2 -> Summary -> segments -> major isoform -> future usage/multi-sample
```

其中：

- `Scan1 -> Scan2 -> Summary` 继续负责 circRNA BSJ 识别、FSJ 计数、stringency 过滤和最终 `.out`。
- CIRI-AS 思路中仍有用的局部规则已经拆解进 segments / isoform reconstruction 模块，而不是以旧脚本输出为直接目标。
- 任何 AS/isoform 扩展默认不改变主 `.out`，也不回写 Scan1/Scan2/Summary 的核心判定结果。

## 2. 上游脚本输入输出

上游 Perl 脚本的命令行接口主要包括：

- `-S/--sam`：输入 SAM，脚本要求 BWA-MEM paired-end 模式生成。
- `-C/--ciri`：CIRI circRNA 结果表。
- `-O/--out`：输出前缀。
- `-F/--ref_file` 或 `-R/--ref_dir`：参考基因组 FASTA。
- `-A/--anno`：可选 GTF/GFF 注释。
- `-D/--output_all`：是否输出中间调试文件。

Rust 版本建议简化为复用主程序已有输入：

```text
ciri \
  -i <input.sam|input.bam> \
  -o <prefix> \
  -r <ref.fa> \
  -a <anno.gtf>
```

当前实现不再使用独立的 `--as` / `--as-out` 入口，而是在主流程写完 `<prefix>.out` 后默认继续执行 segments 后处理，并从完成的 `<prefix>.segments` 重新解析 major isoform。

当前正式输出：

- `<prefix>.segments`：read-level segment chain，当前包含 `type=bsj/backward/outward`。
- `<prefix>.isoforms.gtf`：每个 Summary-confirmed circRNA 的 major isoform 审计记录。
- `<prefix>.isoforms.fa`：高可信 major isoform sequence candidate。
- `<prefix>.bedpe` 与 `<prefix>.segments.bam/.bai`：IGV review sidecar。

当前保留但暂不作为主输出的内容：

- `_splice.list`
- `.list`
- `.isoforms`
- `.isoform_summary`
- `.fa`
- `_AS.list`
- `_coverage.list` / `_jav.list` / `_library_length.list`

这些旧输出对应的设计说明仍保留在本文档后续章节中，作为历史参考，而不是当前默认交付物或后续必须实现的协议。

## 3. 重要兼容问题

当前 CIRI-toolkit `.out` 为 13 列：

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

Perl CIRI-AS 直接把 circRNA 表最后一列当成 junction read 列：

```perl
my @reads = split ',', $line[-1];
```

这与当前 CIRI-toolkit 输出不兼容，因为最后一列是 `Score`，真正的 read 列是 `junction_reads_ID`。Rust 版本必须按 header 解析列名，不能按最后一列读取 read IDs。

这是 CIRI-AS Rust 实现的第一个强约束：

```text
必须按 header 定位 junction_reads_ID，禁止使用“最后一列”作为 read 列。
```

## 4. 核心数据模型

建议先定义一组接近 Perl hash 语义的 Rust 结构，便于 parity 阶段逐步对齐。

```rust
struct CircRecord {
    id: String,
    chr: String,
    start: i32,
    end: i32,
    junction_reads: i32,
    pcc: String,
    non_junction_reads: i32,
    junction_ratio: String,
    circ_type: String,
    gene_id: String,
    strand: Option<char>,
    read_ids: Vec<String>,
}

struct CircCluster {
    id: usize,
    chr: String,
    start: i32,
    end: i32,
    circ_ids: Vec<String>,
}

struct ReadMappingRange {
    chr: String,
    start: i32,
    end: i32,
    read_start: i32,
    read_end: i32,
    mapq: i32,
}

struct InternalSpliceCandidate {
    read_id: String,
    segment_id: usize,
    chr: String,
    site1: i32,
    site2: i32,
    adjust1: i32,
    adjust2: i32,
    cigar_slots: [Option<String>; 3],
    mate_flag: i32,
}

struct SpliceCluster {
    id: usize,
    chr: String,
    site1: i32,
    site2: i32,
    supporting_reads: Vec<(String, usize)>,
    sm_ms_sms: [usize; 3],
}

struct Cirexon {
    circ_id: String,
    start: i32,
    end: i32,
    start_support: i32,
    end_support: i32,
    coverage_median: i32,
    is_icf: bool,
}

struct AsEvent {
    circ_id: String,
    exon: String,
    event_types: Vec<AsType>,
    psi_raw: Option<f64>,
    psi_corrected: Option<f64>,
}
```

命名可以在实现阶段调整，但 parity 阶段建议保留 Perl 变量对应关系，便于 trace 对照。

## 5. 算法阶段拆解

### 5.1 输入校验与 read length 推断

Perl 脚本先扫描 SAM 前若干 read，检查：

- 是否 paired-end；
- R1/R2 是否都有足够 read；
- R1/R2 read length 是否一致；
- read length 是否不小于 40。

Rust 版本建议：

- 继续要求 queryname-sorted 或 unsorted PE 输入；
- SAM/BAM 都通过现有 `sam_bam` 与 `AlignmentRecord` 抽象读取；
- 第一版 CIRI-AS 可明确要求 paired-end，单端输入直接跳过 AS 或返回明确错误；
- read length 推断尽量复用 Scan1 已得到的 `read_len`，避免重复全文件探测。

### 5.2 注释读取与长 exon 抽样

Perl 读取 GTF/GFF exon，建立：

- `start_exon{chr}{locus} = start`
- `end_exon{chr}{locus} = end`
- `length_exon{locus}`
- `gene_exon{chr}{locus} = gene_id`
- 用长 exon 估计 insert length distribution。

当前 `src/annotation.rs` 已有 exon boundary 和 gene span，但 CIRI-AS 还需要：

- 按 chr 排序的 exon 列表；
- exon strand；
- exon length；
- gene id；
- 能快速找出 circ 区间内的 annotated exon。

建议扩展 `Annotation`，新增专用于 AS 的 exon records，不移除现有字段，避免影响 Summary parity。

### 5.3 读取 circRNA 结果

Perl 建立：

- `junction_read{read_id} = circ_id`
- `circ_junc_reads{circ_id}{read_id}`
- `circ_info{circ_id} = chr/start/end/...`
- `circ_chr{chr}{circ_id} = #junction_reads`

Rust 版本实现 `load_circ_records(result_path)`：

- 按 header 读取列；
- 支持 CIRI-toolkit 当前 13 列格式；
- 可兼容旧 CIRI 输出中 `junction_reads_ID` 位于最后一列的格式；
- 对空 read list、空 circ list 给出明确错误。

### 5.4 circRNA 区间聚类

Perl 以 chromosome 内 circ start/end 排序，并按以下条件聚类：

```text
next_circ_start < current_cluster_end + min_intron * 2 + 1
```

其中 `min_intron = 70`。每个 cluster 的范围索引扩展：

```text
cluster_start - min_intron .. cluster_end + min_intron
```

Rust 第一版可以直接保留 per-base range index，便于 parity。后续再替换为 interval tree 或分桶索引。

### 5.5 重扫 SAM/BAM 与覆盖度统计

Perl 重扫 SAM 时按 read name 分组，区分两类 read：

1. 已知 junction read：调用 `record_mapping_detail` 和 `mapping_check1`。
2. 非 junction read 但落在 circ cluster 内：调用 `mapping_check1_add`。

同时统计：

- circ cluster 内 per-base coverage；
- junction read 在 circ 内的 mapping range；
- 用于 insert length distribution 的长 exon paired-end mapping。

Rust 实现建议：

- 新增 `CiriAs::scan_alignments()`；
- 复用 Scan1/Scan2 的 read-group traversal 思路；
- 输出内部状态而不是边扫边写，方便后续单元测试；
- 第一版保持 per-base `HashMap<(chr, pos), depth>`，确认后再优化。

### 5.6 CIGAR 解析与内部 splice 候选识别

Perl 的关键 helper：

- `MSID`
- `MSID_start`
- `mapping_check1`
- `mapping_check1_add`
- `mapping_check2`
- `mapping_check2_add`

这些函数决定 split alignment 如何转为 circ 内部 splice candidate。它们使用 1-based genomic coordinate 和 inclusive end，且会把 `H` 视作 `S`。

Rust 迁移建议：

- 独立放入 `src/as_cigar.rs`；
- 先为 `MSID` / `MSID_start` 建 golden unit tests；
- 注释必须写清楚每个返回槽位的含义；
- 不在第一阶段“重构”成更 Rust 风格的 CIGAR parser，先镜像 Perl 控制流。

### 5.7 splice signal 校验、链向判断与坐标微调

Perl 对候选 splice 使用 `splice_loci_check` 与 `index_compare`：

- 从参考序列取候选边界附近短序列；
- 检查 AC/CT 或 AG/GT 是否在相同 offset 出现；
- 根据 motif offset 调整 `site1/site2/adjust1/adjust2`；
- 校验通过的候选加入 `loci_validated`。

方向含义：

- AC/CT 表示负链线索；
- AG/GT 表示正链线索。

Perl v1.2 的 `index_compare` 有两个重要限制：

1. 原脚本使用 `if AC/CT ... elsif AG/GT ...`，因此当同一短窗口同时存在负链和正链 motif 解释时，会先尝试 AC/CT，只有 AC/CT 不满足时才尝试 AG/GT。
2. 同一 motif 类型内若存在多个合法 offset，原脚本通过 `keys %hash` 取 `$result[0]`，结果受 Perl hash 遍历顺序影响；同一输入在不同 `PERL_HASH_SEED` / `PERL_PERTURB_KEYS` 下可能产生不同 `_splice.list` 坐标和 cluster 数。

Rust 实现不追求复刻这种不稳定行为。当前约定为：

- 同时枚举 AC/CT 与 AG/GT 的合法 offset。
- 如果提供 annotation，则优先选择能让 `site1` 命中 exon start、`site2` 命中 exon end 的解释。
- 若多个解释均有 annotation 支持，同 gene 与 annotation strand 一致的解释优先。
- annotation 无法区分时，才回退到 CIRI-AS v1.2 的链向判断顺序；同一链向内回退到最小 offset，保证输出稳定。
- 无 annotation 时仍使用稳定 fallback，不依赖 Perl hash 顺序。

这是一个有意偏离 Perl CIRI-AS 的地方：目标是得到更稳定、可解释、与 exon boundary 更一致的 internal splice 坐标，而不是字节级复刻 `_splice.list`。

这里要特别注意 Perl `substr` 是 0-based，长度参数与 Rust slice 的 end-exclusive 语义不同。实现应保留 `perl_substr(seq, start, len)` 或专用窗口函数，避免在主逻辑里散落 `+1/-1`。

### 5.8 splice candidate 聚类

Perl 的 `cluster_reads` 两层聚类：

1. 先按 `chr, site2` 排序，`site2` 差值小于等于 3 的候选归入同一 group；
2. 再在 group 内按 `site1` 排序，`site1` 差值小于等于 3 的候选归入 splice cluster；
3. cluster 坐标取中位 read 的 `site1/site2`；
4. 按 3 个 CIGAR slot 统计 distinct CIGAR 数，三槽总数达到 `strigency` 才保留。

Perl 默认 `strigency = 1`。Rust 可先固定该默认值，后续再暴露 CLI 参数。

### 5.9 circ 内 cirexon 预测

对每个 circRNA：

1. 根据 reference motif 或注释推断 circ strand；
2. 找出落在 circ 内的 splice clusters；
3. 汇总 junction read 在 circ 内的 mapping range 和 junction read coverage；
4. 根据 splice cluster 的 start/end 支持构造 candidate exon；
5. 用 coverage validation 判断 exon 是否有效；
6. 判断是否为 ICF。

候选 exon 条件来自 Perl：

```text
end >= start + min_exon_length - 1
end <= start + max_exon_length - 1
start/end 至少一侧有支持，或位于 circ 起止边界
```

默认：

- `min_exon_length = 20`
- `max_exon_length = 2000`
- `min_intron = 70`

### 5.10 coverage validation 与 intron retention

Perl 使用 coverage profile 和 Mann-Whitney U-test 风格统计验证 exon / IR：

- `exon_coverage_validation_single`
- `intron_retention_validation`
- `Z_calculation`

筛选条件包括：

- candidate 内是否有 coverage 0；
- junction-read coverage gap；
- candidate 两侧背景 coverage；
- splice-site 支持 read 数；
- U-test Z 值是否超过 `Z_alpha`。

默认 `U_test_strigency = 1`，对应：

```text
Z_alpha = 1.6449
```

Rust 实现时建议先保留 Perl 的返回码语义，例如 `-1/-2/-8/-9/1/2`，并在注释中记录每个码的含义。

### 5.11 exon path 与 AS 分类

Perl 对 validated exons 构建 circ 内路径：

- 从 `circ_start` 开始；
- 根据 splice junction 或相邻 exon 连接；
- 保留以 `circ_end` 结束的路径；
- 如果存在多条 path，则比较 path 中不共享的 key units。

AS 类型输出规则：

- key unit 被其他 path 整体跳过：`ES`
- 相同 start，不同 end：按 circ strand 判定 `A5SS` 或 `A3SS`
- 相同 end，不同 start：按 circ strand 判定 `A5SS` 或 `A3SS`
- validated intron retained exon：`IR`

输出格式：

```text
circRNA_id
alternatively_spliced_exon
AS_type
psi_estimation_without_correction
psi_estimation_after_correction
```

### 5.12 PSI 与 insert length 校正

Perl 的 PSI 分两类：

- raw PSI：直接用 splice-in / splice-out read count 计算；
- corrected PSI：基于 paired-end insert length distribution 进行 normalization。

校正依赖：

- 从长 exon 中估计 insert length distribution；
- `paired_read_end_mapping`
- `normalize_count4splice`
- `normalize_count4splice_intron_retention`

建议分两步实现：

1. 第一版输出 raw PSI，并在 corrected PSI 无法计算时输出 `n/a`；
2. 第二版再完整迁移 insert length correction。

这样可以先验证 ES/A5SS/A3SS/IR 分类是否正确。

## 6. Rust 模块拆分建议

建议新增：

```text
src/ciri_as.rs
src/as_cigar.rs
src/as_model.rs
```

职责：

- `ciri_as.rs`：CIRI-AS runner，负责阶段调度、I/O、输出。
- `as_cigar.rs`：CIRI-AS 专用 CIGAR/clip 解析与 split candidate 识别。
- `as_model.rs`：AS 阶段数据结构。

也可以先只建 `src/ciri_as.rs`，等逻辑变大后再拆分。但从维护性看，CIGAR 与模型结构最好早拆。

需要扩展：

- `src/main.rs`：新增 `--as` 等参数，在 Summary 后调用。
- `src/lib.rs`：导出 `ciri_as` 模块。
- `src/annotation.rs`：补充 AS 所需 exon records。
- `src/utils.rs` 或新工具模块：新增 AS 输出路径 helper。

## 7. 开发顺序建议

### Phase 0：建立 Perl golden

先用当前测试数据跑 Perl CIRI-AS，得到参考输出。

注意：需要先把 CIRI-toolkit 当前 `.out` 转成 Perl 兼容格式，确保 `junction_reads_ID` 是最后一列，否则 Perl 会误读 `Score`。

建议 golden 文件放在：

```text
tests/chr1/CIRI-AS/
```

或临时验证时放在：

```text
tmp/ciri_as_gold/
```

### Phase 1：输入解析与 circ cluster

目标：

- 能读取当前 `.out`；
- 能建立 `junction_read` / `circ_junc_reads` / `circ_info`；
- 能复现 Perl circ cluster 数量、范围和 circ 归属。

验证：

- 输出 cluster manifest；
- 与 Perl `_cluster.list` 或自制 probe 对齐。

### Phase 2：CIGAR helper 和 read mapping detail

目标：

- 移植 `MSID`；
- 移植 `MSID_start`;
- 移植 `record_mapping_detail`。

验证：

- 针对典型 CIGAR 建单元测试：`M`、`MS`、`SM`、`SMS`、`MDM`、`MIM`、含 `H`。
- 确认 read genomic range、read start/end 与 Perl 一致。

### Phase 3：SAM/BAM 重扫与候选 splice

目标：

- junction read 和 cluster-overlap read 均能进入 mapping check；
- 输出候选 internal splice 数；
- coverage 与 junction read mapping range 可 dump。

验证：

- 与 Perl log 中 candidate splice junction 数对齐；
- 在小 fixture 上比较 coverage 列表。

### Phase 4：splice signal 与 splice cluster

目标：

- 实现 `splice_loci_check` 等价窗口提取；
- 使用 annotation-aware motif / strand / offset tie-break；
- 复现 `cluster_reads`；
- 生成 `<prefix>_splice.list`。

验证：

- CIRI3 主流程仍必须保持 `.out` 的 circ/read/read-assignment/FSJ 四层零差异；
- AS 阶段先比较候选数、motif validated 数、典型 read 的坐标解释；
- `_splice.list` 可与 Perl CIRI-AS 作为参考对照，但不要求 100% parity；若差异来自 annotation-aware tie-break 或 Perl hash 顺序，应记录为有意偏离。

### Phase 5：cirexon 输出（历史方案，当前不再实现）

历史目标：

- 实现 circ strand 推断；
- 实现 candidate exon 构造；
- 实现 exon coverage validation；
- 输出 `<prefix>.list`。

当前 Rust 主线已删除旧 cirexon/list writer。内部结构推断改由 `<prefix>.segments` parser、circ-local graph、annotation-guided block projection 和 read-supported junction/exclusive evidence 完成。若后续需要输出 candidate exon 审计，也应作为新的 isoform usage sidecar 设计，而不是恢复 CIRI-AS `.list` 协议。

验证：

- 比较 cirexon start/end、support read、coverage median、ICF。

### Phase 6：full-length isoform path 输出（已由 major isoform 路线替代）

历史目标：

- 从 validated cirexon 构建 circ-local exon graph；
- 以 Summary confirmed BSJ start/end 为硬锚点；
- 用 BSJ-read-supported internal splice 作为显式边；
- 对没有显式 splice 边的 exon 使用 CIRI-AS 相邻 exon fallback；
- 枚举从 circ start 到 circ end 的 anchored full-length path；
- 输出 `<prefix>.isoforms`、`<prefix>.isoform_summary` 与 `<prefix>.fa`。

当前实现状态：

- 不输出 CIRI-AS `.isoforms/.isoform_summary/.fa`。
- 默认输出 `<prefix>.isoforms.gtf` 和高可信 `<prefix>.isoforms.fa`。
- 每个 circRNA 先输出一个 major isoform；GTF 保留 mature/estimate、coverage、structure_hash 和 FASTA eligibility 审计字段。
- 下一阶段扩展 candidate isoform search space 与 usage 计算，但输出协议另行设计，不复刻 Perl 文件格式。

验证：

- 比较同一 circ 下 isoform 数量、exon chain、junction chain 和 FASTA 长度；
- 使用 simulator truth 时优先比较 `exon_chain` 和 `isoform_len`，暂不比较 AS event type。

### Phase 7：isoform usage / side evidence（当前下一阶段）

目标：

- 在当前 `<prefix>.segments -> circ-local graph` 上建立候选 isoform search space；
- 为候选 isoform 增加 read assignment / fractional usage 估计；
- 输出 per-sample support、usage、confidence tier、structure_hash 和 major/minor 状态；
- 为 multi-sample integration 预留结构合并、major isoform switching 和样本间 usage 比较。

验证：

- 比较 isoform-level read assignment；
- candidate path 不得回写主流程 BSJ 判断；
- 使用 simulator truth 同时评估 major accuracy、minor recovery、usage rank correlation 和 false sequence candidate rate。

## 8. 验证策略

新增 segments/isoform/usage 后应保持两类验证分离：

1. CIRI3 主流程 parity：
   - `scripts/ciri_result_diff.py` 继续用于 `.out`；
   - `circ/read/read-assignment/FSJ` 必须保持零差异。

2. segments / isoform / usage 对照：
   - `scripts/ciri_segments_eval.py` 继续用于 simulator read-level segments 对照；
   - major isoform 继续比较 exon chain、junction chain、FASTA length 和 truth major/exact-any；
   - 后续 usage 模块再增加 isoform-level assignment、minor recovery 和 usage rank correlation；
   - CIRI-AS v1.2 输出只作为历史参考，不作为强制 parity oracle。

建议 AS diff 分层：

```text
splice-level:
  (chr, splice_start, splice_end)

cirexon-level:
  (circRNA_id, cirexon_start, cirexon_end)

AS-event-level:
  (circRNA_id, alternatively_spliced_exon, AS_type)

PSI-level:
  raw PSI / corrected PSI with exact or formatted-string comparison
```

当前优先保证 `<prefix>.segments` 与 simulator `.reads.tsv` 的 read-level 解释稳定，再比较 major isoform 和未来 usage。

## 9. 性能与内存注意事项

历史 CIRI-AS 路线的热点包括：

- 重扫 BAM/SAM；
- per-base coverage 统计；
- read-group 内多 alignment 两两比较；
- circ 内 exon path 构建；
- insert length normalization。

当前主线不再实现旧 CIRI-AS writer。后续 usage/multi-sample 若新增大数据路径，仍应遵守同一原则：先用可验证的结构完成行为验证，再针对全量 BAM/SAM 或多样本规模做内存上界优化。

- 用 interval index 替代 per-base `circ_cluster_range`；
- coverage 改为分 chromosome sparse vector 或 run-length structure；
- read group 重扫复用 Scan1/Scan2 的 BAM shard 框架；
- 对 circ cluster 分 chromosome 并行处理；
- 大样本下按 chromosome spill 中间状态，降低 RSS。

## 10. 开发红线

- segments / isoform / usage 必须是 Summary 后处理或 sidecar 阶段，默认不改变现有 CIRI3 `.out/.bsj` 判定。
- 不允许为了 AS 改动 Scan1/Scan2/Summary 已对齐逻辑，除非有独立回归证明主流程零差异。
- 迁移初期以 Perl 脚本的证据链和坐标系统为参考；遇到 Perl hash 顺序、未利用 annotation boundary 等不稳定或弱启发式时，允许采用更稳定的 Rust 规则，但必须在本文档记录。
- 坐标系统必须文档化：Perl 逻辑使用 1-based genomic coordinate 和 inclusive end。
- CIGAR helper、splice signal 检查、coverage validation、AS classification 都需要 targeted unit tests。
- 临时文件、probe 输出、golden 生成脚本统一放入 `tmp/` 或明确的测试目录，不散落到仓库根目录。

## 11. 当前未决问题

旧 CIRI-AS 输出命名、corrected PSI、AS event taxonomy 和 `.list/.isoform_summary` 兼容性不再是当前路线问题。后续开发前需要进一步确认的是：

- multi-isoform usage sidecar 的文件名、字段和是否需要单独 FASTA eligibility 表；
- backward/outward reads 在多个 circRNA 和多个 isoform 之间的 fractional assignment 规则；
- multi-sample integration 的输入列表格式、sample metadata 字段和 structure merge key；
- 是否在 usage 阶段引入新的 weak evidence，例如 boundary clip 或 per-circ mini-BAM review 结果。

---
最后更新：2026-05-19
