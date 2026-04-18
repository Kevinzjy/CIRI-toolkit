# RO feature 与 circRNA 全长结构识别计划

本文档记录 CIRI-toolkit 后续引入 RO（Reverse Overlap）feature 的设计计划。当前目标不是完整复刻 CIRI-AS / CIRI-full 的全部输出，而是提取其中对 CIRI-toolkit 有价值的思路：先在 `Scan1` 过程中识别 paired-end reads 的 RO 序列，输出 `.ro.fq`，再为后续 RO remap、二次扫描和 circRNA isoform reconstruction 提供证据。

## 1. 当前设计口径

CIRI-AS / CIRI-full 在当前阶段只作为算法参考，不作为 100% parity 目标：

- CIRI3 核心流程 `Scan1 -> Scan2 -> Summary` 仍保持 Java parity 约束。
- CIRI-AS / CIRI-full 相关新功能不要求复刻全部中间文件、历史 header、调试输出和未验证的启发式结果。
- 新功能只保留与目标直接相关的证据：
  - RO read pair 识别。
  - RO merged read。
  - RO read 的后续 mapping 结果。
  - 同一 BSJ 下的 isoform 结构与 usage ratio。

第一阶段先实现：

```text
Scan1 read group -> RO detector -> <prefix>.ro.fq
                                  -> <prefix>.ro.tsv
```

不在第一阶段改变主 `.bsj1/.bsj/.out`。

## 2. 为什么先做 RO

RO feature 适合作为 CIRI-AS / CIRI-full 整合的第一步：

- 它可以直接复用 CIRI-full RO1 的核心思路，边界相对清晰。
- 它发生在 read pair 层面，与当前 `Scan1` 的 read group 组织方式天然匹配。
- `.ro.fq` 可以作为明确的中间产物，后续既能用于 BWA remap，也能用于 full-length isoform reconstruction。
- 在不改主判定结果的前提下，能先评估 RO 对新增 BSJ、低 mappability read 和 full-length structure 的贡献。

## 3. Scan1 集成点

当前 `Scan1` 已经按 read id 聚合 read group，并在 group 内区分 mate：

```text
read group
  pair1 alignments
  pair2 alignments
  stand_map: mate -> representative sequence
```

RO detector 应插在 read group 层：

```text
process one read group:
  1. 按现有逻辑执行 Scan1 BSJ candidate 判断。
  2. 如果 read1/read2 都有可用原始序列和质量值：
       run RO detector
       if RO accepted:
         write merged read to <prefix>.ro.fq
         write metadata to <prefix>.ro.tsv
```

第一阶段的约束：

- 不影响当前 Scan1 `.bsj1` 输出。
- 不改变 `IsBSJHg2` 的判定输入。
- `.ro.fq` 和 `.ro.tsv` 作为 sidecar 输出。
- 默认关闭，通过显式参数开启。

建议 CLI：

```text
--ro-feature
--ro-min-identity 95
--ro-min-overlap 13
--ro-output <prefix.ro.fq>
```

默认输出：

```text
<prefix>.ro.fq
<prefix>.ro.tsv
```

## 4. RO 检测逻辑

第一版复用 CIRI-full RO1 的核心 overlap 规则，但输出面向 CIRI-toolkit：

- 使用 read1 原始测序方向序列。
- 使用 read2 原始测序方向的 reverse-complement。
- 默认 identity 阈值 `95`。
- 默认最小 overlap 长度 `13`。
- seed 长度 `8`。
- seed 至少 `7/8` 个非 N 碱基匹配。
- first-hit 即接受，不寻找全局最优 overlap。

CIRI-full 兼容方向：

```text
seq1 = R1
seq2 = reverse_complement(R2)

for indexn in [0, 1, 2]:
  seed = seq1[indexn..indexn+8]
  scan seed in seq2
  if overlap identity >= 95 and overlap_len >= 13:
    merged = seq2[0..startpoint] + seq1[indexn..]
```

这个方向对应第一阶段的 `5p_ro` 类型。

为了支持用户计划中的类型标注，第一版建议额外实现一个对称方向：

```text
seq1 = reverse_complement(R2)
seq2 = R1
```

该方向对应 `3p_ro` 类型。若同一 read pair 在两个方向都通过阈值，标记为 `bidirectional_ro` 或 `full_length_candidate`。

重要说明：

- `full_length_candidate` 是 sequence-level 标记，只表示 read pair 在 RO 层面具备更完整结构的可能。
- CIRI-full 中真正的 `Full/Part` 判定依赖 RO read remap 后的 exon chain；不能在 Scan1 RO detector 阶段直接等价输出。
- 后续只有经过 RO BAM / RO Scan2 或 isoform reconstruction 验证后，才能升级为 full-length circRNA structure evidence。

## 5. RO 类型定义

第一阶段建议使用以下类型：

| 类型 | 含义 | 是否最终 full-length 证据 |
| --- | --- | --- |
| `5p_ro` | CIRI-full RO1 方向的 reverse overlap | 否，需后续 mapping 验证 |
| `3p_ro` | 对称方向的 reverse overlap | 否，需后续 mapping 验证 |
| `bidirectional_ro` | 同一 read pair 同时满足 `5p_ro` 和 `3p_ro` | 否，标记为 `full_length_candidate` |
| `full_length_candidate` | `.ro.tsv` 中的附加状态，表示 sequence-level 可能覆盖完整 circRNA | 否，需 RO remap 验证 |

不要在第一阶段直接输出 CIRI-full RO2 的 `Full/Part`。该分类应留到 RO read mapping 之后。

## 6. `.ro.fq` 输出设计

FASTQ header 建议保留原始 read id，并附带 RO metadata：

```text
@<read_id> ro_type=5p_ro identity=97 overlap=18 r1=1-18 r2=43-60 source=scan1
ACGT...
+
IIII...
```

字段：

- `read_id`：原始 read id，不加随机后缀，便于后续按 origin read 去重。
- `ro_type`：`5p_ro`、`3p_ro` 或 `bidirectional_ro`。
- `identity`：overlap identity。
- `overlap`：overlap 长度。
- `r1`：R1 overlap 坐标，1-based 闭区间。
- `r2`：R2 reverse-complement 后的 overlap 坐标，1-based 闭区间。
- `source=scan1`：标记来源。

如果同一个 read pair 产生多个 RO 类型，第一版建议只输出优先级最高的一条：

```text
bidirectional_ro > 5p_ro > 3p_ro
```

这样可以避免同一 origin read 在后续 BWA/Scan2 中重复计数。后续如果需要保留多候选，可在 `.ro.tsv` 中记录全部候选，但 `.ro.fq` 仍只保留 selected candidate。

## 7. `.ro.tsv` 输出设计

`.ro.tsv` 用于调试、统计和后续 evidence 追踪：

```text
read_id
ro_type
selected
identity
overlap_len
r1_overlap_start
r1_overlap_end
r2_overlap_start
r2_overlap_end
merged_len
r1_len
r2_len
group_mapped_state
scan1_bsj_hit
notes
```

其中：

- `selected`：是否写入 `.ro.fq`。
- `group_mapped_state`：`both_mapped`、`one_unmapped`、`both_unmapped`、`mixed`。
- `scan1_bsj_hit`：当前 read group 是否已经被 normal Scan1 判定为 BSJ read。
- `notes`：保留异常原因，例如 `missing_quality`、`short_read`、`ambiguous_pair`。

## 8. 质量值与原始方向

当前 `AlignmentRecord` 只保存：

```rust
flag, chrom, pos, mapq, cigar, seq
```

输出 `.ro.fq` 必须有 quality，因此第一阶段需要为 RO detector 补齐质量值来源。建议不要把 RO 所需字段全部塞进现有 BSJ hot path，而是新增 group-level read representation：

```rust
struct MateRead {
    mate_index: i32,
    seq: Vec<u8>,
    qual: Vec<u8>,
    flag: i32,
}
```

构建规则：

- 从 SAM/BAM record 中读取 `SEQ` 和 `QUAL`。
- 对同一 mate 的多条 alignment，选择最长 `SEQ` 作为 representative。
- 如果最长 `SEQ` 并列，保留先出现记录，保持输出稳定。
- 如果 record flag 表示 reverse-complement alignment，则恢复为原始测序方向：
  - `seq = reverse_complement(seq)`
  - `qual = reverse(qual)`

这一点与当前 Scan1 的 alignment-oriented sequence 逻辑不同，需要明确隔离，避免影响 Java parity 的 BSJ 判定路径。

如果输入缺少 quality：

- 第一版建议该 read pair 不写入 `.ro.fq`，只在 `.ro.tsv` 中记录 `missing_quality`。
- 不建议默认填充固定质量值，因为后续 BWA scoring 和可重复性会变得不透明。

## 9. mapped / unmapped read pair 的处理

第一阶段统一规则：

- 只要 read1/read2 能形成 RO，就可以写入 `.ro.fq`。
- `mapped` 状态只写入 `.ro.tsv`，不改变是否输出 merged RO read。
- 不在 Scan1 阶段直接用 mapped RO pair 增加 BSJ read count。

原因：

- RO overlap 本身不是 BSJ 证据。
- mapped pair 的局部位置可作为 evidence，但仍需经过 BSJ validator 或 RO remap 验证。
- 直接把 RO overlap 计入 BSJ 会引入双计数和 false positive 风险。

后续第二阶段可以增加：

```text
<prefix>.ro.fq -> bwa mem -> <prefix>.ro.bam -> RO Scan1/Scan2
```

## 10. 后续 RO remap 计划

RO remap 阶段建议使用外部 `bwa mem` 子进程，而不是在 Rust 内部实现 BWA：

```text
Scan1 RO writer -> <prefix>.ro.fq
bwa mem -T 19 -t <threads> <ref.fa> <prefix>.ro.fq
  -> <prefix>.ro.bam
```

后续可以优化为流式 worker：

```text
Scan1 worker threads
  -> bounded channel<RoFastqRecord>
  -> bwa stdin writer
  -> bwa stdout SAM parser
  -> BAM writer
```

但第一版建议先落地文件型 `.ro.fq`，便于验证和复现。

`<prefix>.ro.bam` 的用途：

- 作为 RO read 的独立 mapping 证据。
- 用于 RO Scan1 生成 `<prefix>.ro.bsj1`。
- 用于 RO Scan2 rescue 已有或新增 candidate。
- 后续用于同一 BSJ 下 isoform chain reconstruction。

## 11. 与 Scan2 / Summary 的关系

第一阶段：

```text
normal Scan1 -> normal .bsj1
RO detector  -> .ro.fq / .ro.tsv
normal Scan2 -> normal .bsj
Summary      -> normal .out
```

主结果不变。

第二阶段：

```text
.ro.fq -> bwa -> .ro.bam
.ro.bam -> RO Scan1 -> .ro.bsj1
normal .bsj1 + .ro.bsj1 -> candidate index
Scan2 original BAM + Scan2 RO BAM -> merged .bsj
Summary 按 origin read id 去重
```

正式合并 RO evidence 前必须实现：

- origin read id 去重。
- evidence source 标记：`original` / `ro`。
- 同一 BSJ 下同一 origin read 只能计一次。
- `.out` 可选增加 RO support 统计，或输出 sidecar evidence 文件。

## 12. 实现阶段

### 阶段 1：RO detector 与 `.ro.fq`

- 新增 `src/ro.rs`。
- 实现 CIRI-full-compatible `5p_ro` 检测。
- 实现对称 `3p_ro` 检测。
- 实现 RO 类型选择。
- 为 Scan1 read group 补齐原始方向 `seq/qual`。
- 输出 `<prefix>.ro.fq` 和 `<prefix>.ro.tsv`。
- 不改变主结果。

阶段 1 进一步拆成以下开发任务：

1. 新增 `src/ro.rs`
   - 定义 `RoDetectConfig`、`RoType`、`RoCandidate`、`RoDetectionResult`。
   - 实现 CIRI-full-compatible first-hit overlap detector。
   - 实现 `5p_ro`、`3p_ro` 与 selected candidate 选择。
   - 在模块级 `//!` 注释中说明：该模块只产生 sequence-level RO sidecar evidence，不参与 BSJ 判定。

2. 增加 RO 单元测试
   - 覆盖 first-hit、整数 identity、`N`、最小 overlap、seed match、quality merge、`5p_ro/3p_ro/bidirectional_ro`。
   - 测试应优先构造短序列直接验证，不依赖 BAM/SAM 文件。

3. 扩展 Scan1 read group 的 RO 输入
   - 新增 RO 专用 mate representation，保存原始测序方向 `seq/qual`。
   - 从 SAM/BAM record 读取 `SEQ` 与 `QUAL`。
   - 对同一 mate 选择最长 representative；长度相同保留先出现记录。
   - 若 alignment flag 为 reverse-complement，则恢复原始方向：`seq = reverse_complement(seq)`，`qual = reverse(qual)`。
   - 该结构必须与现有 `stand_map` 隔离，不能改变 Java parity 的 BSJ hot path。

4. 实现 shard-local RO writer
   - SAM 路径可先写单一 `<prefix>.ro.fq/.ro.tsv`。
   - BAM 分片路径必须写 shard-local part 文件，最后合并，避免并发写锁和输出交错。
   - `.ro.tsv` 只写一次 header；part 文件合并时保留最终文件单 header。

5. 接入 CLI
   - `--ro-feature` 默认关闭。
   - `--ro-min-identity` 默认 `95`。
   - `--ro-min-overlap` 默认 `13`。
   - 默认输出 `<prefix>.ro.fq` 和 `<prefix>.ro.tsv`。

6. 主流程回归
   - `--ro-feature` 关闭时，主结果与当前基线保持一致。
   - `--ro-feature` 开启时，主 `.bsj1/.bsj/.out` 仍保持一致，只额外生成 RO sidecar 文件。

### 阶段 2：RO remap sidecar

- 增加 CLI 或手动步骤运行 `bwa mem`。
- 生成 `<prefix>.ro.bam`。
- 验证 `.ro.bam` 可被现有 Scan1/Scan2 读取。
- 输出 `<prefix>.ro.bsj1`，暂不合并主结果。

### 阶段 3：RO-assisted BSJ

- 合并 normal `.bsj1` 与 `.ro.bsj1` 的候选位点。
- Scan2 支持额外扫描 RO BAM。
- Summary 增加 origin read 去重。
- 输出 RO support sidecar。

### 阶段 4：full-length isoform reconstruction

- 使用 original BAM + RO BAM 的 split junction / RO chain evidence。
- 对同一 BSJ 构建多个 isoform。
- 计算 isoform usage ratio。
- 输出 `<prefix>.isoform.tsv`。

## 13. 测试与验证

第一阶段测试重点：

- RO overlap first-hit 行为。
- identity 整数除法语义。
- `N` 碱基不计入有效 match。
- R2 reverse-complement 和 quality reverse。
- `5p_ro / 3p_ro / bidirectional_ro` 类型选择。
- 同一 read pair 只写一个 selected `.ro.fq` record。
- `--ro-feature` 关闭时，主 `.bsj1/.bsj/.out` 字节级不变。

建议新增测试：

```text
ro_detect_5p_ciri_full_compatible
ro_detect_3p_symmetric
ro_detect_bidirectional_priority
ro_fastq_header_metadata
ro_scan1_disabled_no_output_change
```

验证口径：

- 开启 `--ro-feature` 前后，主结果必须一致。
- `.ro.tsv` 的候选数量和类型可单独统计。
- `.ro.fq` 可以独立用 `bwa mem` 比对，不依赖主流程。

## 14. 代码注释与文档要求

RO1 实现必须遵守仓库 `AGENTS.md` 的 Rust 注释规范，并额外满足以下要求：

- `src/ro.rs` 必须有模块级 `//!` 注释，说明它在 `Scan1 -> Scan2 -> Summary` 主流程之外，只负责 RO sidecar evidence。
- `RoDetectConfig`、`RoType`、`RoCandidate`、RO detector 入口函数必须有 Rust doc comment。
- overlap detector 的 first-hit 行为、整数 identity、`N` 处理、`indexn in [0,1,2]` 这些 CIRI-full RO1 来源的边界必须有注释说明。
- Scan1 中新增的 RO mate representative 结构必须解释为什么和 `stand_map` 分离：一个服务原始测序方向 RO 检测，一个服务 Java parity 的 alignment-oriented BSJ 判定。
- shard-local writer 与合并逻辑必须解释为什么不能多线程共享同一 writer。
- CLI 参数和输出文件若发生变化，必须同步更新：
  - `docs/05-ro-feature-plan.md`
  - `docs/01-development-status.md`
  - 必要时更新 `README.md`，但只在功能正式面向用户时更新。

## 15. 实现前检查清单

进入代码实现前，需要确认以下事项：

- `--ro-feature` 第一版只生成 sidecar，不改变主结果。
- `.ro.fq` 每条 read 只保留一个 selected RO candidate。
- `.ro.tsv` 可以记录调试信息，但不作为主流程输入。
- `full_length_candidate` 只是 sequence-level 状态，不等价于 RO2 的 `Full`。
- 缺少 quality 的 read pair 不写 `.ro.fq`，只在 `.ro.tsv` 记录跳过原因。
- BAM 分片输出使用 part 文件，最终合并。
- 所有新增代码补齐文档注释和关键行内注释。

## 16. 当前结论

下一步最合理的开发切入点是：

```text
Scan1 sidecar RO detector + .ro.fq/.ro.tsv output
```

这一步足够小，不会破坏当前 CIRI3 parity；同时又能为后续 RO remap、second scanning 和 isoform usage 分析打下可验证的中间层。
