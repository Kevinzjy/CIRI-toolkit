# CIRI-toolkit 开发状态

本文档用于描述 Rust 版本 CIRI-toolkit 的当前实现进度、对齐基线、验证口径与长期维护约束。

## 项目概述
CIRI-toolkit 是 CIRI3（Java）的 Rust 复现版本，目标是在保持判定逻辑一致的前提下提升性能与可扩展性。

## 相关文档
- 调试手册：`docs/02-parity-debug-playbook.md`
- 性能总结：`docs/03-performance-optimization.md`
- RO feature 计划：`docs/05-ro-feature-plan.md`
- 模拟数据与 truth 输出设计：`docs/06-simulation-truth-design.md`
- CIRI-AS-style full-length 结构识别设计：`docs/07-full-length-reconstruction.md`
- CIRI-AS 拆解：`docs/CIRI-AS.md`
- CIRI-full 拆解：`docs/CIRI-full.md`
- 文档导航：`docs/00-index.md`

## 当前状态（Functional v1）
- 核心流程 `Scan1 -> Scan2 -> Summary` 已完成并稳定可用。
- 已实现 SAM/BAM 自动识别与对应执行路径。
- BAM 原生支持已并入主流程，不再单独维护历史归档文档。
- `stringency` 过滤、FSJ 统计、注释输出等关键能力均可运行。
- 双端（paired-end）与单端（single-end）输入均已完成验证，当前输出与 Java CIRI3 保持一致。
- `Scan1` 与 `Scan2` 的性能优化阶段已完成，当前版本整体性能已达到并超过 Java 基线，可作为正式功能版本持续使用。

## 对齐基线
- circRNA-level：100%
- read-level：100%
- read-assignment-level：100%
- FSJ-level：100%

## 最新验证结论（2026-05）
- `tests/chr1` BAM / SAM 基线继续保持 circ/read/read-assignment/FSJ 四层 100% 对齐。
- `<prefix>.bsj1` 和 `<prefix>.bsj2` 已升级为 mate-level BSJ 中间文件，包含 `mate_label` 与 `priority`；`Summary` 只消费 `priority=1` 行，最终 `<prefix>.bsj` 保留 `priority=0` 行供后续内部 splice site 识别使用。
- 当前运行路径不再生成额外 `<prefix>.scan1.tmp`、`<prefix>.scan2.tmp` 或 `<prefix>.bsj.raw.tmp`。
- 单端测序数据的专项验证已完成，当前输出与 Java 基线完全一致。
- 大规模 hg38 验证表明：当 Rust 与 Java 使用一致的参考 FASTA / GTF 版本时，`BSJ reads` 与 `FSJ counts` 已完成全量对齐。
- 当前推荐的 hg38 对齐口径为：
  - FASTA：`/data/public/database/gencode/hg38/_BWAindex/hg38.fa`
  - GTF：`/data/public/database/gencode/hg38/gencode.v44.annotation.gtf`
- 最近排查表明：注释版本差异，尤其是 exon 覆盖是否包含对应位点，本身就足以制造表面上的 parity gap；这类差异不应误判为 Scan1/Scan2 核心逻辑偏移。
- 全量 hg38 仍保留极少数 family-level 残余差异，但逐例复核后 Rust 的局部判定更合理；这些 case 作为 Java context-sensitive artifact 归档，不再阻塞当前版本发布。
- 因此，全量 parity 复核时必须先锁定以下环境：
  - 相同 FASTA 版本
  - 相同 GTF 版本
  - 相同 stringency 参数
  - 关闭所有 trace/profile 环境变量

统一校验命令：

```bash
python scripts/ciri_result_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 架构要点（已落地）
- Scan2 使用 Java 语义一致的去重索引与桶遍历顺序。
- `Scan1 -> Scan2 -> Summary` 的算法输入与 mate-level 展示证据通过 `priority` 隔离：`priority=1` 保持 CIRI3-compatible 主流程，`priority=0` 只作为附加 mate evidence 保留。
- Scan2 的 FSJ 统计现已对齐 Java `GetFSJClass.getFSJ(...)` 的 bucket-gating 与 BAM mate-run 语义。
- Summary 在关键路径中对齐 Java 容器迭代行为，避免首命中差异。
- 大文件路径采用 mmap/分片与批量写出策略，兼顾吞吐与内存占用。
- SAM/BAM 两条输入路径最终落到相同的判定语义上，避免格式差异带来的输出漂移。
- `Scan2` 在进入 `Summary` 前会主动释放候选索引工作集，以控制 whole-genome 场景下的 RSS。

## 调试口径
- 先跑 `scripts/ciri_result_diff.py`，再决定是查 `Scan1`、`Scan2` 还是 `Summary`。
- 追踪单条或少量 read 时，统一使用：
  - `CIRI_TRACE_READS`
  - 必要时加 `CIRI_TRACE_ALL_CANDS=1`
  - 必要时加 `CIRI_TRACE_HG2=1`
- 需要性能热点分布时，统一使用：
  - `CIRI_PROFILE_SCAN1=1`
  - `CIRI_PROFILE_SCAN2=1`
- 详细命令、日志标签与 SOP 见 `docs/02-parity-debug-playbook.md`。

## 后续开发计划
- 高优先级功能点：在当前 `priority` 协议基础上继续推进真正的 paired-end R1/R2 BSJ 决策。当前默认路径仍采用 CIRI3 read-pair first-hit 语义；后续更合理的方向是先做 mate-level candidate 收集，再做 pair-level 决策，最后仍按 read-pair-level 计数。R1/R2 支持同一 BSJ 时应记录一致支持，支持不同 BSJ 时应显式标记冲突或 ambiguity。
- 当前优先路线转为 CIRI-AS-style read-level segments reconstruction：以 Summary confirmed BSJ 为锚点，在后处理阶段重扫原 BAM/SAM，先稳定恢复 circRNA 相关 reads 的 segment chain，再把这些结果作为 future full-length reconstruction 的输入层。
- 当前默认后处理不依赖 RO remap，不新增 circ seed，不改变 `priority=1` 主流程证据和 `.out`；当前正式 sidecar 输出先收敛为 `<prefix>.segments`。
- CIRI-AS / CIRI-full 尚未经过本项目同等级别的严格验证和性能优化，后续不以完整复刻其所有历史输出为目标；CIRI-AS 主要提供内部结构识别思路，CIRI-full/RO 暂作为未来 side evidence 风险清单。
- CIRI-AS sidecar 的 splice signal 阶段已明确采用 annotation-aware motif / strand / offset tie-break：当 annotation exon boundary 支持某个解释时，优先于 Perl v1.2 的 `AC/CT elsif AG/GT` 与 hash-order offset 行为；这是有意偏离 Perl 输出、追求稳定和边界可解释性的设计选择。
- CIRI-AS sidecar 的 splice signal 阶段、annotation-aware motif / strand / offset tie-break 与 read-group 识别逻辑继续保留，用于支撑 `<prefix>.segments` 的 best-chain 选择和 segment 边界解释。
- 当前 full-length 目标暂时后移；CLI 不再使用 `--as`，而是在写完 `<prefix>.out` 后默认输出 `<prefix>.segments`。当前阶段优先验证 `circ_id`、`is_circular` 和 `r1/r2_segments` 是否能与 simulator truth 对齐。
- 独立 Rust 模拟器已作为开发 fixture 基础接入，通过 `ciri-simulator` 生成 paired FASTQ、linear-only `<prefix>.annotation.gtf`、circ/isoform-level `<prefix>.isoforms.tsv` 和 read-pair-level `<prefix>.reads.tsv`，用于后续 internal structure 和 full-length path 验证。
- 新增模块仍应遵守小步验证原则：segments 输出优先、主流程 parity 不回退，等 read-level truth 稳定后再恢复 path-level full-length 输出。

## 开发规范
- 后续代码改动应继续遵守已沉淀的 hg38 排障结论、Java parity 约束与关键注释说明，避免在重构或优化中重新引入回归。
- 后续若再进行性能优化或结构调整，仍必须执行完整 parity 校验，出现差异即先回到行为对齐。
- full-length reconstruction 实现前应先阅读 `docs/07-full-length-reconstruction.md`；新增 bundle、graph、FASTA 输出等模块必须补齐 Rust 文档注释，说明职责边界、证据层级和不影响主流程 parity 的原因。
- CIRI-AS / CIRI-full 相关实现当前不以完整 parity 为目标，任何取舍必须写入对应设计文档，避免后续误把历史输出当作强制兼容契约；尤其是 annotation-aware offset/strand 判断这类有意偏离，必须在代码注释和 `docs/CIRI-AS.md` 同步维护。

---
最后更新：2026-05-04
