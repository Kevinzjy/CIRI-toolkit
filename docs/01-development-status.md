# CIRI-toolkit 开发状态

本文档用于描述 Rust 版本 CIRI-toolkit 的当前实现进度、对齐基线、验证口径与长期维护约束。

## 项目概述
CIRI-toolkit 是 CIRI3（Java）的 Rust 复现版本，目标是在保持判定逻辑一致的前提下提升性能与可扩展性。

## 相关文档
- 调试手册：`docs/02-parity-debug-playbook.md`
- 性能总结：`docs/03-performance-optimization.md`
- RO feature 历史设计：`docs/05-ro-feature-plan.md`
- 模拟数据与 truth 输出设计：`docs/06-simulation-truth-design.md`
- segments 与 major isoform 结构识别设计：`docs/07-full-length-reconstruction.md`
- CIRI-AS 拆解（历史参考）：`docs/CIRI-AS.md`
- CIRI-full 拆解（历史参考）：`docs/CIRI-full.md`
- 文档导航：`docs/00-index.md`

## 当前状态（Functional v1）
- 核心流程 `Scan1 -> Scan2 -> Summary` 已完成并稳定可用。
- 已实现 SAM/BAM 自动识别与对应执行路径。
- BAM 原生支持已并入主流程，不再单独维护历史归档文档。
- `stringency` 过滤、FSJ 统计、注释输出等关键能力均可运行。
- 双端（paired-end）与单端（single-end）输入均已完成验证，当前输出与 Java CIRI3 保持一致。
- `Scan1` 与 `Scan2` 的性能优化阶段已完成，当前版本整体性能已达到并超过 Java 基线，可作为正式功能版本持续使用。
- 默认后处理已完成 `<prefix>.segments` 和 major isoform 阶段：主流程写完 `<prefix>.out` 后继续生成 read-level circRNA segment chain，并从稳定的 `<prefix>.segments` 重新解析生成 major isoform GTF/FASTA。

## 对齐基线
- circRNA-level：100%
- read-level：100%
- read-assignment-level：100%
- FSJ-level：100%

## 最新验证结论（2026-05）
- `tests/chr1` BAM / SAM 基线继续保持 circ/read/read-assignment/FSJ 四层 100% 对齐。
- `<prefix>.bsj1` 和 `<prefix>.bsj2` 已升级为 mate-level BSJ 中间文件，包含 `mate_label` 与 `priority`；`Summary` 只消费 `priority=1` 行，最终 `<prefix>.bsj` 保留 `priority=0` 行供后续内部 splice site 识别使用。
- 当前运行路径不再生成额外 `<prefix>.scan1.tmp`、`<prefix>.scan2.tmp` 或 `<prefix>.bsj.raw.tmp`。
- 当前正式用户输出为 `<prefix>.out`、`<prefix>.bsj`、`<prefix>.segments`、`<prefix>.isoforms.gtf`、`<prefix>.isoforms.fa`、`<prefix>.bedpe` 和 `<prefix>.segments.bam/.bai`；`.bsj1/.bsj2/.segments1/.segments2/.segments.non_bsj` 以及 `.part_XXXX.tmp` shard 文件均为内部临时文件，默认成功运行后删除，`--debug` 时保留。
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
- chr1 simulator truth 上当前 read-level segments / major isoform 评估口径为：`tmp/chain_final` 完整流程输出 `10813` 条高可信 FASTA isoform；在有 truth circRNA 的 `13643` 个预测 circRNA 上，major isoform exon-chain 严格匹配率为 `83.89%`，1 bp 容忍为 `84.06%`，2 bp 容忍为 `84.23%`。这些指标用于评估 sequence candidate 可信度，不替代 `.out/.bsj/.segments` 的 parity 验收。

统一校验命令：

```bash
python scripts/ciri_result_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI.result \
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
- segments finalize 采用 shard-local retained streams 与两次流式 pass 收集/校正 junction support，避免把所有 read-level rows 全量常驻内存。
- 所有长时程阶段必须有进度条或清晰 INFO 阶段提示；进度条结束后若还有 finalize / merge / collect / correction 工作，必须立即输出新的用户可见提示。

## 调试口径
- 先跑 `scripts/ciri_result_diff.py`，再决定是查 `Scan1`、`Scan2` 还是 `Summary`。
- 追踪单条或少量 read 时，统一使用：
  - `CIRI_TRACE_READS`
  - 必要时加 `CIRI_TRACE_ALL_CANDS=1`
  - 必要时加 `CIRI_TRACE_HG2=1`
- 需要性能热点分布时，统一使用：
  - CLI `--perf`，用于 Scan1/Scan2 主流程 profiling
  - 兼容旧入口：`CIRI_PROFILE_SCAN1=1` / `CIRI_PROFILE_SCAN2=1`
  - `CIRI_PROFILE_SEGMENTS=1`，用于 segments 内部阶段计时
- 详细命令、日志标签与 SOP 见 `docs/02-parity-debug-playbook.md`。

## 当前性能优化状态（2026-05）
- RNA015434 whole-genome BAM 上，当前总耗时约 `1064s`，目标是在不牺牲内存可扩展性的前提下继续向 `10min` 级别收敛。
- 最近一次真实数据 profiling 显示，Scan2 wall time 约 `526s`，主要热点是 `is_bsj_hg2` validator、mate-level display 路径和未拆细的 `other` 框架成本。
- `non_bsj_sidecar` 在 Scan2 内只占约 `1.5%` shard work；它生成的文件很大，但不是 Scan2 当前主瓶颈。
- 下一步优化应先拆细 display / group_process / other profiling，再优先尝试 display gating 与去重复构造；直接修改 HG2 判定逻辑前必须有更窄的热点证据和 parity 回归。

## 后续开发计划
- 高优先级功能点：在当前 `priority` 协议基础上继续推进真正的 paired-end R1/R2 BSJ 决策。当前默认路径仍采用 CIRI3 read-pair first-hit 语义；后续更合理的方向是先做 mate-level candidate 收集，再做 pair-level 决策，最后仍按 read-pair-level 计数。R1/R2 支持同一 BSJ 时应记录一致支持，支持不同 BSJ 时应显式标记冲突或 ambiguity。
- 当前 read-level segments reconstruction 和 major isoform 输出已进入可用的第一版：以 Summary confirmed BSJ、Scan1/Scan2 sidecar evidence 和 non-BSJ topology sidecar 为输入，默认输出 `<prefix>.segments`、`<prefix>.isoforms.gtf`、可信度更严格的 `<prefix>.isoforms.fa`，以及用于 IGV review 的 `<prefix>.bedpe` / `<prefix>.segments.bam`。Scan1/Scan2 segments sidecar 会携带 short-form `cs` 与原始 BWA `XA:Z` payload；这些信息只用于 post-Summary read-level chain selection，不改变 `.out/.bsj` 主流程判定。
- 当前默认后处理不依赖 RO remap，不新增 circ seed，不改变 `priority=1` 主流程证据和 `.out`；isoform 阶段从已写出的 `<prefix>.segments` 重新解析，方便 `--continue` 调试和后续多样本整合。
- IGV visualization sidecar 当前先提供两个最小 track：`.bedpe` 在 `.out` 写出时展示识别到的 BSJ anchor pair，`.segments.bam/.bai` 在 `.segments` 写出时以合成 alignment 展示 bsj/backward/outward read 的 R1/R2 aligned blocks。BAM 中 internal junction 用 `N` CIGAR，`<bsj>` / `B` marker 或 chain 内坐标回跳会拆成同 read 的多条 alignment；当 `r*_cs` 与拆分后的 part span 一致时同步重建 `SEQ`、part-level `cs:Z` 和含 `I/D/N/M` 的 BAM CIGAR，否则回退为 `N` 序列。`YC` 与 `RG` tags 供 IGV 按固定 RGB 或 read group 着色。后续再补 internal junction arcs、confidence tier tracks 和 per-circ mini-BAM/read-list。
- CIRI-AS / CIRI-full 不再作为当前开发路线；后续不以完整复刻其历史输出为目标。相关文档只保留为历史参考、算法术语解释和风险清单，代码中不再使用的 CIRI-AS/cirexon/full-length path 旧实现应删除。
- CIRI-AS sidecar 的 splice signal 阶段已明确采用 annotation-aware motif / strand / offset tie-break：当 annotation exon boundary 支持某个解释时，优先于 Perl v1.2 的 `AC/CT elsif AG/GT` 与 hash-order offset 行为；这是有意偏离 Perl 输出、追求稳定和边界可解释性的设计选择。
- CIRI-AS sidecar 的 splice signal 阶段、annotation-aware motif / strand / offset tie-break 与 read-group 识别逻辑继续保留，用于支撑 `<prefix>.segments` 的 best-chain 选择和 segment 边界解释。
- 当前 full-length 目标已经收敛为单样本 major isoform 输出；下一阶段补齐多 isoform usage 计算和 multi-sample integration。usage 阶段应复用 `<prefix>.segments` parser / circ-local graph，输出 per-sample support、isoform usage、major isoform switching 和跨样本 structure_hash 合并结果。
- FASTA 输出应继续保持“高可信 sequence candidate”定位：`mature` 一律输出；`estimate` 只有在 reason 与 `segment_coverage_pct` 足够可信时输出。GTF 保留未进入 FASTA 的 unresolved/partial estimate 作为审计记录。
- 独立 Rust 模拟器已作为开发 fixture 基础接入，通过 `ciri-simulator` 生成 paired FASTQ、linear-only `<prefix>.annotation.gtf`、circ/isoform-level `<prefix>.isoforms.tsv` 和 read-pair-level `<prefix>.reads.tsv`，用于后续 internal structure 和 full-length path 验证。
- 新增模块仍应遵守小步验证原则：segments 和 major isoform 输出保持稳定，主流程 parity 不回退；multi-isoform usage 和 multi-sample 先作为后处理/sidecar 推进。

## 开发规范
- 后续代码改动应继续遵守已沉淀的 hg38 排障结论、Java parity 约束与关键注释说明，避免在重构或优化中重新引入回归。
- 后续若再进行性能优化或结构调整，仍必须执行完整 parity 校验，出现差异即先回到行为对齐。
- 大规模数据路径的优化必须先确认峰值 RSS 不随 read/row 全量线性增长，再讨论 CPU 时间；不允许为了速度回退到全量 read-level 结构常驻内存。
- 临时文件命名、清理和 debug 保留策略必须与 `docs/00-index.md`、`docs/03-performance-optimization.md` 保持一致。
- 可复用开发脚本统一放在 `scripts/`，一次性探针和分析输出放在 `tmp/`；迁入 `scripts/` 的脚本必须提供通用 CLI 参数，不依赖硬编码 `tmp/...` 路径。
- isoform usage 或 multi-sample integration 实现前应先阅读 `docs/07-full-length-reconstruction.md`；新增 graph、usage、multi-sample merge、FASTA/GTF 输出等模块必须补齐 Rust 文档注释，说明职责边界、证据层级和不影响主流程 parity 的原因。
- CIRI-AS / CIRI-full 相关文档当前只作为历史参考，不作为强制兼容契约；仍被当前 segments 逻辑复用的 annotation-aware offset/strand 判断，必须在代码注释和 `docs/CIRI-AS.md` 同步维护。

---
最后更新：2026-05-20
