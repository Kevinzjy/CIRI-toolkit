# CIRI-toolkit 开发状态

本文档用于描述 Rust 版本 CIRI-toolkit 的当前实现进度、对齐基线、验证口径与长期维护约束。

## 项目概述
CIRI-toolkit 是 CIRI3（Java）的 Rust 复现版本，目标是在保持判定逻辑一致的前提下提升性能与可扩展性。

## 相关文档
- 调试手册：`docs/02-parity-debug-playbook.md`
- 性能总结：`docs/03-performance-optimization.md`
- RO feature 计划：`docs/05-ro-feature-plan.md`
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

## 最新验证结论（2026-03）
- `tests/chr1` 基线继续保持三层 100% 对齐。
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
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 架构要点（已落地）
- Scan2 使用 Java 语义一致的去重索引与桶遍历顺序。
- Scan2 的 FSJ 统计现已对齐 Java `GetFSJClass.getFSJ(...)` 的 bucket-gating 与 BAM mate-run 语义。
- Summary 在关键路径中对齐 Java 容器迭代行为，避免首命中差异。
- 大文件路径采用 mmap/分片与批量写出策略，兼顾吞吐与内存占用。
- SAM/BAM 两条输入路径最终落到相同的判定语义上，避免格式差异带来的输出漂移。
- `Scan2` 在进入 `Summary` 前会主动释放候选索引工作集，以控制 whole-genome 场景下的 RSS。

## 调试口径
- 先跑 `tests/analyze_diff.py`，再决定是查 `Scan1`、`Scan2` 还是 `Summary`。
- 追踪单条或少量 read 时，统一使用：
  - `CIRI_TRACE_READS`
  - 必要时加 `CIRI_TRACE_ALL_CANDS=1`
  - 必要时加 `CIRI_TRACE_HG2=1`
- 需要性能热点分布时，统一使用：
  - `CIRI_PROFILE_SCAN1=1`
  - `CIRI_PROFILE_SCAN2=1`
- 详细命令、日志标签与 SOP 见 `docs/02-parity-debug-playbook.md`。

## 后续开发计划
- 在当前 CIRI3 核心流程稳定的基础上，后续规划支持 circRNA full-length structure / isoform usage 相关能力。
- CIRI-AS / CIRI-full 尚未经过本项目同等级别的严格验证和性能优化，后续不以完整复刻其所有输出为目标；它们主要作为算法思路来源和风险清单。
- 当前优先实现 RO feature 的第一阶段：在 `Scan1` read group 层识别 paired-end reads 的 RO 序列，输出 `<prefix>.ro.fq` 和 `<prefix>.ro.tsv`，并标注 `5p_ro / 3p_ro / bidirectional_ro / full_length_candidate` 等 sequence-level 类型。
- 第一阶段 RO sidecar 输出不得改变现有 `.bsj1/.bsj/.out`；只有在 RO remap、RO Scan1/Scan2、origin read 去重和回归验证完成后，才考虑把 RO evidence 合并进主 BSJ 判定。
- 后续 full-length isoform reconstruction 的目标输出是同一 BSJ 下的不同 isoform 结构及样本 usage ratio，而不是 CIRI-AS/CIRI-full 的历史中间文件。
- 新增模块仍应遵守小步验证原则：默认关闭、sidecar 输出优先、主流程 parity 不回退。

## 开发规范
- 后续代码改动应继续遵守已沉淀的 hg38 排障结论、Java parity 约束与关键注释说明，避免在重构或优化中重新引入回归。
- 后续若再进行性能优化或结构调整，仍必须执行完整 parity 校验，出现差异即先回到行为对齐。
- RO feature 实现前应先阅读 `docs/05-ro-feature-plan.md`；新增 `src/ro.rs`、Scan1 RO 接入点、原始方向 `seq/qual` 恢复、RO writer 与 part 文件合并逻辑，必须补齐 Rust 文档注释，说明职责边界、输出稳定性和不影响主流程 parity 的原因。
- CIRI-AS / CIRI-full 相关实现当前不以完整 parity 为目标，任何取舍必须写入对应设计文档，避免后续误把历史输出当作强制兼容契约。

---
最后更新：2026-04-18
