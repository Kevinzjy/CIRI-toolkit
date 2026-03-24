# CIRI-toolkit 开发状态

本文档用于描述 Rust 版本 CIRI-toolkit 的当前实现进度、对齐基线、验证口径与下一阶段方向。

## 项目概述
CIRI-toolkit 是 CIRI3（Java）的 Rust 复现版本，目标是在保持判定逻辑一致的前提下提升性能与可扩展性。

## 相关文档
- 调试手册：`docs/parity-debug-playbook.md`
- 性能路线图：`docs/performance-roadmap.md`
- 文档导航：`docs/index.md`

## 当前状态（Functional v1）
- 核心流程 `Scan1 -> Scan2 -> Summary` 已完成并稳定可用。
- 已实现 SAM/BAM 自动识别与对应执行路径。
- `stringency` 过滤、FSJ 统计、注释输出等关键能力均可运行。
- `Scan1` 与 `Scan2` 的核心 BSJ 识别链路已完成性能收敛，当前优化重点转向文档固化、可维护性与剩余边缘 parity case 归档。

## 对齐基线
- circRNA-level：100%
- read-level：100%
- read-assignment-level：100%

## 最新验证结论（2026-03）
- `tests/chr1` 基线继续保持三层 100% 对齐。
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
- 详细命令、日志标签与 SOP 见 `docs/parity-debug-playbook.md`。

## 下一阶段（Post-Parity）
- 在“输出零漂移”前提下继续优化 Scan2。
- 把当前 hg38 排障经验继续沉淀为文档与注释，降低后续回归成本。
- 每次优化必须执行完整 parity 校验，出现差异即回退定位。

---
最后更新：2026-03-24
