# CIRI-toolkit 开发状态

本文档用于描述 Rust 版本 CIRI-toolkit 的当前实现进度、对齐基线与下一阶段方向。

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

## 对齐基线
- circRNA-level：100%
- read-level：100%
- read-assignment-level：100%

统一校验命令：

```bash
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 架构要点（已落地）
- Scan2 使用 Java 语义一致的去重索引与桶遍历顺序。
- Summary 在关键路径中对齐 Java 容器迭代行为，避免首命中差异。
- 大文件路径采用 mmap/分片与批量写出策略，兼顾吞吐与内存占用。

## 下一阶段（Post-Parity）
- 在“输出零漂移”前提下优化 Scan1/Scan2 I/O 性能。
- 每次优化必须执行完整 parity 校验，出现差异即回退定位。

---
最后更新：2026-03-22
