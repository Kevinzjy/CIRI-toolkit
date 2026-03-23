# 文档导航

本目录用于维护项目状态、排查流程与性能优化路线，文档命名统一使用小写 kebab-case。

## 推荐阅读顺序
1. `development-status.md`
   - 当前实现状态、功能基线与阶段目标。
2. `parity-debug-playbook.md`
   - Java 与 Rust 不一致时的标准排查手册（SOP）。
3. `performance-roadmap.md`
   - 在功能对齐完成后的性能优化计划与度量规范。
4. `bam-support-archive.md`
   - BAM 支持历史说明（归档文档，仅供追溯）。

## 维护规则
- 状态类信息写入 `development-status.md`。
- 可执行调试命令与排障步骤写入 `parity-debug-playbook.md`。
- 优化 backlog 与性能实验规范写入 `performance-roadmap.md`。
- 避免跨文档重复维护同一检查清单。
- 所有临时文件统一放在 `tmp/`，不要散落在 `tests/` 或项目根目录。

---
最后更新：2026-03-22
