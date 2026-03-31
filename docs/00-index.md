# 文档导航

本目录用于维护项目状态、排查流程与性能优化路线，文档命名统一使用 `NN-name.md` 的编号前缀加小写 kebab-case。

说明：
- `README.md` 面向最终用户，只保留安装、运行、参数和输出说明。
- 开发状态、排障 SOP、性能总结与专项记录统一维护在 `docs/`。
- 能并入主线文档的历史说明不再单独保留，避免目录里出现只负责“跳转”的文档。

## 推荐阅读顺序
1. `01-development-status.md`
   - 当前实现状态、输入格式支持范围、功能基线、hg38/v44 验证结论与阶段目标。
2. `02-parity-debug-playbook.md`
   - Java 与 Rust 不一致时的标准排查手册（SOP），包含按具体 read trace 的命令模板与 FSJ 对齐检查口径。
3. `03-performance-optimization.md`
   - 性能优化阶段的统一总结、测量规范、维护约束，以及 Scan1 profiling 专项记录。
4. `04-BSJ_scoring.md`
   - BSJ 重评分模型的专题设计笔记，属于探索性方案，不代表当前主流程实现。

## 维护规则
- 状态类信息写入 `01-development-status.md`。
- 可执行调试命令与排障步骤写入 `02-parity-debug-playbook.md`。
- 性能约束、测量规范、优化总结与 profiling 记录统一写入 `03-performance-optimization.md`。
- 避免跨文档重复维护同一检查清单。
- 所有临时文件统一放在 `tmp/`，不要散落在 `tests/` 或项目根目录。

---
最后更新：2026-03-30
