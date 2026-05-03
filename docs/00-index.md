# 文档导航

本目录用于维护项目状态、排查流程与性能优化路线，文档命名统一使用 `NN-name.md` 的编号前缀加小写 kebab-case。

说明：
- `README.md` 面向最终用户，只保留安装、运行、参数和输出说明。
- 开发状态、排障 SOP、性能总结与专项记录统一维护在 `docs/`。
- 能并入主线文档的历史说明不再单独保留，避免目录里出现只负责“跳转”的文档。

## 推荐阅读顺序
1. `01-development-status.md`
   - 当前实现状态、输入格式支持范围、CIRI3 baseline、CIRI-AS-style 内部结构方向、hg38/v44 验证结论与阶段目标。
2. `02-parity-debug-playbook.md`
   - Java 与 Rust 不一致时的标准排查手册（SOP），包含按具体 read trace 的命令模板与 FSJ 对齐检查口径。
3. `03-performance-optimization.md`
   - 性能优化阶段的统一总结、测量规范、维护约束，以及 Scan1 profiling 专项记录。
4. `04-BSJ_scoring.md`
   - BSJ 重评分模型的专题设计笔记，属于探索性方案，不代表当前主流程实现。
5. `05-ro-feature-plan.md`
   - RO feature 和后续 side evidence 方向的历史设计；当前不作为第一优先实现路线。
6. `06-simulation-truth-design.md`
   - 独立 Rust 模拟器、结构化 truth 表、circ/isoform fixture 生成和后续结构验证口径。
7. `07-full-length-reconstruction.md`
   - 基于 confirmed BSJ、second BAM/SAM sweep 和局部 splice graph 的 CIRI-AS-style full-length 结构识别设计。
8. `CIRI-AS.md`
   - CIRI-AS 上游脚本功能拆解；当前作为 circRNA 内部结构识别思路参考，不作为完整 parity 目标，并记录 annotation-aware motif/offset 等有意偏离 Perl 的规则。
9. `CIRI-full.md`
   - CIRI-full 上游 Java 模块拆解；当前主要作为未来 RO/side evidence 风险清单，不复刻全部历史输出。

## 维护规则
- 状态类信息写入 `01-development-status.md`。
- 可执行调试命令与排障步骤写入 `02-parity-debug-playbook.md`。
- 性能约束、测量规范、优化总结与 profiling 记录统一写入 `03-performance-optimization.md`。
- CIRI-AS-style full-length reconstruction 的阶段性设计写入 `07-full-length-reconstruction.md`；RO side evidence 设计保留在 `05-ro-feature-plan.md`；独立模拟器与 truth 输出设计写入 `06-simulation-truth-design.md`。
- CIRI-AS / CIRI-full 相关实现若有意偏离上游历史输出，必须同步写入对应拆解文档，说明偏离原因、验证口径和对主 CIRI3 parity 的隔离方式。
- 避免跨文档重复维护同一检查清单。
- 所有临时文件统一放在 `tmp/`，不要散落在 `tests/` 或项目根目录。

---
最后更新：2026-04-30
