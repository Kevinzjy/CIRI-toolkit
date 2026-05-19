# 文档导航

本目录用于维护项目状态、排查流程与性能优化路线，文档命名统一使用 `NN-name.md` 的编号前缀加小写 kebab-case。

说明：
- `README.md` 面向最终用户，只保留安装、运行、参数和输出说明。
- 开发状态、排障 SOP、性能总结与专项记录统一维护在 `docs/`。
- 能并入主线文档的历史说明不再单独保留，避免目录里出现只负责“跳转”的文档。

## 推荐阅读顺序
1. `01-development-status.md`
   - 当前实现状态、输入格式支持范围、CIRI3 baseline、segments/major isoform 阶段、hg38/v44 验证结论与下一阶段目标。
2. `02-parity-debug-playbook.md`
   - Java 与 Rust 不一致时的标准排查手册（SOP），包含按具体 read trace 的命令模板与 FSJ 对齐检查口径。
3. `03-performance-optimization.md`
   - 性能优化阶段的统一总结、测量规范、维护约束，以及 Scan1 / Scan2 / segments profiling 专项记录。
4. `CIRI3.md`
   - CIRI3 Java 逻辑拆解、Rust parity 经验、`.bsj1/.bsj2/.bsj` mate-level 协议和 `priority` 语义。
5. `04-BSJ_scoring.md`
   - BSJ 重评分模型的专题设计笔记，属于探索性方案，不代表当前主流程实现。
6. `05-ro-feature-plan.md`
   - RO feature 和后续 side evidence 方向的历史设计；当前不作为活跃开发路线。
7. `06-simulation-truth-design.md`
   - 独立 Rust 模拟器、结构化 truth 表、circ/isoform fixture 生成和后续结构验证口径。
8. `07-full-length-reconstruction.md`
   - 默认 `<prefix>.segments` 后处理设计，包含 BSJ/backward/outward read-level chain、内部 junction 校正、major isoform GTF/FASTA sidecar、chr1 FASTA 准确率口径、IGV visualization sidecar 规划和后续 multi-isoform usage / multi-sample integration 入口。
9. `CIRI-AS.md`
   - CIRI-AS 上游脚本功能拆解；当前作为历史参考和 segments 边界解释背景，不作为完整 parity 目标，并记录 annotation-aware motif/offset 等有意偏离 Perl 的规则。
10. `CIRI-full.md`
   - CIRI-full 上游 Java 模块拆解；当前主要作为历史风险清单，不复刻全部历史输出。

## 维护规则
- 状态类信息写入 `01-development-status.md`。
- 可执行调试命令与排障步骤写入 `02-parity-debug-playbook.md`。
- 性能约束、测量规范、优化总结与 profiling 记录统一写入 `03-performance-optimization.md`。
- segments、major isoform、后续 isoform usage 和 multi-sample integration 的阶段性设计写入 `07-full-length-reconstruction.md`；RO side evidence 历史设计保留在 `05-ro-feature-plan.md`；独立模拟器与 truth 输出设计写入 `06-simulation-truth-design.md`。
- CIRI-AS / CIRI-full 相关内容当前只作为历史参考；仍被当前代码复用的局部逻辑若有意偏离上游历史输出，必须同步写入对应拆解文档，说明偏离原因、验证口径和对主 CIRI3 parity 的隔离方式。
- 避免跨文档重复维护同一检查清单。
- 正式用户输出协议当前限定为 `<prefix>.out`、`<prefix>.bsj`、`<prefix>.segments`、`<prefix>.isoforms.gtf`、`<prefix>.isoforms.fa`、`<prefix>.bedpe`、`<prefix>.segments.bam` 和 `<prefix>.segments.bam.bai`；`.bedpe` 在 `.out` 写出时同步生成，`.segments.bam/.bai` 在 `.segments` 写出时同步生成并作为大数据 IGV review 主入口；`.bsj1/.bsj2/.segments1/.segments2/.segments.non_bsj` 及 `.part_XXXX.tmp` 都是内部临时文件。
- 多线程临时文件统一使用 `<merged-path>.part_XXXX.tmp` 命名，默认成功运行后删除；只有 `--debug` 才保留。
- `--continue` 是普通命令的附加执行模式，不改变必填参数；它只从 `<prefix>.segments` 这个合并完成的断点恢复并重建 isoforms，`.out + .bsj`、`.part_XXXX.tmp` 和其他 shard-local 临时文件不作为断点。
- 所有一次性探针脚本、临时分析产物和手工测试输出统一放在 `tmp/`，不要散落在 `tests/` 或项目根目录；如果脚本需要复用，整理到 `scripts/` 并写清 CLI 参数，不保留硬编码 `tmp/...` 输入。

---
最后更新：2026-05-19
