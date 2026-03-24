# AGENTS.md

## 角色与目标
你是 CIRI-toolkit 的开发代理（Agent）。你的目标是实现 CIRI3 的高性能 Rust 版本，并确保行为与 Java 版本一致。

## 核心原则
- **Java 逻辑即规范**：`vendor/CIRI3` 是唯一行为标准，Rust 必须复现相同决策与输出。
- **一致性优先于“聪明修复”**：已知 Java 行为时，不允许引入补偿性启发式逻辑。
- **先对齐再优化**：仅在行为一致后再做性能优化。
- **测试覆盖关键逻辑**：BSJ 识别、CIGAR 分类、Scan2 救援（rescue）、Summary 合并与严格度（stringency）过滤都应有验证。
- **小步快跑、每步验证**：每次改动后都要立刻做差异比较。

## 技术架构约束
- **SAM/BAM**：使用 `noodles-sam` 和 `noodles-bam`。
- **FASTA/GTF 解析**：使用 `needletail` 和 `noodles-gtf`。
- **两遍扫描流程**：遵循 CIRI3 两遍扫描架构（Scan1 -> Scan2 -> Summary）。
- **过滤条件**：严格按 `Summary.java` 的严格度（stringency）逻辑实现。

## 代码注释规范
- **Rust 文档化注释是强约束**：`src/` 下新增或修改的模块、结构体、函数都应补齐规范的 Rust 文档注释。
- **函数注释必须回答两个问题**：这个函数“做什么”，以及“为什么这样设计”；如果实现是为了保持 Java parity、稳定输出顺序、热路径性能、内存复用或兼容某类输入，必须在注释中明确写出。
- **模块注释要说明职责边界**：使用 `//!` 说明模块在 Scan1 / Scan2 / Summary 总流程中的位置，以及与其他模块的关系。
- **重要分支必须有行间/行内注释**：尤其是以下情况不能省略：
  - Java 逻辑照搬、看起来不够“Rust 风格”的代码
  - 经过验证保留的性能优化
  - 容易误改的边界条件、字符串协议、方向判定、排序/遍历顺序
  - 为了减少分配、复用 buffer、保持稳定输出而做的实现细节
- **不要写低价值注释**：禁止只复述字面语义的注释，例如“给变量赋值”“调用函数处理数据”；注释应优先保留设计动机、约束条件、历史背景和性能原因。
- **修改代码时同步更新注释**：如果函数行为、约束或优化策略发生变化，必须同时更新对应文档注释和关键行内注释，不允许代码已变但注释滞后。
- **评审标准**：后续任何代码改动，都应做到“只看函数签名和注释，就能理解职责、关键约束和不能随意重构的原因”。

## 标准工作流
1. **建立参考基线**：必要时先运行 Java CIRI3 生成基线结果。
2. **模块化实现**：先独立实现并验证基础模块（解析、工具函数等）。
3. **按执行顺序对齐**：先 Scan1，再 Scan2，最后 Summary。
4. **阶段性核对**：每个阶段都和 Java 中间/最终结果比较。

## 对齐关键规则（必须遵守）
- **Scan2 索引构建**
  - 必须基于唯一 circ 位点（Java `chrCircSiteMap` 语义），不能直接用原始 BSJ1 行。
  - 候选数据载荷（payload）字段布局必须保持：
    `[site1, site2, strand, signal1, signal2, sum_q]`。
- **Scan2 候选遍历**
  - 必须保持 Java 的桶门控与方向顺序：`num1` 逆序、`num2` 正序。
  - 正常模式保持“首个有效非 `2` tag 即返回”。
- **Scan2 配对序列方向**
  - 必须按“当前比对链方向（alignment strand）”逐条计算，不能按片段（segment）固定方向。
- **默认参数**
  - 命令行（CLI）默认值应与 CIRI3 默认行为一致，除非有明确证据需要调整。

## 调试与验证流程（SOP）
- 使用 `tests/analyze_diff.py` 做三层比较：
  - circRNA 层面
  - 读段 ID（read ID）层面
  - 读段归属（read-assignment）层面
- 可选定向追踪：
  - `CIRI_TRACE_READS`：按读段 ID（read ID）定向追踪
  - `CIRI_TRACE_ALL_CANDS=1`：仅用于调试时查看被追踪读段（traced read）的 Scan2 全候选
- 调试结束后必须关闭追踪（trace）环境变量，再跑最终验证。
- 详细操作手册见：
  - `docs/parity-debug-playbook.md`
  - `docs/index.md`

## 临时文件规范
- 所有临时脚本输出、探针产物、一次性依赖下载统一放在 `tmp/` 目录下。
- 不要在 `tests/`、`src/`、仓库根目录散落临时文件或编译中间产物。
- `tmp/` 为本地工作区目录，不纳入版本控制。

## 验证数据
- **参考基因组**：`tests/chr1.fa`
- **注释文件**：`tests/chr1.gtf`
- **输入数据**：`tests/test.sam`（或对应 BAM）
- **参考输出**：`tests/ref.txt`（以及 `tests/chr1/CIRI3_result.txt`）

## 当前基线状态
当前 `tests/chr1` 基线对齐状态：
- circRNA 层面：100%
- 读段 ID（read ID）层面：100%
- 读段归属（read-assignment）层面：100%
- FSJ 计数层面：100%

## 当前阶段定义（2026-03）
- **阶段结论**：功能对齐已完成，当前版本可作为首个正式功能版本。
- **全量 hg38 口径**：使用 `/data/public/database/gencode/hg38/gencode.v44.annotation.gtf` 与 Java 基线一致地复核 whole-genome parity。
- **版本目标**：在不改变任何判定结果的前提下，继续推进性能优化（优先 Scan1/Scan2 I/O 路径）。
- **变更红线**：
  - 任何优化提交都必须通过 `tests/analyze_diff.py` 三层零差异检查。
  - 若出现差异，先回到行为对齐再谈性能。
  - 不允许以“统计上接近”替代“逐条一致”。

## 交付前检查清单
- 关闭所有 trace 环境变量后执行一次完整流程。
- 对比 Java 与 Rust 结果，确认：
  - `circ_only_java = 0`
  - `circ_only_rust = 0`
  - `read_only_java = 0`
  - `read_only_rust = 0`
  - `read_assignment_only_java = 0`
  - `read_assignment_only_rust = 0`
- 文档与注释同步更新（README / docs / AGENTS），避免“代码已变更但说明未更新”。
