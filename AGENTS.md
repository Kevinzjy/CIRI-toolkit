# AGENTS.md

## 角色与目标
你是 CIRI-toolkit 的开发代理（Agent）。你的目标是实现 CIRI3 的高性能 Rust 版本，并确保行为与 Java 版本一致。

## 核心原则
- **Java 逻辑即规范**：`vendor/CIRI3` 是唯一行为标准，Rust 必须复现相同决策与输出。
- **一致性优先于“聪明修复”**：已知 Java 行为时，不允许引入补偿性启发式逻辑。
- **先对齐再优化**：仅在行为一致后再做性能优化。
- **大规模可扩展性优先于速度**：面向真实全量 BAM/SAM、超大 sidecar 或 whole-genome hg38 数据时，优化目标首先是保证内存有明确上界、可流式/分片处理、不会随 read/row/alignment 全量线性驻留；在 scalability 和内存安全成立后，才进一步优化 wall-clock 时间和 CPU 吞吐。
- **测试覆盖关键逻辑**：BSJ 识别、CIGAR 分类、Scan2 救援（rescue）、Summary 合并与严格度（stringency）过滤都应有验证。
- **小步快跑、每步验证**：每次改动后都要立刻做差异比较。
- **扩展模块不反向污染主流程**：segments、isoform usage、multi-sample 等扩展能力默认作为 sidecar 或后处理推进，未完成验证前不得改变 `Scan1 -> Scan2 -> Summary` 的既有输出。

## 技术架构约束
- **SAM/BAM**：使用 `noodles-sam` 和 `noodles-bam`。
- **FASTA/GTF 解析**：使用 `needletail` 和 `noodles-gtf`。
- **两遍扫描流程**：遵循 CIRI3 两遍扫描架构（Scan1 -> Scan2 -> Summary）。
- **过滤条件**：严格按 `Summary.java` 的严格度（stringency）逻辑实现。

## 扩展模块策略（segments / isoforms / multi-sample）
- **CIRI3 parity 与扩展功能分层**：`vendor/CIRI3` 仍是核心 BSJ 检测行为标准；`vendor/CIRI-AS` 与 `vendor/CIRI-full` 只作为历史算法参考和风险清单，不再作为当前开发路线或完整 parity 目标。
- **当前已落地能力**：默认后处理已实现 `<prefix>.segments`、`<prefix>.isoforms.gtf`、高可信 `<prefix>.isoforms.fa`、`<prefix>.bedpe` 和 `<prefix>.segments.bam/.bai`；`<prefix>.segments` 包含 chain-level CIGAR 与 short-form `cs`，`.segments.bam` 从这些字段重建 synthetic sequence / `cs:Z` 用于 IGV review；isoform 阶段必须从已写出的 `<prefix>.segments` 重新解析，保证断点调试和后续多样本整合共用同一输入边界。
- **当前 multi-sample 状态**：two-pass shared-catalog workflow 已有首版闭环，包括 `--1st-pass`、`ciri-merge`、`--2nd-pass --circ` 和 `ciri-assemble`。`ciri-assemble` 当前使用串行 isoform reconstruction / usage assignment；由于现阶段速度足够快，暂不实现并行化，后续等样本数增加并拿到大 multi-sample benchmark 后再基于 profiling 决定是否并行加速。
- **下一阶段目标**：继续围绕同一 circRNA 的多个候选 isoform 完善 usage 计算、置信度分层和 multi-sample integration；重点是比较 major isoform switching、结构稳定性和样本间 usage 变化，而不是沿 CIRI-AS / CIRI-full / RO remap 路线继续复刻历史输出。
- **证据边界**：BSJ/backward/outward segments 是当前 isoform 图的主要证据层。后续新增 evidence 必须先进入审计字段或 sidecar，不能反向改变 `.out/.bsj/.segments` 的既有判定。
- **文档优先**：isoform usage、多样本整合、输出协议或证据分层发生变化时，必须同步更新 `docs/07-full-length-reconstruction.md`、`docs/01-development-status.md` 和必要的用户文档。
- **用户定稿文案优先**：用户已经调整过的 clap help 文案和 README 用法说明视为用户界面契约；除非用户明确要求，不要为了同步实现细节而主动改写。输出协议需要记录时优先更新内部设计文档，确实需要改 README 或 clap 文案时先单独确认。
- **项目本地 skills**：CIRI 专用重复工作流保存在 `.agents/skills/`。涉及 segments / isoform / FASTA 可信度评估时优先读取 `.agents/skills/ciri-segments-isoform-eval/SKILL.md`；涉及 release、README、AGENTS 或开发文档同步时优先读取 `.agents/skills/ciri-release-and-doc-sync/SKILL.md`。

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
5. **长时程任务必须可见**：任何可能处理大文件、全量 BAM/SAM、超大 sidecar 或耗时超过数十秒的阶段，都必须提供进度条或明确的阶段性日志；如果进度条结束后仍有 finalize/merge/load/correction 等重计算步骤，也必须立即输出用户可见提示，避免终端表现为“卡住”。

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
  - 已确认例外：`ciri -s/--stringency` 默认值为 `0`，用于保留候选 circRNA 供 segments、isoform 和多样本规则后续过滤；Java CIRI3 `-S/--strigency` 默认仍记录为 `2`，stringency 判定公式本身必须保持 Java parity。
  - 已确认例外：`ciri --min-span` 默认值为 `50`，用于保留较短 circRNA 候选供后处理过滤；Java CIRI3 `-Min` 默认仍记录为 `140`，显式传参时仍按用户指定值执行。

## 调试与验证流程（SOP）
- 使用 `scripts/ciri_result_diff.py` 做三层比较：
  - circRNA 层面
  - 读段 ID（read ID）层面
  - 读段归属（read-assignment）层面
- 除非正在做必须依赖 debug build 的单元级调试，否则涉及真实 SAM/BAM fixture、CIRI3/CIRI-AS parity 或性能观察的流程应使用 `--release` 模式运行，以避免 debug build 的额外耗时干扰迭代。
- 可选定向追踪：
  - CLI `--trace <READS>`：按读段 ID（read ID）定向追踪，写入 `<prefix>.trace.log`
  - `CIRI_TRACE_READS`：兼容旧脚本的按 read ID 定向追踪入口
  - `CIRI_TRACE_ALL_CANDS=1`：仅用于调试时查看被追踪读段（traced read）的 Scan2 全候选
- 性能 profiling：
  - CLI `--perf`：主流程 profiling 写入 `<prefix>.perf.log`
  - 兼容旧入口：`CIRI_PROFILE_SCAN1=1` / `CIRI_PROFILE_SCAN2=1`
  - `CIRI_PROFILE_SEGMENTS=1`：segments 内部阶段计时
- `--debug` 只表示保留内部临时文件，不等同于 trace 或 profiling。
- 调试结束后必须关闭追踪（trace）环境变量，再跑最终验证。
- 常用开发脚本：
  - `scripts/ciri_result_diff.py`：CIRI3 vs Rust `.out` 四层差异检查。
  - `scripts/ciri_segments_eval.py`：segments 与 simulator truth 的 read-level 评估。
  - `scripts/ciri_read_subset.py`：按 read list 生成调试子集。
  - `scripts/ciri_extract_interval_read_names.py`：按 genomic interval 从 BAM/SAM 收集 read-name list，便于生成 focused subset。
- 详细操作手册见：
  - `docs/02-parity-debug-playbook.md`
  - `docs/00-index.md`

## 临时文件规范
- 当前正式用户输出限定为 `<prefix>.out`、`<prefix>.bsj`、`<prefix>.segments`、`<prefix>.isoforms.gtf`、`<prefix>.isoforms.fa`、`<prefix>.bedpe`、`<prefix>.segments.bam` 和 `<prefix>.segments.bam.bai`；`.bedpe` 在 `.out` 写出时同步生成，`.segments.bam/.bai` 在 `.segments` 写出时同步生成并作为大数据 IGV review 主入口；`.bsj1/.bsj2/.segments1/.segments2/.segments.non_bsj` 以及 `.part_XXXX.tmp` shard 文件都是内部临时文件。
- 所有多进程/多线程 shard 临时文件统一使用合并后文件名加 `.part_XXXX.tmp` 后缀，例如 `<prefix>.segments.non_bsj.part_0001.tmp`。
- 默认成功运行后删除内部临时文件；只有 CLI `--debug` 才保留。
- 所有临时脚本输出、探针产物、一次性依赖下载统一放在 `tmp/` 目录下；需要长期复用的脚本整理到 `scripts/`，并提供通用 CLI 参数，不保留硬编码 `tmp/...` 输入。
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
- **性能优先级**：所有真实大数据路径的优化应先控制峰值 RSS 和临时 I/O 上界，优先采用 Scan1/Scan2 式 shard-local spill、bounded merge、按需回读和流式写出；只有在确认不会 OOM 或长时间持有全量 read-level 结构后，才允许用更多内存换取计算速度。
- **变更红线**：
  - 任何优化提交都必须通过 `scripts/ciri_result_diff.py` 三层零差异检查。
  - 若出现差异，先回到行为对齐再谈性能。
  - 不允许以“统计上接近”替代“逐条一致”。

## 当前扩展阶段定义（2026-05）
- **阶段结论**：read-level segments 与单样本 major isoform 全长识别已实现，当前默认输出可作为 first usable version。严格 segment 边界仍存在可解释误差，尤其是 `<10bp` terminal fragment 被 mapper 吸收到相邻 block 或残留 clip 低于 rescue 阈值的场景；这类弱证据后续若要修复，必须作为 sub-10bp rescue 独立评估，不能放宽普通 junction correction 规则。
- **实现范围**：主流程默认输出 `<prefix>.segments`、`<prefix>.isoforms.gtf`、高可信 `<prefix>.isoforms.fa`、`<prefix>.bedpe` 和 `<prefix>.segments.bam/.bai`；`--continue` 只从已完成的 `<prefix>.segments` 断点重建 isoforms。
- **下一阶段目标**：实现多 isoform usage 计算与 multi-sample 整合，包括候选 isoform search space、per-sample support/usage、major isoform switching 和跨样本结构合并。
- **暂不推进内容**：不再把 CIRI-AS / CIRI-full / RO remap 作为主开发路线；相关旧设计只保留为参考，废弃或未使用代码应优先删除，除非仍直接支撑 segments 或 major isoform 输出。

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
