# 性能优化总结与维护约束

本文档用于汇总当前版本已经完成的性能优化结论、统一后续性能改动的测量规范，并保留 `Scan1`、`Scan2` 与 segments 后处理阶段的 profiling 记录。

## 当前结论
- CIRI3-compatible `Scan1 -> Scan2 -> Summary` 主流程已经完成 Java parity 对齐。
- 当前活跃优化对象是 whole-genome BAM 上的 `Scan2` 与默认 `<prefix>.segments` 后处理阶段。
- 优化优先级是：先保证大规模数据下峰值 RSS 与临时 I/O 有明确上界，再优化 wall-clock 时间。
- 因此，本文件同时作为性能约束、测量规范、历史档案和当前 Scan2/segments 优化 backlog。

## 已完成的优化范围
- 不改变 Java 对齐结果（行为与输出保持一致）。
- 已完成 Scan1/Scan2 的 I/O 路径与热路径分配开销优化。
- 已完成 BAM 解码、候选构造与关键数据结构上的低风险性能收敛。

## 强制门槛
每次性能相关改动后都要通过零差异校验：
- `circ_only_java = 0`
- `circ_only_rust = 0`
- `read_only_java = 0`
- `read_only_rust = 0`
- `read_assignment_only_java = 0`
- `read_assignment_only_rust = 0`

参考命令：

```bash
python scripts/ciri_result_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 历史优化重点

### Scan1/Scan2 I/O 吞吐
- 降低热循环中的临时分配与字符串构造。
- 强化批量写出策略，减少小写入次数。
- 结合数据规模调优分片大小与 flush 策略。

### BAM 解码与解析路径
- 避免不必要的 CIGAR/SEQ 中间字符串物化。
- 评估虚拟位置更新与进度显示开销。
- 保留诊断路径，同时隔离调试分支的运行时成本。

### 数据结构局部性
- 检查高频 HashMap/HashSet key 组织方式。
- 在安全前提下以可复用缓冲替代重复 `format!`。
- 通过大 BAM 数据评估缓存命中变化。

## 后续维护要求
- 若后续没有明确收益和验证计划，不建议为“可能更快”而重启大规模性能重构。
- 若后续再次进行性能优化，仍应遵守“先 parity、后提速”的原则。
- 任何真实大数据优化都必须先评估峰值 RSS、是否存在全量 read/row 常驻内存、以及临时文件是否可流式合并。
- 可以接受为了降低 OOM 风险牺牲一部分 CPU 时间；只有内存可扩展性确认后，才允许用更多缓存或预计算换速度。
- 固定输入数据、线程数与执行参数后再比较。
- 至少记录以下指标：
  - 总耗时
  - Scan1 耗时
  - Scan2 耗时
  - 峰值 RSS
  - parity 差异结果
- 任一优化若引入输出差异，直接判定为无效。

## 单步迭代模板
1. 记录基线指标。
2. 只做一个最小改动。
3. 重新构建并复测性能。
4. 运行 parity 校验。
5. 按“速度 + 对齐”双指标决定保留或回退。

## profiling 与进度条规范

### 入口

- 真实 BAM/SAM、CIRI3 parity、CIRI-AS/segments 性能观察统一使用 `--release` 构建。
- 推荐用 CLI `--perf` 打开主流程 profiling，输出固定写入 `<prefix>.perf.log`。
- 兼容环境变量入口：
  - `CIRI_PROFILE_SCAN1=1`
  - `CIRI_PROFILE_SCAN2=1`
  - `CIRI_PROFILE_SEGMENTS=1`，用于 segments 细分阶段计时，当前输出到 stderr

### 日志解释

- `wall_ms` 是用户实际等待时间。
- `shard_work_ms` 是所有 worker 线程累计工作时间，通常大于 `wall_ms`。
- 多线程阶段判断优化潜力时，应同时看 `wall_ms` 和 shard work 占比；只看单个线程累计时间容易高估收益。
- 长时程阶段必须有进度条或明确的 INFO 阶段提示。若进度条结束后还有 finalize / merge / load / correction 等重计算步骤，必须立刻输出新的 INFO 行并在可行时加进度条。

### 临时文件与内存边界

- 多进程 / 多线程临时文件统一使用 `<merged-path>.part_XXXX.tmp` 命名。
- 临时文件只能作为 shard-local spill、bounded merge 或 debug artifact 使用，不应成为正式输出协议。
- 默认成功运行后删除内部临时文件；用户显式传入 `--debug` 时保留。
- 正式用户输出当前限定为 `<prefix>.out`、`<prefix>.bsj`、`<prefix>.segments`、`<prefix>.isoforms.gtf`、`<prefix>.isoforms.fa`、`<prefix>.bedpe`、`<prefix>.segments.bam` 和 `<prefix>.segments.bam.bai`；`.bedpe` 跟随 `.out` 写出，`.segments.bam/.bai` 跟随 `.segments` 写出并作为大数据 IGV review 主入口。

## Scan2 专项记录（2026-05）

### RNA015434 whole-genome baseline

测试命令口径：

```bash
target/release/ciri \
  -i ./bam/RNA015434_S1.bam \
  -o ./CIRI-rs.ciri \
  -r /data/public/database/gencode/hg38/_BWAindex/hg38.fa \
  -a /data/public/database/gencode/hg38/gencode.v44.annotation.gtf \
  -t 16 -s 0 --perf
```

真实数据运行中，几个主要阶段的 wall time 约为：

- Scan1：`1min27s`
- Scan2：`8min32s`
- Scan2 non-BSJ topology sidecar loading：`3min49s`
- segments junction support collection：`1min46s`
- segments ambiguous row correction：`1min57s`

最近一次 Scan2 profiling：

```text
[PROFILE_SCAN2] wall_ms=525631.514 shard_work_ms=8125676.793 merge_ms=182.006 records=796992560 groups=314123425 candidate_checks=134202627 candidate_hits=265427
[PROFILE_SCAN2] shard_breakdown_ms group_process=3101549.327 (38.2%) validator=2778147.879 (89.6% of group) display=2706837.458 (33.3%) segment_sidecar=764943.619 (9.4%) non_bsj_sidecar=125373.266 (1.5%) write=7637.016 (0.1%) other=1419336.107 (17.5%)
[PROFILE_SCAN2_SIDECAR] bsj_rows=2191949 bsj_bytes=221655182 non_bsj_rows=82828148 non_bsj_bytes=19668758102
[PROFILE_SCAN2_HG2] calls=269693240 sw_calls=70537486 linear11_calls=2556514 linear12_calls=2486426 circ2_calls=155528 circ3_calls=560148
[PROFILE_SCAN2_HG2] total_ms=5149988.812 sw_ms=2958217.035 (57.4%) linear11_ms=1189990.592 (23.1%) linear12_ms=620147.359 (12.0%) circ2_ms=5289.751 (0.1%) circ3_ms=276921.395 (5.4%) other_ms=99422.680 (1.9%)
```

### Scan2 当前结论

- `non_bsj_sidecar` 在 Scan2 内只占 `1.5%` shard work。它生成的文件很大，但不是 Scan2 的主要 wall-clock 瓶颈。
- `group_process` 占 `38.2%` shard work，其中 `validator` 占 group 内 `89.6%`，说明 `is_bsj_hg2` 仍是主流程 CPU 热点。
- `display` 占 `33.3%` shard work，已经和主 validator 同量级。它属于 mate-level `.bsj` display / segments 支撑路径，优化风险低于直接改 Java parity 判定逻辑。
- `other` 占 `17.5%`，需要继续拆分 BAM decode、read group assembly、alignment record 构造和未覆盖的 side work。
- `segment_sidecar` 占 `9.4%`，更可能是证据行构造与字符串格式化成本，而不是纯 I/O。

### Scan2 优化优先级

1. 先拆细 `display` profiling：统计 Scan1 claim 后跳过的 read/mate、实际进入 display validator 的 read/mate、display HG2 调用数和命中数。
2. 再拆细 `group_process`：区分 candidate enumeration、`circ_c` 构造、HG2 调用和 FSJ 更新。
3. 再拆细 `other`：确认未覆盖时间来自 BAM 解码、group assembly、字符串分配还是外围框架。
4. 优先尝试低风险 display gating / 去重复构造；只有定位清楚后再进入 `is_bsj_hg2` 内部逻辑优化。
5. 若优化 `is_bsj_hg2`，必须保持 Java branch order、early return 和 tag 语义，不允许用近似判断替代 Smith-Waterman / linear competition 结果。

## Segments 后处理专项记录（2026-05）

### 当前结构

- Scan1/Scan2 在主流程中写 `<prefix>.segments1`、`<prefix>.segments2` 和 non-BSJ topology sidecar。
- segments finalize 使用 shard-local retained streams，避免把所有 read-level rows、junction support 和 ambiguous rows 全量放入内存。
- `Collecting junction support...` 与 `Correcting ambiguous segment rows...` 是两个独立 pass，必须分别打印 INFO 和进度条。

### 当前结论

- streaming finalize 后，真实 RNA015434 的峰值 RSS 已从不可接受的全量常驻下降到可运行范围。
- 当前内存策略可以继续保留；后续速度优化不应回退到全量 HashMap / Vec 常驻所有 segment row 的设计。
- 后续 segments 速度优化应优先减少重复解析 retained shard、减少最终 row materialization 前的字符串拆装，而不是扩大内存缓存。

## Scan1 专项记录

### 历史基线
- Java 版本：约 14s
  - First scan：约 5s
  - Second scan：约 9s
- Rust 版本：约 47s
  - First scan：约 22s
  - Second scan：约 22s
- 当前性能验证口径：`cargo run --release -- ... -t 16`
- 当前一致性验证口径：`python scripts/ciri_result_diff.py tests/chr1/CIRI3_result.txt <rust_result>`

### 已完成的 Scan1 优化

#### 1. BAM 并行分片处理已恢复
- `Scan1` 的 BAM 路径已从单线程顺序读取切回并行分片执行。
- shard 输出会在最后 merge 成统一 First Scan 中间文件 `<prefix>.bsj1`。
  Second Scan 的 shard 输出会 merge 成 `<prefix>.bsj2`。
- 最终 `<prefix>.bsj` 由 `<prefix>.bsj1 + <prefix>.bsj2` 在 Summary 之后排序生成；当前不再额外重扫输入，也不再生成 `<prefix>.scan1.tmp`、`<prefix>.scan2.tmp` 或 `<prefix>.bsj.raw.tmp`。
- 这一步已经通过 `tests/chr1` 的 circ/read/read-assignment/FSJ 四层零差异验证。

#### 2. BAM 热循环的低风险去分配优化
位置：[src/scan1.rs](/data/zhangjy/workspace/CIRI-toolkit/src/scan1.rs)

已落地的优化包括：
- 复用 `cigar_buf` / `seq_buf`，减少每条 BAM record 的临时 `String` 分配。
- `stand_map.entry(...).or_insert_with(...)`，避免 eager clone。
- shard merge 改成复用 `read_line` buffer，不再每行都新建 `String`。
- `process_group_view` 中预计算 `misd(cigar, seq_len)`，避免在双重循环里重复计算。
- `seq_oriented` 在方向一致时直接借用 `read_seq`，不再无条件整段复制。
- `line_arr` 从 `Vec<String>` 改成定长数组切片，减少热路径堆分配。
- 返回结果行的嵌套 `format!` 已扁平化。

#### 3. `is_bsj_hg1` 内部的低风险字符串优化
位置：[src/is_bsj_hg2.rs](/data/zhangjy/workspace/CIRI-toolkit/src/is_bsj_hg2.rs)

已落地的优化包括：
- `chr\tpos` key 构造改为复用 buffer，而不是反复 `format!`。
- 去掉 `split(...).collect::<Vec<_>>()` 这类纯中间分配，改为迭代器逐段取值。
- 多处 `format!("{}{}", ...)` 改成预分配 `String` + `push_str`。

### Scan1 profiling 结果
为了确认 `Scan1` 的真实瓶颈，已经在 [src/scan1.rs](/data/zhangjy/workspace/CIRI-toolkit/src/scan1.rs) 增加了仅在 `CIRI_PROFILE_SCAN1=1` 时启用的分段计时。

使用命令：

```bash
CIRI_PROFILE_SCAN1=1 cargo run --release --   -i tests/chr1/test.bam   -o tmp/scan1_profile_release   -r tests/chr1/chr1.fa   -a tests/chr1/chr1.gtf   -t 16
```

profiling 输出：

```text
[PROFILE_SCAN1] wall_ms=20888.792 shard_work_ms=234356.861 merge_ms=6.756 records=12608979 groups=1834141 hg1_calls=89330 hg1_hits=69480
[PROFILE_SCAN1] shard_breakdown_ms group_process=223592.201 (95.4%) bsj_judge=223001.374 (99.7% of group) write=2.413 (0.0%) other=10762.247 (4.6%)
```

解释：
- `wall_ms` 是 `Scan1` 实际墙钟时间，约 20.9s。
- `shard_work_ms` 是所有线程累计工作时间，所以会大于墙钟时间。
- `group_process` 已经占到 shard 工作时间的 `95.4%`。
- `bsj_judge` 又占到 `group_process` 的 `99.7%`。
- `merge` 和写文件几乎可以忽略。

### Scan1 阶段结论
当前 `Scan1` 的主耗时已经明确：
- 不是 shard merge
- 不是结果写出
- 不是外围 BAM 路径框架
- 核心瓶颈就是 `process_group_view -> is_bsj_hg1` 这条 BSJ 判断链

这些 profiling 结果已经完成了当前阶段的定位任务。随着后续优化落地，项目整体性能已达到并超过 Java 基线，因此 `Scan1` 不再作为活跃性能攻关项单独推进。

### Scan1 后续维护要求
1. 若未来再次修改 `Scan1` 热路径，仍应先固定基线，再做最小改动。
2. 性能结论统一基于 `--release`。
3. 任何结构调整或性能优化都必须重新执行 parity 校验。
4. 必要时可重新开启 `CIRI_PROFILE_SCAN1=1` 复用本文件中的 profiling 口径。

---
最后更新：2026-05-13
