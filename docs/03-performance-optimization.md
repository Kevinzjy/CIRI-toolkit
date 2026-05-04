# 性能优化总结与维护约束

本文档用于汇总当前版本已经完成的性能优化结论、统一后续性能改动的测量规范，并保留 `Scan1` 优化阶段的 profiling 记录。

## 当前结论
- 当前版本的性能优化阶段已完成。
- 在保持 Java parity 的前提下，整体性能已达到并超过 Java 基线。
- 因此，本文件不再作为活跃优化 backlog 使用，而是作为后续维护时的性能约束、测量规范与历史档案。

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
最后更新：2026-03-30
