# Scan1 优化状态（2026-03-23）

本文档记录当前 `Scan1` 的性能现状、已经完成的低风险优化、release profiling 结果，以及下一阶段的优化重点。

## 当前基线
- Java 版本：约 14s
  - First scan：约 5s
  - Second scan：约 9s
- Rust 版本：约 47s
  - First scan：约 22s
  - Second scan：约 22s
- 当前性能验证口径：`cargo run --release -- ... -t 16`
- 当前一致性验证口径：`python tests/analyze_diff.py tests/chr1/CIRI3_result.txt <rust_result>`

## 强制约束
- 所有优化都必须保持与 Java 三层零差异：
  - `circ_only_java = 0`
  - `circ_only_rust = 0`
  - `read_only_java = 0`
  - `read_only_rust = 0`
  - `read_assignment_only_java = 0`
  - `read_assignment_only_rust = 0`
- 性能结论只看 `--release`，不再使用 debug 结果判断是否提速。
- 当前默认性能验证线程数为 `-t 16`。

## 已完成的优化

### 1. BAM 并行分片处理已恢复
- `Scan1` 的 BAM 路径已从单线程顺序读取切回并行分片执行。
- shard 输出会在最后 merge 成统一 `.BSJ1`。
- 这一步已经通过 `tests/chr1` 的三层零差异验证。

### 2. BAM 热循环的低风险去分配优化
位置：[src/scan1.rs](/data/zhangjy/workspace/CIRI-toolkit/src/scan1.rs)

已落地的优化包括：
- 复用 `cigar_buf` / `seq_buf`，减少每条 BAM record 的临时 `String` 分配。
- `stand_map.entry(...).or_insert_with(...)`，避免 eager clone。
- shard merge 改成复用 `read_line` buffer，不再每行都新建 `String`。
- `process_group_view` 中预计算 `misd(cigar, seq_len)`，避免在双重循环里重复计算。
- `seq_oriented` 在方向一致时直接借用 `read_seq`，不再无条件整段复制。
- `line_arr` 从 `Vec<String>` 改成定长数组切片，减少热路径堆分配。
- 返回结果行的嵌套 `format!` 已扁平化。

### 3. `is_bsj_hg1` 内部的低风险字符串优化
位置：[src/is_bsj_hg2.rs](/data/zhangjy/workspace/CIRI-toolkit/src/is_bsj_hg2.rs)

已落地的优化包括：
- `chr\tpos` key 构造改为复用 buffer，而不是反复 `format!`。
- 去掉 `split(...).collect::<Vec<_>>()` 这类纯中间分配，改为迭代器逐段取值。
- 多处 `format!("{}{}", ...)` 改成预分配 `String` + `push_str`。

## 当前 profiling 结果
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

## 结论
当前 `Scan1` 的主耗时已经明确：
- 不是 shard merge
- 不是结果写出
- 不是外围 BAM 路径框架
- 核心瓶颈就是 `process_group_view -> is_bsj_hg1` 这条 BSJ 判断链

也就是说，继续在外围 I/O、merge、progress bar 上做优化，收益会非常有限。后续要把 `Scan1` 明显压到接近 Java 水平，必须继续深入 `is_bsj_hg1`。

## 下一阶段优化重点
按优先级排序：

1. 降低 `is_bsj_hg1` 中重复 `java_substring(...).to_string()` 的次数。
2. 减少 `circ_range_seq`、`linear_range`、`str_new` 等候选级字符串的重复物化。
3. 继续消除只为中间判断服务的短生命周期 `String`。
4. 在不改变 Java 分支顺序的前提下，检查能否把部分字符串判断延迟到真正需要时再构造。

## 当前建议工作方式
1. 只做一类最小优化。
2. `cargo check`。
3. `cargo run --release -- ... -t 16`。
4. `python tests/analyze_diff.py ...` 做三层零差异检查。
5. 必要时重新打开 `CIRI_PROFILE_SCAN1=1` 对比分段耗时变化。

## 当前状态结论
- 功能 parity：稳定保持 100%
- BAM 并行分片：已恢复并验证
- 低风险外围优化：已基本榨干
- 真正下一步：直接优化 `is_bsj_hg1` 热路径

---
最后更新：2026-03-23
