# CIRI-toolkit 对齐排障手册

本手册用于在 **Java CIRI3 vs Rust CIRI-toolkit** 出现差异时，快速、可重复地定位根因并回归验证。

目标：避免临时“拍脑袋”排查，统一团队调试流程与判定口径。

## 1. 调试原则

- 以 `vendor/CIRI3` 为唯一行为规范（spec）。
- 先保证行为一致，再考虑性能优化。
- 每次只改一个小点，改后立即复测。
- 优先做“可证伪”的对比：先数据定位，再代码定位。
- 严禁用“补偿逻辑”掩盖与 Java 的真实分歧。

## 2. 标准输入与基线

- 参考基因组：`tests/chr1/chr1.fa`
- 注释文件：`tests/chr1/chr1.gtf`
- 输入 BAM：`tests/chr1/test.bam`
- Java 结果：`tests/chr1/CIRI3_result.txt`

建议每次调试都固定上述数据，以减少噪音变量。

whole-genome / hg38 复核时，还需要先锁定：

- 参考 FASTA 版本一致
- 注释 GTF 版本一致
- exon 覆盖与 Java 基线一致

最近的 hg38 排查已经证明：注释版本差异本身就足以制造表面上的 parity gap。

## 3. 三层差异检查

使用统一脚本一次性看三层指标：

```bash
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

三层定义：

- circRNA-level：circ 位点集合是否一致
- read-level：read ID 集合是否一致
- read-assignment-level：`(circRNA_ID, read_ID)` 是否一致

排查优先级建议：

1. 先把 circRNA 层面对齐到 100%
2. 再把 read ID 层面对齐到 100%
3. 最后收敛 read-assignment 层面

## 4. 定向追踪开关

### 4.1 追踪指定 read

```bash
CIRI_TRACE_READS="simulate:7037,simulate:7050" \
./target/release/ciri-toolkit \
  -i tests/chr1/test.bam \
  -o tests/chr1/CIRI-rs.trace \
  -r tests/chr1/chr1.fa \
  -a tests/chr1/chr1.gtf \
  -t 16
```

输出到 `stderr`，关键标签：

- `TRACE_SCAN1_CAND`
- `TRACE_SCAN1_HG1`
- `TRACE_SCAN2_CAND`
- `TRACE_SCAN1_SHARD`（BAM 分片边界）

### 4.2 扫描 Scan2 全候选（仅诊断）

```bash
CIRI_TRACE_READS="simulate:398600" \
CIRI_TRACE_ALL_CANDS=1 \
./target/release/ciri-toolkit \
  -i tests/chr1/test.bam \
  -o tests/chr1/CIRI-rs.trace_all \
  -r tests/chr1/chr1.fa \
  -a tests/chr1/chr1.gtf \
  -t 16
```

注意：

- `CIRI_TRACE_ALL_CANDS=1` 会改变 traced read 的候选遍历行为（用于看全量分支）。
- 该模式只用于定位问题，不用于最终 parity 统计。

### 4.3 追踪 Scan2 validator 细分分支

当需要确认某条 read 在 `is_bsj_hg2` 里到底卡在哪个阶段时，打开 `CIRI_TRACE_HG2=1`：

```bash
CIRI_TRACE_READS="A00785:126:HJFMGDRXX:1:1153:23086:11350" \
CIRI_TRACE_HG2=1 \
./target/release/ciri-toolkit \
  -i tests/hg38/diff.subset.bam \
  -o tmp/hg38.trace_hg2 \
  -r tests/hg38/hg38.fa \
  -a tests/hg38/gencode.v29.annotation.gtf \
  -s 0 -t 4 \
  2> tmp/hg38.trace_hg2.log
```

关键标签：

- `TRACE_SCAN2_HG2`
- `stage=enter`
- `stage=sm_*` / `stage=ms_*`
- `stage=circ2_fail`
- `stage=pass`

这类日志适合判断：是 seed 检查失败、linear competition 失败，还是 circ 序列验证失败。

### 4.4 打开 release profiling

用于定位热点，不用于比较结果正确性：

```bash
CIRI_PROFILE_SCAN1=1 ./target/release/ciri-toolkit ...
CIRI_PROFILE_SCAN2=1 ./target/release/ciri-toolkit ...
```

关键输出：

- `PROFILE_SCAN1`
- `PROFILE_SCAN1_HG1`
- `PROFILE_SCAN2`
- `PROFILE_SCAN2_HG2`

## 5. 推荐排查流程（SOP）

1. 运行正常模式，拿到 Rust `.out`。
2. 用 `analyze_diff.py` 看三层差异。
3. 从 `ONLY_*` 集合中挑 2~4 条典型 read。
4. 开 `CIRI_TRACE_READS` 看 Scan1/Scan2 分支与 tag。
5. 如怀疑候选顺序问题，再开 `CIRI_TRACE_ALL_CANDS=1`。
6. 如怀疑 `is_bsj_hg2` 内部分支，再开 `CIRI_TRACE_HG2=1`。
7. 在 Java 源码中对应位置做逐分支比对（变量级）。
8. 只做一个最小改动，立即复测三层指标。
9. 关闭所有 trace/profile 环境变量，做最终验证。

## 6. 单条 read 定位建议

针对具体 read，建议按这个顺序缩小范围：

1. 先在 Rust 最终 `.bsj` 里找 read 是否出现，并看尾列来源是 `scan1` 还是 `scan2`。
2. 如果不在 `.bsj`，先查 `Scan1`：
   - 看是否有 `TRACE_SCAN1_CAND`
   - 看是否进入 `TRACE_SCAN1_HG1`
3. 如果在 `.bsj` 里有，但不在最终 `.out`，再查 `Summary`：
   - 看 `TRACE_SCAN2_CAND`
   - 必要时开 `CIRI_TRACE_ALL_CANDS=1`
   - 必要时看 `TRACE_SCAN2_HG2`
4. 如果 subset 与 full 结论不一致，优先怀疑“缺了中间竞争上下文”，不要直接判定算法错误。

## 7. 常见根因清单（本项目已验证）

- Scan2 索引来源不一致：
  - Java：`BSJ1 -> chrCircSiteMap(HashSet) -> siteArray/siteMap`
  - 若 Rust 直接用原始 BSJ1 行建索引，可能引入重复候选，影响 early-return 路径。
- Scan2 payload 字段布局不一致：
  - 应为 `[site1, site2, strand, signal1, signal2, sum_q]`。
- 配对序列方向计算不一致：
  - 应按“当前 alignment strand”逐条判断，而非按 segment 固定方向。
- 候选遍历顺序不一致：
  - Java 是桶门控 + 方向遍历：`num1` 逆序，`num2` 正序。
- 注释版本不一致：
  - exon 覆盖变化会直接改变 `circRNA_type`、`gene_id`，并可能影响“看起来是否支持”某些 case 的判断。
- subset 缺少中间上下文：
  - 若 full 与 subset 结论不同，说明该 case 依赖未进入最终 `.out` 的竞争候选或 supporting context。

## 8. 变更验收门槛

每次提交前至少满足：

- circ-level：100%
- read-level：100%
- read-assignment-level：100%（若当前任务目标覆盖到 assignment）

并附上命令与摘要指标，确保可复现。

## 9. 输出与记录建议

- 调试产物命名建议：
  - `tmp/<dataset>_<short_name>.out`
  - `tmp/<dataset>_<short_name>.log`
- 每次实验记录四件事：
  - 改了什么
  - 为什么改
  - 指标变化
  - 是否保留

这样可以快速回放“哪一步真正有效”。
