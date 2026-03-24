# CIRI-toolkit（Rust）

`CIRI-toolkit` 是对 Java CIRI3 的 Rust 高性能复现实现，目标是在保持识别逻辑一致的前提下，提供更好的可维护性与执行性能。

## 当前发布状态
- **Functional v1 可用**：核心流程 `Scan1 -> Scan2 -> Summary` 已稳定运行。
- **对齐基线已达成**（`tests/chr1`）：
  - circRNA-level：100%
  - read-level：100%
  - read-assignment-level：100%
- **大规模 hg38 验证已收敛**：在使用与 Java 基线一致的参考与注释版本时，`BSJ reads` 已达到 100% 对齐；剩余差异优先检查注释版本、exon 覆盖与输出口径，而不是先怀疑 Scan1/Scan2 判定链。

## 主要特性
- **Java 一致性优先**：关键判断路径按 Java CIRI3 行为对齐。
- **Scan1 分裂比对识别**：识别 BSJ 候选信号。
- **Scan2 PEM/SMS Rescue**：对歧义读段进行救援判定并统计 FSJ。
- **Summary 合并与过滤**：执行位点聚合、stringency 过滤与最终输出。
- **SAM/BAM 双格式支持**：自动识别输入格式并走对应处理路径。

## 构建

```bash
cargo build --release
```

## 运行示例

```bash
./target/release/ciri-toolkit \
  -i input.sam \
  -o sample.ciri \
  -r reference.fa \
  -a annotation.gtf
```

## 参数说明
- `-i, --in`：输入文件（SAM 或 BAM）
- `-o, --out`：输出前缀
- `-r, --ref`：参考基因组 FASTA
- `-a, --anno`：注释文件 GTF（可选）
- `-m, --mapq`：最小 MAPQ（默认 `10`，与 Java `-U` 默认一致）
- `-s, --stringency`：过滤等级（`0/1/2`，默认 `2`，与 Java `-S` 默认一致）
- `--min-span`：最小环长（默认 `140`，与 Java `-Min` 默认一致）
- `--max-span`：最大环长（默认 `200000`，与 Java `-Max` 默认一致）
- `--linear-range-size-min`：线性竞争区间（默认 `50000`）
- `-t, --threads`：线程数（默认自动）
- `-M, --mem-per-thread`：每线程内存预算（如 `2G`、`512M`，默认 `512M`）

> 默认参数已对齐 CIRI3，一般无需手动设置上述核心阈值。

输出命名规则：

- 若 `-o sample.ciri`
  - 最终结果：`sample.ciri.out`
  - First Scan 合并结果：`sample.ciri.bsj1`
  - Final BSJ 合并结果：`sample.ciri.bsj`
  - First Scan 临时文件：`sample.ciri.bsj1.part_0001.tmp`
  - Second Scan 临时文件：`sample.ciri.bsj2.part_0001.tmp`
  - FSJ 临时文件：`sample.ciri.fsj.part_0001.tmp`
  - `Scan2` 合并完成后会自动删除中间的 `sample.ciri.bsj1`

最终 `.bsj` 会比内部中间格式多一个尾列，标记该行来自 `scan1` 还是 `scan2`。

## 一键差异统计

快速比较 Java 与 Rust 输出差异：

```bash
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result --show-read-ids
```

完整三层对齐检查：

```bash
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 可选调试追踪

按 read ID 定向追踪：

```bash
CIRI_TRACE_READS="simulate:7037,simulate:7050" \
./target/release/ciri-toolkit \
  -i tests/chr1/test.bam \
  -o tests/chr1/CIRI-rs.trace \
  -r tests/chr1/chr1.fa \
  -a tests/chr1/chr1.gtf \
  -t 16
```

仅用于诊断的 Scan2 全候选追踪：

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

> `CIRI_TRACE_ALL_CANDS=1` 仅用于排障，不用于最终 parity 指标统计。

更细的 Scan2 validator 分支追踪：

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

性能 profiling 开关：

```bash
CIRI_PROFILE_SCAN1=1 ./target/release/ciri-toolkit ...
CIRI_PROFILE_SCAN2=1 ./target/release/ciri-toolkit ...
```

> 实际排障时优先把 `stderr` 重定向到 `tmp/*.log`，并在最终验证前关闭所有 trace/profile 环境变量。

## 输出格式

最终 `.out` 为 13 列，兼容 CIRI3 常见下游分析：
1. `circRNA_ID`
2. `chr`
3. `circRNA_start`
4. `circRNA_end`
5. `#junction_reads`
6. `SM_MS_SMS`
7. `#non_junction_reads`
8. `junction_reads_ratio`
9. `circRNA_type`
10. `gene_id`
11. `strand`
12. `junction_reads_ID`
13. `Score`

## 文档索引
- 项目状态：`docs/development-status.md`
- 对齐排障手册：`docs/parity-debug-playbook.md`
- 性能优化路线：`docs/performance-roadmap.md`
- 文档导航：`docs/index.md`
