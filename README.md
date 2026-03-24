# CIRI-toolkit（Rust）

`CIRI-toolkit` 是 CIRI3 的高性能 Rust 升级版，实现了高速BAM/SAM处理，低内存的同时保证结果的高度一致性。

实测 ~80G BAM (~300G SAM) 文件，CIRI3 总用时 80min (scan1 67min + scan2 14min)。CIRI-toolkit 用时 6 min (1.25min + 4.5min)；提速>10x

后续将加入对CIRI-AS/CIRI-full功能的整合，实现环形RNA的高效分析

## 安装与构建

```bash
cargo build --release
```

编译后的可执行文件位于：

```bash
./target/release/ciri-toolkit
```

## 快速开始

```bash
./target/release/ciri-toolkit \
  -i input.bam \
  -o sample.ciri \
  -r reference.fa \
  -a annotation.gtf
```

常见全转录组数据示例：

```bash
./target/release/ciri-toolkit \
  -i RNA015434_S1.bam \
  -o RNA015434 \
  -r /data/public/database/gencode/hg38/_BWAindex/hg38.fa \
  -a /data/public/database/gencode/hg38/gencode.v44.annotation.gtf \
  -s 0 \
  -t 16
```

## 参数说明
- `-i, --in`：输入文件（SAM 或 BAM）
- `-o, --out`：输出前缀
- `-r, --ref`：参考基因组 FASTA
- `-a, --anno`：注释文件 GTF（可选）
- `-m, --mapq`：最小 MAPQ（默认 `10`，与 CIRI3 `-U` 默认一致）
- `-s, --stringency`：过滤等级（`0/1/2`，默认 `2`，与 CIRI3 `-S` 默认一致）
- `--min-span`：最小环长（默认 `140`，与 CIRI3 `-Min` 默认一致）
- `--max-span`：最大环长（默认 `200000`，与 CIRI3 `-Max` 默认一致）
- `--linear-range-size-min`：线性竞争区间（默认 `50000`）
- `-t, --threads`：线程数（默认为CPU可用核心数）
- `-M, --mem-per-thread`：每线程内存预算（如 `2G`、`512M`，默认 `512M`）

> 默认参数已对齐 CIRI3，推荐使用 `-s 0` 输出所有潜在 circRNA 后手动过滤，其他参数一般无需手动设置。

## 输出文件

当 `-o <prefix>` 时，程序会生成：

- 环形RNA识别结果：`<prefix>.out`
- 运行日志：`<prefix>.log`
- BSJ reads 比对情况：`<prefix>.bsj`

## 结果格式

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

## 开发文档

开发、对齐排障、性能优化与验证说明统一放在 `docs/`：

- 文档导航：`docs/index.md`
- 项目状态：`docs/development-status.md`
- 对齐排障手册：`docs/parity-debug-playbook.md`
- 性能优化路线：`docs/performance-roadmap.md`

---
最后更新：2026-03-24
