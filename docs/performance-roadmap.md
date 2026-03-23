# 性能优化路线图（对齐完成后）

本文档定义在功能对齐完成后的性能优化策略。

## 目标范围
- 不改变 Java 对齐结果（行为与输出保持一致）。
- 优先优化 Scan1/Scan2 的 I/O 与热路径分配开销。

## 强制门槛（必须满足）
每次优化后都要通过零差异校验：
- `circ_only_java = 0`
- `circ_only_rust = 0`
- `read_only_java = 0`
- `read_only_rust = 0`
- `read_assignment_only_java = 0`
- `read_assignment_only_rust = 0`

参考命令：

```bash
python tests/analyze_diff.py tests/chr1/CIRI3_result.txt tests/chr1/CIRI-rs.result \
  --show-read-ids --show-read-assignments
```

## 优先级队列

### P0：Scan1/Scan2 I/O 吞吐
- 降低热循环中的临时分配与字符串构造。
- 强化批量写出策略，减少小写入次数。
- 结合数据规模调优分片大小与 flush 策略。

### P1：BAM 解码与解析路径
- 避免不必要的 CIGAR/SEQ 中间字符串物化。
- 评估虚拟位置更新与进度显示开销。
- 保留诊断路径，同时隔离调试分支的运行时成本。

### P2：数据结构局部性
- 检查高频 HashMap/HashSet key 组织方式。
- 在安全前提下以可复用缓冲替代重复 `format!`。
- 通过大 BAM 数据评估缓存命中变化。

## 统一测量规范
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

---
最后更新：2026-03-22
