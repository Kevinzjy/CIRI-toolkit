# circRNA 候选 BSJ 位点重评分模型设计思路

本文档是探索性设计笔记，不代表当前主流程实现。当前 CIRI-toolkit 主流程仍保持 CIRI3-compatible `Scan1 -> Scan2 -> Summary` parity；已落地的默认扩展是 `<prefix>.segments`、major isoform GTF/FASTA 和 IGV review sidecar。后续活跃路线是 multi-isoform usage 与 multi-sample integration，而不是先引入新的 BSJ scoring 模型。

## 1. 背景与核心问题

现有 circRNA 识别方法在实际应用中往往依赖以下两类先验：

1. **canonical 剪接信号**，尤其是 GT/AG；
2. **现有 GTF/GFF 注释中的已知外显子边界**。

这种策略虽然能够显著降低假阳性，但也带来两个明显局限：

- **未注释 BSJ 位点容易被漏检**；
- **非典型但真实的 back-splicing 位点可能被过早过滤**。

这一问题在植物中尤其突出，因为植物剪接位点识别通常比动物更依赖整体序列环境，而不是高度固定的单一 motif；同时植物 circRNA 的成环规则和动物并不完全一致，直接套用哺乳动物规则容易损失召回率。

因此，更合理的策略不是把 GT/AG 或 GTF 作为硬过滤条件，而是：

> **先高召回地收集所有 candidate BSJ reads，再针对每个候选区域内可能的 donor/acceptor 组合进行重评分与重排序，从而判断哪个 junction site 最可能是真实 circRNA 的 back-splice 位点。**

---

## 2. 总体方法框架

建议采用一个 **两阶段（high-recall detection + BSJ rescoring）** 的框架。

### Stage I：高召回候选检测

目标是尽可能保留所有潜在 BSJ，而不是在这一阶段追求极高特异性。

可选输入包括：

- STAR chimeric / split alignment 输出；
- BWA/PCC 类 partial mapping 结果；
- paired-end support 信息；
- long-read split alignment 结果（如纳米孔数据）；
- 现有 circRNA 软件的原始候选输出（如 CIRI、CIRCexplorer、find_circ、CircPlant 等）。

这一阶段建议：

- **不强依赖 GTF 注释**；
- **不把 GT/AG 作为硬门槛**；
- 仅做最基础的 mapping quality、anchor length、重复比对过滤；
- 将所有 candidate BSJ reads 聚类为候选 BSJ 区域（candidate cluster）。

### Stage II：局部断点枚举与重评分

对于每一个 candidate cluster，在候选 read 支持的局部窗口内枚举可能的：

- 下游 donor site；
- 上游 acceptor site；
- donor × acceptor 组合对应的候选 BSJ。

然后对每个候选 BSJ 计算综合得分，并输出：

- 最优 donor/acceptor 组合；
- 该组合的 posterior probability / confidence score；
- 候选 BSJ 的排序列表。

---

## 3. 需要解决的关键问题

这个模型的目标不应定义为：

> “仅根据 splice motif 判断某个位点是不是 circRNA。”

而更应定义为：

> **“对于比对给出的多个局部可行断点解释，利用剪接信号、比对一致性和聚合证据，选出最可能的真实 back-splice donor/acceptor 组合。”**

这是一个 **局部排序 / 断点解析** 问题，而不是简单的全局二分类问题。

---

## 4. 特征设计

建议将特征分为四大类。

### 4.1 剪接位点强度特征（splice-site prior）

用于回答：这个 donor/acceptor 组合从“剪接语法”角度是否合理。

可包括：

#### 5' donor 相关
- 是否为 GT；
- 是否为 GC 等次典型位点；
- donor 周围序列窗口（如 -3~+6）的 PWM / MaxEnt 风格打分；
- 物种特异 donor motif 概率。

#### 3' acceptor 相关
- 是否为 AG；
- acceptor 周围序列窗口（如 -20~+3）的 PWM / MaxEnt 风格打分；
- 上游局部 U-rich / AU-rich 程度；
- 弱 branch point-like motif（如植物中可弱参考 YURAY）信号；
- 物种特异 acceptor motif 概率。

#### 组合特征
- donor 与 acceptor 的联合分数；
- 两侧 flanking intron/exon 的碱基组成差异；
- 是否与植物已知短内含子、U-rich 背景相一致。

> 注意：这里不要只做 canonical / non-canonical 的二元判定，而应输出连续得分。

---

### 4.2 比对与断点解释特征（alignment likelihood）

用于回答：给定某一 donor/acceptor 组合，当前 reads 是否更支持这一断点解释。

可包括：

- 两侧 split anchor length；
- 左右 anchor 是否平衡；
- 断点附近 mismatch 数；
- indel 数；
- soft-clipping 模式；
- microhomology 长度；
- 是否落在低复杂度区域；
- 是否位于重复序列 / 多重比对区域；
- mapping quality；
- 同一 read 在多个候选断点解释下的相对比对得分差值。

这一类特征非常关键，因为很多假阳性并不是“motif 不合理”，而是“比对解释不唯一”。

---

### 4.3 junction 聚合证据特征（cross-read consistency）

用于回答：多个 reads 是否一致支持同一个 BSJ，而不是各自支持相邻但不一致的伪断点。

可包括：

- 支持该 BSJ 的独立 reads 数；
- 去 PCR duplicates 后的支持 read 数；
- 不同 read 的断点是否高度收敛；
- 正负链一致性；
- paired-end mate 是否落在合理 circRNA 区域；
- 跨样本复现性；
- 不同测序批次的一致性；
- junction 支持是否集中于 read 末端（防模板切换伪影）；
- 是否存在 full-length circular transcript 证据。

真实 circRNA 的一个重要特征是：**支持 reads 往往收敛到同一个 donor/acceptor 组合**。

---

### 4.4 外部软先验特征（soft prior）

用于提高排序能力，但不作为硬性限制。

可包括：

- 是否位于已注释 exon boundary；
- 是否位于已知基因区域内；
- donor/acceptor 是否与同一线性转录本兼容；
- 是否存在对应 forward splice junction (FSJ) 支持；
- 是否被多个 circRNA 工具共同检测；
- 是否在 RNase R 样本中富集；
- flanking intron 是否存在局部反向互补序列；
- flanking 区域是否存在 TE/repeat pairing 特征。

> 在植物中，反向互补内含子并不一定像动物那样普遍，因此这类特征只适合作为辅助信息。

---

## 5. 模型形式建议

### 5.1 局部排序模型（首选）

对每一个 candidate cluster：

- 枚举局部 donor × acceptor 候选；
- 为每个组合提取上述特征；
- 训练一个 ranker，对候选组合排序；
- 输出 top-1 断点及其置信度。

适合的模型包括：

- XGBoost / LightGBM ranking；
- pairwise ranking loss；
- listwise ranking；
- small neural scorer。

这种形式最贴合实际问题，因为你真正要解决的是：

> “在多个可能断点中，哪个断点最合理？”

而不是抽象地问某个位点是不是 circRNA。

---

### 5.2 半监督 / PU-learning

如果缺乏完整金标准，可采用正例-未标注学习。

#### 正例来源
- RNase R 富集且多工具共识的高可信 BSJ；
- 有 full-length circRNA 支持的位点；
- 多个样本复现、断点高度一致的位点；
- orthogonal validation 支持的已知 circRNA。

#### 负例来源
- scrambled junction；
- 重复区伪断点；
- template-switching 高风险区域；
- read-through / tandem duplication 形成的假 backsplice；
- 仅单条 read、断点发散严重的低可信候选。

未标注样本可作为混合背景，通过 PU-learning 或 semi-supervised learning 学习决策边界。

---

### 5.3 EM / latent variable 模型

如果希望更严格处理 read-level 不确定性，可将“真实断点”视为潜变量。

基本思路：

1. 初始化 splice-site prior；
2. 对每条 read 在多个候选 donor/acceptor 组合之间分配 posterior；
3. 聚合高 posterior 断点更新模型参数；
4. 迭代优化直到收敛。

这种做法尤其适合：

- 未注释位点很多；
- reads 噪声较大；
- long-read 数据断点误差较复杂。

---

## 6. 训练数据构建策略

### 不建议的做法

- 直接把所有 candidate BSJ reads 当正样本；
- 直接把 canonical GT/AG 且在注释边界上的位点定义为真值；
- 直接把未注释位点定义为负例。

这样会使模型学到严重偏差：

- 偏向高表达、已注释、canonical 的动物式 circRNA；
- 把未注释真实位点错学成假阳性；
- 把 mapping artifact 学成 motif 规律。

### 推荐做法

构建一个分层训练集：

#### 高可信正例
- RNase R enriched；
- 多工具交集；
- 多 reads 支持；
- full-length 支持；
- 手工验证或数据库验证。

#### 困难负例（hard negatives）
- 与正例位点距离很近但比对稍差的竞争断点；
- 同一 read 可支持的局部替代 donor/acceptor；
- 重复序列附近伪断点；
- 仅因微同源产生的虚假 backsplice。

#### 未标注样本
- 其他 candidate BSJ，用于半监督建模。

> 最佳策略不是“随机负采样”，而是**局部竞争位点负采样**，这样模型会真正学会断点判别。

---

## 7. 针对植物 circRNA 的特别建议

植物系统中需要避免直接套用动物规则。建议特别加入以下植物相关特征：

- flanking intron 的 AU-rich / U-rich 程度；
- intron 长度分布与植物背景模型的偏离程度；
- donor / acceptor 周围植物特异 PWM；
- 弱 branch-point-like motif 的位置分布；
- TE 相关配对序列；
- 物种特异 repeat family 特征；
- exon/intron GC contrast。

同时需要注意：

- 不要过度依赖反向互补内含子作为成环必要条件；
- 不要把非 GT/AG 位点直接视为假阳性；
- 不要把“未注释”直接等同于“错误”。

---

## 8. 一个可实现的打分函数示意

可以将每个候选 BSJ 的最终分数写成：

```text
Score(BSJ) = 
    w1 * AlignmentLikelihood
  + w2 * SpliceSitePrior
  + w3 * CrossReadConsistency
  + w4 * AnnotationSoftPrior
  + w5 * CircularityEvidence
```

其中：

- **AlignmentLikelihood**：该 donor/acceptor 组合对所有支持 reads 的解释能力；
- **SpliceSitePrior**：物种特异 splice motif、U-rich/AU-rich 背景等；
- **CrossReadConsistency**：多 reads 是否一致收敛到该断点；
- **AnnotationSoftPrior**：是否靠近注释边界、是否兼容已知转录本；
- **CircularityEvidence**：RNase R、full-length、paired-end、跨样本复现等证据。

最终输出应包括：

- top-ranked BSJ；
- posterior probability / calibrated confidence；
- alternative junction candidates；
- 各特征分解得分，便于解释。

---

## 9. 与现有软件的关系

这个框架不一定要替代现有 circRNA 软件，更适合作为一个 **rescoring / refinement module**：

- 上游使用 STAR、BWA、CIRI、CircPlant 等生成高召回候选；
- 下游使用本模型对 BSJ 精细断点进行重评分；
- 最终输出更稳健、对未注释位点更友好的高置信 circRNA 集合。

因此它的定位可以是：

> **A species-aware back-splice junction refinement model for de novo circRNA discovery.**

---

## 10. 核心创新点总结

这个思路的创新不在于“再做一个 motif 过滤器”，而在于以下几点：

1. **把 circRNA 鉴定问题重定义为局部断点排序问题**；
2. **弱化对 GTF 注释和 GT/AG 硬过滤的依赖**；
3. **将 splice motif 从硬规则改为 soft prior**；
4. **显式整合 read-level alignment likelihood 与 junction-level 聚合证据**；
5. **允许未注释和非典型位点在后验评分中被保留**；
6. **可构建植物特异的 species-aware 模型，而非直接套用动物规则。**

---

## 11. 最终一句话概括

> **最佳策略不是“先用 motif 定义 circRNA”，而是“先最大化召回 candidate BSJ，再利用物种特异剪接先验、比对一致性和聚合证据，对局部多个 donor/acceptor 组合进行重评分，从而解析最可能的真实 back-splice junction。”**
