# CCP-DP 参数敏感性结果使用说明

## 1. 数据是否完整

本次结果包含两个数据集、18种配置和3次重复：

```text
2 datasets × 18 configurations × 3 repeats = 108 runs
```

Success、Mean ED和Reconstruction在重复运行中完全一致，因此相应标准差为0；三次重复
主要用于估计运行时间波动。

默认配置结果为：

| Dataset | Success (%) | Mean ED | Reconstruction (%) | Runtime (s) | Structured decoding (%) |
|---|---:|---:|---:|---:|---:|
| Srinivas-110 | 97.9968 | 0.0491 | 99.9554 | 1.961 ± 0.076 | 10.01 |
| Chandak-108 | 95.3615 | 0.1221 | 99.8869 | 2.931 ± 0.113 | 100.00 |

## 2. 主要结论

在全部扫描范围内，最大Success变化为0.614个百分点，出现在Chandak-108的
Cross-weight scale=0.75。不同参数的最大绝对Success变化为：

| Parameter | Maximum absolute change (pp) |
|---|---:|
| Cross-weight scale | 0.614 |
| Coverage switch | 0.477 |
| Minimum route fraction | 0.080 |
| Reuse threshold | 0.000 |
| Phase-load thresholds | 0.000 |
| Phase band | 0.000 |

这说明默认参数附近总体稳定，但不能笼统写成“所有参数完全不敏感”：coverage switch和
Cross-weight scale对Chandak-108仍有可见影响。

### Coverage switch

- Srinivas-110在20–60 reads之间最大波动仅0.060个百分点；
- Chandak-108在阈值20时下降0.477个百分点，在30时下降0.136个百分点；
- 阈值40–60形成相对稳定区间。

### Cross-weight scale

- Chandak-108在0.75时下降0.614个百分点，在1.25时下降0.273个百分点；
- Srinivas-110最大下降0.100个百分点；
- 默认scale=1.0在两个数据集上均给出最优或并列最优结果。

### Minimum route fraction

- Srinivas-110将下限从10%降至5%后，实际解码7.28%的簇，Success下降0.080个百分点；
- 提高至15%和20%仅增加0.010个百分点，但运行时间分别增加至2.947 s和3.340 s；
- 默认10%以1.961 s达到97.9968%，是更合理的精度—效率折中；
- Chandak-108由high phase load触发100%解码，因此route floor不影响结果。

### Reuse、phase-load thresholds和phase band

这些扫描没有改变Success。正确解释是：在测试范围内，它们没有改变当前数据集的活跃
路由，或没有改变DP的最终最优路径；不能据此声称这些参数在所有数据上都没有作用。

## 3. 正文图注

> **Parameter sensitivity and efficiency trade-offs of CCP-DP.** Changes in
> exact reconstruction rate relative to the default configuration are shown
> for (a) the read-count threshold separating the two local-guide procedures,
> (b) the scale applied to the automatically calibrated cross-cluster weight,
> and (c) the minimum fraction of clusters routed to structured decoding.
> (d) Runtime under different routing floors. Dashed lines indicate the default
> values. Runtime shows mean ± s.d. over three runs; reconstruction outputs were
> deterministic across repeated runs.

## 4. 补充图注

### Supplementary Fig. S1

> **Sensitivity of exact reconstruction to all tested parameters.** Each panel
> reports the change in exact reconstruction rate relative to the default
> configuration. Dashed lines denote the default values. Parameters that did
> not alter the active routing regime or final optimal path produced identical
> reconstructions over the tested ranges.

### Supplementary Fig. S2

> **Sensitivity of residual edit distance to all tested parameters.** Values
> are reported relative to the default configuration; negative values indicate
> lower residual edit distance.

### Supplementary Fig. S3

> **Effect of parameter variation on structured-decoding frequency.**
> Chandak-108 remained in the dense-decoding regime, whereas the routing floor
> directly controlled the decoded fraction on Srinivas-110.

### Supplementary Fig. S4

> **Runtime sensitivity.** Runtime change relative to the default configuration
> is reported for all parameters. Each point summarizes three runs.

## 5. 英文正文结果段落

> We evaluated the sensitivity of CCP-DP to six routing and decoding parameters
> on Srinivas-110 and Chandak-108 while changing one parameter at a time. Across
> all tested settings, the largest absolute change in exact reconstruction was
> 0.61 percentage points. The coverage switch was stable over 40–60 reads on
> Chandak-108 and varied by at most 0.06 points on Srinivas-110. Scaling the
> cross-cluster contribution away from its default value caused the most visible
> degradation, but the maximum decrease remained 0.61 points. On Srinivas-110,
> reducing the routing floor from 10% to 5% decreased success by 0.08 points,
> whereas increasing it to 15–20% improved success by only 0.01 points while
> increasing runtime from 1.96 s to 2.95–3.34 s. Chandak-108 remained in the
> high-phase-load dense-decoding regime and was therefore unaffected by the
> routing floor. The reuse threshold, phase-load thresholds and phase-band width
> did not change exact reconstruction within the tested ranges, because these
> perturbations did not alter the active routing regime or the final optimal
> path. Overall, the default configuration lies in a stable region and provides
> a favorable accuracy–runtime trade-off without dataset-specific retuning.

## 6. 中文翻译

> 我们在Srinivas-110和Chandak-108上采用单因素方式，评估了CCP-DP对六个路由与解码
> 参数的敏感性。在全部测试设置中，精确重建率的最大绝对变化为0.61个百分点。
> Chandak-108的coverage switch在40–60 reads范围内保持稳定，而Srinivas-110上的最大
> 波动仅为0.06个百分点。将跨簇贡献偏离默认尺度会造成最明显的性能下降，但最大降幅仍
> 为0.61个百分点。在Srinivas-110上，将最低路由比例从10%降低至5%使成功率下降0.08个
> 百分点；将其提高至15%–20%仅带来0.01个百分点的提升，却使运行时间从1.96 s增加至
> 2.95–3.34 s。Chandak-108始终处于高phase-load的全簇解码状态，因此不受最低路由比例
> 影响。在测试范围内，reuse threshold、phase-load thresholds和phase-band width没有改变
> 精确重建结果，因为这些扰动没有改变活跃路由或最终最优路径。总体而言，默认配置位于
> 稳定区域，并在不针对特定数据集重新调参的情况下取得了较好的精度—时间平衡。

## 7. 推荐放置位置

- 正文：`Fig_parameter_sensitivity_main`和第5节英文段落；
- 补充材料：Fig. S1–S4、`table_full_sensitivity.tex`和完整source data；
- `table_default_results.tex`可放补充材料，若正文已有主结果表则无需重复。
