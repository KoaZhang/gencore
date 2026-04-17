# gencore 第二轮性能优化 Spec

## Why

第一轮优化已实现：unordered_map 替换、字符串/内存优化、按 contig 并行。当前性能：串行 146s，并行 -t 4 约 114s（提速 22%）。但并行效果受限于数据集中在少数大 contig（chr1 占大部分时间），且并行路径与串行路径存在约 0.15% 的结果差异（全量聚类 vs 增量聚类）。第二轮优化目标：消除并行/串行结果差异、在单个大 contig 内实现并行、优化核心算法的 O(N²) 热点。

## What Changes

### 改动 1：统一串行与并行路径的聚类语义
- 将 `consensusSerial()` 的增量聚类改为全量聚类（与 `processContig` 一致），消除 `-t 1` 与 `-t 4` 的结果差异
- 具体做法：串行路径也先全部读入再按 contig 处理，复用 `processContig` 方法
- **BREAKING**：串行路径的输出将与原版有微小差异（与当前并行路径一致），但语义更准确

### 改动 2：单 contig 内按 position range 分片并行
- 对大 contig（reads 数量超过阈值，如 100000），按 position range 切分为多个 chunk
- 每个 chunk 独立执行 clustering + consensus（同一 chunk 内的 reads 按 mapping position 聚类）
- chunk 之间可能有边界 pair（跨 chunk 的 pair），需要特殊处理
- 最终按 chunk 顺序合并输出

### 改动 3：`Group::consensusMergeBam` 中 O(N²) `isPartOf` 优化
- 当前：对每个 pair，与所有其他 pair 做 `isPartOf` 比较，O(N²)
- 优化：先按 CIGAR 分组，只有相同 CIGAR 模式的 reads 才需要比较，大幅减少比较次数
- 或：先按 read length 排序，利用 `isPartOf` 的单调性剪枝

### 改动 4：`Cluster::clusterByUMI` 中 UMI 聚类优化
- 当前：贪心策略每次遍历所有 UMI 找 top，O(N²)
- 优化：先按 UMI 排序/分组，相同 UMI 直接合并，相近 UMI 再做 diff 比较
- 或：使用 `unordered_map<string, int>` 的 `umiCount` 找 top 时，维护一个 max-heap

### 改动 5：`Pair::computeScore` 内存分配优化
- 当前：每次 `computeScore` 都 `new char[l_qseq]`
- 优化：使用 `unique_ptr<char[]>` 或预分配的 buffer pool

### 改动 6：`BamUtil::getQName` / `getUMI` 减少临时 string 构造
- `getQName` 在 `Cluster::addRead`、`Group::addRead`、`Pair::getQName` 中频繁调用
- 优化：使用 `string_view` 或直接在 bam1_t 上比较 qname，避免构造临时 string

## Impact

- Affected code:
  - [gencore.cpp](file:///home/zhangrui/gencore/src/gencore.cpp) - 统一串行/并行路径、单 contig 分片并行
  - [gencore.h](file:///home/zhangrui/gencore/src/gencore.h) - 新增分片并行方法声明
  - [cluster.cpp](file:///home/zhangrui/gencore/src/cluster.cpp) - UMI 聚类优化
  - [group.cpp](file:///home/zhangrui/gencore/src/group.cpp) - isPartOf 优化
  - [pair.cpp](file:///home/zhangrui/gencore/src/pair.cpp) - 内存分配优化
  - [bamutil.h](file:///home/zhangrui/gencore/src/bamutil.h) / [bamutil.cpp](file:///home/zhangrui/gencore/src/bamutil.cpp) - string_view 优化

- Affected behavior:
  - 串行路径输出将与原版有微小差异（与当前并行路径一致），但语义更准确
  - 并行路径结果不变
  - UMI 聚类、consensus、duplex 的核心规则不变

## ADDED Requirements

### Requirement: 统一串行与并行路径
系统 SHALL 在 `-t 1` 和 `-t 4` 时产生完全一致的输出。

#### Scenario: 串行并行一致性
- **WHEN** 分别以 `-t 1` 和 `-t 4` 运行 gencore
- **THEN** 输出 BAM 的 read 数量、SSCS/DCS 计数完全一致

### Requirement: 单 contig 内分片并行
系统 SHALL 支持在单个大 contig 内按 position range 分片并行处理。

#### Scenario: 大 contig 分片
- **WHEN** 输入 BAM 中某 contig 的 reads 数量超过阈值（默认 100000）
- **THEN** 该 contig 被切分为多个 chunk 并行处理，各 chunk 结果按序合并

#### Scenario: 跨 chunk pair 处理
- **WHEN** 存在跨 chunk 边界的 read pair
- **THEN** 该 pair 被正确处理（归入左 chunk 或右 chunk），输出与不分片时一致

### Requirement: isPartOf O(N²) 优化
系统 SHALL 优化 `Group::consensusMergeBam` 中的 `isPartOf` 检查，减少不必要的比较。

#### Scenario: CIGAR 分组优化
- **WHEN** 一个 Group 中有 N 个 pairs
- **THEN** `isPartOf` 比较次数从 O(N²) 降低到 O(N * K)，其中 K 是相同 CIGAR 模式的平均 pair 数

### Requirement: UMI 聚类优化
系统 SHALL 优化 `Cluster::clusterByUMI` 中的 top UMI 查找。

#### Scenario: top UMI 查找优化
- **WHEN** 一个 Cluster 中有 N 个 pairs
- **THEN** top UMI 查找从 O(N) 降到 O(logN) 或 O(1)（使用 max-heap 或排序）

### Requirement: Pair 内存分配优化
系统 SHALL 减少 `Pair::computeScore` 中的动态内存分配。

#### Scenario: unique_ptr 替换
- **WHEN** `computeScore` 被调用
- **THEN** 使用 `unique_ptr<char[]>` 而非 `new char[]`，避免内存泄漏风险

### Requirement: 字符串构造优化
系统 SHALL 减少 `BamUtil::getQName` 等函数的临时 string 构造。

#### Scenario: string_view 优化
- **WHEN** `Cluster::addRead` 或 `Group::addRead` 查找 qname
- **THEN** 使用 `string_view` 比较，避免构造临时 string

## MODIFIED Requirements

### Requirement: Gencore::consensus() 串行路径
串行路径不再使用增量聚类（`addToProperCluster` 的每 10000 条处理逻辑），改为与并行路径一致的全量聚类。`consensusSerial()` 将被简化为调用 `processContig`。

### Requirement: Gencore::consensusParallel() 分片逻辑
`consensusParallel()` 新增对大 contig 的分片并行处理。

## REMOVED Requirements

### Requirement: addToProperCluster 增量处理
**Reason**: 增量聚类导致串行/并行结果不一致，且全量聚类语义更准确
**Migration**: 串行路径改用 `processContig` 全量处理
