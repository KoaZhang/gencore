# gencore 性能优化与并行化 Spec

## Why

gencore 当前以单线程串行方式处理整个 BAM 文件，680 万条 reads 的测试数据需要 147 秒。核心瓶颈在于：整个 consensus 流程（读入 → 按 mapping position 聚类 → UMI 聚类 → consensus 合并 → 输出）全部在单线程中完成，且大量使用 `map<string, Pair*>` 等红黑树结构做查找，内存分配/释放频繁。在保持算法语义和输出一致的前提下，通过分阶段优化可以显著提升运行速度。

## What Changes

### 阶段 0：Profiling 与基线建立
- 建立可复现的性能基线（时间、内存、CPU 利用率）
- 使用 `perf` / `gprof` / `valgrind callgrind` 定位热点函数
- 确认构建系统可用，建立自动化测试脚本

### 阶段 1：数据结构与算法微优化（不改并行，不改语义）
- `Cluster::mPairs` / `Group::mPairs`：将 `map<string, Pair*>` 替换为 `unordered_map<string, Pair*>`，降低查找复杂度从 O(logN) 到 O(1)
- `Cluster::clusterByUMI` 中 `umiCount` 同理替换为 `unordered_map`
- `Gencore::mProperClusters` / `mUnProperClusters`：三层嵌套 `map<int, map<int, map<long, Cluster*>>>` 的最内层替换为 `unordered_map`
- `BamUtil::getQName()` / `BamUtil::getUMI()` / `BamUtil::getCigar()` 等频繁调用函数：减少临时 `string` 构造，使用 `string_view` 或缓存
- `Pair::computeScore()` 中的 `new char[]` 分配：预分配或使用对象池
- `Group::makeConsensus()` 中的 `new char[]` / `delete[]`：改为栈分配或复用缓冲区
- `Gencore::mOutSet`（`set<bam1_t*, bamComp>`）：评估是否可用更高效的结构

### 阶段 2：按染色体（contig）级并行
- **核心思路**：BAM 按 coordinate 排序，不同 contig 之间数据完全独立，可按 contig 并行处理
- 将 `Gencore::consensus()` 主循环改为两阶段：
  1. **Phase 1 - 分发**：单线程读 BAM，按 contig tid 分发 reads 到各 contig 的缓冲区
  2. **Phase 2 - 并行处理**：每个 contig 独立执行 `addToCluster → clusterByUMI → consensusMerge → outputPair` 逻辑
  3. **Phase 3 - 有序输出**：按 contig 顺序收集各线程结果，保证输出 BAM 的 coordinate 排序
- Stats 统计需要线程安全：每个 contig 线程持有独立的 `Stats` 对象，最后合并
- Reference 单例需要线程安全读取（只读，天然安全）
- 跨 contig 的 pair（`CrossRefMapped`）需要特殊处理：归入较小 tid 的 contig 组，或单独串行处理

### 阶段 3：I/O 优化
- BAM 读取：评估 htslib 的多线程解码（`hts_set_threads`）
- BAM 写入：各 contig 线程写入临时文件，最后 `sam_cat` 合并；或使用线程安全写入队列
- 输出顺序：保证与原版一致的 coordinate sort

### 阶段 4：细粒度优化（可选，视阶段 1-3 效果）
- `Group::consensusMergeBam()` 中 O(N²) 的 `isPartOf` 检查：优化为 O(N) 或提前剪枝
- `Cluster::clusterByUMI()` 中 UMI 比较的贪心策略：评估是否有更高效的聚类方法
- 内存池：为 `bam1_t` / `Pair` / `Cluster` / `Group` 引入对象池，减少 malloc/free 开销

## Impact

- Affected code:
  - [gencore.cpp](file:///home/zhangrui/gencore/src/gencore.cpp) - 主流程重构（并行化核心）
  - [gencore.h](file:///home/zhangrui/gencore/src/gencore.h) - 类定义变更
  - [cluster.h](file:///home/zhangrui/gencore/src/cluster.h) / [cluster.cpp](file:///home/zhangrui/gencore/src/cluster.cpp) - 数据结构替换
  - [group.h](file:///home/zhangrui/gencore/src/group.h) / [group.cpp](file:///home/zhangrui/gencore/src/group.cpp) - 数据结构替换
  - [pair.h](file:///home/zhangrui/gencore/src/pair.h) / [pair.cpp](file:///home/zhangrui/gencore/src/pair.cpp) - 内存分配优化
  - [bamutil.h](file:///home/zhangrui/gencore/src/bamutil.h) / [bamutil.cpp](file:///home/zhangrui/gencore/src/bamutil.cpp) - 字符串优化
  - [stats.h](file:///home/zhangrui/gencore/src/stats.h) / [stats.cpp](file:///home/zhangrui/gencore/src/stats.cpp) - 线程安全改造 + 合并
  - [Makefile](file:///home/zhangrui/gencore/Makefile) - 添加 pthread 链接、优化编译选项
  - [options.h](file:///home/zhangrui/gencore/src/options.h) / [options.cpp](file:///home/zhangrui/gencore/src/options.cpp) - 新增线程数参数

- Affected behavior:
  - 输出 BAM 的 reads 内容和顺序必须与原版一致（允许因浮点精度导致的 stats 微小差异）
  - gencore.json / gencore.html 的数值必须一致
  - UMI 聚类、consensus、duplex 的核心规则不变

## ADDED Requirements

### Requirement: 性能基线与 Profiling
系统 SHALL 在任何优化前建立可复现的性能基线，包括运行时间、峰值内存、CPU 利用率。

#### Scenario: 建立基线
- **WHEN** 执行原版 gencore 处理测试数据
- **THEN** 记录运行时间、内存占用、输出 BAM 的 md5 / read count / stats 数值作为基线

#### Scenario: Profiling 热点定位
- **WHEN** 使用 perf/callgrind 对原版做 profiling
- **THEN** 识别出 CPU 时间占比最高的前 5 个函数，记录到文档

### Requirement: 数据结构微优化
系统 SHALL 将核心路径中的 `map` 替换为 `unordered_map`，减少临时字符串分配，优化内存分配模式，且不改变算法语义。

#### Scenario: unordered_map 替换
- **WHEN** 将 `Cluster::mPairs`、`Group::mPairs`、`umiCount` 等从 `map` 改为 `unordered_map`
- **THEN** 输出 BAM 与原版完全一致，运行时间缩短

#### Scenario: 字符串与内存优化
- **WHEN** 优化 `BamUtil::getQName()` 等频繁调用函数的临时 string 分配
- **THEN** 输出不变，内存分配次数减少

### Requirement: 按染色体并行
系统 SHALL 支持按 contig 级别的并行处理，通过命令行参数 `-t/--threads` 控制线程数（默认 1，兼容原行为）。

#### Scenario: 多线程运行
- **WHEN** 用户指定 `-t 4` 且输入 BAM 包含多个 contig
- **THEN** 各 contig 的 consensus 处理并行执行，输出 BAM 按 coordinate 排序与原版一致

#### Scenario: 单线程兼容
- **WHEN** 用户不指定 `-t` 或指定 `-t 1`
- **THEN** 行为与原版完全一致

#### Scenario: 跨 contig pair 处理
- **WHEN** 存在跨 contig mapping 的 read pair
- **THEN** 该 pair 被正确处理（归入较小 tid 的 contig 或单独串行处理），输出与原版一致

### Requirement: 线程安全 Stats
系统 SHALL 保证 Stats 统计在多线程环境下的正确性。

#### Scenario: Stats 合并
- **WHEN** 多个 contig 线程各自持有独立 Stats 对象
- **THEN** 最终合并的 Stats 数值与单线程结果一致

### Requirement: 输出一致性验证
系统 SHALL 在每次优化后验证输出与基线的一致性。

#### Scenario: BAM 一致性
- **WHEN** 优化后运行 gencore
- **THEN** 输出 BAM 的 read 数量、总 base 数、mismatch 数、SSCS/DCS 计数与基线一致

#### Scenario: 报告一致性
- **WHEN** 优化后运行 gencore
- **THEN** gencore.json 中各字段数值与基线一致

## MODIFIED Requirements

### Requirement: Gencore::consensus() 主流程
原版单线程串行读取 → 聚类 → consensus → 输出。修改为支持两阶段并行：分发阶段单线程读 BAM 按 contig 分发，处理阶段多线程并行 consensus，输出阶段按序合并。

### Requirement: Options 类
新增 `int threads` 字段（默认 1），新增 `-t/--threads` 命令行参数。

### Requirement: Makefile
新增 `-lpthread` 链接（已有但未实际使用），CFLAGS 添加 `-pthread`，确认 `-O3` 已启用。

## REMOVED Requirements

无。所有原有功能保持不变。
