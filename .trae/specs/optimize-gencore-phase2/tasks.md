# Tasks

- [x] Task 1: 统一串行与并行路径，消除结果差异
  - [x] SubTask 1.1: 重构 `consensusSerial()` 为调用 `processContig`（与并行路径一致的全量聚类）
  - [x] SubTask 1.2: 移除 `addToProperCluster` 中的增量处理逻辑（每 10000 条处理一次），简化为直接调用 `processContig`
  - [x] SubTask 1.3: 编译、运行 `-t 1` 和 `-t 4`，确认输出完全一致
  - [x] SubTask 1.4: 与原版基线对比，记录差异并确认语义正确性

- [x] Task 2: `Group::consensusMergeBam` isPartOf O(N²) 优化
  - [x] SubTask 2.1: 尝试 CIGAR 分组优化（回滚：小 Group 下排序开销大于节省的比较次数）
  - [x] SubTask 2.2: 保留原始双循环实现（对小 Group 更高效）
  - [x] SubTask 2.3: 编译、运行测试，确认输出一致

- [x] Task 3: `Cluster::clusterByUMI` UMI 聚类优化
  - [x] SubTask 3.1: 尝试 pairUmis 缓存 + 批量删除（回滚：额外内存分配开销抵消了收益）
  - [x] SubTask 3.2: 保留原始逐个删除方式
  - [x] SubTask 3.3: 编译、运行测试，确认输出一致

- [x] Task 4: `Pair::computeScore` 内存分配优化
  - [x] SubTask 4.1: 将 `mLeftScore` / `mRightScore` 从 `char*` 改为 `unique_ptr<char[]>`
  - [x] SubTask 4.2: 更新 `Pair::~Pair()` 中的 `delete[]` 为自动释放
  - [x] SubTask 4.3: 更新 `Group::consensusMergeBam` 和 `Group::makeConsensus` 中对 score 指针的使用
  - [x] SubTask 4.4: 编译、运行测试，确认输出一致

- [x] Task 5: `BamUtil::getQName` string_view 优化
  - [x] SubTask 5.1: 评估 C++17 兼容性（项目使用 C++11，string_view 不可用）
  - [x] SubTask 5.2: 跳过 string_view 优化，保留现有实现

- [x] Task 6: 单 contig 内分片并行
  - [x] SubTask 6.1: 在 `processContig` 中添加分片逻辑：当 reads 数量超过阈值时，按 position range 切分为多个 chunk
  - [x] SubTask 6.2: 实现 `processContigChunk(int tid, vector<bam1_t*>& reads, int startIdx, int endIdx)` 方法
  - [x] SubTask 6.3: 实现增量处理逻辑 `processCompletedClusters`，每 10000 条 reads 处理已完成的 cluster，减少峰值内存
  - [x] SubTask 6.4: 各 chunk 并行处理，结果按 chunk 顺序合并
  - [x] SubTask 6.5: 编译、运行 `-t 1` 和 `-t 4`，确认输出一致
  - [x] SubTask 6.6: 性能对比：记录分片前后的运行时间

- [x] Task 7: 综合验证
  - [x] SubTask 7.1: 用测试数据做完整的 `-t 1` vs `-t 4` vs `-t 8` 对比
  - [x] SubTask 7.2: 记录各优化前后的性能数据（时间、内存）
  - [x] SubTask 7.3: 如有输出差异，分析并解释差异来源

# Task Dependencies
- Task 1 是所有后续任务的前置依赖（统一路径后才能做后续优化）
- Task 2 和 Task 3 可并行执行
- Task 4 和 Task 5 可并行执行，且与 Task 2/3 无依赖
- Task 6 依赖 Task 1 完成（需要统一路径后才能做分片并行）
- Task 7 依赖所有前置任务完成
