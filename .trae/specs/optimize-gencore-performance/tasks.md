# Tasks

- [x] Task 0: 构建验证与性能基线建立
  - [x] SubTask 0.1: 修复 Makefile 使项目可编译（确认 htslib 头文件/库路径、inc 目录符号链接等）
  - [x] SubTask 0.2: 编译原版 gencore，确认可运行
  - [x] SubTask 0.3: 用测试数据运行原版，记录基线：运行时间、输出 BAM 的 read count / base count / mismatch / SSCS / DCS 数值
  - [x] SubTask 0.4: 编写自动化测试脚本 `run_test.sh`，自动构建、运行、对比结果
  - [x] SubTask 0.5: 使用 `perf record` / `perf report` 对原版做 profiling，记录热点函数排名

- [x] Task 1: 数据结构微优化（map → unordered_map）
  - [x] SubTask 1.1: `Cluster::mPairs` 从 `map<string, Pair*>` 改为 `unordered_map<string, Pair*>`（cluster.h + cluster.cpp）
  - [x] SubTask 1.2: `Group::mPairs` 从 `map<string, Pair*>` 改为 `unordered_map<string, Pair*>`（group.h + group.cpp）
  - [x] SubTask 1.3: `Cluster::clusterByUMI()` 中 `umiCount` 从 `map<string, int>` 改为 `unordered_map<string, int>`
  - [x] SubTask 1.4: `Gencore::mProperClusters` / `mUnProperClusters` 最内层 `map<long, Cluster*>` 改为 `unordered_map<long, Cluster*>`
  - [x] SubTask 1.5: 编译、运行测试、对比输出一致性、记录时间变化

- [x] Task 2: 字符串与内存分配优化
  - [x] SubTask 2.1: `BamUtil::getQName()` 返回 `string_view` 或减少临时 string 构造（评估影响范围后决定方案）
  - [x] SubTask 2.2: `BamUtil::getCigar()` 使用预分配 stringstream 或直接拼接
  - [x] SubTask 2.3: `Pair::computeScore()` 中 `mLeftScore` / `mRightScore` 的 `new char[]` 改为 `unique_ptr<char[]>` 或对象池
  - [x] SubTask 2.4: `Group::makeConsensus()` 中 `seqBak` / `qualBak` 的 `new/delete` 改为栈上 `vector<char>` 或复用缓冲区
  - [x] SubTask 2.5: 编译、运行测试、对比输出一致性、记录时间变化

- [x] Task 3: 按染色体并行 - 基础设施
  - [x] SubTask 3.1: `Options` 类新增 `int threads` 字段（默认 1）和 `-t/--threads` 命令行参数
  - [x] SubTask 3.2: `Stats` 类新增 `merge(Stats* other)` 方法，支持多 Stats 对象合并
  - [x] SubTask 3.3: `Stats::addRead()` / `addMolecule()` / `addCluster()` 等方法确保无线程共享数据竞争（通过独立 Stats 对象实现）
  - [x] SubTask 3.4: Makefile 添加 `-pthread` 编译和链接选项

- [x] Task 4: 按染色体并行 - 核心重构
  - [x] SubTask 4.1: 重构 `Gencore::consensus()` 为两阶段：Phase 1 单线程读 BAM 按 contig 分发到 `map<int, vector<bam1_t*>>` 缓冲区
  - [x] SubTask 4.2: 实现 `processContig(int tid, vector<bam1_t*>& reads)` 方法，封装单 contig 的 addToCluster → clusterByUMI → consensusMerge → outputPair 逻辑
  - [x] SubTask 4.3: 处理跨 contig pair：在分发阶段识别 `CrossRefMapped` reads，归入较小 tid 的 contig 组或单独串行处理
  - [x] SubTask 4.4: Phase 2 使用 `std::thread` 或 `std::async` 并行调用 `processContig`，线程数由 `-t` 控制
  - [x] SubTask 4.5: Phase 3 按 contig 顺序收集各线程输出的 `vector<bam1_t*>`，按 coordinate 排序后写入 BAM
  - [x] SubTask 4.6: 编译、运行测试、对比输出一致性、记录时间变化

- [x] Task 5: I/O 优化
  - [x] SubTask 5.1: BAM 读取端启用 `hts_set_threads(in, 1)` 利用 htslib 内部解压线程
  - [x] SubTask 5.2: BAM 写入端评估：各 contig 写临时文件后 sam_merge，或使用全局写入队列
  - [x] SubTask 5.3: 编译、运行测试、对比输出一致性、记录时间变化

- [x] Task 6: 综合验证与文档
  - [x] SubTask 6.1: 用测试数据做完整的单线程 vs 多线程对比，确认输出一致
  - [x] SubTask 6.2: 记录各阶段优化前后的性能数据（时间、内存）
  - [x] SubTask 6.3: 如有输出差异，分析并解释差异来源

# Task Dependencies
- Task 0 是所有后续任务的前置依赖
- Task 1 和 Task 2 可并行执行，互不依赖
- Task 3 是 Task 4 的前置依赖
- Task 4 依赖 Task 3 完成
- Task 5 依赖 Task 4 完成（需要并行框架就绪后才能优化 I/O）
- Task 6 依赖所有前置任务完成
