# ALLCools 链方向修改和UMI矫正机制详解

## 概述

ALLCools 是一个用于单细胞甲基化数据分析的工具包，其中包含了两个重要的数据处理机制：
1. **链方向修改** - 统一不同比对工具的链方向表示（第一步）
2. **UMI矫正** - 用于去除PCR重复并矫正测序错误（第二步）

本文档详细解释这两个机制的原理、实现和应用。**注意：ALLCools的处理顺序是先进行链方向转换，再进行UMI矫正。**

## 1. 链方向修改机制（第一步处理）

### 1.1 背景问题

不同的BS-seq比对工具对链方向的处理方式不同：

#### Bismark 方式
- C→T转换的reads强制映射到正链 (+)
- G→A转换的reads强制映射到负链 (-)
- 使用XG标签标记转换类型 (XG="CT" 或 XG="GA")

#### Hisat-3n 方式
- 保持R1/R2的原始方向
- 同一转换类型可能有不同的链方向
- 使用YZ标签标记转换类型 (YZ="+" 或 YZ="-")

### 1.2 ALLCools 解决方案

ALLCools 通过 `_convert_bam_strandness` 函数统一链方向：

```python
def _convert_bam_strandness(in_bam_path, out_bam_path):
    for read in in_bam:
        # 自动检测标签类型
        if read.has_tag("YZ"):
            is_ct_func = _is_read_ct_conversion_hisat3n
        elif read.has_tag("XG"):
            is_ct_func = _is_read_ct_conversion_bismark
        
        # 统一链方向
        if is_ct_func(read):  # C→T conversion
            read.is_forward = True
        else:  # G→A conversion
            read.is_forward = False
```

### 1.3 标签检测函数

```python
def _is_read_ct_conversion_hisat3n(read):
    """Hisat-3n: YZ="+" 表示 C→T 转换"""
    return read.get_tag("YZ") == "+"

def _is_read_ct_conversion_bismark(read):
    """Bismark: XG="CT" 表示 C→T 转换"""
    return read.get_tag("XG") == "CT"
```

### 1.4 对后续分析的影响

统一链方向后，mpileup分析变得一致：
- **C位点**: 只统计正链reads
  - `.` 表示未甲基化 (C保持不变)
  - `T` 表示甲基化 (C→T转换)
- **G位点**: 只统计负链reads
  - `,` 表示未甲基化 (G保持不变)  
  - `a` 表示甲基化 (G→A转换)

## 2. UMI矫正机制（第二步处理）

### 2.1 背景

在单细胞甲基化测序中，UMI (Unique Molecular Identifier) 用于：
- 标识原始分子，区分PCR重复
- 通过多个reads的一致性检查来矫正测序错误
- 提高数据质量和定量准确性

**重要**：UMI矫正是在链方向统一之后进行的，这样可以确保mpileup数据的一致性。

### 2.2 输入数据格式

UMI矫正处理的是 samtools mpileup 的输出，格式如下：
```
chr1    46930353    C    5    .T.T,    IIHII    GAAGGTGTGTAT,GAAGGTGTGTAT,GAAGGTGTGTAC,GAAGGTGTGTAT,GAAGGTGTGTAT
```

各列含义：
- 列1: 染色体
- 列2: 位置
- 列3: 参考碱基
- 列4: 覆盖度
- 列5: 碱基序列 (. 表示与参考相同，字母表示不同)
- 列6: 质量分数
- 列7: UMI序列 (逗号分隔)

### 2.3 矫正流程

#### 步骤0: Indel处理 (预处理步骤)

在进行UMI矫正之前，ALLCools首先处理mpileup输出中的插入缺失(indel)信息：

**Indel标记格式**:
- `+N[bases]`: 表示插入，N为插入长度，后跟插入的碱基序列
- `-N[bases]`: 表示缺失，N为缺失长度，后跟缺失的碱基序列

**处理逻辑**:
```python
# 检测indel标记
incons_basecalls = read_bases.count("+") + read_bases.count("-")
if incons_basecalls > 0:
    # 逐字符解析，移除indel信息
    while index < len(read_bases):
        if read_bases[index] == "+" or read_bases[index] == "-":
            # 解析indel长度
            indel_size = ""
            while read_bases[ind].isdigit():
                indel_size += read_bases[ind]
                ind += 1
            # 跳过indel序列
            index = ind + int(indel_size)
```

**处理示例**:
```
原始: .+2AG.T-1C,
处理后: ..T,
```

**处理原因**:
- Indel信息会干扰碱基计数
- 甲基化分析关注的是碱基替换而非indel
- 确保后续UMI矫正的准确性

#### 步骤1: 过滤相关碱基
- **C位点**: 只保留 `.` (未甲基化) 和 `T` (C→T转换) 的reads
- **G位点**: 只保留 `,` (未甲基化) 和 `a` (G→A转换) 的reads

#### 步骤2: UMI合并 (基于汉明距离)
```python
def hamming_distance(str1, str2):
    """计算两个字符串的汉明距离"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))
```

合并规则：
- 汉明距离 ≤ 1 的UMI被认为是同一原始分子
- 合并时保留数量多的UMI
- 数量相同时保留两个UMI不合并

#### 步骤3: UMI内碱基矫正
对每个UMI内的多个reads进行碱基矫正：

1. **按数量优先**：选择出现次数最多的碱基
2. **质量分数决胜**：数量相同时，选择平均质量分数最高的碱基
3. **统一矫正**：将该UMI的所有reads统一为矫正后的碱基

```python
# 质量分数转换
def phred33_to_quality(ascii_char):
    return ord(ascii_char) - 33

# 碱基矫正逻辑
if len(tied_bases) == 1:
    corrected_base = most_common_base
else:
    # 计算平均质量分数
    for base, qualities in base_quality_map.items():
        avg_quality = sum(qualities) / len(qualities)
        if avg_quality > best_avg_quality:
            best_base = base
```

#### 步骤4: 输出结果
- 每个UMI只保留一个代表性read
- 覆盖度等于unique UMI的数量
- 所有碱基和质量分数都是矫正后的结果

### 2.4 矫正示例

**输入**:
```
chr1  46930353  C  5  .T.T,  IIHII  GAAGGTGTGTAT,GAAGGTGTGTAT,GAAGGTGTGTAC,GAAGGTGTGTAT,GAAGGTGTGTAT
```

**处理过程**:
0. Indel处理: 移除插入缺失标记 (如有)
1. UMI计数: GAAGGTGTGTAT(4), GAAGGTGTGTAC(1)
2. 汉明距离: GAAGGTGTGTAT vs GAAGGTGTGTAC = 1
3. 合并UMI: 保留GAAGGTGTGTAT (数量更多)
4. 碱基矫正: . (3次) > T (2次)，选择 `.`
5. 统一矫正: 所有碱基改为 `.`

**输出**:
```
chr1  46930353  C  1  .  I  GAAGGTGTGTAT
```

## 3. 实际应用

### 3.1 使用场景

1. **单细胞甲基化分析**: 需要链方向统一和UMI去重矫正
2. **多比对工具兼容**: 统一Bismark和Hisat-3n的输出格式
3. **数据质量控制**: 通过UMI矫正提高数据准确性

### 3.2 处理流程顺序

ALLCools的完整处理流程：
1. **第一步：链方向转换** - 使用 `convert_bam_strandness=True`
2. **第二步：UMI矫正** - 使用 `tag="UMI"`

### 3.3 参数配置

在 `bam_to_allc` 函数中：
```python
bam_to_allc(
    bam_path="input.bam",
    reference_fasta="ref.fa", 
    convert_bam_strandness=True,  # 第一步：启用链方向修改
    tag="UMI"  # 第二步：启用UMI矫正
)
```

### 3.4 输出文件

- `output.allc.tsv.gz`: 主要的ALLC格式文件
- `output_mpl_old.txt`: 原始mpileup输出
- `output_mpl_correction.txt`: UMI矫正后的mpileup输出
- `output.allc.tsv.gz.count.csv`: 统计信息

## 4. 优势和局限性

### 4.1 优势

**链方向修改**:
- 兼容多种比对工具
- 统一分析流程
- 自动检测和处理
- 保持数据完整性
- 为后续UMI矫正提供一致的数据基础

**UMI矫正**:
- 有效去除PCR重复
- 矫正测序错误
- 提高定量准确性
- 支持质量分数加权
- 基于统一链方向的可靠矫正

### 4.2 局限性

**链方向修改**:
- 依赖特定标签 (XG/YZ)
- 不支持其他比对工具
- 需要额外的处理步骤

**UMI矫正**:
- 需要UMI信息
- 计算复杂度较高
- 可能过度矫正低覆盖区域
- 依赖链方向统一的前置步骤

## 5. 最佳实践

1. **数据预处理**: 确保BAM文件包含正确的转换标签
2. **处理顺序**: 严格按照链方向转换→UMI矫正的顺序进行
3. **参数调优**: 根据数据特点调整质量阈值和覆盖度要求
4. **质量控制**: 检查每个步骤前后的统计信息
5. **验证结果**: 比较不同处理方式的结果一致性

## 6. 总结

ALLCools的链方向修改和UMI矫正机制为单细胞甲基化数据分析提供了强大的数据预处理能力。**关键在于正确的处理顺序**：首先进行链方向统一，确保mpileup数据的一致性，然后进行UMI矫正，通过智能的UMI合并和碱基矫正显著提高数据质量和分析的准确性。这种两步处理策略充分考虑了实际应用中的各种情况，为研究人员提供了可靠的分析工具。