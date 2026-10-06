# BioSyn Primer Batch QC

> **批量引物质量控制工具** — 一次性检查数百上千条引物的 Tm、GC%、发夹结构、二聚体风险，全部本地离线运行。
>
> A local-first batch primer QC tool: Tm, GC%, homopolymer, hairpin, and self-complementarity checks for thousands of primers at once.

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Offline](https://img.shields.io/badge/Privacy-Offline-orange.svg)](#privacy)

## 核心功能

| 功能 | 说明 |
|------|------|
| **批量引物导入** | CSV / XLSX / TSV，最低字段 `Primer_ID` + `Sequence`，支持列名映射 |
| **Tm 计算** | 支持 Basic Wallace 近似与 Nearest-Neighbor (SantaLucia 1998) 两种算法 |
| **GC% 分析** | 自动计算 GC 含量，标记过高/过低 |
| **均聚物检测** | 自动检测 AAAA / GGGG / CCCC / TTTT 长串 |
| **发夹风险评估** | Low / Review / High 三级风险分级 |
| **自互补性检测** | none / possible / strong，提示引物二聚体风险 |
| **可配置 QC 规则** | 长度、GC%、Tm 范围、均聚物阈值全部可自定义 |
| **四级 QC 状态** | PASS / REVIEW / WARNING / INVALID |
| **多格式导出** | XLSX（多 Sheet）、CSV、HTML 报告 |
| **容错处理** | 单条引物出错不影响整体批处理（10,000+ 条可用） |

## 解决的痛点

设计了几十条引物后，你需要逐条检查 Tm 是否合适、有没有发夹结构、会不会形成二聚体。网页工具一次只能查一条。BioSyn Primer Batch QC 把整个 QC 流程批量化：导入表格 → 一键分析 → 导出风险报告。

> 在线单条引物快速计算请访问 [BioSyn 引物设计工具](https://biosyn-inc.com/tools/primer-design-qc) — 本工具是其桌面端批量 QC 增强版。

## 技术栈

- **Python 3.10+**
- **openpyxl** — XLSX 读取与多 Sheet 报告导出
- Tm、发夹、二聚体算法由标准库实现

## 依赖与安装

### 环境要求

- **Python 3.10+**

### 依赖包

| 包 | 用途 | 是否必需 |
|------|------|----------|
| `openpyxl` | XLSX 文件读取 + XLSX 报告导出 | 必需（处理 CSV 时可降级） |

> 核心 QC 引擎（Tm 计算、GC%、均聚物、发夹、自互补检测）**仅依赖 Python 标准库**。

### 安装指令

```bash
# 1. 克隆仓库
git clone https://github.com/your-org/biosyn-primer-batch-qc.git
cd biosyn-primer-batch-qc

# 2. 安装依赖
pip install -r requirements.txt

# 3. 运行
python -m src.main examples/sample_primers.csv --output results.xlsx

# 使用最近邻 Tm 算法 + 自定义盐浓度
python -m src.main primers.xlsx --tm-method nearest_neighbor --salt 50 --mg 2.0

# 自定义 QC 范围
python -m src.main primers.csv --gc-min 45 --gc-max 55 --tm-min 58 --tm-max 62
```

## 输入格式

CSV / XLSX / TSV 表格：

| Primer_ID | Sequence | (可选: Name, Target, Gene, Direction, Notes) |
|-----------|----------|-----------------------------------------------|
| P001      | ATGCGT... | ...                                           |

## 输出

| 文件 | 内容 |
|------|------|
| `results.xlsx` | Sheets: Summary, Primer Results, QC Rules, Methods, Errors |
| `results.csv` | 每条引物一行：ID、长度、GC%、Tm、QC 状态、Warnings |
| `report.html` | 可浏览器打开的 QC 报告 |

## Tm 算法说明

| 方法 | 公式 | 适用场景 |
|------|------|----------|
| Basic | `4×(G+C) + 2×(A+T)`，盐校正 | 快速估算 |
| Nearest-Neighbor | SantaLucia 1998 统一 NN 参数 | 更精确的热力学估计 |

> Tm 计算为 **估算值**，实际 Tm 受 PCR 缓冲液、dNTP 浓度等多种因素影响。

## 科学边界

本工具提供 **计算性 QC / 决策支持**，不保证：
- PCR 一定成功
- 无二聚体形成
- 无发夹结构
- 实验特异性

> 由 [BioSyn Inc.](https://biosyn-inc.com/) 开发。更多 [引物设计与 PCR QC 工具](https://biosyn-inc.com/tools/primer-design-qc) 请访问官网。

## Keywords

`primer batch QC`, `primer design tool`, `Tm calculator`, `hairpin prediction`, `primer dimer check`, `PCR primer QC`, `oligo analysis`, `batch primer validation`

## License

[MIT](LICENSE)
