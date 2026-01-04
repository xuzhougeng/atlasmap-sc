# SOMA-Tiles Preprocess

将 H5AD 单细胞数据转换为多分辨率 Zarr bins，为 SOMA-Tiles 可视化提供高性能数据支持。

## 功能特性

- **多分辨率分箱**：基于四叉树的空间分箱算法，支持 1-12 级缩放
- **智能基因选择**：自动选择高变异基因、用户指定的标记基因和高表达基因
- **表达量聚合**：为每个 bin 预计算平均表达量和最大表达量
- **TileDBSOMA 存储**：可选存储完整表达矩阵，支持查询任意基因
- **类别数据处理**：支持多个 obs 列的类别统计（如 cell_type、leiden 等）
- **高效压缩**：使用 Zarr 格式存储，支持 zstd 压缩
- **批量处理**：支持大规模数据集的内存高效处理

## 安装

```bash
cd preprocessing

# 使用 uv（推荐）
uv venv
source .venv/bin/activate
uv pip install -e .

# 如果出现硬链接警告（跨文件系统时），使用：
uv pip install -e . --link-mode=copy

# 重新安装（强制重装）
uv pip install -e . --reinstall

# 或使用 pip
pip install -e .
```

## 快速开始

### 命令行使用

基本用法

```bash
wget -O blood.h5ad https://datasets.cellxgene.cziscience.com/c3d1a5e6-780b-4fe9-a39b-1864f927e87b.h5ad
soma-preprocess run -i blood.h5ad -o ./data/blood
```

使用所有的表达基因

```bash
wget -O retina.h5ad https://datasets.cellxgene.cziscience.com/fc4b5c6c-26ae-4c65-9798-51f604b69ba6.h5ad
soma-preprocess run -i retina.h5ad -o ./data/retina --all-expressed --category cell_type
```

完整参数

```bash
soma-preprocess run \
    -i data.h5ad \
    -o ./output \
    --umap-key X_umap \
    --zoom-levels 8 \
    --n-genes 500 \
    --category cell_type \
    --category leiden \
    --chunk-size 256 \
    --batch-size 100000 \
    --name "My Dataset" \
    --verbose
```


### 使用配置文件

```bash
# 生成示例配置文件
soma-preprocess init-config -o config.yaml

# 使用配置文件运行
soma-preprocess from-config -c config.yaml
```

## 命令行参数

### `soma-preprocess run`

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--input` | `-i` | 必填 | 输入 H5AD 文件路径 |
| `--output` | `-o` | 必填 | 输出目录路径 |
| `--umap-key` | | `X_umap` | UMAP 坐标在 `adata.obsm` 中的键名 |
| `--zoom-levels` | `-z` | `8` | 缩放级别数量（1-12） |
| `--n-genes` | `-g` | `500` | 预聚合的基因数量 |
| `--all-expressed` | `-a` | | 使用所有表达的基因而非 top N 基因 |
| `--min-cells` | | `3` | 基因至少在多少个细胞中表达（仅与 `--all-expressed` 一起使用） |
| `--no-soma` | | | 禁用 TileDBSOMA 存储（默认启用） |
| `--category` | `-c` | | 要包含的类别列（可多次指定） |
| `--chunk-size` | | `256` | Zarr chunk 大小 |
| `--batch-size` | | `100000` | 处理批次大小 |
| `--name` | | | 数据集名称 |
| `--verbose` | `-v` | | 启用详细日志 |

### `soma-preprocess visualize`

验证预处理输出，生成静态可视化图片。

> 提示：这里的 `zoom` 表示“分箱分辨率”（越大越细）。以默认 `--zoom-levels 8` 为例：
> - `zoom=0`：1×1 个 bin，通常只会看到 1 个点/块（非常粗）
> - `zoom=7`：128×128 个 bin，能看到完整的 UMAP 结构（建议作为默认检查/展示层）

```bash
# 基本用法（默认会画 3,5,7 三个 zoom，便于对比；随机选择 3 个基因）
soma-preprocess visualize -i ./output/zarr -o ./figures

# 推荐：只看最高分辨率（最接近前端打开时的全局效果）
soma-preprocess visualize -i ./output/zarr -o ./figures -z 7

# 指定 zoom 级别和基因
soma-preprocess visualize -i ./output/zarr -o ./figures -z 3,5,7 -g CD3D -g CD8A

# 指定格式和分辨率
soma-preprocess visualize -i ./output/zarr -o ./figures -f svg --dpi 300
```

| 参数 | 简写 | 默认值 | 说明 |
|------|------|--------|------|
| `--input` | `-i` | 必填 | Zarr 输出目录路径 |
| `--output` | `-o` | 必填 | 图片输出目录 |
| `--zoom-levels` | `-z` | `3,5,7` | 要可视化的 zoom 级别（逗号分隔）；建议至少包含最高 zoom（如 zoom_levels=8 则为 7） |
| `--category` | `-c` | 全部 | 要可视化的类别列 |
| `--gene` | `-g` | 随机选择 | 要可视化的基因（可多次指定） |
| `--n-genes` | | `3` | 随机选择的基因数量 |
| `--format` | `-f` | `png` | 输出格式（png/svg/pdf） |
| `--dpi` | | `150` | 图片分辨率 |
| `--verbose` | `-v` | | 启用详细日志 |

输出结构：
```
figures/
├── category/
│   ├── cell_type_zoom_3.png
│   ├── cell_type_zoom_5.png
│   ├── cell_type_zoom_7.png
│   └── cell_type_multi_zoom.png
└── expression/
    ├── GENE1_zoom_3.png
    ├── GENE1_zoom_5.png
    ├── GENE1_zoom_7.png
    └── GENE1_multi_zoom.png
```

## 配置文件格式

```yaml
# 输入输出
input_path: data.h5ad
output_dir: ./preprocessed

# 坐标设置
umap_key: X_umap
coordinate_range: 256.0

# 分箱设置
zoom_levels: 8
tile_size: 256

# 基因选择
n_genes: 500
hvg_n_top: 300
marker_genes:
  - CD3D
  - CD8A
  - MS4A1
  - CD14
use_all_expressed: false  # 设为 true 使用所有表达的基因
min_cells_expressed: 3    # 基因至少在多少个细胞中表达

# 类别设置
category_columns:
  - cell_type
  - leiden
default_category: cell_type

# Zarr 设置
zarr_compressor: zstd
zarr_compression_level: 3
chunk_size: 256

# 内存设置
batch_size: 100000

# 元数据
dataset_name: My Dataset
dataset_description: Single-cell dataset description
```

## 输出结构

```
output/
├── zarr/
│   ├── bins.zarr/           # 多分辨率 bin 数据
│   │   ├── zoom_0/          # 缩放级别 0
│   │   ├── zoom_1/          # 缩放级别 1
│   │   └── ...
│   ├── metadata.json        # 数据集元数据
│   └── gene_index.json      # 基因索引映射
│
└── soma/                    # TileDBSOMA 存储（默认启用）
    └── experiment.soma/     # 完整单细胞数据
        ├── obs              # 细胞元数据（含归一化坐标和 bin 索引）
        ├── ms/RNA/var       # 基因元数据
        └── ms/RNA/X/data    # 完整表达矩阵（稀疏）
```

### metadata.json 结构

```json
{
  "dataset_name": "My Dataset",
  "n_cells": 100000,
  "n_genes_total": 20000,
  "n_genes_preaggregated": 500,
  "zoom_levels": 8,
  "tile_size": 256,
  "coordinate_range": 256.0,
  "preaggregated_genes": ["gene1", "gene2", ...],
  "categories": {
    "cell_type": {
      "values": ["T cell", "B cell", ...],
      "mapping": {"T cell": 0, "B cell": 1, ...}
    }
  },
  "bounds": {
    "min_x": 0.0,
    "max_x": 256.0,
    "min_y": 0.0,
    "max_y": 256.0
  },
  "soma_enabled": true,
  "soma_path": "soma/experiment.soma",
  "all_genes_queryable": true
}
```

## 处理流程

1. **加载数据**：读取 H5AD 文件
2. **坐标处理**：从 `obsm` 提取 UMAP 坐标并归一化到 `[0, 256)` 范围
3. **基因选择**：
   - 默认模式：选择 top N 基因（HVG + 标记基因 + 高表达基因）
   - `--all-expressed` 模式：使用所有在至少 `--min-cells`（默认 3）个细胞中表达的基因
4. **类别映射**：为指定的 obs 列构建类别到索引的映射
5. **多分辨率分箱**：
   - 对每个缩放级别，使用四叉树算法分配 bin
   - 聚合每个 bin 内的表达量（平均值、最大值）
   - 统计每个 bin 内的类别分布
6. **写入 Zarr**：将分箱数据写入压缩的 Zarr 格式
7. **生成元数据**：输出 `metadata.json` 和 `gene_index.json`

## 技术细节

### 四叉树分箱

采用 2 的幂次分箱大小，确保高效的四叉树细分：

| 缩放级别 | 每轴 bin 数 | bin 大小（坐标单位） |
|----------|-------------|----------------------|
| 0 | 1 | 256 |
| 1 | 2 | 128 |
| 2 | 4 | 64 |
| 3 | 8 | 32 |
| 4 | 16 | 16 |
| 5 | 32 | 8 |
| 6 | 64 | 4 |
| 7 | 128 | 2 |

### 表达量聚合

每个 bin 存储：
- `expression_mean`：bin 内所有细胞的平均表达量
- `expression_max`：bin 内最大表达量
- `cell_count`：bin 内细胞数量
- `cell_ids`：bin 内细胞索引列表
- `category_counts`：各类别的细胞计数

## 依赖项

- Python >= 3.9
- scanpy >= 1.9.0
- anndata >= 0.10.0
- zarr >= 2.16.0
- numpy >= 1.24.0
- numba >= 0.58.0
- scipy >= 1.11.0
- pandas >= 2.0.0
- click >= 8.1.0
- scikit-misc >= 0.1.4（HVG 计算需要）

## 许可证

MIT License
