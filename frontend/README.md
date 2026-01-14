# AtlasMap 前端

基于 TypeScript + Vite + Leaflet 的单细胞空间转录组可视化前端。

## 项目结构

```
frontend/
├── src/
│   ├── api/
│   │   └── client.ts           # API 客户端和类型定义
│   ├── components/
│   │   ├── GeneSelector.ts     # 基因搜索选择器
│   │   ├── CategoryFilter.ts   # 分类过滤器
│   │   ├── CategoryLegend.ts   # 分类图例
│   │   ├── CategoryColumnSelector.ts  # 分类列选择器
│   │   ├── TabPanel.ts         # 标签面板 (Category/Expression)
│   │   └── CellQueryPanel.ts   # Bin 查询面板
│   ├── map/
│   │   └── MapController.ts    # Leaflet 地图控制器
│   ├── state/
│   │   └── StateManager.ts     # 应用状态管理
│   ├── styles/
│   │   └── main.css
│   └── main.ts                 # 入口文件
├── index.html
├── package.json
├── tsconfig.json
└── vite.config.ts
```

## 数据类型

所有类型定义在 `src/api/client.ts` 中：

### 元数据

```typescript
interface Metadata {
    dataset_name: string;       // 数据集名称
    n_cells: number;            // 细胞总数
    n_genes_total: number;      // 基因总数
    n_genes_preaggregated: number;  // 预聚合基因数
    zoom_levels: number;        // 缩放层级数
    tile_size: number;          // 瓦片大小
    preaggregated_genes: string[];  // 预聚合基因列表
    categories: Record<string, CategoryInfo>;  // 分类信息
    bounds: {                   // 空间边界
        min_x: number;
        max_x: number;
        min_y: number;
        max_y: number;
    };
}

interface CategoryInfo {
    values: string[];                  // 分类值列表
    mapping: Record<string, number>;   // 值到索引的映射
}
```

### 基因相关

```typescript
interface GeneInfo {
    name: string;           // 基因名称
    index: number;          // 基因索引
    preaggregated: boolean; // 是否预聚合
}

interface GeneStats {
    gene: string;           // 基因名称
    index: number;          // 基因索引
    expressing_bins: number;// 表达 bin 数量
    total_bins: number;     // 总 bin 数量
    total_cells: number;    // 表达细胞数
    mean_expression: number;// 平均表达量
    max_expression: number; // 最大表达量
}
```

### Bin 数据

```typescript
interface BinExpressionInfo {
    bin_index: number;      // Bin 索引
    bin_x: number;          // X 坐标
    bin_y: number;          // Y 坐标
    cell_count: number;     // 细胞数量
    expression: number;     // 表达量
    categories: Record<string, string>;  // 分类信息
}

interface BinQueryResult {
    bins: BinExpressionInfo[];  // Bin 列表
    total_count: number;        // 总数
    offset: number;             // 偏移量
    limit: number;              // 每页数量
    gene: string;               // 基因名称
    threshold: number;          // 阈值
}

interface BinQueryParams {
    threshold?: number;     // 表达阈值
    offset?: number;        // 偏移量
    limit?: number;         // 限制数量
}
```

### 分类图例

```typescript
interface CategoryLegendItem {
    value: string;      // 分类值
    color: string;      // 颜色 (十六进制)
    index: number;      // 索引
}
```

### 应用状态

定义在 `src/state/StateManager.ts`：

```typescript
interface AppState {
    // 视图状态
    zoom: number;
    center: [number, number];

    // 过滤状态
    selectedCategories: string[];
    expressionFilter: {
        gene: string | null;
        min: number;
        max: number;
    };

    // 显示状态
    colorMode: 'category' | 'expression';
    colorGene: string | null;
    colorScale: string;
}
```

## API 端点

所有 API 和 Tiles 端点均需在路径前加 `/d/{dataset}`。可通过 `GET /api/datasets` 获取数据集列表。

### 数据集（全局）

| 方法 | 端点 | 返回类型 | 描述 |
|------|------|----------|------|
| GET | `/api/datasets` | `DatasetsResponse` | 获取可用数据集列表 |

### 元数据与统计

| 方法 | 端点 | 返回类型 | 描述 |
|------|------|----------|------|
| GET | `/d/{dataset}/api/metadata` | `Metadata` | 获取数据集元信息 |
| GET | `/d/{dataset}/api/stats` | `Record<string, any>` | 获取统计信息 |

### 基因相关

| 方法 | 端点 | 返回类型 | 描述 |
|------|------|----------|------|
| GET | `/d/{dataset}/api/genes` | `{ genes: string[], total: number }` | 获取所有基因列表 |
| GET | `/d/{dataset}/api/genes/{gene}` | `GeneInfo` | 获取基因详情 |
| GET | `/d/{dataset}/api/genes/{gene}/stats` | `GeneStats` | 获取基因表达统计 |
| GET | `/d/{dataset}/api/genes/{gene}/bins?threshold=X&offset=Y&limit=Z` | `BinQueryResult` | 查询表达 Bin |

### 分类相关

| 方法 | 端点 | 返回类型 | 描述 |
|------|------|----------|------|
| GET | `/d/{dataset}/api/categories` | `Record<string, CategoryInfo>` | 获取所有分类 |
| GET | `/d/{dataset}/api/categories/{column}/colors` | `Record<string, string>` | 获取颜色映射 |
| GET | `/d/{dataset}/api/categories/{column}/legend` | `CategoryLegendItem[]` | 获取图例数据 |

### 瓦片端点

| 端点 | 描述 |
|------|------|
| `/d/{dataset}/tiles/{z}/{x}/{y}.png` | 基础瓦片 |
| `/d/{dataset}/tiles/{z}/{x}/{y}/category/{column}.png` | 分类着色瓦片 |
| `/d/{dataset}/tiles/{z}/{x}/{y}/expression/{gene}.png?colormap={scale}` | 表达着色瓦片 |

## URL 参数

### 基因表达直接跳转

主页面支持通过 URL 参数 `gene` 直接跳转到指定基因的表达可视化视图。

**用法：**

```
/?dataset={dataset_id}&gene={gene_name}
```

**示例：**

```
http://localhost:3000/?dataset=retina&gene=AT5G43870
```

**行为：**
- 自动切换到 Expression 标签页
- 加载指定基因的表达瓦片
- 在基因搜索框中显示当前基因名
- 更新 Cell Query Panel 中的基因信息
- 跳过默认的 Category 初始化

**应用场景：**
- 从差异表达分析结果页面（DE Results）点击基因名直接查看其空间表达模式
- 书签保存特定基因的可视化链接
- 在报告或文档中引用特定基因的可视化

## 组件说明

components
- GeneSelector：基因搜索和选择组件。加载完整基因列表后在客户端进行过滤，支持防抖输入，返回前 20 个匹配结果。
- CategoryColumnSelector：分类列下拉选择器，用于切换当前显示的分类列。
- CategoryFilter：分类值过滤器，为每个分类列渲染复选框组，支持多选过滤。
- CategoryLegend：分类图例组件，显示当前分类列的颜色图例，支持折叠/展开。
- CellQueryPanel：Bin 查询面板，包含表达阈值滑块、查询按钮、统计信息显示和分页结果表格。
- TabPanel：标签面板，在 Category 和 Expression 两个标签页之间切换。
- MapController：Leaflet 地图控制器，使用 CRS.Simple 坐标系统显示 UMAP 可视化，管理基础瓦片层、分类瓦片层和表达瓦片层。
- StateManager：应用状态管理器，使用发布-订阅模式，支持按键订阅和全局订阅。

## 开发指南

### 安装依赖

```bash
npm install
```

### 启动开发服务器

```bash
npm run dev
```

开发服务器运行在 http://localhost:3000，自动代理 `/api` 和 `/d` 请求到后端 (http://localhost:8080)。

### 构建生产版本

```bash
npm run build
```

### 类型检查

```bash
npm run type-check
```

### 代码检查

```bash
npm run lint
```
