# Fomo.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.9%20|%201.10%20|%201.11-blue)](https://julialang.org/)

[中文文档](README_zh.md) | [English](README.md)

**Fomo** - **Fo**rward **Mo**deling：高性能二维各向同性弹性波数值模拟器。

## ✨ 特性

- 🚀 **后端调度架构** - 一套代码，CPU/GPU 自动切换
- 📐 **高阶交错网格有限差分** - 支持 2 至 10 阶空间精度
- 🛡️ **混合吸收边界 (HABC)** - 有效抑制边界反射
- 🌊 **自由地表建模** - 准确模拟 Rayleigh 面波
- ⚡ **多 GPU 并行** - 自动负载均衡
- 📁 **多格式支持** - SEG-Y、Binary、MAT、NPY、HDF5、JLD2
- 🎬 **视频录制** - 实时波场可视化，支持多属性波场

## 📋 系统要求

- **Julia 1.9、1.10 或 1.11**（暂不支持 1.12，CairoMakie 兼容性问题）
- CUDA 显卡（可选，用于 GPU 加速）

## 🔧 安装

### 从 GitHub 安装

```julia
using Pkg
Pkg.add(url="https://github.com/Wuheng10086/Fomo.jl")
```

### 本地开发

```bash
git clone https://github.com/Wuheng10086/Fomo.jl.git
cd Fomo.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## 🚀 快速开始

### 高层 API（推荐）

```julia
using Fomo

# 创建速度模型
nx, nz = 400, 200
dx = 10.0f0

vp = fill(3000.0f0, nz, nx)
vs = fill(1800.0f0, nz, nx)
rho = fill(2200.0f0, nz, nx)

vp[100:end, :] .= 4000.0f0
vs[100:end, :] .= 2400.0f0

model = VelocityModel(vp, vs, rho, dx, dx; name="双层模型")

# 运行模拟（不录制视频）
result = simulate!(
    model,
    2000.0f0, 50.0f0,                    # 震源位置 (x, z)，单位：米
    Float32.(100:20:3900),               # 检波器 x 坐标
    fill(10.0f0, 190);                   # 检波器 z 坐标
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs")
)

# 运行模拟（录制视频）- VideoConfig 作为独立参数
result = simulate!(
    model,
    2000.0f0, 50.0f0,
    Float32.(100:20:3900),
    fill(10.0f0, 190);
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs"),
    video_config = VideoConfig(fields=[:vz], skip=5, fps=30)
)
```

## 📂 示例

### `examples/run_demo.jl` - 综合演示

演示三种场景：

```julia
using Fomo

# 演示 1：快速测试（均匀模型，无视频）
result1 = simulate!(model, src_x, src_z, rec_x, rec_z;
    config = SimulationConfig(nt=1000, output_dir="outputs/demo1"))

# 演示 2：面波可视化（带视频）
result2 = simulate!(model2, src_x, src_z, rec_x, rec_z;
    config = SimulationConfig(nt=4000, f0=20.0f0, output_dir="outputs/demo2"),
    video_config = VideoConfig(fields=[:vz], skip=5, fps=30))
```

运行：
```bash
julia --project=. examples/run_demo.jl
```

## 📂 项目结构

```
Fomo.jl/
├── src/
│   ├── Fomo.jl                 # 主模块
│   ├── backends/               # CPU/CUDA 抽象层
│   │   └── backend.jl
│   ├── types/                  # 数据结构
│   │   ├── model.jl            # 速度模型定义
│   │   └── structures.jl       # 模拟结构体
│   ├── kernels/                # 有限差分核函数
│   │   ├── velocity.jl         # 速度更新核函数
│   │   ├── stress.jl           # 应力更新核函数
│   │   ├── boundary.jl         # 边界条件核函数
│   │   ├── source_receiver.jl  # 震源和检波器核函数
│   │   └── vacuum.jl           # 真空公式核函数
│   ├── simulation/             # 模拟逻辑
│   │   ├── api.jl              # 高层 API
│   │   ├── simple_api.jl       # 简化 API
│   │   ├── init.jl             # 初始化例程
│   │   ├── init_vacuum.jl      # 真空初始化
│   │   ├── time_stepper.jl     # 时间步进例程
│   │   ├── shots.jl            # 炮点处理
│   │   ├── batch.jl            # 批处理
│   │   └── parallel.jl         # 并行处理
│   ├── io/                     # 输入/输出操作
│   │   ├── model_io.jl         # 模型 I/O
│   │   ├── gather_io.jl        # 道集数据 I/O
│   │   └── geometry_io.jl      # 几何 I/O
│   └── visualization/          # 绘图与视频
│       ├── plots.jl            # 静态绘图
│       └── video.jl            # 视频生成
├── examples/
│   ├── run_demo.jl
│   ├── run_irregular_with_video.jl
│   └── run_vacuum_topography.jl
├── README.md
├── README_zh.md
└── Project.toml
```

## 📚 API 参考

### 配置结构体

```julia
SimulationConfig(
    nt = 3000,              # 时间步数
    f0 = 15.0f0,            # 震源频率 (Hz)
    nbc = 50,               # 吸收边界层数
    fd_order = 8,           # 有限差分阶数
    free_surface = true,    # 启用自由地表
    output_dir = "outputs"  # 输出目录
)

VideoConfig(
    fields = [:vz],         # 录制的场 (:vx, :vz, :vel, :p)
    skip = 10,              # 每 N 步录制一帧
    fps = 30,               # 视频帧率
    colormap = :seismic     # 色图
)
```

### 高层函数

| 函数 | 描述 |
|------|------|
| `simulate!(model, src_x, src_z, rec_x, rec_z; config, video_config)` | 规则地表模拟 |

### 仿真控制

| 配置 | 描述 |
|------|------|
| `SimulationConfig(free_surface=true)` | 启用面波（瑞利波） |
| `SimulationConfig(free_surface=false)` | 禁用面波 |
| `VideoConfig(fields=[:vz, :vx, :vel])` | 多属性波场视频录制 |

### 模型创建函数

| 函数 | 描述 |
|------|------|
| `create_homogeneous_model(vp, vs, rho, (nz, nx), dx)` | 创建均匀模型 |
| `create_layered_model(layers, dx)` | 创建分层模型 |
| `create_gradient_model(vp_func, vs_func, rho_func, (nz, nx), dx)` | 创建梯度模型 |

### 地形辅助函数

| 函数 | 描述 |
|------|------|
| `flat_surface(nx, dx, depth)` | 平坦地表 |
| `sinusoidal_surface(nx, dx; amplitude, wavelength)` | 正弦地表 |
| `gaussian_valley(nx, dx; valley_depth, width)` | 高斯谷地 |
| `gaussian_hill(nx, dx; hill_height, width)` | 高斯山丘 |
| `tilted_surface(nx, dx; depth_left, depth_right)` | 倾斜地表 |
| `step_surface(nx, dx; depth_left, depth_right)` | 阶梯/悬崖 |
| `random_surface(nx, dx; amplitude, smoothness)` | 随机粗糙地表 |
| `combine_surfaces(s1, s2, ...)` | 组合多种形状 |

## 📖 参考文献

1. Luo, Y., & Schuster, G. (1990). Parsimonious staggered grid finite-differencing of the wave equation. *Geophysical Research Letters*, 17(2), 155-158.

2. 任志明, 刘洋. (2014). 一阶弹性波方程数值模拟中的混合吸收边界条件. *地球物理学报*, 57(2), 595-606. doi:10.6038/cjg20140223

3. Li, X., Yao, G., Niu, F., Wu, D., & Liu, N. (2023). Waveform inversion of seismic first arrivals acquired on irregular surface. *Geophysics*, 88(3), R289-R302.

## 📄 许可证

MIT License - 详见 [LICENSE](LICENSE)

## 👤 作者

Wuheng - 2025