# Fomo.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.9%20|%201.10%20|%201.11-blue)](https://julialang.org/)

[中文文档](README_zh.md) | [English](README.md)

**Fomo** - **Fo**rward **Mo**deling：高性能二维各向同性弹性波数值模拟器。

```ibm_method=:mirror```  还不稳定，请勿使用。

## ✨ 特性

- 🚀 **后端调度架构** - 一套代码，CPU/GPU 自动切换
- 📐 **高阶交错网格有限差分** - 支持 2 至 10 阶空间精度
- 🛡️ **混合吸收边界 (HABC)** - 有效抑制边界反射
- 🌊 **自由地表建模** - 准确模拟 Rayleigh 面波
- 🏔️ **不规则地形 (IBM)** - 浸入边界法处理复杂地表
- ⚡ **多 GPU 并行** - 自动负载均衡
- 📁 **多格式支持** - SEG-Y、Binary、MAT、NPY、HDF5、JLD2
- 🎬 **视频录制** - 实时波场可视化

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

### 不规则自由表面

```julia
using Fomo

model = VelocityModel(vp, vs, rho, dx, dx)

# 使用辅助函数定义地表形状
z_surface = sinusoidal_surface(nx, dx; base_depth=50, amplitude=30, wavelength=1000)

# 或组合多种形状
z_surface = combine_surfaces(
    sinusoidal_surface(nx, dx; amplitude=20),
    gaussian_valley(nx, dx; valley_depth=25, width=300)
)

# 运行模拟
result = simulate_irregular!(
    model,
    z_surface,                           # 你定义的地表形状
    2000.0f0,                            # 震源 x 坐标
    Float32.(100:20:3900);               # 检波器 x 坐标
    config = IrregularSurfaceConfig(
        nt = 3000,
        ibm_method = :direct_zero,       # 或 :mirror 更高精度
        src_depth = 30.0f0,              # 地表以下深度
        output_dir = "outputs_irregular"
    ),
    video_config = VideoConfig(fields=[:vz], skip=10)
)
```

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

# 演示 3：不规则地表 IBM
z_surface = combine_surfaces(
    sinusoidal_surface(nx, dx; base_depth=50, amplitude=25),
    gaussian_valley(nx, dx; valley_depth=20, width=250)
)
result3 = simulate_irregular!(model3, z_surface, src_x, rec_x;
    config = IrregularSurfaceConfig(nt=3000, ibm_method=:direct_zero),
    video_config = VideoConfig(fields=[:vz], skip=10))
```

运行：
```bash
julia --project=. examples/run_demo.jl
```

### `examples/run_regular_surface_video.jl` - 面波可视化

双层模型展示 P 波、S 波和 Rayleigh 面波：

```julia
using Fomo

# 双层模型
vp[1:160, :] .= 2500.0f0    # 上层
vp[161:end, :] .= 4000.0f0  # 下层

model = VelocityModel(vp, vs, rho, 5.0f0, 5.0f0)

video_cfg = VideoConfig(fields=[:vz], skip=5, fps=30, colormap=:seismic)

result = simulate!(model, nx*dx/2, 50.0f0, rec_x, rec_z;
    config = SimulationConfig(nt=4000, f0=20.0f0, output_dir="outputs_regular"),
    video_config = video_cfg)
```

视频中可观察到的波：
- **P 波**：~2500 m/s（最快）
- **S 波**：~1500 m/s
- **Rayleigh 面波**：~1380 m/s（沿地表传播）

运行：
```bash
julia --project=. examples/run_regular_surface_video.jl
```

### `examples/run_irregular_with_video.jl` - 不规则地形

演示自定义地表形状：

```julia
using Fomo

# 示例 1：使用辅助函数
z_surface = sinusoidal_surface(nx, dx; base_depth=50, amplitude=25, wavelength=1500)

# 示例 2：完全自定义形状
x = Float32.((0:nx-1) .* dx)
z_custom = Float32.(60.0 .+ 20.0 .* sin.(2π .* x ./ 1500.0) .+
                    10.0 .* sin.(2π .* x ./ 300.0))

# 添加峡谷
for i in 1:nx
    xi = (i-1) * dx
    if 1300 < xi < 1700
        z_custom[i] += 30.0f0 * (1 - abs(xi - 1500) / 200)
    end
end

result = simulate_irregular!(model, z_custom, src_x, rec_x;
    config = IrregularSurfaceConfig(nt=3000, src_depth=40.0f0),
    video_config = VideoConfig(fields=[:vz], skip=10))
```

运行：
```bash
julia --project=. examples/run_irregular_with_video.jl
```

## 🏔️ 浸入边界法 (IBM)

IBM 无需细化网格即可准确模拟复杂地形：

| 方法 | 网格规模 | 时间步数 | 总成本 |
|------|----------|----------|--------|
| 细网格 + 阶梯逼近 | 4N | 2T | **8×** |
| **IBM（本方法）** | N | T | **~1.08×** |

两种 IBM 方法可选：
- `:direct_zero` - 稳定，推荐大多数情况使用
- `:mirror` - 更高精度，可能需要更小时间步

## 📂 项目结构

```
Fomo.jl/
├── src/
│   ├── Fomo.jl                 # 主模块
│   ├── backends/               # CPU/CUDA 抽象层
│   ├── types/                  # 数据结构
│   ├── kernels/                # 有限差分核函数
│   │   ├── velocity.jl
│   │   ├── stress.jl
│   │   ├── boundary.jl
│   │   └── ibm.jl
│   ├── surface/                # 不规则地表
│   ├── simulation/             # 时间步进
│   │   ├── api.jl              # 高层 API
│   │   ├── time_stepper.jl
│   │   └── time_stepper_ibm.jl
│   ├── io/                     # 模型/数据读写
│   └── visualization/          # 绑图与视频
├── examples/
│   ├── run_demo.jl
│   ├── run_regular_surface_video.jl
│   └── run_irregular_with_video.jl
└── README.md
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

IrregularSurfaceConfig(
    nt = 3000,
    f0 = 15.0f0,
    ibm_method = :direct_zero,  # :direct_zero 或 :mirror
    ibm_iterations = 3,
    src_depth = 30.0f0,         # 震源位于地表以下深度
    rec_depth = 0.0f0,          # 检波器深度（0 = 地表）
    output_dir = "outputs"
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
| `simulate_irregular!(model, z_surface, src_x, rec_x; config, video_config)` | 不规则地表模拟 |

### 地形辅助函数

| 函数 | 描述 |
|------|------|
| `flat_surface`, `sinusoidal_surface`, `gaussian_valley`, `gaussian_hill` | 基础形状 |
| `tilted_surface`, `step_surface`, `random_surface` | 更多形状 |
| `combine_surfaces(s1, s2, ...; method=:add)` | 组合形状 |

## 📖 参考文献

1. Luo, Y., & Schuster, G. (1990). Parsimonious staggered grid finite-differencing of the wave equation. *Geophysical Research Letters*, 17(2), 155-158.

2. 任志明, 刘洋. (2014). 一阶弹性波方程数值模拟中的混合吸收边界条件. *地球物理学报*, 57(2), 595-606. doi:10.6038/cjg20140223

3. Li, X., Yao, G., Niu, F., Wu, D., & Liu, N. (2023). Waveform inversion of seismic first arrivals acquired on irregular surface. *Geophysics*, 88(3), R289-R302.

## 📄 许可证

MIT License - 详见 [LICENSE](LICENSE)

## 👤 作者

Wuheng - 2025
