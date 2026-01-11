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
- 🏔️ **不规则地形 (IBM)** - 浸入边界法处理复杂地表
- ⚡ **多 GPU 并行** - 自动负载均衡，榨干显卡性能
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

### 可选依赖

读取不同格式的模型文件：

```julia
using Pkg
Pkg.add("SegyIO")  # SEG-Y 文件
Pkg.add("MAT")     # MATLAB 文件  
Pkg.add("NPZ")     # NumPy 文件
Pkg.add("HDF5")    # HDF5 文件
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

# 不录制视频运行模拟
result = simulate!(
    model,
    2000.0f0, 50.0f0,                    # 震源 (x, z)
    Float32.(100:20:3900),               # 检波器 x 坐标
    fill(10.0f0, 190);                   # 检波器 z 坐标
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs")
)

# 录制视频运行模拟 - VideoConfig 是独立参数
result = simulate!(
    model,
    2000.0f0, 50.0f0,
    Float32.(100:20:3900),
    fill(10.0f0, 190);
    config = SimulationConfig(nt=3000, f0=15.0f0, output_dir="outputs"),
    video_config = VideoConfig(fields=[:vz], skip=5, fps=30)
)

# 自动保存结果：
# - outputs/gather.bin
# - outputs/gather.png  
# - outputs/wavefield_vz.mp4（如果提供了 video_config）
```

### 不规则自由表面

```julia
using Fomo

# 创建模型
model = VelocityModel(vp, vs, rho, dx, dx)

# 使用辅助函数定义地表形状
z_surface = sinusoidal_surface(nx, dx; 
    base_depth=50, amplitude=30, wavelength=1000)

# 或组合多种形状
z_surface = combine_surfaces(
    sinusoidal_surface(nx, dx; amplitude=20),
    gaussian_valley(nx, dx; valley_depth=25, width=300)
)

# 或完全自定义
x = Float32.((0:nx-1) .* dx)
z_surface = Float32.(50.0 .+ 20.0 .* sin.(2π .* x ./ 1000.0))

# 运行模拟 - video_config 是独立参数
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

## 🔧 底层 API

对于需要完全控制的高级用户：

```julia
using Fomo

# 后端选择
be = is_cuda_available() ? backend(:cuda) : backend(:cpu)

# 初始化组件
medium = init_medium(model, nbc, fd_order, be; free_surface=true)
habc = init_habc(medium.nx, medium.nz, nbc, dt, dx, dz, vp_max, be)
fd_coeffs = to_device(get_fd_coefficients(fd_order), be)
wavefield = Wavefield(medium.nx, medium.nz, be)
params = SimParams(dt, nt, dx, dz, fd_order)

# 设置震源和检波器
src = Source(src_i, src_j, to_device(wavelet, be))
rec = Receivers(to_device(rec_i, be), to_device(rec_j, be), 
                to_device(zeros(Float32, nt, n_rec), be), :vz)

# 带回调函数运行时间循环
run_time_loop!(be, wavefield, medium, habc, fd_coeffs, src, rec, params;
    progress = true,
    on_step = (W, info) -> begin
        # 自定义每步操作
        return true
    end
)
```

## ⚡ 多 GPU 并行

```julia
using Fomo

model = load_model("marmousi.jld2")

# 定义观测系统
src_x = Float32.(100:200:16900)
src_z = fill(10.0f0, length(src_x))
rec_x = Float32.(0:15:17000)
rec_z = fill(20.0f0, length(rec_x))

wavelet = ricker_wavelet(25.0f0, dt, nt)
params = SimParams(dt, nt, model.dx, model.dz, 8)

# 自动使用所有可用 GPU！
results = run_shots_auto!(
    model, rec_x, rec_z, src_x, src_z, wavelet, params;
    nbc=50, fd_order=8, output_dir="outputs/"
)
```

## 🎬 视频录制

```julia
using Fomo

# 配置视频
config = VideoConfig(
    fields = [:vz],     # 录制 vz 分量（可选 :vx, :vz, :vel, :p）
    skip = 10,          # 每 10 步保存一帧
    downsample = 1      # 空间降采样
)

recorder = MultiFieldRecorder(medium.nx, medium.nz, dt, config)

# 带录制回调运行
run_time_loop!(be, wavefield, medium, habc, fd_coeffs, src, rec, params;
    on_step = (W, info) -> begin
        Fomo.record!(recorder.recorder, W, info.k, dt)
        return true
    end
)

# 生成 MP4
generate_video(recorder.recorder, "wavefield.mp4"; fps=30, colormap=:seismic)
```

## 📁 模型读写

```julia
using Fomo

# 从 JLD2 加载（推荐）
model = load_model("marmousi.jld2")

# 从分离的 SEG-Y 文件加载（需要 SegyIO）
using SegyIO
model = load_model_files(
    vp = "vp.segy",
    vs = "vs.segy", 
    rho = "rho.segy",
    dx = 12.5
)

# 查看模型信息
model_info(model)

# 重采样模型
model_sim = resample_model(model, 10.0, 10.0; method=:linear)

# 保存为 JLD2
save_model("model.jld2", model)
```

## 📂 项目结构

```
Fomo.jl/
├── src/
│   ├── Fomo.jl                 # 主模块
│   ├── backends/               # CPU/CUDA 抽象层
│   │   └── backend.jl
│   ├── types/                  # 数据结构
│   │   ├── structures.jl
│   │   └── model.jl
│   ├── kernels/                # 有限差分核函数
│   │   ├── velocity.jl
│   │   ├── stress.jl
│   │   ├── boundary.jl
│   │   ├── source_receiver.jl
│   │   └── ibm.jl
│   ├── surface/                # 不规则地表
│   │   └── irregular.jl
│   ├── simulation/             # 时间步进与炮管理
│   │   ├── init.jl
│   │   ├── time_stepper.jl
│   │   ├── time_stepper_ibm.jl
│   │   ├── shots.jl
│   │   ├── parallel.jl
│   │   └── api.jl              # 高层 API
│   ├── io/                     # 模型/数据读写
│   │   ├── model_io.jl
│   │   ├── gather_io.jl
│   │   └── geometry_io.jl
│   └── visualization/          # 绑图与视频
│       ├── video.jl
│       └── plots.jl
├── examples/
│   ├── run_demo.jl                    # 综合演示
│   ├── run_regular_surface_video.jl  # 面波可视化
│   └── run_irregular_with_video.jl   # 不规则地形示例
├── test/
└── docs/
```

## 🏔️ IBM：不规则自由表面

浸入边界法 (IBM) 无需细化网格即可准确模拟复杂地形：

| 方法 | 网格规模 | 时间步数 | 总成本 |
|------|----------|----------|--------|
| 细网格 + 阶梯逼近 | 4N | 2T | **8×** |
| **IBM（本方法）** | N | T | **1.08×** |

两种 IBM 方法可选：
- `:direct_zero` - 稳定，推荐大多数情况使用
- `:mirror` - 更高精度，可能需要更小时间步

## 🚀 性能优化

| 优化措施 | 加速比 | 描述 |
|----------|--------|------|
| 预计算浮力 (1/ρ) | 15-25% | 消除速度更新中的除法 |
| 预计算 λ+2μ | 5-10% | 减少应力更新计算量 |
| 展开 FD 模板 | 10-15% | 更好的 SIMD 向量化 |
| 优化 GPU 块 (32×8) | 10-20% | 更好的内存合并访问 |

## 🧪 运行测试

```bash
cd Fomo.jl
julia --project=. -e "using Pkg; Pkg.test()"
```

## 📚 API 参考

### 配置结构体
- `SimulationConfig` - 规则地表模拟配置
- `IrregularSurfaceConfig` - 不规则地表模拟配置
- `SimulationResult` - 结果容器，包含数据和文件路径

### 高层函数
- `simulate!()` - 运行规则地表模拟
- `simulate_irregular!()` - 运行不规则地表模拟

### 地形辅助函数
- `flat_surface()`, `sinusoidal_surface()`, `gaussian_valley()` 等
- `combine_surfaces()` - 组合多种地形

### 核心类型
- `VelocityModel` - 速度模型容器
- `Medium` - 计算网格与材料属性
- `Wavefield` - 波场数组 (vx, vz, txx, tzz, txz)
- `SimParams` - 模拟参数

### 底层函数
- `init_medium()` - 初始化计算介质
- `init_habc()` - 初始化吸收边界
- `run_time_loop!()` - 带回调运行时间循环
- `run_shots!()` - 顺序执行多炮
- `run_shots_auto!()` - 自动多 GPU 并行

## 📖 参考文献

1. Luo, Y., & Schuster, G. (1990). *Parsimonious staggered grid finite-differencing of the wave equation*. Geophysical Research Letters.

2. Liu, Y., & Sen, M. K. (2012). *A hybrid absorbing boundary condition for elastic staggered-grid modelling*. Geophysical Prospecting.

## 📄 许可证

MIT License - 详见 [LICENSE](LICENSE)

## 👤 作者

zswh - 2025
