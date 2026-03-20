# MATLAB 颗粒–气泡液膜排液研究框架（起步版）

这个目录是一个**从零开始、便于后续继续扩展**的 MATLAB 研究框架。

## 当前结构

- `pb_run.m`：主入口
- `pb_default_parameters.m`：参数初始化
- `pb_film_thickness.m`：液膜厚度模型
- `pb_marker_density.m`：拉格朗日点投影密度
- `pb_refinement_indicator.m`：局部加密判据
- `pb_build_adaptive_grid.m`：非均匀径向网格生成
- `pb_solve_reynolds_pressure.m`：轴对称 Reynolds 方程求解
- `pb_taylor_lubrication_force.m`：Taylor 理论参考阻力
- `pb_disjoining_force.m`：近接触弱排斥正则项
- `pb_run_validation_case.m`：验证算例
- `pb_run_dynamic_case.m`：准静态接近过程
- `pb_plot_validation_results.m` / `pb_plot_dynamic_results.m`：后处理绘图

## 设计目的

相比旧版单文件原型，这一版首先解决的是**结构可扩展性**：

1. 先把“几何、判据、网格、求解器、动力学、绘图”拆开；
2. 让后续加入更复杂物理（变形、表面张力、DLVO、真实 IBM/LBM 外流场）更容易；
3. 让你能逐个函数独立改动，而不是每次都在一个超长脚本里改。

## 运行方式

### 方式 A：最简单的调用方式（推荐）

在 MATLAB 当前目录切到仓库根目录后，直接执行：

```matlab
results = particle_bubble_collision_demo();
results = particle_bubble_collision_demo('validation');
results = particle_bubble_collision_demo('dynamic');
```

其中 `particle_bubble_collision_demo.m` 是**兼容入口包装器**，会自动把路径切到 `matlab_particle_bubble/`，然后转发到 `pb_run.m`。

### 方式 B：直接运行快速上手脚本

如果你不想手动记函数名，可以直接运行：

```matlab
particle_bubble_quickstart
```

这个脚本里已经给出了：

- 整体调用示例；
- 只跑验证算例的示例；
- 只跑动力学算例的示例；
- 自定义参数后调用模块函数的示例；
- 最底层 `pb_build_adaptive_grid` + `pb_solve_reynolds_pressure` 的示例。

### 方式 C：直接调用模块函数（进阶）

如果你要自己改参数、改判据、改求解器，可以直接这样写：

```matlab
addpath('matlab_particle_bubble');
par = pb_default_parameters();
par.initialGap = 15e-6;
par.validationUrel = 2e-3;

validation = pb_run_validation_case(par, true);
dynamic    = pb_run_dynamic_case(par, true);
```

再往下，你也可以手动调用最底层模块：

```matlab
addpath('matlab_particle_bubble');
par  = pb_default_parameters();
grid = pb_build_adaptive_grid(1e-6, par);
sol  = pb_solve_reynolds_pressure(grid, 1e-6, 1e-3, par);
```

## 下一步建议

下一阶段最值得继续做的三件事：

1. 把当前轴对称润滑模型换成“薄膜区 + 外流场”耦合模型；
2. 将 `pb_refinement_indicator.m` 升级为真正可发表的方法判据；
3. 在该模块化框架基础上，再考虑迁移到更高保真的 IBM-LBM 或 C++ 实现。
