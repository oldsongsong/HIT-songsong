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

在 MATLAB 当前目录切到仓库根目录后，可直接执行：

```matlab
results = particle_bubble_collision_demo();
results = particle_bubble_collision_demo('validation');
results = particle_bubble_collision_demo('dynamic');
```

其中 `particle_bubble_collision_demo.m` 现在只是一个**兼容入口包装器**，真正的实现已迁移到 `matlab_particle_bubble/` 目录下。

## 下一步建议

下一阶段最值得继续做的三件事：

1. 把当前轴对称润滑模型换成“薄膜区 + 外流场”耦合模型；
2. 将 `pb_refinement_indicator.m` 升级为真正可发表的方法判据；
3. 在该模块化框架基础上，再考虑迁移到更高保真的 IBM-LBM 或 C++ 实现。
