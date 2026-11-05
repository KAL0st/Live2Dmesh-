# Live2D 自动布点算法

自研 ArtMesh 自动布点算法的独立源码快照，供阅读、讨论和移植。包含透明轮廓、镂空、锁边、纹理结构线、材质与对称采样、约束三角化和网格质量优化。

这是本地实现，不是官方 Cubism 5.3 源码，也不宣称与其数值结果等价。

## 源码范围

| 文件 | 内容 |
|---|---|
| `src/core/MeshAutoGen.{h,cpp}` | 完整生成流水线，包括本轮小镂空降密、重复描边弱化、交点合并，以及文档接入位置 |
| `src/core/MeshTextureFeatures.{h,cpp}` | RGB 结构候选、邻近色块线拟合、高光弱参考、柔化曲线、约束交点 |
| `src/core/MeshStripSampling.h` | 长条/细线中心线和点距换算 |
| `src/core/MeshSymmetry.h` | 透明轮廓对称检测和成对采样 |
| `src/edit/DelaunayTri.{h,cpp}` | 固定边约束下的三角化 |
| `src/edit/MeshTopology.{h,cpp}` | 自动布点所需的边界环、几何查询和轮廓裁剪；仅保留相关部分 |
| `tools/delaunay_test.cpp` | 不依赖模型、Qt 的几何回归测试 |

不包含编辑器界面、渲染、物理、参数系统、文件格式解析器、模型纹理、可执行文件、SDK/JAR、主项目 Git 历史或其他工作区改动。`SOURCE_MANIFEST.json` 记录实际导出范围及文件哈希。

## 构建边界

下面命令编译**纯几何和纹理分析库**，并运行约束三角化测试，仅需 C++20 编译器和 CMake：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

`MeshAutoGen.cpp` 保留完整流水线，便于查看镂空、锁边、覆盖检查等实现，但**不属于上述独立构建目标**。其入口仍使用宿主项目的 `Cmo3Document`、`XsParser`、`XsPrologue`、`TextureUvMapping` 和 `PngCodec`。这些文档/图片适配器没有被复制进算法仓库。独立应用要先把纹理读取、UV/关键形映射和文档写回替换为自己的数据接口；不能将纯几何库构建成功理解为完整自动布点接口已脱离宿主。

## 算法流程

1. 读取原生 RGBA 和 UV 映射，提取非零 alpha 外轮廓与完全透明的封闭孔洞。
2. 拟合外圈，区分软质/硬质支撑需求；可用对称检测进入镜像采样。
3. 按目标点距采样镂空，安全简化像素阶梯，再朝实体纹理方向生成反向锁边。
4. 识别纹理结构线。高光作弱参考；镂空部件中紧贴并沿轮廓走向的重复描边可降为弱参考。
5. 合并锁边与结构线的约束图，再补内部点并约束三角化。保护原始纹理线的端点、转折及手工线。
6. 合并安全的近点，优化面片形状；检查孔洞留空、纹理覆盖、UV 与关键形映射，最后写回文档。

## 主要控制项

| 参数 | 含义 |
|---|---|
| `polygonOuterDensity` | 外圈点距，越大越稀疏 |
| `polygonBorderDensity` | 外侧锁边点距，0 跟随外圈 |
| `polygonInnerDensity` | 内部点距，0 关闭自动补点 |
| `interiorTriangleTarget` | 内部三角形近似目标数量，结构约束可能使其超出 |
| `material` / `autoSymmetry` | 材质策略与对称检测 |
| `generateHoles` / `holeMinArea` | 保留镂空及最小面积过滤 |
| `holePointSpacing` / `holeBoundaryMargin` | 镂空轮廓点距及朝透明侧的安全退让 |
| `generateHoleBorder` / `holeBorderSpacing` / `holeBorderDistance` | 反向锁边开关、点距及向实体侧距离 |
| `MeshFeatureSettings::weakBoundaryDuplicates` | 重复描边弱参考开关；仅作用于有镂空的自动候选 |
| `MeshFeatureSettings::curveSmoothing` | 自动结构线柔化容差 |

点距采用 `MeshStripSampling.h` 中的工作尺度换算，不等于直接输入源图像素距离。轮廓退让、锁边距离等参数按源图像素解释。

## 本轮镂空优化

- 尖角判断结合一段轮廓的整体转向，不把所有像素阶梯直角都锁死。
- 相邻近点允许在更大的目标点距范围内简化，同时限制累计偏差、面积变化、最长合并边、自交和非透明像素覆盖。
- 原始 CV 端点和转折固定；锁边切割产生的共线冗余交点，可沿原直线合并。新的锁边必须仍包围对应孔洞。
- 重复描边弱化同时要求接近已有支撑轮廓和走向接近平行；不使用“色差小就删除”的规则。
- `mergedHoleVertices`、`mergedBorderVertices`、`mergedDerivedFeatureVertices` 和 `MeshFeatureLine::boundaryReference` 分别用于归因，避免把总点数减少误报成单一阶段的收益。

## 验证摘要

2026-11-06 宿主集成测试为 **509 passed / 0 failed**。同一模型四个真实纹理图层经生成、保存读回及独立覆盖检查通过；所有非零 alpha 像素中心和四角保持覆盖，孔洞留空。七孔袖子的点数从 694 降到 550，短于 1 源像素的边从 105 降到 47，仍存在需要进一步优化的窄带面片。

这些宿主测试和模型不在本仓库中；本仓库可直接运行的是 `delaunay_test` 纯几何测试。两类验证不能混为一谈。

## 限制

识别是几何与色差启发式，不是真正理解衣褶语义。等边/等腰三角形是优化目标，固定约束、尖端、窄带仍可能产生细长三角形。空间不足以放置独立锁边的极窄桥会拒绝生成，可减小距离或关闭镂空锁边；不会强行填孔。

镂空模式当前不强制对称采样。静态几何、保存读回和编译检查不等于完成所有动态变形或官方编辑器行为验证。本仓库未新增或改变项目自有源码的许可证。
