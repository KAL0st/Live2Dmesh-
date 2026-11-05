#pragma once
//
// MeshTopology —— ArtMesh 网格拓扑工具（P0-2，docs/33 §2.10-§2.17）
//
// 这里放三类东西，全部是**纯几何**、不碰 Qt / 文档 / 渲染：
//   ① 从三角形集合反推边界环（官方靠「边类型」知道外轮廓；我们的 cmo3 只存 indices，
//      所以轮廓只能从三角形自己反推 —— 这正是官方 calcOutlinesAsCW 的等价物）；
//   ② 点是否在多边形内（射线法，含边上的判据），给轮廓裁剪用；
//   ③ **三角化入口**：约束 Delaunay + 按轮廓裁掉凸包多出来的那一圈三角形。
//
// 为什么必须裁剪：官方 setup_indices_Delaunay 的输入带 outlineIndexCW，最后一步
// removeOutsideTriangles 把重心落在轮廓外的三角形扔掉（AGTriangulationBase:35）。
// 一个手臂/发梢这种**凹形**网格，纯 Delaunay 会把凹口也用三角形填平 —— 直接照抄
// 索引就会让 mesh 长出一块原本不存在的皮。
//
// ⚠️ 三角化要按 **mesh 自身空间的顶点** 做（cmo3 的 CArtMeshSource.positions），
//   不能按画布坐标：画布经过非均匀缩放，Delaunay 的对角线选择会变（不是相似变换）。
//
#include <string>
#include <array>
#include <utility>
#include <vector>

#include "edit/DelaunayTri.h"

namespace l2d {

/// 一个顶点环（闭合，首尾不重复）
using VertexLoop = std::vector<int>;

struct MeshTopologyInput {
    std::vector<float> pts;                     ///< 2n 个浮点（x,y 交错）
    std::vector<std::pair<int,int>> fixedEdges; ///< 必须保留的边（官方 = 非 AUTO_TRIANGULATION 边）
    std::vector<VertexLoop> outlines;           ///< 裁剪轮廓（空 = 不裁剪，只做 Delaunay）
};

struct MeshTopologyResult {
    std::vector<int> tris;              ///< 3 个索引一个三角形，逆时针
    std::vector<VertexLoop> boundary;   ///< 三角化结果自身的边界环（裁剪后）
    int outsideRemoved = 0;             ///< 被轮廓裁掉的三角形数
    int flips = 0;
    bool ok = true;
    const char* note = "";
};

/// 从三角形集合反推边界环（沿有向边计数、成对抵消，再按顶点串成环）。
///   · nonManifold 非空时累加：同向边重复 / 顶点出度 > 1 / 环走不通 的次数。
///   · 有洞的网格（多个环）会全部返回；顺序不保证。
std::vector<VertexLoop> boundaryLoops(const std::vector<int>& tris, int vertexCount, int* nonManifold = nullptr);

/// 点是否在环内（射线法 + 边上判据）；坐标与 pts 同空间。
bool pointInLoop(const std::vector<float>& pts, const VertexLoop& loop, float x, float y);

/// 环的有向面积（正 = 逆时针；y 向上的坐标系）
double loopArea(const std::vector<float>& pts, const VertexLoop& loop);

/// 三角形集合覆盖的有向面积之和（有符号；正常三角化应为正 = 逆时针）
double trianglesArea(const std::vector<float>& pts, const std::vector<int>& tris);

/// 三角化 + 轮廓裁剪。失败时 ok=false 且 note 说明原因（tris 可能是空）。
MeshTopologyResult rebuildMeshTriangles(const MeshTopologyInput& in, const DelaunayOptions& opt = DelaunayOptions{});

} // namespace l2d
