#pragma once
//
// DelaunayTri —— 约束 Delaunay 三角化（网格编辑的拓扑内核，docs/33 §2.10-§2.14）
// ---------------------------------------------------------------------------
// 为什么需要它：本仓库原来的三角形只来自 .cmo3 的 `indices`（一份已烘焙的三角化），
// 一旦允许拓扑编辑（加点 / 加边 / 删边 / 翻转对角线），就必须自己重算三角化。
// 官方的口径（docs/33 卡 3/4/5/6）：
//   · 三角化器是「带固定边约束的 Delaunay」（GTriangulation_Delaunay3），不是耳切；
//   · 固定边 = 非自动边（NORMAL / LOCKED / USER_TRIANGULATION），AUTO_TRIANGULATION 全部重算；
//   · 合法性判据 = 空外接圆：`radius - dist > EPS` 就要翻；EPS = 1.0E-6f；
//   · 翻转前先判「新对角线是否还与固定边相交」，再判空圆（官方 checkDelaunay 第 4/6 步）；
//   · 输出必须归一化绕序（正有向面积），否则面剔除 / 法线全反（docs/33 §2.14）。
//
// 本实现的算法结构（与官方一致的部分 / 有意不同的部分）：
//   1. 超大三角形 + 逐点插入（Bowyer-Watson 空圆驱逐）：等价于官方「先插固定边端点、
//      再插其余点」的增量插入，但不需要把固定边端点单独排序——点集一样，结果一样；
//   2. **固定边恢复**：对每条固定边，若不在边集里，就反复 flip 与它严格相交的内部边，
//      直到它成为一条边（官方用 pending 队列 + 500 次硬上限；我们换成「有进展就继续、
//      无进展才失败」的收敛判据 —— docs/33 §4.4 明确允许超越这一点）；
//   3. 合法化（LOP）：固定边不动，其余内部边按空圆判据翻转直到稳定；
//   4. 去掉与超大三角形相连的三角形，绕序归一化为 CCW，丢弃退化三角形。
//
// 输入输出契约：
//   · pts：2 float / 顶点，**已经换算到同一坐标空间**（调用方负责本地/画布变换）；
//   · fixedEdges：顶点下标对（0..n-1），顺序无关、允许重复；
//   · 输出 tris：3 int / 三角形，绕序统一为正面积（CCW），可直接当 GL indices 用；
//   · 输出不包含退化（|面积| <= eps）三角形；顶点若重复坐标会跳过（affectedNote 里说明）。
//
#include <cstddef>
#include <utility>
#include <vector>

namespace l2d {

struct DelaunayOptions {
    /// 浮点容差。官方 `L.a.f()` = 1.0E-6f（docs/33 §0.1）；坐标系是画布 px 时够用。
    double eps = 1.0e-6;
    /// 翻转总次数上限（防御性；正常网格远用不到）。到上限仍不收敛 ⇒ ok=false。
    int maxFlips = 2000000;
};

struct DelaunayResult {
    std::vector<int> tris;            ///< 3 int / 三角形，CCW
    int flips = 0;                    ///< 翻转次数（诊断 / 性能观测）
    int recovered = 0;                ///< 需要恢复才成为边的固定边条数
    int degenerateDropped = 0;        ///< 因为 |面积| <= eps 被丢掉的三角形数
    int duplicatePoints = 0;          ///< 与已插点重合（< eps）而被跳过的点数
    bool ok = true;
    const char* note = "";            ///< ok=false 时的原因（静态字符串）
};

/// 约束 Delaunay 三角化。失败时 ok=false，tris 内容不应使用。
DelaunayResult delaunayTriangulate(const std::vector<float>& pts,
                                   const std::vector<std::pair<int, int>>& fixedEdges,
                                   const DelaunayOptions& opt = DelaunayOptions {});

/// 多边形内三角化（polygon-first）：先把轮廓环耳切成三角形，再把内部点逐点插进去。
///
/// 为什么光有 delaunayTriangulate 不够：它是「对全部点做 DT，再按三角形重心裁掉轮廓外的」，
/// 这个裁剪对**近乎共线**的边界顶点是不成立的 —— 弦把发丝宽的凹口整块填平。
/// 实测梦寝兔兔 rabbit.saved.cmo3 的 548 个 ArtMesh 里有 103 个踩到这条：
/// 面积只差 3e-9（相对），但那个边界顶点被推进内部、轮廓少一条边，
/// 用户以后拖它时轮廓不再跟着动，等于悄悄改了拓扑。
/// 官方的做法就是直接对轮廓三角化（docs/33 §2.10 的 outlineIndexCW），不是事后裁剪，
/// 所以这里照做：耳切保证每条轮廓边都在三角化里，内点插入再把它拉回 Delaunay。
///
/// outline：轮廓的顶点下标环，按顺序排列，**不自交、无重复顶点**；只支持单环（无洞）。
/// fixedEdges：额外的固定边（轮廓边会自动成为固定边，不必重复传）。
DelaunayResult delaunayTriangulatePolygon(const std::vector<float>& pts,
                                          const std::vector<int>& outline,
                                          const std::vector<std::pair<int, int>>& fixedEdges = {},
                                          const DelaunayOptions& opt = DelaunayOptions {});

/// 由三角形列表推边集（a<b 规范化的 pair，去重、升序）。拓扑层与测试共用。
std::vector<std::pair<int, int>> edgesOfTriangles(const std::vector<int>& tris);

}  // namespace l2d
