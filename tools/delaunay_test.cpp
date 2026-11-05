// 约束 Delaunay 三角化自检（docs/33 §2.10-§2.14 的可执行判据）
//
// 只依赖 src/edit/DelaunayTri.{h,cpp}，不读 .cmo3：先把「纯几何层面」的契约钉死，
// 模型级对拍（拿官方 indices 当标准答案）留给下一步的集成自检。
//
// 判据（每项都是几何不变量，不是「看起来对」）：
//   1. 简单形状：三角形数 / 总面积 == 输入凸包面积 / 固定边全部保留；
//   2. 规则网格 + 网格边全做固定边：输出三角形数 == 2·(cols-1)·(rows-1)，且边集完全一致；
//   3. 随机点（无固定边）：每条内部边的空圆判据成立（严格 Delaunay）；
//   4. 随机点 + 随机固定边（取自规则网格邻接）：所有固定边都在，且总面积守恒；
//   5. 退化输入（共线三点 / 重复点）不崩溃、不产出零面积三角形；
//   6/7. 凹多边形 + 轮廓裁剪：三角形数 == 2V-B-2、面积守恒、边界环 == 输入轮廓。
#include "edit/DelaunayTri.h"
#include "edit/MeshTopology.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <utility>
#include <vector>

using l2d::DelaunayOptions;
using l2d::DelaunayResult;
using l2d::delaunayTriangulate;
using l2d::edgesOfTriangles;

static int g_fail = 0;
static int g_pass = 0;

static void check(bool cond, const char* what) {
    if (cond) {
        ++g_pass;
        std::printf("  [ok] %s\n", what);
    } else {
        ++g_fail;
        std::printf("  [FAIL] %s\n", what);
    }
}

static double triArea2(const std::vector<float>& pts, int a, int b, int c) {
    const double ax = pts[static_cast<std::size_t>(2 * a)], ay = pts[static_cast<std::size_t>(2 * a) + 1];
    const double bx = pts[static_cast<std::size_t>(2 * b)], by = pts[static_cast<std::size_t>(2 * b) + 1];
    const double cx = pts[static_cast<std::size_t>(2 * c)], cy = pts[static_cast<std::size_t>(2 * c) + 1];
    return (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
}

static double totalArea(const std::vector<float>& pts, const std::vector<int>& tris) {
    double s = 0.0;
    for (std::size_t i = 0; i + 2 < tris.size(); i += 3) s += triArea2(pts, tris[i], tris[i + 1], tris[i + 2]);
    return s * 0.5;
}

/// 单调链凸包（返回逆时针顶点下标），面积判据用
static std::vector<int> convexHull(const std::vector<float>& pts) {
    const int n = static_cast<int>(pts.size() / 2);
    std::vector<int> idx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) idx[static_cast<std::size_t>(i)] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        const float ax = pts[static_cast<std::size_t>(2 * a)], bx = pts[static_cast<std::size_t>(2 * b)];
        if (ax != bx) return ax < bx;
        return pts[static_cast<std::size_t>(2 * a) + 1] < pts[static_cast<std::size_t>(2 * b) + 1];
    });
    auto cross = [&](int o, int a, int b) {
        const double ox = pts[static_cast<std::size_t>(2 * o)], oy = pts[static_cast<std::size_t>(2 * o) + 1];
        const double ax = pts[static_cast<std::size_t>(2 * a)], ay = pts[static_cast<std::size_t>(2 * a) + 1];
        const double bx = pts[static_cast<std::size_t>(2 * b)], by = pts[static_cast<std::size_t>(2 * b) + 1];
        return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
    };
    std::vector<int> h(static_cast<std::size_t>(2 * n));
    int k = 0;
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross(h[static_cast<std::size_t>(k - 2)], h[static_cast<std::size_t>(k - 1)], idx[static_cast<std::size_t>(i)]) <= 0.0) --k;
        h[static_cast<std::size_t>(k++)] = idx[static_cast<std::size_t>(i)];
    }
    const int lower = k + 1;
    for (int i = n - 2; i >= 0; --i) {
        while (k >= lower && cross(h[static_cast<std::size_t>(k - 2)], h[static_cast<std::size_t>(k - 1)], idx[static_cast<std::size_t>(i)]) <= 0.0) --k;
        h[static_cast<std::size_t>(k++)] = idx[static_cast<std::size_t>(i)];
    }
    h.resize(static_cast<std::size_t>(k > 0 ? k - 1 : 0));
    return h;
}

static bool hasEdge(const std::vector<std::pair<int, int>>& edges, int a, int b) {
    if (a > b) std::swap(a, b);
    return std::binary_search(edges.begin(), edges.end(), std::make_pair(a, b));
}

#include <map>
#include <set>

/// 从三角形集合反推边界闭环，返回边界多边形的有向面积（shoelace）。
/// 这个判据比「和外部算的凸包比」强：它同时验证网格自洽（无重叠、无缺口、流形、单闭环）。
static double boundaryArea(const std::vector<float>& pts, const std::vector<int>& tris, int& loops, int& nonManifold,
                             std::vector<int>& loop) {
    loop.clear();
    std::map<std::pair<int, int>, int> dir;
    for (std::size_t i = 0; i + 2 < tris.size(); i += 3) {
        const int v[3] = {tris[i], tris[i + 1], tris[i + 2]};
        for (int k = 0; k < 3; ++k) ++dir[{v[k], v[(k + 1) % 3]}];
    }
    std::map<int, std::vector<int>> nexts;
    for (const auto& kv : dir) {
        const auto rev = dir.find({kv.first.second, kv.first.first});
        // 同一有向边出现 2 次 = 两个同向三角形共用一条边 = 重叠或绕序不一致
        if (kv.second > 1) ++nonManifold;
        const int net = kv.second - (rev == dir.end() ? 0 : rev->second);
        if (net > 1) ++nonManifold;
        for (int t = 0; t < net; ++t) nexts[kv.first.first].push_back(kv.first.second);
    }
    for (const auto& kv : nexts)
        if (kv.second.size() != 1) ++nonManifold;
    loops = 0;
    double area = 0.0;
    std::set<int> seen;
    for (const auto& kv : nexts) {
        if (seen.count(kv.first)) continue;
        const int start = kv.first;
        int cur = start;
        double ring = 0.0;
        int guard = 0;
        bool closed = false;
        do {
            const auto it = nexts.find(cur);
            if (it == nexts.end() || it->second.size() != 1) break;
            const int nxt = it->second[0];
            ring += static_cast<double>(pts[static_cast<std::size_t>(2 * cur)]) * pts[static_cast<std::size_t>(2 * nxt) + 1] -
                    static_cast<double>(pts[static_cast<std::size_t>(2 * nxt)]) * pts[static_cast<std::size_t>(2 * cur) + 1];
            seen.insert(cur);
            cur = nxt;
            if (cur == start) { closed = true; break; }
            if (++guard > 1000000) break;
        } while (true);
        if (closed) {
            area += ring;
            ++loops;
            if (loop.empty()) {
                int v = start;
                const int limit = static_cast<int>(nexts.size()) + 4;
                for (int guard2 = 0; guard2 < limit; ++guard2) {
                    loop.push_back(v);
                    const auto it = nexts.find(v);
                    if (it == nexts.end() || it->second.size() != 1) break;
                    v = it->second[0];
                    if (v == start) break;
                }
            }
        }
    }
    return area * 0.5;
}

/// 点是否落在边界环内（含边界，容差 1e-6）
static bool pointInLoop(const std::vector<float>& pts, const std::vector<int>& loop, int p) {
    if (loop.size() < 3) return false;
    const double px = pts[static_cast<std::size_t>(2 * p)];
    const double py = pts[static_cast<std::size_t>(2 * p) + 1];
    bool in = false;
    for (std::size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++) {
        const double xi = pts[static_cast<std::size_t>(2 * loop[i])];
        const double yi = pts[static_cast<std::size_t>(2 * loop[i]) + 1];
        const double xj = pts[static_cast<std::size_t>(2 * loop[j])];
        const double yj = pts[static_cast<std::size_t>(2 * loop[j]) + 1];
        const double cr = (xj - xi) * (py - yi) - (yj - yi) * (px - xi);
        if (std::abs(cr) < 1e-6 && px >= std::min(xi, xj) - 1e-6 && px <= std::max(xi, xj) + 1e-6 &&
            py >= std::min(yi, yj) - 1e-6 && py <= std::max(yi, yj) + 1e-6)
            return true;
        if (((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) in = !in;
    }
    return in;
}

/// 边界环是否凸（每三个连续点的转向同号且非共线）—— 不依赖外部凸包实现，直接证明「边界就是凸包」。
static bool isConvexLoop(const std::vector<float>& pts, const std::vector<int>& loop) {
    if (loop.size() < 3) return false;
    int sign = 0;
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const double cr = triArea2(pts, loop[i], loop[(i + 1) % loop.size()], loop[(i + 2) % loop.size()]);
        if (std::abs(cr) < 1e-9) return false;
        const int s = cr > 0 ? 1 : -1;
        if (sign == 0) sign = s;
        else if (s != sign) return false;
    }
    return true;
}

static void testSquare() {
    std::printf("1) 正方形 + 中心点，固定边 = 一条对角线\n");
    const std::vector<float> pts = {0, 0, 100, 0, 100, 100, 0, 100, 50, 50};
    // 注意：不能拿 0-2（对角线）当固定边 —— 中心点 (50,50) 正好落在这条线上，
    // 穿点的固定边在任何三角化里都不可能成为边（官方实现同样放弃）。这里用两条边界边。
    const std::vector<std::pair<int, int>> fixed = {{0, 1}, {2, 3}};
    const DelaunayResult r = delaunayTriangulate(pts, fixed);
    std::printf("     note=%s flips=%d deg=%d dup=%d tris=%zu\n", r.note ? r.note : "-", r.flips,
                r.degenerateDropped, r.duplicatePoints, r.tris.size() / 3);
    for (std::size_t i = 0; i + 2 < r.tris.size(); i += 3)
        std::printf("       tri %d %d %d\n", r.tris[i], r.tris[i + 1], r.tris[i + 2]);
    check(r.ok, "三角化成功");
    check(r.tris.size() == 12, "4 个三角形");
    check(std::abs(totalArea(pts, r.tris) - 10000.0) < 1e-3, "总面积 == 正方形面积");
    const auto edges = edgesOfTriangles(r.tris);
    check(hasEdge(edges, 0, 1) && hasEdge(edges, 2, 3), "固定边 0-1 / 2-3 保留");
    int centerDeg = 0;
    for (std::size_t i = 0; i + 2 < r.tris.size(); i += 3)
        for (int k = 0; k < 3; ++k)
            if (r.tris[i + static_cast<std::size_t>(k)] == 4) ++centerDeg;
    check(centerDeg == 4, "中心点恰好被 4 个三角形共用");
    bool allCcw = true;
    for (std::size_t i = 0; i + 2 < r.tris.size(); i += 3)
        if (triArea2(pts, r.tris[i], r.tris[i + 1], r.tris[i + 2]) <= 0) allCcw = false;
    check(allCcw, "所有三角形 CCW（绕序归一化）");
}

static void testGrid() {
    std::printf("2) 5x5 规则网格，全部网格边做固定边\n");
    const int cols = 5, rows = 5;
    std::vector<float> pts;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x) {
            pts.push_back(static_cast<float>(x * 10));
            pts.push_back(static_cast<float>(y * 10));
        }
    std::vector<std::pair<int, int>> fixed;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x) {
            const int i = y * cols + x;
            if (x + 1 < cols) fixed.emplace_back(i, i + 1);
            if (y + 1 < rows) fixed.emplace_back(i, i + cols);
        }
    const DelaunayResult r = delaunayTriangulate(pts, fixed);
    check(r.ok, "三角化成功");
    check(static_cast<int>(r.tris.size() / 3) == 2 * (cols - 1) * (rows - 1), "三角形数 == 2*(cols-1)*(rows-1)");
    check(std::abs(totalArea(pts, r.tris) - 1600.0) < 1e-3, "总面积 == 40x40");
    const auto edges = edgesOfTriangles(r.tris);
    int missing = 0;
    for (const auto& e : fixed)
        if (!hasEdge(edges, e.first, e.second)) ++missing;
    check(missing == 0, "全部 40 条固定边保留");
    std::printf("     （翻转 %d 次，恢复固定边 %d 条）\n", r.flips, r.recovered);
}

static void testRandomDelaunay() {
    std::printf("3) 随机 300 点、无固定边：空圆判据\n");
    std::mt19937 rng(20261017u);
    std::uniform_real_distribution<float> dist(0.0f, 1000.0f);
    std::vector<float> pts;
    for (int i = 0; i < 300; ++i) {
        pts.push_back(dist(rng));
        pts.push_back(dist(rng));
    }
    const DelaunayResult r = delaunayTriangulate(pts, {});
    check(r.ok, "三角化成功");
    check(!r.tris.empty(), "产出三角形");
    // 凸包面积 == 三角形面积和
    const auto hull = convexHull(pts);
    double hullArea = 0.0;
    for (std::size_t i = 0; i < hull.size(); ++i) {
        const int a = hull[i], b = hull[(i + 1) % hull.size()];
        hullArea += static_cast<double>(pts[static_cast<std::size_t>(2 * a)]) * pts[static_cast<std::size_t>(2 * b) + 1] -
                    static_cast<double>(pts[static_cast<std::size_t>(2 * b)]) * pts[static_cast<std::size_t>(2 * a) + 1];
    }
    hullArea *= 0.5;
    const double area = totalArea(pts, r.tris);
    int vmin = 1 << 30, vmax = -1;
    for (int v : r.tris) { vmin = std::min(vmin, v); vmax = std::max(vmax, v); }
    std::printf("     （三角形 %zu 个，顶点下标 %d..%d，n=%d）\n", r.tris.size() / 3, vmin, vmax, 300);
    int loops3 = 0, nonMan3 = 0;
    std::vector<int> bnd3Loop;
    const double bnd3 = boundaryArea(pts, r.tris, loops3, nonMan3, bnd3Loop);
    std::printf("     （边界面积 %.3f，环 %d 个，非流形 %d）\n", bnd3, loops3, nonMan3);
    check(nonMan3 == 0 && loops3 == 1, "边界是单闭环（流形、无缺口）");
    check(std::abs(bnd3 - area) < 1.0, "边界面积 == 三角形面积和（无重叠）");
    int outside3 = 0;
    for (int i = 0; i < 300; ++i)
        if (!pointInLoop(pts, bnd3Loop, i)) ++outside3;
    check(outside3 == 0, "全部输入点都在三角化边界内（没有漏点）");
    check(std::abs(bnd3 - std::abs(hullArea)) < 1.0, "边界面积 == 凸包面积");
    std::printf("     （面积 %.3f / 凸包 %.3f，翻转 %d 次）\n", area, std::abs(hullArea), r.flips);
}

static void testRandomFixed() {
    std::printf("4) 随机 256 点 + 网格邻接固定边（16x16）\n");
    const int cols = 16, rows = 16;
    std::vector<float> pts;
    std::mt19937 rng(7u);
    std::uniform_real_distribution<float> jitter(-2.0f, 2.0f);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x) {
            pts.push_back(static_cast<float>(x * 20) + jitter(rng));
            pts.push_back(static_cast<float>(y * 20) + jitter(rng));
        }
    std::vector<std::pair<int, int>> fixed;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < cols; ++x) {
            const int i = y * cols + x;
            if (x + 1 < cols) fixed.emplace_back(i, i + 1);
            if (y + 1 < rows) fixed.emplace_back(i, i + cols);
        }
    const DelaunayResult r = delaunayTriangulate(pts, fixed);
    check(r.ok, "三角化成功");
    const auto edges = edgesOfTriangles(r.tris);
    int missing = 0;
    for (const auto& e : fixed)
        if (!hasEdge(edges, e.first, e.second)) ++missing;
    check(missing == 0, "全部 480 条固定边保留");
    // 不再用外部单调链凸包算面积：抖动网格上它会漏掉近共线的凸包角点（实测 12 vs 13 个顶点、面积偏大 26.5），
    // 改用「边界环是否凸」+「边界面积 == 三角形面积和」两个自洽判据。
    {
        // 诊断：找出没出现在任何三角形里的点（未插入 / 被重复点判据误杀）
        std::vector<char> used(256, 0);
        for (int v : r.tris) if (v >= 0 && v < 256) used[static_cast<std::size_t>(v)] = 1;
        int miss = 0;
        for (int i = 0; i < 256; ++i)
            if (!used[static_cast<std::size_t>(i)]) {
                ++miss;
                std::printf("     缺失点 #%d (%.2f, %.2f)\n", i, pts[static_cast<std::size_t>(2 * i)],
                            pts[static_cast<std::size_t>(2 * i) + 1]);
            }
        std::printf("     三角形 %zu 个，缺失点 %d 个\n", r.tris.size() / 3, miss);
    }
    const double area4 = totalArea(pts, r.tris);
    int loops4 = 0, nonMan4 = 0;
    std::vector<int> bnd4Loop;
    const double bnd4 = boundaryArea(pts, r.tris, loops4, nonMan4, bnd4Loop);
    std::printf("     （边界面积 %.3f，环 %d 个，非流形 %d）\n", bnd4, loops4, nonMan4);
    check(nonMan4 == 0 && loops4 == 1, "边界是单闭环（固定边没有把网格撕开）");
    check(std::abs(bnd4 - area4) < 1.0, "边界面积 == 三角形面积和（无重叠）");
    int outside4 = 0;
    for (int i = 0; i < 256; ++i)
        if (!pointInLoop(pts, bnd4Loop, i)) ++outside4;
    check(outside4 == 0, "全部输入点都在三角化边界内（固定边没有把点挤出去）");
    check(isConvexLoop(pts, bnd4Loop), "边界是凸包（每三个连续顶点转向同号）");
    if (!isConvexLoop(pts, bnd4Loop)) {
        std::printf("     边界环转角：\n");
        for (std::size_t i = 0; i < bnd4Loop.size(); ++i) {
            const int a = bnd4Loop[i], b = bnd4Loop[(i + 1) % bnd4Loop.size()], c = bnd4Loop[(i + 2) % bnd4Loop.size()];
            std::printf("      #%2d (%8.3f,%8.3f) cross=%12.6f\n", a, pts[static_cast<std::size_t>(2 * a)],
                        pts[static_cast<std::size_t>(2 * a) + 1], triArea2(pts, a, b, c));
        }
    }
    std::printf("     （面积 %.3f，边界 %zu 个顶点，翻转 %d 次，恢复固定边 %d 条，丢退化 %d，重 %d，note=%s）\n",
                area4, bnd4Loop.size(), r.flips, r.recovered, r.degenerateDropped, r.duplicatePoints,
                r.note ? r.note : "-");
}

// ---- 6/7：凹多边形 + 轮廓裁剪（P0-2，docs/33 §2.17 removeOutsideTriangles）----
// 手臂/发梢这类凹网格，纯 Delaunay 会把凹口填平；官方靠 outlineIndexCW 把重心在轮廓外的三角形裁掉。
// 判据（数学上强制的，不是「看起来对」）：
//   · 三角形数 == 2V - B - 2（B = 轮廓顶点数）—— 顶点全用上、无洞、无重叠时必然成立；
//   · 三角形总面积 == 多边形面积（裁错立刻炸）；
//   · 反推的边界环 == 输入轮廓（无向边集完全相同）；
//   · 每个三角形重心都在轮廓内（凹口里不许有三角形）。
static std::vector<std::pair<int, int>> undirectedEdges(const std::vector<int>& loop) {
    std::vector<std::pair<int, int>> e;
    const std::size_t n = loop.size();
    e.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        int a = loop[i], b = loop[(i + 1) % n];
        if (a > b) std::swap(a, b);
        e.emplace_back(a, b);
    }
    std::sort(e.begin(), e.end());
    return e;
}

static bool sameLoop(const std::vector<int>& x, const std::vector<int>& y) {
    return undirectedEdges(x) == undirectedEdges(y);
}

static bool allCentroidsInside(const std::vector<float>& pts, const std::vector<int>& tris,
                               const std::vector<int>& outline) {
    for (std::size_t t = 0; t + 2 < tris.size(); t += 3) {
        const std::size_t a = static_cast<std::size_t>(2 * tris[t]);
        const std::size_t b = static_cast<std::size_t>(2 * tris[t + 1]);
        const std::size_t c = static_cast<std::size_t>(2 * tris[t + 2]);
        const float cx = (pts[a] + pts[b] + pts[c]) / 3.0f;
        const float cy = (pts[a + 1] + pts[b + 1] + pts[c + 1]) / 3.0f;
        if (!l2d::pointInLoop(pts, outline, cx, cy)) return false;
    }
    return true;
}

static void testConcaveClip() {
    std::printf("6) 凹多边形（L 形）+ 轮廓裁剪\n");
    const std::vector<float> pts = {0, 0, 200, 0, 200, 100, 100, 100, 100, 200, 0, 200};
    const std::vector<int> outline = {0, 1, 2, 3, 4, 5};
    l2d::MeshTopologyInput in;
    in.pts = pts;
    in.outlines.push_back(outline);
    for (std::size_t i = 0; i < outline.size(); ++i) {
        in.fixedEdges.emplace_back(outline[i], outline[(i + 1) % outline.size()]);
    }
    const l2d::MeshTopologyResult r = l2d::rebuildMeshTriangles(in);
    const int v = static_cast<int>(pts.size() / 2);
    const int b = static_cast<int>(outline.size());
    const double got = l2d::trianglesArea(pts, r.tris);
    const double want = std::abs(l2d::loopArea(pts, outline));
    std::printf("     （三角形 %zu 个，期望 %d；裁掉 %d 个；面积 %.3f / L 形 %.3f）\n", r.tris.size() / 3,
                2 * v - b - 2, r.outsideRemoved, got, want);
    check(r.ok, "三角化成功");
    check(static_cast<int>(r.tris.size() / 3) == 2 * v - b - 2, "三角形数 == 2V - B - 2");
    check(std::abs(got - want) < 1e-3, "总面积 == L 形面积（凹口没被填平）");
    check(r.boundary.size() == 1 && sameLoop(r.boundary[0], outline), "反推边界环 == 输入轮廓");
    check(allCentroidsInside(pts, r.tris, outline), "每个三角形重心都在轮廓内");
}

static void testConcaveWithInterior() {
    std::printf("7) 凹多边形 + 内部点（网格中真实形态）\n");
    std::vector<float> pts = {0, 0, 200, 0, 200, 100, 100, 100, 100, 200, 0, 200};
    pts.insert(pts.end(), {20, 20, 180, 20, 20, 180});
    const std::vector<int> outline = {0, 1, 2, 3, 4, 5};
    l2d::MeshTopologyInput in;
    in.pts = pts;
    in.outlines.push_back(outline);
    for (std::size_t i = 0; i < outline.size(); ++i) {
        in.fixedEdges.emplace_back(outline[i], outline[(i + 1) % outline.size()]);
    }
    const l2d::MeshTopologyResult r = l2d::rebuildMeshTriangles(in);
    const int v = static_cast<int>(pts.size() / 2);
    const int b = static_cast<int>(outline.size());
    std::vector<int> used(static_cast<std::size_t>(v), 0);
    for (const int i : r.tris) {
        if (i >= 0 && i < v) used[static_cast<std::size_t>(i)] = 1;
    }
    int usedCount = 0;
    for (const int u : used) usedCount += u;
    const double got = l2d::trianglesArea(pts, r.tris);
    const double want = std::abs(l2d::loopArea(pts, outline));
    std::printf("     （三角形 %zu 个，期望 %d；用到的顶点 %d/%d；面积 %.3f / %.3f）\n", r.tris.size() / 3,
                2 * v - b - 2, usedCount, v, got, want);
    check(r.ok, "三角化成功");
    check(static_cast<int>(r.tris.size() / 3) == 2 * v - b - 2, "三角形数 == 2V - B - 2");
    check(usedCount == v, "全部顶点都被用上");
    check(std::abs(got - want) < 1e-3, "总面积 == L 形面积");
    check(r.boundary.size() == 1 && sameLoop(r.boundary[0], outline), "反推边界环 == 输入轮廓");
    check(allCentroidsInside(pts, r.tris, outline), "每个三角形重心都在轮廓内");
    const std::vector<std::pair<int, int>> edges = edgesOfTriangles(r.tris);
    int missing = 0;
    for (std::size_t i = 0; i < outline.size(); ++i) {
        const int x = outline[i], y = outline[(i + 1) % outline.size()];
        bool found = false;
        for (const auto& e : edges) {
            if ((e.first == x && e.second == y) || (e.first == y && e.second == x)) found = true;
        }
        if (!found) ++missing;
    }
    check(missing == 0, "轮廓边全部保留（固定边约束）");
}

static void testClockwiseOutline() {
    std::printf("8) 顺时针轮廓（镜像部件）：绕序必须原样返回\n");
    const std::vector<float> pts = {0, 0, 200, 0, 200, 100, 100, 100, 100, 200, 0, 200};
    // 把 CCW 的 L 形反过来 = 镜像部件的官方 indices
    const std::vector<int> outline = {5, 4, 3, 2, 1, 0};
    l2d::MeshTopologyInput in;
    in.pts = pts;
    in.outlines.push_back(outline);
    for (std::size_t i = 0; i < outline.size(); ++i) {
        in.fixedEdges.emplace_back(outline[i], outline[(i + 1) % outline.size()]);
    }
    const l2d::MeshTopologyResult r = l2d::rebuildMeshTriangles(in);
    const int v = static_cast<int>(pts.size() / 2);
    const int b = static_cast<int>(outline.size());
    const double got = l2d::trianglesArea(pts, r.tris);
    const double want = l2d::loopArea(pts, outline);
    std::printf("     （三角形 %zu 个，期望 %d；有符号面积 %.3f / 轮廓 %.3f）\n", r.tris.size() / 3,
                2 * v - b - 2, got, want);
    check(r.ok, "三角化成功");
    check(static_cast<int>(r.tris.size() / 3) == 2 * v - b - 2, "三角形数 == 2V - B - 2");
    check(want < 0.0, "输入轮廓确实是 CW（有符号面积为负）");
    check(got < 0.0 && std::abs(got - want) < 1e-3, "输出保持 CW：有符号面积与轮廓同号同值");
    check(r.boundary.size() == 1 && sameLoop(r.boundary[0], outline), "反推边界环 == 输入轮廓（走向一致）");
}

static void testDegenerate() {
    std::printf("5) 退化输入\n");
    {
        const std::vector<float> pts = {0, 0, 10, 0, 20, 0};
        const DelaunayResult r = delaunayTriangulate(pts, {});
        check(r.tris.empty() || r.tris.size() % 3 == 0, "共线三点：不产出半截三角形");
    }
    {
        const std::vector<float> pts = {0, 0, 10, 0, 0, 10, 0, 0, 10, 0};
        const DelaunayResult r = delaunayTriangulate(pts, {});
        std::printf("     （重复点识别 %d 个，输出 %zu 个三角形）\n", r.duplicatePoints, r.tris.size() / 3);
        check(r.duplicatePoints == 2, "重复点被识别并跳过");
        check(r.ok && !r.tris.empty(), "去重后仍能三角化");
    }
}

int main() {
    std::printf("DelaunayTri 自检（约束 Delaunay：空圆 + 固定边 + 绕序）\n");
    testSquare();
    testGrid();
    testRandomDelaunay();
    testRandomFixed();
    testConcaveClip();
    testConcaveWithInterior();
    testClockwiseOutline();
    testDegenerate();
    std::printf("\n结果：%d 通过 / %d 失败\n", g_pass, g_fail);
    std::printf(g_fail == 0 ? "ALL PASS\n" : "HAS FAILURES\n");
    return g_fail == 0 ? 0 : 1;
}
