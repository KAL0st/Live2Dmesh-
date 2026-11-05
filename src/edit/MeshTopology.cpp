// MeshTopology.cpp —— 见 MeshTopology.h 的文件头说明
#include "edit/MeshTopology.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

namespace l2d {
namespace {

std::uint64_t packEdge(int a, int b) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
}

int edgeFirst(std::uint64_t k) { return static_cast<int>(k >> 32); }
int edgeSecond(std::uint64_t k) { return static_cast<int>(k & 0xffffffffu); }

double px(const std::vector<float>& pts, int i) { return pts[static_cast<std::size_t>(2 * i)]; }
double py(const std::vector<float>& pts, int i) { return pts[static_cast<std::size_t>(2 * i) + 1]; }

/// 有向边 (a,b) 与点 p 的叉积（×2 面积）；>0 表示 p 在 a->b 左侧
double crossOf(const std::vector<float>& pts, int a, int b, double x, double y) {
    return (px(pts, b) - px(pts, a)) * (y - py(pts, a)) - (py(pts, b) - py(pts, a)) * (x - px(pts, a));
}

}  // namespace

std::vector<VertexLoop> boundaryLoops(const std::vector<int>& tris, int vertexCount, int* nonManifold) {
    std::vector<VertexLoop> loops;
    if (tris.size() < 3 || vertexCount <= 0) return loops;

    std::unordered_map<std::uint64_t, int> dirCount;
    dirCount.reserve(tris.size());
    for (std::size_t t = 0; t + 2 < tris.size(); t += 3) {
        const int v[3] = {tris[t], tris[t + 1], tris[t + 2]};
        for (int k = 0; k < 3; ++k) {
            if (v[k] < 0 || v[k] >= vertexCount || v[(k + 1) % 3] < 0 || v[(k + 1) % 3] >= vertexCount) continue;
            ++dirCount[packEdge(v[k], v[(k + 1) % 3])];
        }
    }

    // 同向边出现多次 = 重叠三角形 / 绕序不一致
    int nm = 0;
    for (const auto& kv : dirCount) {
        if (kv.second > 1) ++nm;
    }

    // 边界边：该方向有、反方向没有（成对抵消）
    std::unordered_map<int, std::vector<int>> nextOf;
    std::vector<std::uint64_t> boundary;
    boundary.reserve(dirCount.size());
    for (const auto& kv : dirCount) {
        if (kv.second <= 0) continue;
        const int a = edgeFirst(kv.first), b = edgeSecond(kv.first);
        const auto rev = dirCount.find(packEdge(b, a));
        const int revCount = (rev == dirCount.end()) ? 0 : rev->second;
        if (kv.second - revCount <= 0) continue;
        boundary.push_back(kv.first);
        nextOf[a].push_back(b);
    }
    for (auto& kv : nextOf) {
        if (kv.second.size() > 1) ++nm;   // 顶点有两条出边 = 边界自交或非流形
    }

    std::unordered_set<std::uint64_t> visited;
    visited.reserve(boundary.size() * 2);
    for (const std::uint64_t start : boundary) {
        if (visited.count(start)) continue;
        VertexLoop loop;
        std::uint64_t cur = start;
        const int loopStart = edgeFirst(start);
        for (std::size_t guard = 0; guard <= boundary.size(); ++guard) {
            visited.insert(cur);
            const int a = edgeFirst(cur), b = edgeSecond(cur);
            loop.push_back(a);
            const auto it = nextOf.find(b);
            if (it == nextOf.end() || it->second.empty()) {   // 断了：不是闭合环
                ++nm;
                break;
            }
            const int nxt = it->second.front();
            if (nxt == loopStart) {                            // 闭合
                visited.insert(packEdge(b, nxt));               // 最后闭合边同样属于已走完的环
                loop.push_back(b);                             // ★ 收尾顶点 b 必须补进来
                break;
            }
            const std::uint64_t nextEdge = packEdge(b, nxt);
            if (visited.count(nextEdge)) {                     // 环里绕进了已访问边
                loop.push_back(b);
                ++nm;
                break;
            }
            cur = nextEdge;
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }

    if (nonManifold) *nonManifold += nm;
    return loops;
}

bool pointInLoop(const std::vector<float>& pts, const VertexLoop& loop, float x, float y) {
    const std::size_t n = loop.size();
    if (n < 3) return false;
    bool inside = false;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        const int a = loop[j], b = loop[i];
        const double ax = px(pts, a), ay = py(pts, a);
        const double bx = px(pts, b), by = py(pts, b);
        // 边上（含端点）：判据按边长的相对容差，和 DelaunayTri::inside 同口径
        const double ex = bx - ax, ey = by - ay;
        const double len = std::hypot(ex, ey);
        if (len > 0.0) {
            const double cr = ex * (static_cast<double>(y) - ay) - ey * (static_cast<double>(x) - ax);
            if (std::abs(cr) <= 1.0e-6 * len) {
                const double dot = (static_cast<double>(x) - ax) * ex + (static_cast<double>(y) - ay) * ey;
                if (dot >= -1.0e-9 * len * len && dot <= len * len * (1.0 + 1.0e-9)) return true;
            }
        }
        // 射线法
        if (((ay > y) != (by > y))) {
            const double xi = ax + (static_cast<double>(y) - ay) * ex / ey;
            if (xi > static_cast<double>(x)) inside = !inside;
        }
    }
    return inside;
}

double loopArea(const std::vector<float>& pts, const VertexLoop& loop) {
    double sum = 0.0;
    const std::size_t n = loop.size();
    for (std::size_t i = 0; i < n; ++i) {
        const int a = loop[i], b = loop[(i + 1) % n];
        sum += px(pts, a) * py(pts, b) - px(pts, b) * py(pts, a);
    }
    return sum * 0.5;
}

double trianglesArea(const std::vector<float>& pts, const std::vector<int>& tris) {
    double sum = 0.0;
    for (std::size_t t = 0; t + 2 < tris.size(); t += 3) {
        const int a = tris[t], b = tris[t + 1], c = tris[t + 2];
        sum += (px(pts, b) - px(pts, a)) * (py(pts, c) - py(pts, a)) -
               (px(pts, c) - px(pts, a)) * (py(pts, b) - py(pts, a));
    }
    return sum * 0.5;
}

MeshTopologyResult rebuildMeshTriangles(const MeshTopologyInput& in, const DelaunayOptions& opt) {
    MeshTopologyResult out;
    const int n = static_cast<int>(in.pts.size() / 2);
    if (n < 3 || in.pts.size() % 2 != 0) {
        out.ok = false;
        out.note = "顶点不足 3 个";
        return out;
    }

    // ★ 单轮廓走 polygon-first（耳切轮廓 + 内点插入），多轮廓暂时仍走「全点 DT + 重心裁剪」。
    //   单轮廓是真实 ArtMesh 的绝对主流；先对全部点做 DT 再裁剪在**近乎共线**的边界顶点上
    //   会失效（弦把发丝宽的凹口填平，轮廓少一条边），耳切从构造上就保证轮廓边全在。
    std::vector<std::pair<int, int>> polygonEdges = in.fixedEdges;
    const DelaunayResult dr =
        in.outlines.size() == 1
            ? delaunayTriangulatePolygon(in.pts, in.outlines[0], polygonEdges, opt)
            : delaunayTriangulate(in.pts, polygonEdges, opt);
    out.flips = dr.flips;
    if (!dr.ok) {
        out.ok = false;
        out.note = dr.note ? dr.note : "Delaunay 失败";
        return out;
    }

    std::vector<int> tris = dr.tris;
    if (!in.outlines.empty() && in.outlines.size() > 1) {
        // 轮廓边一条都不能丢：Delaunay 已把它们当固定边，这里再做一次保险
        std::unordered_set<std::uint64_t> outlineEdges;
        for (const VertexLoop& loop : in.outlines) {
            const std::size_t m = loop.size();
            for (std::size_t i = 0; i < m; ++i) {
                const int a = loop[i], b = loop[(i + 1) % m];
                outlineEdges.insert(packEdge(a, b));
                outlineEdges.insert(packEdge(b, a));
            }
        }
        // ★ 判据是「重心在轮廓内」，没有别的豁免。
        //   曾经给「有任意一条边是轮廓边」的三角形开过保留口子，本意是救薄片三角形，
        //   结果凹口三角形恰好有两条边是轮廓边（L 形的 (100,100)-(200,100)-(100,200)
        //   就是），于是整块凹口被留下、面积多出 5000（实测 35000 vs 30000）。
        //   几何上这个判据是严格的：完全落在简单多边形内部的三角形，重心必在内部；
        //   反之跨出凹口的三角形重心必在外部。留口子只会放走应该裁的。
        (void)outlineEdges;
        std::vector<int> kept;
        kept.reserve(tris.size());
        for (std::size_t t = 0; t + 2 < tris.size(); t += 3) {
            const int a = tris[t], b = tris[t + 1], c = tris[t + 2];
            const double cx = (px(in.pts, a) + px(in.pts, b) + px(in.pts, c)) / 3.0;
            const double cy = (py(in.pts, a) + py(in.pts, b) + py(in.pts, c)) / 3.0;
            bool keep = false;
            for (const VertexLoop& loop : in.outlines) {
                if (pointInLoop(in.pts, loop, static_cast<float>(cx), static_cast<float>(cy))) {
                    keep = true;
                    break;
                }
            }
            if (keep) {
                kept.push_back(a);
                kept.push_back(b);
                kept.push_back(c);
            } else {
                ++out.outsideRemoved;
            }
        }
        tris.swap(kept);
    }

    // ★ 绕序必须跟输入轮廓一致：cmo3 里被镜像过的部件（左右眼、左右眉）官方 indices 是 CW 的，
    //   实测 rabbit 511 个可对拍网格里有 10 个（mesh#82/#87/#89/#145/#162/#172/#183/#228/#423…）
    //   面积符号为负。我们的三角化内核一律产出 CCW，如果不在这里整体翻回来，
    //   重建结果就与原网格绕序相反 —— 对有背面剔除的渲染路径是肉眼可见的错，
    //   对「重算三角形」这种原地替换操作更是把数据改坏了。
    //   判据取全部轮廓环的有符号面积之和：轮廓环的走向继承自原 indices，符号就是原网格的绕序。
    double loopSign = 0.0;
    for (const VertexLoop& loop : in.outlines) loopSign += loopArea(in.pts, loop);
    if (loopSign < 0.0) {
        for (std::size_t t = 0; t + 2 < tris.size(); t += 3) std::swap(tris[t + 1], tris[t + 2]);
    }

    int nm = 0;
    out.boundary = boundaryLoops(tris, n, &nm);
    out.tris = std::move(tris);
    return out;
}

} // namespace l2d
