#include "edit/DelaunayTri.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace l2d {
namespace {

struct Tri {
    int v[3];
    bool alive = true;
};

/// 坐标量化键：这里**不能**复用 edgeKey —— 它会把两个参数排序，把 (x,y) 和 (y,x) 折成同一个键，
/// 于是 (100,0) 与 (0,100) 被判成重复点（实测 5 点正方形直接少一个顶点、固定边恢复失败）。
inline uint64_t packXY(int x, int y) {
    return (static_cast<uint64_t>(static_cast<uint32_t>(x)) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(y));
}

inline uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return packXY(a, b);
}

/// 几何内核。点云一律用 double 存：super 三角形顶点离点云有三个数量级，
/// 用 float 的话 (ax-dx) 这类差值的舍入会被 4 次方判据放大到和结论同量级。
struct Geom {
    const std::vector<double>& p;
    double eps;

    const double* pt(int i) const { return &p[static_cast<std::size_t>(i) * 2]; }

    /// 有向面积 ×2：>0 = CCW
    double orient(int a, int b, int c) const {
        const double* A = pt(a);
        const double* B = pt(b);
        const double* C = pt(c);
        return (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    }

    /// 「d 是否在 abc 的外接圆内」，返回 半径 - 距离（> eps 就是在圆内）。
    /// 与官方同口径：`c.radius - dist > EPS`（EPS = com.live2d.util.L.a.f() = 1.0E-6f）。
    /// 不用 4×4 行列式判据的原因见文件头注释。
    double circleTest(int a, int b, int c, int d) const {
        const double* A = pt(a);
        const double* B = pt(b);
        const double* C = pt(c);
        const double* D = pt(d);
        const double den = 2.0 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]));
        if (std::abs(den) < 1e-300) return -1.0e300;   // 三点共线：外接圆在无穷远，永远不算「圆内」
        const double a2 = A[0] * A[0] + A[1] * A[1];
        const double b2 = B[0] * B[0] + B[1] * B[1];
        const double c2 = C[0] * C[0] + C[1] * C[1];
        const double ux = (a2 * (B[1] - C[1]) + b2 * (C[1] - A[1]) + c2 * (A[1] - B[1])) / den;
        const double uy = (a2 * (C[0] - B[0]) + b2 * (A[0] - C[0]) + c2 * (B[0] - A[0])) / den;
        return std::hypot(A[0] - ux, A[1] - uy) - std::hypot(D[0] - ux, D[1] - uy);
    }

    /// 三角形尺度（长度^2）：给 orient / inside 的容差做无量纲归一
    double areaScale(int a, int b, int c) const {
        const double* A = pt(a);
        const double* B = pt(b);
        const double* C = pt(c);
        const double ux = B[0] - A[0], uy = B[1] - A[1];
        const double vx = C[0] - A[0], vy = C[1] - A[1];
        return std::sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy));
    }

    /// 两点距离（把「面积容差」换算成「点到边的距离容差」）
    double dist(int a, int b) const {
        const double* A = pt(a);
        const double* B = pt(b);
        return std::hypot(B[0] - A[0], B[1] - A[1]);
    }

    /// d 是否落在三角形 abc 内（含边）。
    ///   ★ 容差必须按边长归一，不能按三角形面积：orient(A,B,pi) 等于 |AB| x 点到直线的距离，
    ///     抖动网格边长可达 1e2~1e5，用 eps*|AB|*|AC| 当容差等于允许点离线 eps*|AC| 单位还算「在内」，
    ///     于是同一个点同时命中两个不相邻的薄片三角形，洞分裂成两环（实测 300 随机点 p13）。
    bool inside(int a, int b, int c, int d) const {
        const double area = orient(a, b, c);
        if (std::abs(area) <= 1e-300) return false;
        double d1 = orient(a, b, d);
        double d2 = orient(b, c, d);
        double d3 = orient(c, a, d);
        if (area < 0) {
            d1 = -d1;
            d2 = -d2;
            d3 = -d3;
        }
        return d1 >= -eps * dist(a, b) && d2 >= -eps * dist(b, c) && d3 >= -eps * dist(c, a);
    }

    /// 两线段严格相交（端点相触不算）
    bool crossStrict(int a, int b, int c, int d) const {
        const double o1 = orient(c, d, a);
        const double o2 = orient(c, d, b);
        const double o3 = orient(a, b, c);
        const double o4 = orient(a, b, d);
        return ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) && ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0));
    }
};

inline int otherVertex(const Tri& t, int a, int b) {
    for (int k = 0; k < 3; ++k)
        if (t.v[k] != a && t.v[k] != b) return t.v[k];
    return -1;
}

/// 三角形是否含无向边 a-b
inline bool hasUndirectedEdge(const Tri& t, int a, int b) {
    for (int k = 0; k < 3; ++k) {
        const int x = t.v[k], y = t.v[(k + 1) % 3];
        if ((x == a && y == b) || (x == b && y == a)) return true;
    }
    return false;
}

inline Tri makeTriCcw(const Geom& g, int a, int b, int c) {
    Tri t;
    if (g.orient(a, b, c) >= 0) {
        t.v[0] = a;
        t.v[1] = b;
        t.v[2] = c;
    } else {
        t.v[0] = a;
        t.v[1] = c;
        t.v[2] = b;
    }
    return t;
}

/// 边 → 相邻三角形（最多 2 个；-1 空）
using EdgeMap = std::unordered_map<uint64_t, std::array<int, 2>>;

void buildEdgeMap(const std::vector<Tri>& tris, EdgeMap& map) {
    map.clear();
    for (int i = 0; i < static_cast<int>(tris.size()); ++i) {
        const Tri& t = tris[static_cast<std::size_t>(i)];
        if (!t.alive) continue;
        for (int k = 0; k < 3; ++k) {
            const uint64_t key = edgeKey(t.v[k], t.v[(k + 1) % 3]);
            auto it = map.find(key);
            if (it == map.end()) {
                map.emplace(key, std::array<int, 2> {i, -1});
            } else if (it->second[1] < 0) {
                it->second[1] = i;
            }
            // >2：网格已非流形，保留前两个即可（后续翻转会把它修正回来）
        }
    }
}

/// 全局翻转一次 a-b 边（两个相邻三角形）。成功返回 true。
bool flipEdge(std::vector<Tri>& tris, const Geom& g, int ta, int tb, int a, int b, double eps,
              int& flips, int maxFlips) {
    if (ta < 0 || tb < 0 || ta == tb) return false;
    if (flips >= maxFlips) return false;
    Tri& A = tris[static_cast<std::size_t>(ta)];
    Tri& B = tris[static_cast<std::size_t>(tb)];
    if (!A.alive || !B.alive) return false;
    // ★ 防御：两个三角形必须真的含边 a-b。调用方可能拿着「翻转前建好的边表」来翻，
    //   此时下标指向的三角形已经被改写过 —— otherVertex 会返回一个毫不相干的顶点，
    //   翻出来的三角形与原四边形毫无关系（实测 mesh#79 面积凭空多 5.8%、边界裂成 4 环）。
    if (!hasUndirectedEdge(A, a, b) || !hasUndirectedEdge(B, a, b)) return false;
    const int c = otherVertex(A, a, b);
    const int d = otherVertex(B, a, b);
    if (c < 0 || d < 0 || c == d) return false;
    // 四边形必须是凸的（c-d 与 a-b 严格相交），否则翻出来会重叠。
    //   ★ 这里只要「严格相交」，**不能**再叠加绝对面积门槛（原来是
    //     |orient| <= eps*areaScale 就拒绝）。真实网格里固定边经常「几乎穿过」某个顶点：
    //     梦寝兔兔 mesh#0 的固定边 16-20 离顶点 0 只有 0.0098，而门槛是 0.094 ⇒
    //     每一次能恢复它的翻转都被挡回去，最终 96 次翻转后报「找不到可翻转的相交边」。
    //     crossStrict 已经保证 a、b 分居 c-d 两侧（符号相反），面积再小也是合法翻转，
    //     翻出来的薄片三角形会在最终输出时按退化丢掉。
    if (!g.crossStrict(a, b, c, d)) return false;
    A = makeTriCcw(g, a, c, d);
    B = makeTriCcw(g, b, d, c);
    ++flips;
    return true;
}

/// 网格中是否已有边 a-b
bool hasEdge(const std::vector<Tri>& tris, int a, int b) {
    for (const Tri& t : tris) {
        if (!t.alive) continue;
        for (int k = 0; k < 3; ++k) {
            const int x = t.v[k], y = t.v[(k + 1) % 3];
            if ((x == a && y == b) || (x == b && y == a)) return true;
        }
    }
    return false;
}


/// 耳切：把轮廓环三角化（polygon-first 的第一步）。
///   · 环的绕序先归一化成 CCW，切出来的三角形一律 makeTriCcw；
///   · 凹点的耳朵（orient <= 0）不切，除非剩下的环里根本没有凸耳 —— 那时按「归一化凸性」最强的顶点强切一刀，
///     全是凹的说明环自交或面积为 0，返回 false 让调用方报失败；
///   · 共线的顶点会被相邻的耳朵一路带着走，最后落进某个三角形里，不会凭空丢掉。
bool earClipLoop(const Geom& g, int n, const std::vector<int>& loop, double eps, std::vector<Tri>& out) {
    const int m = static_cast<int>(loop.size());
    if (m < 3) return false;
    std::vector<int> check = loop;
    for (const int v : check) {
        if (v < 0 || v >= n) return false;
    }
    std::sort(check.begin(), check.end());
    if (std::adjacent_find(check.begin(), check.end()) != check.end()) return false;   // 重复顶点

    std::vector<int> poly = loop;
    double a2 = 0.0;
    for (int i = 0; i < m; ++i) {
        const int a = poly[static_cast<std::size_t>(i)];
        const int b = poly[static_cast<std::size_t>((i + 1) % m)];
        a2 += g.pt(a)[0] * g.pt(b)[1] - g.pt(b)[0] * g.pt(a)[1];
    }
    if (a2 < 0.0) std::reverse(poly.begin(), poly.end());

    int guard = 0;
    const int guardMax = m * m + 64;
    while (poly.size() > 3) {
        if (++guard > guardMax) return false;
        const int sz = static_cast<int>(poly.size());
        int clipped = -1;
        for (int i = 0; i < sz && clipped < 0; ++i) {
            const int a = poly[static_cast<std::size_t>((i + sz - 1) % sz)];
            const int c = poly[static_cast<std::size_t>(i)];
            const int b = poly[static_cast<std::size_t>((i + 1) % sz)];
            if (g.orient(a, c, b) <= 0.0) continue;   // 凹点或共线点不是耳朵
            bool blocked = false;
            for (int k = 0; k < sz && !blocked; ++k) {
                const int v = poly[static_cast<std::size_t>(k)];
                if (v == a || v == c || v == b) continue;
                if (g.inside(a, c, b, v)) blocked = true;   // 耳朵里有别的顶点
            }
            if (!blocked) clipped = i;
        }
        if (clipped < 0) {
            // 没有凸耳：环里全是凹点/共线点 ⇒ 自交或退化。挑「归一化凸性」最强的切一刀，
            // 若它仍然是负的（真凹），这个环不可能被耳切成简单多边形，直接失败。
            double best = -1.0e300;
            int bestI = -1;
            for (int i = 0; i < sz; ++i) {
                const int a = poly[static_cast<std::size_t>((i + sz - 1) % sz)];
                const int c = poly[static_cast<std::size_t>(i)];
                const int b = poly[static_cast<std::size_t>((i + 1) % sz)];
                const double sc = g.areaScale(a, c, b);
                const double val = sc > 0.0 ? g.orient(a, c, b) / sc : -1.0e300;
                if (val > best) {
                    best = val;
                    bestI = i;
                }
            }
            if (bestI < 0 || best <= -eps) return false;
            clipped = bestI;
        }
        const int sz2 = static_cast<int>(poly.size());
        const int a = poly[static_cast<std::size_t>((clipped + sz2 - 1) % sz2)];
        const int c = poly[static_cast<std::size_t>(clipped)];
        const int b = poly[static_cast<std::size_t>((clipped + 1) % sz2)];
        out.push_back(makeTriCcw(g, a, c, b));
        poly.erase(poly.begin() + clipped);
    }
    out.push_back(makeTriCcw(g, poly[0], poly[1], poly[2]));
    return true;
}

}  // namespace

static DelaunayResult triangulateImpl(const std::vector<float>& pts,
                                      const std::vector<std::pair<int, int>>& fixedEdges,
                                      const std::vector<int>* outline,
                                      const DelaunayOptions& opt) {
    DelaunayResult res;
    const double eps = opt.eps > 0 ? opt.eps : 1.0e-6;
    const bool polygonMode = outline != nullptr;   // true：先耳切轮廓，再插内点（见头文件）
    const int n = static_cast<int>(pts.size() / 2);
    if (n < 3) {
        res.ok = false;
        res.note = "顶点少于 3 个，无法三角化";
        return res;
    }

    // ---- 点云 / 去重 / 种子三角形 ----
    //   ★ 这里**不用**「超级三角形」。抖动网格里近共线的三点外接圆半径可以达到 span^2/抖动量级，
    //     超级顶点一旦落进这些巨圆，S∪{A,B,C} 的 Delaunay 限制回 S 就不再是 DT(S)：
    //     删掉含超级顶点的三角形后，凸包会被啃掉一块（实测抖动网格 #127 处凹进 26.49）。
    //     把超级三角形放得更大只是把问题推给 double 精度，治不了根。
    //     所以改成「只对真实点插入」：内部点走空圆驱逐，凸包外的点直接连到可视凸包边，
    //     三角化覆盖的区域始终严格等于「已插入点集的凸包」。
    std::vector<double> P(static_cast<std::size_t>(n) * 2, 0.0);
    for (int i = 0; i < n; ++i) {
        P[static_cast<std::size_t>(2 * i)] = pts[static_cast<std::size_t>(2 * i)];
        P[static_cast<std::size_t>(2 * i) + 1] = pts[static_cast<std::size_t>(2 * i) + 1];
    }
    Geom g {P, eps};

    std::unordered_set<uint64_t> seen;
    seen.reserve(static_cast<std::size_t>(n) * 2);
    std::vector<int> uniq;
    uniq.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const uint64_t qk = packXY(static_cast<int>(std::floor(pts[static_cast<std::size_t>(2 * i)] / (eps * 1.0e3))),
                                   static_cast<int>(std::floor(pts[static_cast<std::size_t>(2 * i) + 1] / (eps * 1.0e3))));
        if (seen.insert(qk).second) {
            uniq.push_back(i);
        } else {
            ++res.duplicatePoints;
        }
    }
    std::vector<Tri> tris;
    tris.reserve(static_cast<std::size_t>(n) * 2 + 8);
    std::vector<char> inserted(static_cast<std::size_t>(n), 0);
    // polygon 模式：轮廓顶点掩码。输出阶段要保住「三个顶点都在轮廓上」的退化薄片，
    //   否则共线的边界顶点会因为三角形按退化被丢掉而从网格里消失（轮廓少一个顶点）。
    std::vector<char> outlineMask;

    if (polygonMode) {
        // ---- polygon-first：先耳切轮廓，再把内点插进去 ----
        //   ★ 不做「全点 DT + 按重心裁剪」。实测梦寝兔兔 rabbit.saved.cmo3 里 103/511 个
        //     ArtMesh 有近乎共线的边界顶点，弦会把发丝宽的凹口整块填平：面积只差 3e-9（相对），
        //     但那个顶点被推进内部、轮廓少一条边，用户以后拖它轮廓不再跟着动。
        //     官方本来就是直接对轮廓三角化（docs/33 §2.10 的 outlineIndexCW），这里照做。
        if (!earClipLoop(g, n, *outline, eps, tris)) {
            res.ok = false;
            res.note = "轮廓环无法耳切（重复顶点 / 自交 / 面积为 0）";
            return res;
        }
        outlineMask.assign(static_cast<std::size_t>(n), 0);
        for (const int v : *outline) {
            inserted[static_cast<std::size_t>(v)] = 1;
            outlineMask[static_cast<std::size_t>(v)] = 1;
        }
    } else {
        // ---- 通用模式：种子三角形 + 逐点插入（三角化区域随插入扩张到凸包）----
        int ia = uniq.empty() ? -1 : uniq[0];
        int ib = -1, ic = -1;
        for (std::size_t k = 1; k < uniq.size() && ia >= 0; ++k)
            if (g.dist(ia, uniq[k]) > eps) {
                ib = uniq[k];
                break;
            }
        if (ib >= 0) {
            for (std::size_t k = 1; k < uniq.size(); ++k) {
                const int c = uniq[k];
                if (c == ia || c == ib) continue;
                if (std::abs(g.orient(ia, ib, c)) > eps * g.dist(ia, ib) * g.dist(ia, c)) {
                    ic = c;
                    break;
                }
            }
        }
        if (ic < 0) {
            res.ok = false;
            res.note = "顶点全部重合或共线，构不成三角形";
            return res;
        }

        tris.push_back(makeTriCcw(g, ia, ib, ic));
        inserted[static_cast<std::size_t>(ia)] = 1;
        inserted[static_cast<std::size_t>(ib)] = 1;
        inserted[static_cast<std::size_t>(ic)] = 1;

    }
    // ---- 逐点插入（Bowyer-Watson：空圆驱逐 + 连通扩张 + 边界重连）----
    //   ★ 坏三角形集合必须「连通、且以 pi 为核」。若把外接圆含 pi 的三角形一股脑删掉，
    //     共圆点与含 super 顶点的薄壳三角形会形成多个洞，重连出的三角形互相重叠
    //     （实测 300 随机点膨胀到 26523 个三角形、总面积 1740 倍）。
    //     所以：种子 = 真正包含 pi 的三角形，之后沿边 BFS 只收外接圆严格含 pi 的邻居。
    // 以 p 为极点的局部合法化（Lawson）。插入后凡是「以 p 为顶点」的三角形，
    // 其底边都可能是非法的；一翻就会翻到新的三角形，所以循环到稳定为止。
    // 凸包外的点走扇形，结果本来就不是 Delaunay，这一步是把它拉回 Delaunay 的关键。
    auto legalizeAround = [&](int p) {
        for (int guard = 0; guard < 200000; ++guard) {
            EdgeMap em;
            buildEdgeMap(tris, em);
            bool flipped = false;
            for (int ti = 0; ti < static_cast<int>(tris.size()) && !flipped; ++ti) {
                const Tri& t0 = tris[static_cast<std::size_t>(ti)];
                if (!t0.alive) continue;
                if (t0.v[0] != p && t0.v[1] != p && t0.v[2] != p) continue;
                const int vv[3] = {t0.v[0], t0.v[1], t0.v[2]};
                for (int k = 0; k < 3; ++k) {
                    const int x = vv[(k + 1) % 3], y = vv[(k + 2) % 3];
                    if (x == p || y == p) continue;
                    const auto it = em.find(edgeKey(x, y));
                    if (it == em.end() || it->second[1] < 0) continue;   // 凸包边：外侧没有三角形
                    const int ta = it->second[0], tb = it->second[1];
                    const int z = otherVertex(tris[static_cast<std::size_t>(ta)], x, y) == p
                                      ? otherVertex(tris[static_cast<std::size_t>(tb)], x, y)
                                      : otherVertex(tris[static_cast<std::size_t>(ta)], x, y);
                    if (z < 0 || z == p) continue;
                    if (g.circleTest(p, x, y, z) <= eps) continue;
                    if (flipEdge(tris, g, ta, tb, x, y, eps, res.flips, opt.maxFlips)) flipped = true;
                    break;
                }
            }
            if (!flipped) break;
        }
    };

    std::vector<char> mark;
    std::vector<int> queue;
    std::unordered_map<uint64_t, int> boundaryCount;
    for (int pi = 0; pi < n; ++pi) {
        if (inserted[static_cast<std::size_t>(pi)]) continue;

        //   ★ 种子只能用「真正包含 pi 的三角形」，**不能**顺手把外接圆含 pi 的一起扫进来：
        //     细长三角形的外接圆半径可以远大于其尺度，半径差为正却不包含 pi，
        //     于是被标成种子；它和真正包含 pi 的那个三角形在邻接图上不连通，洞就成了多个环，
        //     重连出来的三角形互相重叠（实测 300 随机点：p13 起 border 开始大于 bad+2，最终 24117 个三角形）。
        mark.assign(tris.size(), 0);
        queue.clear();
        for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
            const Tri& t = tris[static_cast<std::size_t>(ti)];
            if (!t.alive) continue;
            if (g.inside(t.v[0], t.v[1], t.v[2], pi)) {
                mark[static_cast<std::size_t>(ti)] = 1;
                queue.push_back(ti);
            }
        }
        if (queue.empty() && polygonMode) {
            // polygon 模式：点一定在轮廓内，「找不到包含它的三角形」只可能是容差问题。
            //   ★ 绝不能走下面的「可见凸包边扇形」：那条路会把点连到轮廓边的外侧，
            //     直接把三角形造到多边形外面。这里退化成「拆分最近的三角形」，几何误差有界。
            int best = -1;
            double bestD = 1.0e300;
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
                const Tri& t = tris[static_cast<std::size_t>(ti)];
                if (!t.alive) continue;
                double d = 0.0;
                for (int k = 0; k < 3; ++k)
                    d += std::hypot(g.pt(t.v[k])[0] - g.pt(pi)[0], g.pt(t.v[k])[1] - g.pt(pi)[1]);
                if (d < bestD) {
                    bestD = d;
                    best = ti;
                }
            }
            if (best >= 0) {
                const Tri t = tris[static_cast<std::size_t>(best)];
                tris[static_cast<std::size_t>(best)].alive = false;
                tris.push_back(makeTriCcw(g, t.v[0], t.v[1], pi));
                tris.push_back(makeTriCcw(g, t.v[1], t.v[2], pi));
                tris.push_back(makeTriCcw(g, t.v[2], t.v[0], pi));
                inserted[static_cast<std::size_t>(pi)] = 1;
                legalizeAround(pi);
            }
            continue;
        }
        if (queue.empty()) {
            // pi 不在任何三角形里 ⇒ 它在当前凸包外（或数值上只落在边界薄三角形边上）。
            //   ★ 凸包外的点**不能**用空圆驱逐：它的坏集合未必连通、甚至可能为空，
            //     BFS 从任意种子出发都会把它连到内部而漏掉真正可视的那圈凸包边。
            //   正确做法：把 pi 逐个连到「可视凸包边」（a→b 是有向凸包边，p 在其右侧）。
            std::unordered_set<uint64_t> dirEdges;
            dirEdges.reserve(tris.size() * 3);
            for (const Tri& t : tris) {
                if (!t.alive) continue;
                for (int k = 0; k < 3; ++k) dirEdges.insert(packXY(t.v[k], t.v[(k + 1) % 3]));
            }
            std::vector<std::pair<int, int>> vis;
            for (const Tri& t : tris) {
                if (!t.alive) continue;
                for (int k = 0; k < 3; ++k) {
                    const int a = t.v[k], b = t.v[(k + 1) % 3];
                    if (dirEdges.count(packXY(b, a))) continue;   // 内部边
                    //   ★ 只收「严格可视」的边（orient < 0）。把共线边也算进来会造出零面积扇形：
                    //     网格点按行插入时 p 常落在某条凸包边的延长线上，那一块区域其实由相邻的
                    //     可视边负责覆盖，多算一条只会留下一条假的凸包边。
                    if (g.orient(a, b, pi) < -eps * g.dist(a, b)) vis.emplace_back(a, b);
                }
            }
            std::sort(vis.begin(), vis.end());
            vis.erase(std::unique(vis.begin(), vis.end()), vis.end());
            if (!vis.empty()) {
                for (const auto& e : vis) tris.push_back(makeTriCcw(g, e.first, e.second, pi));
                legalizeAround(pi);
                continue;
            }
            // 数值兜底：pi 落在所有三角形外（只可能出现在边界极薄三角形上）。
            // 取「半径差最大」的那个当种子，保证 BFS 有起点，且坏集合仍然连通。
            int best = -1;
            double bestVal = -1.0e300;
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
                const Tri& t = tris[static_cast<std::size_t>(ti)];
                if (!t.alive) continue;
                const double val = g.circleTest(t.v[0], t.v[1], t.v[2], pi);
                if (val > bestVal) {
                    bestVal = val;
                    best = ti;
                }
            }
            if (best < 0) continue;
            mark[static_cast<std::size_t>(best)] = 1;
            queue.push_back(best);
        }

        // 邻接表：全量路径只用于首次三角化；真正的大网格走「局部重建」，不该走这里。
        EdgeMap em;
        buildEdgeMap(tris, em);
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            const Tri& t = tris[static_cast<std::size_t>(queue[qi])];
            for (int k = 0; k < 3; ++k) {
                const auto it = em.find(edgeKey(t.v[k], t.v[(k + 1) % 3]));
                if (it == em.end()) continue;
                for (int s = 0; s < 2; ++s) {
                    const int nb = it->second[static_cast<std::size_t>(s)];
                    if (nb < 0 || mark[static_cast<std::size_t>(nb)]) continue;
                    const Tri& nbt = tris[static_cast<std::size_t>(nb)];
                    if (!nbt.alive) continue;
                    // 严格圆内（半径差 > eps，与官方 `radius - dist > EPS` 同口径）：
                    // 连通 BFS 保证坏集合是一整块，边界才可能是单环。
                    if (g.circleTest(nbt.v[0], nbt.v[1], nbt.v[2], pi) > eps) {
                        mark[static_cast<std::size_t>(nb)] = 1;
                        queue.push_back(nb);
                    }
                }
            }
        }

        //   ★ 洞「填实」：坏集合连通 ≠ 洞是单环。共圆点或薄壳三角形会让洞中间夹着一个
        //     未标记的三角形（内岛），边界就成了两个环，直接重连必自重叠。
        //     判据：某未标记三角形的三条边都能在坏三角形的边里找到 ⇒ 它被洞完全包住 ⇒ 强制并进洞里。
        //     （实测 300 随机点：p13 bad=4 border=8 就是这种情形；不填实会一路烂到 1232 个三角形。）
        for (int fill = 0; fill < 64; ++fill) {
            std::unordered_set<uint64_t> badEdges;
            badEdges.reserve(queue.size() * 3);
            for (const int ti : queue) {
                const Tri& t = tris[static_cast<std::size_t>(ti)];
                for (int k = 0; k < 3; ++k) badEdges.insert(edgeKey(t.v[k], t.v[(k + 1) % 3]));
            }
            bool grew = false;
            for (int ti = 0; ti < static_cast<int>(tris.size()); ++ti) {
                if (mark[static_cast<std::size_t>(ti)]) continue;
                const Tri& t = tris[static_cast<std::size_t>(ti)];
                if (!t.alive) continue;
                if (badEdges.count(edgeKey(t.v[0], t.v[1])) && badEdges.count(edgeKey(t.v[1], t.v[2])) &&
                    badEdges.count(edgeKey(t.v[2], t.v[0]))) {
                    mark[static_cast<std::size_t>(ti)] = 1;
                    queue.push_back(ti);
                    grew = true;
                }
            }
            if (!grew) break;
        }

        // 边计数：只统计坏三角形的有向边；成对抵消后剩下的就是洞边界
        boundaryCount.clear();
        for (const int ti : queue) {
            const Tri& t = tris[static_cast<std::size_t>(ti)];
            for (int k = 0; k < 3; ++k) {
                const int a = t.v[k], b = t.v[(k + 1) % 3];
                ++boundaryCount[packXY(a, b)];
            }
        }
        for (const int ti : queue) tris[static_cast<std::size_t>(ti)].alive = false;

        std::vector<std::pair<int, int>> border;
        border.reserve(boundaryCount.size());
        for (const auto& kv : boundaryCount) {
            const int a = static_cast<int>(kv.first >> 32);
            const int b = static_cast<int>(kv.first & 0xffffffffu);
            if (boundaryCount.find(packXY(b, a)) == boundaryCount.end()) border.emplace_back(a, b);
        }
        //   ★ 退化三角形（面积≈0）也必须加进网格，**不能**在这里丢。
        //     丢一个就等于在洞里留个缺口：网格从此非流形，后面的 inside 会同时命中重
        //     叠区域，洞分裂成两个环（实测 p13：bad=4 却 border=8，一条五边形环 +
        //     一条薄片三角形环）。退化三角形面积贡献为 0，只在最终输出时过滤。
        for (const auto& e : border) tris.push_back(makeTriCcw(g, e.first, e.second, pi));
        legalizeAround(pi);

    }

    // ---- 固定边恢复：与固定边严格相交的内部边逐条翻转，直到固定边成为网格边 ----
    //   （对应官方 setup_indices_Delaunay 里「固定边必须是三角形边」这一步；
    //     官方用 500 次硬上限，我们改成「无进展即失败」，且失败时明确报出来而不是静默降级。）
    std::unordered_set<uint64_t> recoveredEdges;
    for (const auto& fe : fixedEdges) {
        const int u = fe.first, v = fe.second;
        if (u == v || u < 0 || v < 0 || u >= n || v >= n) continue;
        if (hasEdge(tris, u, v)) { recoveredEdges.insert(edgeKey(u,v)); continue; }
        ++res.recovered;
        bool done = false;
        EdgeMap initialEdges;
        buildEdgeMap(tris,initialEdges);
        std::vector<uint64_t> crossing;
        for(const auto& kv:initialEdges) {
            const int a=int(kv.first>>32),b=int(kv.first&0xffffffffu);
            if(kv.second[1]>=0&&g.crossStrict(u,v,a,b))crossing.push_back(kv.first);
        }
        std::sort(crossing.begin(),crossing.end());
        std::deque<uint64_t> pending(crossing.begin(),crossing.end());
        size_t stalled=0;
        // FIFO is essential: a new diagonal can still cross the requested
        // segment. Re-scanning an unordered edge map could immediately flip
        // that diagonal back, cycling forever on a closed internal contour.
        for (int guard = 0; guard < 100000 && !done && !pending.empty(); ++guard) {
            const auto key=pending.front();pending.pop_front();
            EdgeMap em;
            buildEdgeMap(tris, em);
            const auto found=em.find(key);
            if(found==em.end()||found->second[1]<0)continue;
            const int a=int(key>>32),b=int(key&0xffffffffu);
            if(!g.crossStrict(u,v,a,b))continue;
            if(recoveredEdges.count(key)) {
                res.ok=false;res.note="固定边互相交叉，需先拆分交点";return res;
            }
            const int ta=found->second[0],tb=found->second[1];
            const int c=otherVertex(tris[size_t(ta)],a,b),d=otherVertex(tris[size_t(tb)],a,b);
            if(flipEdge(tris,g,ta,tb,a,b,eps,res.flips,opt.maxFlips)) {
                stalled=0;
                if(g.crossStrict(u,v,c,d))pending.push_back(edgeKey(c,d));
            } else {
                pending.push_back(key);
                if(++stalled>=pending.size())break;
            }
            done = hasEdge(tris, u, v);
        }
        recoveredEdges.insert(edgeKey(u,v));
        if (!done) {
            res.ok = false;
            res.note = "固定边恢复失败：超过翻转上限";
            return res;
        }
    }

    // ---- 合法化（LOP）：固定边不动，其余内部边按空圆判据翻到稳定 ----
    std::unordered_set<uint64_t> fixedSet;
    fixedSet.reserve(fixedEdges.size() * 2);
    for (const auto& fe : fixedEdges) fixedSet.insert(edgeKey(fe.first, fe.second));

    //   ★ 翻转会当场改写三角形，而边表是「扫描开始时」的快照。所以一次扫描最多翻一条，
    //     翻完立刻重建边表再扫 —— 否则后面的条目指向的三角形已经变了（实测 mesh#79
    //     第 2 条翻转用了第 1 条翻转前的下标，四边形根本不存在，翻出重叠三角形）。
    bool changed = true;
    int rounds = 0;
    while (changed && rounds < 200000) {
        changed = false;
        ++rounds;
        EdgeMap em;
        buildEdgeMap(tris, em);
        for (const auto& kv : em) {
            if (kv.second[1] < 0) continue;
            const int a = static_cast<int>(kv.first >> 32);
            const int b = static_cast<int>(kv.first & 0xffffffffu);
            if (fixedSet.count(kv.first)) continue;   // 固定边永不翻
            const int ta = kv.second[0], tb = kv.second[1];
            const Tri& taT = tris[static_cast<std::size_t>(ta)];
            const int c = otherVertex(taT, a, b);
            const int d = otherVertex(tris[static_cast<std::size_t>(tb)], a, b);
            if (c < 0 || d < 0 || c == d) continue;
            // circleTest 要求前三点的顺序与三角形一致（a<b 只是 edgeKey 的排序）；
            // 顺着相反方向判会把「不满足空圆」的边当成可翻，翻到最后网格自重叠。
            bool forward = false;
            for (int k = 0; k < 3; ++k)
                if (taT.v[k] == a && taT.v[(k + 1) % 3] == b) forward = true;
            if ((forward ? g.circleTest(a, b, c, d) : g.circleTest(b, a, c, d)) <= eps) continue;
            if (flipEdge(tris, g, ta, tb, a, b, eps, res.flips, opt.maxFlips)) {
                changed = true;
                break;   // 边表已过期，回外层重建
            }
        }
    }

    // ---- 输出：统一 CCW、丢退化 ----
    res.tris.clear();
    res.tris.reserve(tris.size() * 3);
    for (const Tri& t : tris) {
        if (!t.alive) continue;
        if (t.v[0] >= n || t.v[1] >= n || t.v[2] >= n) continue;
        if (t.v[0] < 0 || t.v[1] < 0 || t.v[2] < 0) continue;
        const double a2 = g.orient(t.v[0], t.v[1], t.v[2]);
        if (std::abs(a2) <= eps * g.areaScale(t.v[0], t.v[1], t.v[2])) {
            // polygon 模式：三个顶点都在轮廓上的退化薄片要留住，否则共线边界顶点会消失。
            const bool onOutline = !outlineMask.empty() &&
                                   outlineMask[static_cast<std::size_t>(t.v[0])] != 0 &&
                                   outlineMask[static_cast<std::size_t>(t.v[1])] != 0 &&
                                   outlineMask[static_cast<std::size_t>(t.v[2])] != 0;
            if (!onOutline) {
                ++res.degenerateDropped;
                continue;
            }
        }
        res.tris.push_back(t.v[0]);
        res.tris.push_back(t.v[1]);
        res.tris.push_back(t.v[2]);
    }
    return res;
}

DelaunayResult delaunayTriangulate(const std::vector<float>& pts,
                                   const std::vector<std::pair<int, int>>& fixedEdges,
                                   const DelaunayOptions& opt) {
    return triangulateImpl(pts, fixedEdges, nullptr, opt);
}

DelaunayResult delaunayTriangulatePolygon(const std::vector<float>& pts,
                                          const std::vector<int>& outline,
                                          const std::vector<std::pair<int, int>>& fixedEdges,
                                          const DelaunayOptions& opt) {
    return triangulateImpl(pts, fixedEdges, &outline, opt);
}

std::vector<std::pair<int, int>> edgesOfTriangles(const std::vector<int>& tris) {
    std::vector<std::pair<int, int>> edges;
    edges.reserve(tris.size());
    for (std::size_t i = 0; i + 2 < tris.size(); i += 3) {
        int v[3] = {tris[i], tris[i + 1], tris[i + 2]};
        for (int k = 0; k < 3; ++k) {
            int a = v[k], b = v[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            edges.emplace_back(a, b);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

}  // namespace l2d
