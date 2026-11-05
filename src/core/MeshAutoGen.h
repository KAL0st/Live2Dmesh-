#pragma once
//
// MeshAutoGen —— 自动网格的本地实现：贴图 alpha 轮廓、边界采样和约束 Delaunay。
// 参照官方点间距与细条中心线分支；轮廓/覆盖和连边是本地实现，并非完整等价。
//
#include "core/Cmo3Document.h"
#include "core/XsParser.h"
#include "core/MeshTextureFeatures.h"
#include <string>
#include <array>
#include <vector>

namespace l2d {

enum class MeshAutoGenStyle { Automatic, Soft, Rigid, FineLine };
enum class MeshAutoGenMaterial { Unspecified, Soft, Rigid, ByName };

struct MeshNameMode {
    MeshAutoGenStyle style = MeshAutoGenStyle::Automatic;
    std::string reason;
};

struct MeshAutoGenParams {
    MeshAutoGenMaterial material = MeshAutoGenMaterial::Unspecified; ///< explicit intent overrides style; never inferred from RGB straightness
    bool autoSymmetry = true; ///< local alpha-silhouette detection and paired sampling
    int polygonOuterDensity = 80; ///< 外侧点间距；值越大越稀疏（不是点数）
    int polygonInnerDensity = 80; ///< 锁边以内点间距；0 只关闭内部自动补点
    int polygonBorderDensity = 0; ///< local support-row interval; 0 follows the outer setting
    bool generateBorder = true; ///< local support-row switch, independent of interior fill
    bool generateHoles = false; ///< retain enclosed fully transparent regions as mesh holes
    int holeMinArea = 16; ///< minimum enclosed transparent area in source pixels squared
    int holePointSpacing = 80;
    int holeBoundaryMargin = 1; ///< retreat into transparency to avoid pixel stair-step anchors; 0 exact
    bool generateHoleBorder = true;
    int holeBorderSpacing = 0; ///< 0 follows holePointSpacing
    int holeBorderDistance = 2; ///< offset away from transparent hole, into opaque material
    bool staggerBorder = true; ///< local validated midpoint-offset alignment; favors isosceles band triangles
    bool regularInterior = true; ///< fixed crease/ring constraints, then near-equilateral interior refinement
    int interiorTriangleTarget = 0; ///< approximate interior face budget; 0 uses spacing controls
    int interiorBorderMargin = 0; ///< local interior-to-support clearance in source pixels; 0 automatic
    int polygonFeatureDensity = 80; ///< local automatic crease simplification interval, independent of both rows
    int polygonMargin = 14; ///< 本地 alpha 圆形扩边半径（源图像素），可超出 PNG 矩形
    int polygonInnerMargin = 14; ///< 纹理内部锁边环的优先内缩距离，窄处自适应缩小
    int polygonMinMargin = 5; ///< minimum free-fill clearance from the outer contour
    int polygonMinBoundsPt = 15;
    int thresholdAlpha = 0;
    MeshAutoGenStyle style = MeshAutoGenStyle::Automatic; ///< local generation choice, not an official CMO3 field
};

struct MeshAutoGenResult {
    MeshAutoGenStyle effectiveStyle = MeshAutoGenStyle::Automatic;
    std::string nameModeReason;
    bool symmetryDetected = false, symmetrySampling = false;
    double symmetryScore = 0;
    std::array<double,4> symmetryAxis{0,0,0,1}; ///< source pixels: center x/y, axis dx/dy
    std::string symmetryFallback;
    int oldVertices = 0;
    int newVertices = 0;
    int oldTriangles = 0;
    int newTriangles = 0;
    bool changed = false;
    int featureEdges = 0;
    int boundaryVertices = 0;
    int supportVertices = 0;
    int borderVertices = 0, borderEdges = 0, borderRings = 0; ///< 连续内侧支撑环；复用的褶皱点不计新增 supportVertices
    int extrapolatedVertices = 0;
    int components = 0;
    int holeContours = 0, holeVertices = 0, holeBorderRings = 0;
    int mergedHoleVertices = 0; ///< adjacent hole stations removed before support construction
    size_t opaquePixels = 0, coveredPixels = 0;
    bool conservativeFallback = false;
    bool stripSampling = false;
    int flowVertices = 0, flowEdges = 0, crossSections = 0; ///< fixed longitudinal / transverse structure, before free gap fill
    int sharpTips = 0, mergedConstraintVertices = 0; ///< local simplification diagnostics
    int rigidBorderRemoved = 0; ///< redundant near-collinear support stations removed for rigid material
    int staggeredBorderVertices = 0; ///< validated support stations shifted toward outer-edge midpoint offsets
    int priorityMergedVertices = 0; ///< non-outline/non-feature vertices collapsed by role
    int mergedBorderVertices = 0; ///< support-ring stations removed by validated short-edge collapse
    int mergedDerivedFeatureVertices = 0; ///< collinear ring/crease intersections removed without moving original crease anchors
    int refinedInteriorVertices = 0;
    int interiorTriangles = 0;
    int longEdgeVertices = 0;
    int qualityRemovedVertices = 0, qualityRefinedVertices = 0; ///< accepted sliver repairs, excluding the outer lock band
};

struct MeshFeaturePreview {
    std::vector<uint8_t> rgba;
    int width = 0, height = 0;
    float u0 = 0.f, v0 = 0.f, du = 0.f, dv = 0.f;
    std::array<float,6> uvToPixel{1,0,0,0,1,0};
    bool premultiplied = false;
    bool limited = false;
    int mergedNearbyLines = 0;
    bool conservativeFallback = false;
    std::vector<std::vector<float>> contours; // closed outer boundaries, logical UVs
    std::vector<std::vector<float>> holes, holeBorders; // logical UVs
    std::vector<MeshFeatureLine> lines; // logical texture UVs, matching the mesh
};

class MeshAutoGen {
public:
    static MeshNameMode modeFromName(const std::string& name);
    static MeshAutoGenParams readParams(const Cmo3Document& doc, int meshIdx);
    static bool textureFeatures(const Cmo3Document& doc, int meshIdx,
                               const MeshFeatureSettings& settings,
                               MeshFeaturePreview& preview, std::string& error,
                               const MeshAutoGenParams* params = nullptr);
    /// 只改 CMeshGeneratorExtension/MeshGenerateSetting 的 7 字段；其它字段不动。
    static bool writeParams(XsDocument& work, const Cmo3Document& doc, int meshIdx,
                            const MeshAutoGenParams& p, std::string& error);

    /// Fitted cubic alpha contour with bounded-error curvature-adaptive sampling,
    /// a constrained inner border ring and crease anchors. Single-stroke strips
    /// use width-aware cross sections and longitudinal rails before sparse gap
    /// fill. Rigid / fine-line styles suppress unnecessary 2-D fill. Both border
    /// rows and the structure graph remain fixed during edge relaxation.
    /// UV/all forms use the old triangle field, with nearest-triangle signed
    /// extrapolation across the full source image and its outer margin. The
    /// contour is not clipped to the source PNG rectangle. Unused old points are removed;
    /// glue-bound UIDs remain protected. Every nonzero-alpha pixel is verified.
    /// 写回 source positions/uvs/indices + GEditableMesh2 + 所有 CArtMeshForm.positions。
    static bool generate(XsDocument& work, const Cmo3Document& doc, int meshIdx,
                         const MeshAutoGenParams& p, MeshAutoGenResult& result,
                         std::string& error,
                         const std::vector<MeshFeatureLine>* features = nullptr);
};

}  // namespace l2d
