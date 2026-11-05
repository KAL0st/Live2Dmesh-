#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace l2d {

struct MeshFeatureLine {
    // Pixel coordinates in detector output; logical UVs in MeshFeaturePreview.
    std::vector<float> points;
    bool enabled = true;
    float length = 0.f;
    bool centerline = false; // paired narrow crease edges, rather than a color boundary
    bool automatic = false; // detector output may be simplified within raster tolerance; explicit paths stay exact
    bool weakReference = false; // bright-region candidate: sparse guidance, not a fixed crease
    bool curved = false; // bounded curve fit, sampled for straight mesh edges
    bool boundaryReference = false; // automatic outline-following paint edge, distinct from a highlight
};

struct MeshFeatureSettings {
    int contrast = 24;
    int minLength = 16;  // stored-image pixels
    int cornerTolerance = 2; // maximum deviation in stored-image pixels; straight spans keep endpoints
    int maxFoldWidth = 12;   // merge the two sides of a narrow crease into a centerline
    int nearbyMergeDistance = 12; // source pixels; fit adjacent similarly directed color boundaries, 0 disables
    int alphaThreshold = 0;
    bool premultiplied = false;
    bool weakHighlights = true;
    float curveSmoothing = 2.f; // source-pixel fitting corridor; 0 retains detected polylines
    float curveSpacing = 16.f; // arc-length sampling interval, set by mesh feature spacing
    float curveBoundaryMargin = 0.f; // keep fitted structures away from the protected alpha/lock band
    bool weakBoundaryDuplicates = true; // hole-mode painted outlines become sparse references; manual paths stay exact
};

struct MeshFeatureDetection {
    std::vector<MeshFeatureLine> lines;
    bool limited = false;
    int mergedNearbyLines = 0;
    std::string error;
};

// RGB edges, not semantic recognition. Alpha/support boundaries are excluded.
MeshFeatureDetection detectMeshTextureFeatures(const uint8_t* rgba, int width, int height,
    const MeshFeatureSettings& settings, const uint8_t* support = nullptr);

// Build a planar constraint graph: share crossings and split edges through
// existing points. Transactional; invalid/oversized input leaves points unchanged.
bool appendMeshFeatureConstraints(std::vector<float>& points,
    const std::vector<std::vector<float>>& paths,
    std::vector<std::pair<int, int>>& edges, std::string& error);

} // namespace l2d
