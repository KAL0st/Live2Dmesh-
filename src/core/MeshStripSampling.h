#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace l2d::meshSampling {

// Cubism 5.3 prepares a half-size working image, reducing its longest
// content side by 0.8 until <=256 (minimum 20). Point spacing is scaled by
// 0.5 in that image. Keep these units separate from a vertex-count budget.
inline double workScale(int width,int height) {
    const int longest=std::max(width,height);
    int side=longest/2;
    while(side>256)side=int(side*.8);
    return double(std::max(20,side))/std::max(1,longest);
}

inline double pointInterval(int width,int height,int setting) {
    return std::clamp(setting,1,2000)*.5/workScale(width,height);
}

// Candidate locations in an inset alpha region. Unlike a filled lattice,
// these are only searched by the sparse inner-ring sampler. A work-cell OR
// and four-neighbor distance keep faint pixels represented in this guide;
// MeshAutoGen still validates the final triangles against the native image.
inline std::vector<float> insetRegion(const std::vector<uint8_t>& rgba,int width,int height,double inset) {
    const double scale=workScale(width,height);
    const int w=int(std::ceil(width*scale))+2,h=int(std::ceil(height*scale))+2;
    std::vector<int> distance(size_t(w)*h);
    for(int y=0;y<height;++y)for(int x=0;x<width;++x)
        if(rgba[(size_t(y)*width+x)*4+3])distance[size_t(1+int(y*scale))*w+1+int(x*scale)]=w+h;
    for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x) {
        const int p=y*w+x;distance[p]=std::min({distance[p],distance[p-1]+1,distance[p-w]+1});
    }
    for(int y=h-2;y>0;--y)for(int x=w-2;x>0;--x) {
        const int p=y*w+x;distance[p]=std::min({distance[p],distance[p+1]+1,distance[p+w]+1});
    }
    std::vector<float> result;
    const int minimum=std::max(1,int(std::ceil(inset*scale)));
    for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x)if(distance[y*w+x]>=minimum) {
        result.push_back(float((x-.5)/scale));result.push_back(float((y-.5)/scale));
    }
    return result;
}

// Local implementation of the single-stroke thin-region test and medial
// path sampling. It does not duplicate PolygonMaker2's rail construction.
// Work pixels conservatively OR source alpha; final coverage is checked at
// source resolution by MeshAutoGen, independently of this guide path.
struct StrokeGuideInfo {int pixels=0,ends=0,junctions=0,radius=0;double length=0,widthRatio=0;};
inline std::vector<float> strokeCenterline(const std::vector<uint8_t>& rgba,int width,int height,bool narrowOnly=false,StrokeGuideInfo* info=nullptr) {
    if(info)*info={};
    if(width<1||height<1||rgba.size()!=size_t(width)*height*4)return {};
    const double scale=workScale(width,height);
    const int w=int(std::ceil(width*scale))+2,h=int(std::ceil(height*scale))+2;
    std::vector<uint8_t> mask(size_t(w)*h);
    for(int y=0;y<height;++y)for(int x=0;x<width;++x)
        if(rgba[(size_t(y)*width+x)*4+3])mask[size_t(1+int(y*scale))*w+1+int(x*scale)]=1;
    const int dx[8]={0,1,1,1,0,-1,-1,-1},dy[8]={-1,-1,0,1,1,1,0,-1};
    // Reject separate components, matching the official narrow-strip gate.
    std::vector<uint8_t> seen(mask.size());std::vector<int> queue;
    int components=0;
    for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x) {
        const int p=y*w+x;if(!mask[p]||seen[p])continue;
        if(++components>1)return {};
        queue={p};seen[p]=1;
        for(size_t i=0;i<queue.size();++i)for(int k=0;k<8;++k) {
            const int q=queue[i]+dy[k]*w+dx[k];
            if(mask[q]&&!seen[q]){seen[q]=1;queue.push_back(q);}
        }
    }
    if(components!=1)return {};
    // Four-neighbor distance propagation, rather than a Euclidean transform.
    std::vector<int> distance(mask.size(),w+h);
    for(size_t p=0;p<mask.size();++p)if(!mask[p])distance[p]=0;
    for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x) {
        const int p=y*w+x;distance[p]=std::min({distance[p],distance[p-1]+1,distance[p-w]+1});
    }
    int radius=0;
    for(int y=h-2;y>0;--y)for(int x=w-2;x>0;--x) {
        const int p=y*w+x;distance[p]=std::min({distance[p],distance[p+1]+1,distance[p+w]+1});radius=std::max(radius,distance[p]);
    }
    if(narrowOnly&&radius>10)return {};
    if(info)info->radius=radius;
    // Two-phase thinning preserves connected paths and their endpoints.
    auto skeleton=mask;
    for(int pass=0;pass<std::max(w,h);++pass) {
        bool changed=false;
        for(int phase=0;phase<2;++phase) {
            std::vector<int> erase;
            for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x) {
                const int p=y*w+x;if(!skeleton[p])continue;
                int n[8],count=0,turns=0;
                for(int k=0;k<8;++k){n[k]=skeleton[p+dy[k]*w+dx[k]];count+=n[k];}
                for(int k=0;k<8;++k)turns+=!n[k]&&n[(k+1)%8];
                if(count<2||count>6||turns!=1)continue;
                const bool preserve=phase==0?(n[0]*n[2]*n[4]||n[2]*n[4]*n[6]):(n[0]*n[2]*n[6]||n[0]*n[4]*n[6]);
                if(!preserve)erase.push_back(p);
            }
            for(int p:erase)skeleton[p]=0;
            changed|=!erase.empty();
        }
        if(!changed)break;
    }
    auto neighbors=[&](int p) {
        std::vector<int> result;
        for(int k=0;k<8;++k) {
            const int q=p+dy[k]*w+dx[k];if(!skeleton[q])continue;
            // A diagonal with an orthogonal alternative is redundant, and
            // would otherwise turn an ordinary raster bend into a junction.
            if(dx[k]&&dy[k]&&(skeleton[p+dx[k]]||skeleton[p+dy[k]*w]))continue;
            result.push_back(q);
        }
        return result;
    };
    // A wide rectangle has short diagonal skeleton spurs at its end caps.
    // Remove only leaf-to-junction branches shorter than that junction's
    // local width. Long branches (e.g. a forked lock) reject the single-stroke
    // model below instead of being replaced with an arbitrary principal axis.
    if(!narrowOnly)for(int pass=0;pass<8;++pass) {
        std::vector<int> erase;
        for(int y=1;y<h-1;++y)for(int x=1;x<w-1;++x) {
            const int p=y*w+x;if(!skeleton[p]||neighbors(p).size()!=1)continue;
            std::vector<int> branch{p};int previous=-1,at=p;double length=0;
            for(;;) {
                auto adjacent=neighbors(at);if(adjacent.size()!=2&&at!=p)break;
                int next=-1;for(int q:adjacent)if(q!=previous){next=q;break;}
                if(next<0)break;
                length+=std::hypot(double(next%w-at%w),double(next/w-at/w));
                previous=at;at=next;branch.push_back(at);
            }
            if(neighbors(at).size()>2&&length<distance[at]*2.25)
                erase.insert(erase.end(),branch.begin(),branch.end()-1);
        }
        if(erase.empty())break;
        for(int p:erase)skeleton[p]=0;
    }
    std::vector<int> ends;int count=0,junctions=0;
    for(size_t p=0;p<skeleton.size();++p)if(skeleton[p]) {
        ++count;const auto adjacent=neighbors(int(p));
        if(adjacent.size()==1)ends.push_back(int(p));
        else if(adjacent.size()!=2)++junctions;
    }
    if(info){info->pixels=count;info->ends=int(ends.size());info->junctions=junctions;}
    if(ends.size()!=2||junctions||count<6||(narrowOnly&&count<radius*3))return {};
    std::vector<int> path;int previous=-1,at=ends[0];
    while(path.size()<=size_t(count)) {
        path.push_back(at);if(at==ends[1])break;
        const auto adjacent=neighbors(at);int next=-1;
        for(int q:adjacent)if(q!=previous){next=q;break;}
        if(next<0)return {};previous=at;at=next;
    }
    if(path.size()!=size_t(count))return {};
    if(!narrowOnly) {
        double length=0;std::vector<int> radii;
        for(size_t i=0;i<path.size();++i) {
            radii.push_back(distance[path[i]]);
            if(i)length+=std::hypot(double(path[i]%w-path[i-1]%w),double(path[i]/w-path[i-1]/w));
        }
        std::sort(radii.begin(),radii.end());
        // Width-relative elongation also recognizes bent strips whose image
        // bounding box is nearly square; disks and compact panels stay 2-D.
        const double typicalWidth=2*radii[radii.size()/2];
        if(info){info->length=length/scale;info->widthRatio=(length+typicalWidth)/typicalWidth;}
        if(length+typicalWidth<2.2*typicalWidth)return {};
    }
    std::vector<float> result;result.reserve(path.size()*2);
    for(size_t i=0;i<path.size();++i) {
        double x=0,y=0;int n=0;
        const int window=i==0||i+1==path.size()?0:2;
        for(int j=std::max(0,int(i)-window);j<=std::min(int(path.size())-1,int(i)+window);++j){x+=path[j]%w-0.5;y+=path[j]/w-0.5;++n;}
        result.push_back(float(x/n/scale));result.push_back(float(y/n/scale));
    }
    if(!narrowOnly&&result.size()>=8) {
        // Thinning/pruning stops at the centers of rounded end caps. Extend
        // those tangents through the remaining alpha so cross sections also
        // support the caps instead of leaving a large triangle fan at each end.
        auto cap=[&](bool front) {
            std::vector<float> extra;const size_t count=result.size()/2;
            const size_t end=front?0:count-1,inside=front?std::min(size_t(8),count-1):count-1-std::min(size_t(8),count-1);
            const double x=result[2*end],y=result[2*end+1],dx=x-result[2*inside],dy=y-result[2*inside+1],length=std::hypot(dx,dy);
            if(length<1e-5)return extra;
            const double step=std::max(1.,1/scale);
            for(double at=step;at<=2*radius/scale;at+=step) {
                const double xx=x+dx/length*at,yy=y+dy/length*at;const int ix=int(std::floor(xx)),iy=int(std::floor(yy));
                if(ix<0||ix>=width||iy<0||iy>=height||!rgba[(size_t(iy)*width+ix)*4+3])break;
                extra.push_back(float(xx));extra.push_back(float(yy));
            }
            return extra;
        };
        auto front=cap(true),back=cap(false);std::vector<float> extended;
        for(size_t i=front.size()/2;i>0;--i){extended.push_back(front[2*i-2]);extended.push_back(front[2*i-1]);}
        extended.insert(extended.end(),result.begin(),result.end());extended.insert(extended.end(),back.begin(),back.end());result=std::move(extended);
    }
    return result;
}

inline std::vector<float> thinCenterline(const std::vector<uint8_t>& rgba,int width,int height) {
    return strokeCenterline(rgba,width,height,true);
}

} // namespace l2d::meshSampling
