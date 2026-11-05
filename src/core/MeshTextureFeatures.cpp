#include "core/MeshTextureFeatures.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

namespace l2d {
namespace {
double lengthOf(const std::vector<float>& p) {
    double sum = 0;
    for (size_t i = 2; i + 1 < p.size(); i += 2)
        sum += std::hypot(p[i] - p[i - 2], p[i + 1] - p[i - 1]);
    return sum;
}
// Keep endpoints and bends. A long straight crease still has only two anchors.
template<class SegmentAllowed>
std::vector<float> samplePath(const std::vector<float>& p, float tolerance, const SegmentAllowed& allowed) {
    const int n = int(p.size() / 2);
    if (n < 2) return {};
    std::vector<uint8_t> keep(n);
    keep.front() = keep.back() = 1;
    std::vector<std::pair<int,int>> stack{{0, n - 1}};
    while (!stack.empty()) {
        auto [a,b] = stack.back(); stack.pop_back();
        const double dx = p[2*b] - p[2*a], dy = p[2*b+1] - p[2*a+1];
        const double d2 = dx*dx + dy*dy;
        double best = double(tolerance)*tolerance; int at = -1;
        for (int i = a + 1; i < b; ++i) {
            const double t = d2 > 0 ? std::clamp(((p[2*i]-p[2*a])*dx + (p[2*i+1]-p[2*a+1])*dy)/d2, 0., 1.) : 0.;
            const double distance = std::pow(p[2*i]-p[2*a]-t*dx, 2) + std::pow(p[2*i+1]-p[2*a+1]-t*dy, 2);
            if (distance > best) { best = distance; at = i; }
        }
        // Coarse simplification must not bridge a transparent notch. Recover
        // traced stations until each chord remains within the source texture.
        if(at<0&&b>a+1&&!allowed(p[2*a],p[2*a+1],p[2*b],p[2*b+1]))at=(a+b)/2;
        if (at >= 0) { keep[at] = 1; stack.emplace_back(a,at); stack.emplace_back(at,b); }
    }
    std::vector<float> result{p[0], p[1]};
    for (int i = 1; i < n; ++i) if (keep[i]) {
        result.push_back(p[2*i]); result.push_back(p[2*i+1]);
    }
    return result;
}

double pathDistance(const std::vector<float>& path,double x,double y) {
    double best=1e30;
    for(size_t k=2;k<path.size();k+=2){const double dx=path[k]-path[k-2],dy=path[k+1]-path[k-1],sq=dx*dx+dy*dy;
        const double t=sq?std::clamp(((x-path[k-2])*dx+(y-path[k-1])*dy)/sq,0.,1.):0;
        best=std::min(best,std::hypot(x-path[k-2]-t*dx,y-path[k-1]-t*dy));}
    return best;
}

// Round gentle turns with tangent quadratic arcs, then sample by arc length.
// This is a bounded approximation of CV evidence, not a semantic spline. Keep
// sharp corners/junction endpoints, reject alpha excursions and new crossings.
template<class SegmentAllowed>
std::vector<float> fitCurve(const std::vector<float>& original,double tolerance,double spacing,const SegmentAllowed& allowed) {
    if(original.size()<6||tolerance<=0)return original;
    auto path=original;
    const bool closed=std::hypot(path[0]-path[path.size()-2],path[1]-path.back())<1e-5;
    if(closed){ // put the seam on a straight span, not on an arbitrary corner
        path={(original[0]+original[2])*.5f,(original[1]+original[3])*.5f};
        path.insert(path.end(),original.begin()+2,original.end());path.push_back(path[0]);path.push_back(path[1]);
    }
    std::vector<float> dense{path[0],path[1]};std::vector<size_t> stops{0};bool changed=false;
    auto append=[&](double x,double y){if(std::hypot(x-dense[dense.size()-2],y-dense.back())>1e-5){dense.push_back(float(x));dense.push_back(float(y));}};
    for(size_t k=2;k+2<path.size();k+=2){
        const double ux=path[k]-path[k-2],uy=path[k+1]-path[k-1],vx=path[k+2]-path[k],vy=path[k+3]-path[k+1];
        const double a=std::hypot(ux,uy),b=std::hypot(vx,vy),dot=a*b>1e-9?(ux*vx+uy*vy)/(a*b):-1;
        if(dot<.4226182617||a<1e-4||b<1e-4){append(path[k],path[k+1]);stops.push_back(dense.size()/2-1);continue;}
        if(dot>.99999){append(path[k],path[k+1]);continue;}
        const double radius=std::min({a*.45,b*.45,2*tolerance/std::sqrt((1-dot)*.5)});
        const double ax=path[k]-ux/a*radius,ay=path[k+1]-uy/a*radius,bx=path[k]+vx/b*radius,by=path[k+1]+vy/b*radius;
        append(ax,ay);const int steps=std::clamp(int(std::ceil(radius*2)),4,64);
        for(int i=1;i<=steps;++i){const double t=double(i)/steps,s=1-t;append(s*s*ax+2*s*t*path[k]+t*t*bx,s*s*ay+2*s*t*path[k+1]+t*t*by);}
        changed=true;
    }
    if(!changed)return original;
    append(path[path.size()-2],path.back());stops.push_back(dense.size()/2-1);
    auto corridor=[&](const std::vector<float>& from,const std::vector<float>& to,bool checkAlpha){
        for(size_t k=2;k<from.size();k+=2){const double dx=from[k]-from[k-2],dy=from[k+1]-from[k-1];
            if(checkAlpha&&!allowed(from[k-2],from[k-1],from[k],from[k+1]))return false;
            const int n=std::max(1,int(std::ceil(std::hypot(dx,dy)*2)));
            for(int i=0;i<=n;++i)if(pathDistance(to,from[k-2]+dx*i/n,from[k-1]+dy*i/n)>tolerance+1e-4)return false;
        }return true;};
    // Refine only when needed to stay within the fitting corridor. Bound point
    // growth; a difficult fit falls back to the original sparse polyline.
    for(int attempt=0;attempt<3;++attempt){
        std::vector<float> sampled{dense[0],dense[1]};
        for(size_t part=1;part<stops.size();++part){const size_t start=stops[part-1],end=stops[part];
            if(start==end)continue;
            std::vector<double> arc{0};for(size_t i=start+1;i<=end;++i)arc.push_back(arc.back()+std::hypot(dense[2*i]-dense[2*i-2],dense[2*i+1]-dense[2*i-1]));
            const double length=arc.back();int count=std::max(1,int(std::ceil(length/std::max(1.,spacing/std::pow(1.5,attempt)))));
            const std::vector<float> chord{dense[2*start],dense[2*start+1],dense[2*end],dense[2*end+1]};
            bool straight=true;for(size_t i=start+1;i<end;++i)if(pathDistance(chord,dense[2*i],dense[2*i+1])>.05){straight=false;break;}
            if(straight)count=1;else count=std::max(2,count);
            for(int j=1;j<=count;++j){const double at=length*j/count;auto it=std::lower_bound(arc.begin()+1,arc.end(),at);const size_t index=std::min(size_t(it-arc.begin()),arc.size()-1);
                const double t=(at-arc[index-1])/std::max(1e-9,arc[index]-arc[index-1]);const size_t a=start+index-1,b=a+1;
                sampled.push_back(float(dense[2*a]+t*(dense[2*b]-dense[2*a])));sampled.push_back(float(dense[2*a+1]+t*(dense[2*b+1]-dense[2*a+1])));}
        }
        if(sampled.size()>std::min(size_t(512),original.size()*3))break;
        if(corridor(sampled,original,true)&&corridor(original,sampled,false))return sampled;
    }
    return original;
}

int crossingCount(const std::vector<float>& a,const std::vector<float>& b,bool same=false) {
    auto cross=[](double ax,double ay,double bx,double by,double x,double y){return (bx-ax)*(y-ay)-(by-ay)*(x-ax);};
    int count=0;
    for(size_t i=2;i<a.size();i+=2)for(size_t j=same?i+4:2;j<b.size();j+=2){
        const double p=cross(a[i-2],a[i-1],a[i],a[i+1],b[j-2],b[j-1]),q=cross(a[i-2],a[i-1],a[i],a[i+1],b[j],b[j+1]);
        const double r=cross(b[j-2],b[j-1],b[j],b[j+1],a[i-2],a[i-1]),s=cross(b[j-2],b[j-1],b[j],b[j+1],a[i],a[i+1]);
        if(p*q<-1e-8&&r*s<-1e-8)++count;
    }return count;
}

// A fit with the same intersection count can still squeeze a neighboring
// feature into an almost tangent sliver. Compare clearance at matching source
// locations, not just intersections. Existing shared junctions remain exact.
bool keepsClearance(const std::vector<float>& fit,const std::vector<float>& source,const std::vector<float>& neighbor) {
    for(size_t k=2;k<fit.size();k+=2){const double dx=fit[k]-fit[k-2],dy=fit[k+1]-fit[k-1];
        const int count=std::max(1,int(std::ceil(std::hypot(dx,dy))));
        for(int i=0;i<=count;++i){const double x=fit[k-2]+dx*i/count,y=fit[k-1]+dy*i/count;
            double best=1e30,ox=x,oy=y;
            for(size_t j=2;j<source.size();j+=2){const double sx=source[j]-source[j-2],sy=source[j+1]-source[j-1],sq=sx*sx+sy*sy;
                const double t=sq?std::clamp(((x-source[j-2])*sx+(y-source[j-1])*sy)/sq,0.,1.):0;
                const double px=source[j-2]+t*sx,py=source[j-1]+t*sy,d=std::hypot(px-x,py-y);
                if(d<best){best=d;ox=px;oy=py;}}
            const double before=pathDistance(neighbor,ox,oy),after=pathDistance(neighbor,x,y);
            if(after+.02<before*.8)return false;
        }
    }return true;
}

// Zhang-Suen thinning of the small centerline graph produced by paired edges.
// Endpoints and junctions survive; the iteration cap bounds interactive work.
void thinCenters(std::vector<uint8_t>& mask, int w, int h) {
    for (int iteration=0; iteration<64; ++iteration) {
        bool changed=false;
        for (int pass=0; pass<2; ++pass) {
            std::vector<int> erase;
            for (int y=1; y<h-1; ++y) for (int x=1; x<w-1; ++x) {
                const int k=y*w+x;
                if (!mask[k]) continue;
                const int p[8]={mask[k-w],mask[k-w+1],mask[k+1],mask[k+w+1],
                    mask[k+w],mask[k+w-1],mask[k-1],mask[k-w-1]};
                const int count=std::accumulate(p,p+8,0);
                int transitions=0;
                for(int i=0;i<8;++i) transitions+=!p[i]&&p[(i+1)%8];
                if(count<2||count>6||transitions!=1)continue;
                if(pass==0 ? (p[0]*p[2]*p[4]||p[2]*p[4]*p[6])
                           : (p[0]*p[2]*p[6]||p[0]*p[4]*p[6]))continue;
                erase.push_back(k);
            }
            for(int k:erase)mask[k]=0;
            changed|=!erase.empty();
        }
        if(!changed)break;
    }
}
}

MeshFeatureDetection detectMeshTextureFeatures(const uint8_t* rgba, int width, int height,
    const MeshFeatureSettings& s, const uint8_t* support) {
    MeshFeatureDetection out;
    if (!rgba || width < 7 || height < 7 || int64_t(width)*height > 100000000) {
        out.error = "纹理尺寸无效或过大"; return out;
    }
    // Bound interactive detection cost; output stays in source pixel coordinates.
    const double scale = std::min(1., 1024. / std::max(width,height));
    const int w = std::max(7,int(width*scale)), h = std::max(7,int(height*scale)), n = w*h;
    const float sx = float(width)/w, sy = float(height)/h;
    std::vector<std::array<float,3>> rgb(n), smooth(n);
    std::vector<uint8_t> valid(n), interior(n);
    for (int y=0; y<h; ++y) for (int x=0; x<w; ++x) {
        const size_t src = size_t(std::min(height-1,int((y+.5f)*sy)))*width + std::min(width-1,int((x+.5f)*sx));
        const int alpha = rgba[4*src+3];
        valid[y*w+x] = alpha > std::max(12,s.alphaThreshold) && (!support || support[src]);
        for (int c=0;c<3;++c) rgb[y*w+x][c] = s.premultiplied && alpha > 0
            ? std::min(255.f, rgba[4*src+c]*255.f/alpha) : float(rgba[4*src+c]);
    }
    for (int y=2;y<h-2;++y) for (int x=2;x<w-2;++x) {
        bool ok = true;
        for (int j=-2;j<=2 && ok;++j) for (int i=-2;i<=2;++i) if (!valid[(y+j)*w+x+i]) { ok=false; break; }
        interior[y*w+x] = ok;
    }
    const int kernel[3] = {1,2,1};
    for (int y=1;y<h-1;++y) for (int x=1;x<w-1;++x)
        for (int j=-1;j<=1;++j) for(int i=-1;i<=1;++i) for (int c=0;c<3;++c)
            smooth[y*w+x][c] += rgb[(y+j)*w+x+i][c]*kernel[j+1]*kernel[i+1]/16.f;
    std::vector<float> magnitude(n), gradientX(n), gradientY(n);
    std::vector<uint8_t> direction(n), weak(n), edge(n);
    for (int y=2;y<h-2;++y) for (int x=2;x<w-2;++x) {
        const int k=y*w+x; if(!interior[k]) continue;
        float best=0, gx=0, gy=0;
        for(int c=0;c<3;++c) {
            const float dx = (smooth[k-w+1][c]+2*smooth[k+1][c]+smooth[k+w+1][c]
                -smooth[k-w-1][c]-2*smooth[k-1][c]-smooth[k+w-1][c])/4.f;
            const float dy = (smooth[k+w-1][c]+2*smooth[k+w][c]+smooth[k+w+1][c]
                -smooth[k-w-1][c]-2*smooth[k-w][c]-smooth[k-w+1][c])/4.f;
            const float m=std::hypot(dx,dy);
            if(m>best) { best=m; gx=dx; gy=dy; }
        }
        magnitude[k]=best;
        gradientX[k]=gx;gradientY[k]=gy;
        double a=std::atan2(gy,gx)*180/3.141592653589793; if(a<0)a+=180;
        direction[k] = uint8_t(int((a+22.5)/45)%4);
    }
    const float high=float(std::clamp(s.contrast,1,255)), low=high*.45f;
    const int offset[4]={1,w+1,w,w-1};
    std::vector<int> queue;
    for(int k=2*w+2;k<n-2*w-2;++k) if(interior[k]) {
        const int d=offset[direction[k]];
        if(magnitude[k]>=low && magnitude[k]>=magnitude[k+d] && magnitude[k]>magnitude[k-d]) {
            weak[k]=1;
            if(magnitude[k]>=high) { edge[k]=1; queue.push_back(k); }
        }
    }
    const int dx[8]={1,1,0,-1,-1,-1,0,1}, dy[8]={0,1,1,1,0,-1,-1,-1};
    for(size_t q=0;q<queue.size();++q) for(int d=0;d<8;++d) {
        const int k=queue[q]+dy[d]*w+dx[d];
        if(weak[k]&&!edge[k]) { edge[k]=1;queue.push_back(k); }
    }
    // A fold is often a thin dark/light stripe: Canny sees its two sides.
    // Pair opposing gradients only when the colors outside both sides agree.
    // Broad color boundaries remain separate; alpha/support edges never enter.
    std::vector<uint8_t> centers(n),paired(n);
    std::vector<int> centerOf(n,-1);
    const int foldWidth=std::clamp(int(s.maxFoldWidth/std::max(sx,sy)),0,32);
    for(int k:queue) {
        if(magnitude[k]<=0)continue;
        const float nx=gradientX[k]/magnitude[k],ny=gradientY[k]/magnitude[k];
        bool found=false;
        for(int distance=2;distance<=foldWidth&&!found;++distance)for(int sign:{-1,1}) {
            const int x=int(std::lround(k%w+sign*distance*nx)),y=int(std::lround(k/w+sign*distance*ny));
            if(x<2||x>=w-2||y<2||y>=h-2)continue;
            const int v=y*w+x;
            if(!edge[v]||magnitude[v]<=0)continue;
            const float dot=(gradientX[k]*gradientX[v]+gradientY[k]*gradientY[v])/(magnitude[k]*magnitude[v]);
            if(dot>-.8f)continue;
            const int ax=std::clamp(int(std::lround(k%w-sign*2*nx)),0,w-1);
            const int ay=std::clamp(int(std::lround(k/w-sign*2*ny)),0,h-1);
            const int bx=std::clamp(int(std::lround(x+sign*2*nx)),0,w-1);
            const int by=std::clamp(int(std::lround(y+sign*2*ny)),0,h-1);
            if(!interior[ay*w+ax]||!interior[by*w+bx])continue;
            float difference=0;
            for(int c=0;c<3;++c)difference=std::max(difference,std::abs(rgb[ay*w+ax][c]-rgb[by*w+bx][c]));
            if(difference>std::max(12.f,high*.75f))continue;
            const int center=((k/w+y+1)/2)*w+(k%w+x+1)/2;
            centers[center]=1;paired[k]=paired[v]=1;
            centerOf[k]=center;centerOf[v]=center;found=true;break;
        }
    }
    // Join centers of adjacent samples before thinning; this avoids gaps on diagonals.
    for(int k:queue)if(centerOf[k]>=0)for(int d=0;d<8;++d) {
        const int v=k+dy[d]*w+dx[d];
        if(centerOf[v]<0)continue;
        const int a=centerOf[k],b=centerOf[v],dxC=b%w-a%w,dyC=b/w-a/w;
        const int steps=std::max(std::abs(dxC),std::abs(dyC));
        if(steps>3)continue;
        for(int j=0;j<=steps;++j) {
            const float t=steps?float(j)/steps:0;
            centers[int(std::lround(a/w+t*dyC))*w+int(std::lround(a%w+t*dxC))]=1;
        }
    }
    thinCenters(centers,w,h);
    for(int k:queue)if(paired[k])edge[k]=0;
    auto traceMask=[&](const std::vector<uint8_t>& mask,bool centerline) {
    auto neighbors=[&](int k) {
        std::vector<int> r;
        for(int d=0;d<8;++d) {
            const int v=k+dy[d]*w+dx[d];
            // A diagonal next to an axial connection makes spurious junctions.
            if(!mask[v] || (dx[d]&&dy[d]&&(mask[k+dx[d]]||mask[k+dy[d]*w])))continue;
            r.push_back(v);
        }
        return r;
    };
    std::set<std::pair<int,int>> visited;
    auto trace=[&](int start,int next) {
        std::vector<float> p{(start%w+.5f)*sx,(start/w+.5f)*sy};
        int previous=start, current=next;
        for(int guard=0;guard<n;++guard) {
            if(!visited.emplace(std::min(previous,current),std::max(previous,current)).second)break;
            p.push_back((current%w+.5f)*sx);p.push_back((current/w+.5f)*sy);
            const auto ns=neighbors(current);
            if(current==start || ns.size()!=2)break;
            const int v=ns[0]==previous?ns[1]:ns[0];previous=current;current=v;
        }
        const float length=float(lengthOf(p));
        if(length>=std::max(2,s.minLength)) {
            auto onTexture=[&](float ax,float ay,float bx,float by) {
                const int probes=std::max(1,int(std::ceil(std::hypot((bx-ax)/sx,(by-ay)/sy)*2)));
                for(int j=0;j<=probes;++j) {
                    const int x=int((ax+(bx-ax)*j/probes)/sx),y=int((ay+(by-ay)*j/probes)/sy);
                    if(x<0||x>=w||y<0||y>=h||!valid[y*w+x])return false;
                }
                return true;
            };
            auto anchors=samplePath(p,float(std::max(1,s.cornerTolerance)),onTexture);
            if(anchors.size()>=4)out.lines.push_back({std::move(anchors),true,length,centerline,true});
        }
    };
    std::vector<int> pixels;
    for(int y=2;y<h-2;++y)for(int x=2;x<w-2;++x)if(mask[y*w+x])pixels.push_back(y*w+x);
    for(int k:pixels) { const auto ns=neighbors(k); if(ns.size()!=2)for(int v:ns)trace(k,v); }
    for(int k:pixels)for(int v:neighbors(k)) if(!visited.count({std::min(k,v),std::max(k,v)}))trace(k,v);
    };
    traceMask(centers,true);traceMask(edge,false);
    // Heuristic bright-island detection, not semantic material recognition.
    // A bright background is not a highlight; only relatively small, bright,
    // low-chroma components qualify. Preserve dark creases and explicit paths.
    if(s.weakHighlights){
        std::vector<float> luminance;int opaque=0;
        for(int k=0;k<n;++k)if(valid[k]){luminance.push_back(.2126f*rgb[k][0]+.7152f*rgb[k][1]+.0722f*rgb[k][2]);++opaque;}
        std::sort(luminance.begin(),luminance.end());
        const float threshold=luminance.empty()?255.f:std::max(190.f,luminance[luminance.size()/2]+20.f);
        std::vector<int> labels(n,-1),sizes;
        auto bright=[&](int k){const auto& c=rgb[k];return valid[k]&&.2126f*c[0]+.7152f*c[1]+.0722f*c[2]>=threshold&&*std::max_element(c.begin(),c.end())-*std::min_element(c.begin(),c.end())<=60;};
        for(int seed=0;seed<n;++seed)if(labels[seed]<0&&bright(seed)){
            const int label=int(sizes.size());std::vector<int> work{seed};labels[seed]=label;
            for(size_t j=0;j<work.size();++j){const int k=work[j];for(auto [xx,yy]:std::initializer_list<std::pair<int,int>>{{k%w-1,k/w},{k%w+1,k/w},{k%w,k/w-1},{k%w,k/w+1}})
                if(xx>=0&&xx<w&&yy>=0&&yy<h){const int v=yy*w+xx;if(labels[v]<0&&bright(v)){labels[v]=label;work.push_back(v);}}}
            sizes.push_back(int(work.size()));
        }
        for(auto& line:out.lines){int probes=0,hits=0;
            for(size_t i=2;i<line.points.size();i+=2){const double ax=line.points[i-2]/sx,ay=line.points[i-1]/sy,bx=line.points[i]/sx,by=line.points[i+1]/sy;
                const int steps=std::max(1,int(std::ceil(std::hypot(bx-ax,by-ay)/3)));
                for(int j=0;j<steps;++j){const int x=int(ax+(bx-ax)*j/steps),y=int(ay+(by-ay)*j/steps);bool hit=false;
                    for(int yy=std::max(0,y-3);yy<=std::min(h-1,y+3);++yy)for(int xx=std::max(0,x-3);xx<=std::min(w-1,x+3);++xx){const int label=labels[yy*w+xx];if(label>=0&&sizes[label]>=std::max(24,opaque/1000)&&sizes[label]<opaque*.35)hit=true;}
                    ++probes;hits+=hit;}}
            if(probes>=3&&hits>=probes*.7){line.weakReference=true;line.centerline=false;}
        }
    }
    std::stable_sort(out.lines.begin(),out.lines.end(),[](const auto& a,const auto& b){return a.length>b.length;});
    // Fit a nearby group to its longest coherent path. Bound every member
    // against the ORIGINAL representative, so a chain of close lines cannot
    // swallow an arbitrarily wide color region. No alpha/support gaps or
    // crossing branches may be bridged. Explicit user paths never enter here.
    if(s.nearbyMergeDistance>0) {
        const double radius=std::clamp(s.nearbyMergeDistance,1,64);
        struct Probe {double x=0,y=0,tx=0,ty=0;};
        auto samples=[&](const std::vector<float>& path) {
            std::vector<Probe> result;
            for(size_t k=2;k<path.size();k+=2) {
                const double dx=path[k]-path[k-2],dy=path[k+1]-path[k-1],length=std::hypot(dx,dy);
                if(length<1e-6)continue;
                const int count=std::clamp(int(std::ceil(length/std::max(1.,radius*.5))),1,256);
                for(int j=0;j<count;++j)result.push_back({path[k-2]+dx*j/count,path[k-1]+dy*j/count,dx/length,dy/length});
            }
            if(!result.empty()){auto end=result.back();end.x=path[path.size()-2];end.y=path.back();result.push_back(end);}
            return result;
        };
        struct Nearest {double distance=1e30,x=0,y=0,tx=0,ty=0,arc=0;};
        auto nearest=[](const std::vector<float>& path,double x,double y) {
            Nearest best;double arc=0;
            for(size_t k=2;k<path.size();k+=2) {
                const double dx=path[k]-path[k-2],dy=path[k+1]-path[k-1],length=std::hypot(dx,dy);
                if(length<1e-6)continue;
                const double t=std::clamp(((x-path[k-2])*dx+(y-path[k-1])*dy)/(length*length),0.,1.);
                const double px=path[k-2]+t*dx,py=path[k-1]+t*dy,d=std::hypot(x-px,y-py);
                if(d<best.distance)best={d,px,py,dx/length,dy/length,arc+t*length};
                arc+=length;
            }
            return best;
        };
        auto allowed=[&](double ax,double ay,double bx,double by) {
            const int count=std::max(1,int(std::ceil(std::hypot((bx-ax)/sx,(by-ay)/sy)*2)));
            for(int k=0;k<=count;++k){const int x=int((ax+(bx-ax)*k/count)/sx),y=int((ay+(by-ay)*k/count)/sy);
                if(x<0||x>=w||y<0||y>=h||!valid[y*w+x])return false;}
            return true;
        };
        const size_t budget=std::min(size_t(256),out.lines.size());
        std::vector<uint8_t> removed(out.lines.size());
        for(size_t i=0;i<budget;++i)if(!removed[i]) {
            const auto base=out.lines[i].points;const auto stations=samples(base);
            if(stations.size()<2||stations.size()>2048)continue;
            std::vector<std::vector<float>> members;
            for(size_t j=i+1;j<budget&&members.size()<16;++j)if(!removed[j]) {
                const auto& candidate=out.lines[j];const auto tests=samples(candidate.points);
                if(candidate.weakReference!=out.lines[i].weakReference)continue;
                if(tests.size()<2||tests.size()>2048)continue;
                bool fits=true;int side=0,orientation=0;double first=0,last=0;
                for(size_t k=0;k<tests.size();++k) {
                    const auto& q=tests[k];const auto near=nearest(base,q.x,q.y);
                    const double dot=q.tx*near.tx+q.ty*near.ty,cross=near.tx*(q.y-near.y)-near.ty*(q.x-near.x);
                    const int sign=dot>=0?1:-1,currentSide=cross>.35?1:cross<-.35?-1:0;
                    if(near.distance>radius||std::abs(dot)<.94||!allowed(q.x,q.y,near.x,near.y)||
                        (orientation&&orientation!=sign)||(side&&currentSide&&side!=currentSide)||
                        (k&&sign*(near.arc-last)<-.5)){fits=false;break;}
                    if(!k)first=near.arc;last=near.arc;orientation=sign;if(currentSide)side=currentSide;
                }
                if(!fits||std::abs(last-first)<candidate.length*.75)continue;
                // Do not erase an explicit intersection represented in separate
                // graph paths merely because their tangents are almost parallel.
                if(side==0&&candidate.length<out.lines[i].length*.5)continue;
                members.push_back(candidate.points);removed[j]=1;++out.mergedNearbyLines;
            }
            if(members.empty())continue;
            std::vector<float> fitted;
            for(const auto& q:stations) {
                double x=q.x,y=q.y;int count=1;
                for(const auto& member:members){const auto near=nearest(member,q.x,q.y);
                    if(near.distance<=radius&&std::abs(q.tx*near.tx+q.ty*near.ty)>.94&&allowed(q.x,q.y,near.x,near.y)){x+=near.x;y+=near.y;++count;}}
                fitted.push_back(float(x/count));fitted.push_back(float(y/count));
            }
            auto fit=samplePath(fitted,float(std::max(1,s.cornerTolerance)),allowed);
            bool safe=fit.size()>=4;
            for(size_t k=2;safe&&k<fit.size();k+=2)safe=allowed(fit[k-2],fit[k-1],fit[k],fit[k+1]);
            // Conservative fallback retains the original representative when
            // averaging bends would cut a hole; group stays within its corridor.
            out.lines[i].points=safe?std::move(fit):base;
            out.lines[i].length=float(lengthOf(out.lines[i].points));
        }
        std::vector<MeshFeatureLine> kept;for(size_t i=0;i<out.lines.size();++i)if(!removed[i])kept.push_back(std::move(out.lines[i]));
        out.lines=std::move(kept);
        std::stable_sort(out.lines.begin(),out.lines.end(),[](const auto& a,const auto& b){return a.length>b.length;});
    }
    for(auto& line:out.lines)if(line.weakReference){
        auto allowed=[&](float ax,float ay,float bx,float by){const int steps=std::max(1,int(std::ceil(std::hypot((bx-ax)/sx,(by-ay)/sy))));
            for(int i=0;i<=steps;++i){const int x=int((ax+(bx-ax)*i/steps)/sx),y=int((ay+(by-ay)*i/steps)/sy);if(x<0||x>=w||y<0||y>=h||!valid[y*w+x])return false;}return true;};
        line.points=samplePath(line.points,float(std::max(4,s.cornerTolerance*3)),allowed);
        line.length=float(lengthOf(line.points));
    }
    size_t points=0, count=0;
    for(const auto& line:out.lines) {
        if(count>=64 || points+line.points.size()/2>2048)break;
        points+=line.points.size()/2;++count;
    }
    out.limited=count<out.lines.size();out.lines.resize(count);
    std::vector<int> clearance;
    if(s.curveSmoothing>0&&s.curveBoundaryMargin>0){
        // Conservative chessboard distance to excluded alpha/support pixels.
        clearance.resize(n);
        for(int y=0;y<h;++y)for(int x=0;x<w;++x)clearance[y*w+x]=valid[y*w+x]?std::min({x+1,y+1,w-x,h-y}):0;
        for(int y=1;y<h;++y)for(int x=1;x<w-1;++x){const int k=y*w+x;clearance[k]=std::min(clearance[k],1+std::min({clearance[k-1],clearance[k-w-1],clearance[k-w],clearance[k-w+1]}));}
        for(int y=h-2;y>=0;--y)for(int x=w-2;x>0;--x){const int k=y*w+x;clearance[k]=std::min(clearance[k],1+std::min({clearance[k+1],clearance[k+w-1],clearance[k+w],clearance[k+w+1]}));}
    }
    if(s.curveSmoothing>0)for(size_t i=0;i<out.lines.size();++i){auto& line=out.lines[i];
        // Weak highlights remain sparse suggestions, not another dense row.
        if(line.weakReference)continue;
        auto allowed=[&](double ax,double ay,double bx,double by){const int steps=std::max(1,int(std::ceil(std::hypot((bx-ax)/sx,(by-ay)/sy)*2)));
            for(int k=0;k<=steps;++k){const int x=int((ax+(bx-ax)*k/steps)/sx),y=int((ay+(by-ay)*k/steps)/sy);
                if(x<0||x>=w||y<0||y>=h||!valid[y*w+x]||(!clearance.empty()&&clearance[y*w+x]<=s.curveBoundaryMargin/std::min(sx,sy)))return false;}return true;};
        auto fit=fitCurve(line.points,std::clamp(s.curveSmoothing,0.f,12.f),std::max(1.f,s.curveSpacing),allowed);
        if(fit==line.points||points-line.points.size()/2+fit.size()/2>2048)continue;
        bool safe=crossingCount(fit,fit,true)==crossingCount(line.points,line.points,true);
        for(size_t j=0;safe&&j<out.lines.size();++j)if(j!=i&&!out.lines[j].weakReference)
            safe=crossingCount(fit,out.lines[j].points)==crossingCount(line.points,out.lines[j].points)&&keepsClearance(fit,line.points,out.lines[j].points);
        if(safe){points=points-line.points.size()/2+fit.size()/2;line.points=std::move(fit);line.curved=true;line.length=float(lengthOf(line.points));}
    }
    return out;
}

bool appendMeshFeatureConstraints(std::vector<float>& points,
    const std::vector<std::vector<float>>& paths,
    std::vector<std::pair<int,int>>& edges, std::string& error) {
    error.clear();
    std::vector<float> p=points;
    if(p.size()%2 || p.empty()) {error="结构线输入网格无效";return false;}
    double range=1;
    for(float v:p) {if(!std::isfinite(v)){error="结构线坐标不是有限数";return false;}range=std::max(range,std::abs(double(v)));}
    const double eps=std::max(1e-6,range*1e-7);
    auto vertex=[&](double x,double y) {
        for(size_t i=0;i<p.size()/2;++i)if(std::hypot(p[2*i]-x,p[2*i+1]-y)<=eps)return int(i);
        const int id=int(p.size()/2);p.push_back(float(x));p.push_back(float(y));return id;
    };
    std::vector<std::pair<int,int>> segments;
    for(const auto& path:paths) {
        if(path.size()%2) {error="结构线坐标不完整";return false;}
        for(float v:path)if(!std::isfinite(v)){error="结构线坐标不是有限数";return false;}
        for(size_t i=2;i+1<path.size();i+=2) {
            int a=vertex(path[i-2],path[i-1]),b=vertex(path[i],path[i+1]);
            if(a!=b)segments.emplace_back(a,b);
            if(segments.size()>2048 || p.size()/2>16000){error="结构线过多，请提高点距或忽略部分线";return false;}
        }
    }
    // Split every proper crossing before recovering any constraint.
    for(size_t i=0;i<segments.size();++i)for(size_t j=0;j<i;++j) {
        const auto [a,b]=segments[i]; const auto [c,d]=segments[j];
        const double ax=p[2*a],ay=p[2*a+1],bx=p[2*b],by=p[2*b+1];
        const double cx=p[2*c],cy=p[2*c+1],rx=bx-ax,ry=by-ay,sx=p[2*d]-cx,sy=p[2*d+1]-cy;
        const double denom=rx*sy-ry*sx;
        if(std::abs(denom)<=eps*std::max(std::hypot(rx,ry),std::hypot(sx,sy)))continue;
        const double t=((cx-ax)*sy-(cy-ay)*sx)/denom,u=((cx-ax)*ry-(cy-ay)*rx)/denom;
        if(t>0 && t<1 && u>0 && u<1)vertex(ax+t*rx,ay+t*ry);
    }
    if(p.size()/2>16000){error="结构线交点过多";return false;}
    std::set<std::pair<int,int>> result;
    for(auto [a,b]:segments) {
        const double ax=p[2*a],ay=p[2*a+1],dx=p[2*b]-ax,dy=p[2*b+1]-ay,len=std::hypot(dx,dy);
        std::vector<std::pair<double,int>> cuts;
        for(int k=0;k<int(p.size()/2);++k) {
            const double x=p[2*k]-ax,y=p[2*k+1]-ay,t=(x*dx+y*dy)/(len*len);
            if(t>=-eps/len && t<=1+eps/len && std::abs(x*dy-y*dx)<=eps*len)cuts.emplace_back(t,k);
        }
        std::sort(cuts.begin(),cuts.end());
        for(size_t i=1;i<cuts.size();++i)if(cuts[i-1].second!=cuts[i].second)
            result.emplace(std::min(cuts[i-1].second,cuts[i].second),std::max(cuts[i-1].second,cuts[i].second));
    }
    points=std::move(p);edges.assign(result.begin(),result.end());return true;
}
} // namespace l2d
