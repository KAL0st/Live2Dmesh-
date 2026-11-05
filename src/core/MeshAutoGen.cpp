#include "core/MeshAutoGen.h"

#include "edit/MeshTopology.h"
#include "core/XsPrologue.h"
#include "core/TextureUvMapping.h"
#include "core/MeshStripSampling.h"
#include "core/MeshSymmetry.h"
#include "io/PngCodec.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace l2d {
namespace {

// This local generator uses alpha silhouettes and explicit sampling controls.
// Cubism's PolygonMaker2 density/offset algorithm has not been reconstructed;
// these controls must not be described as numerically equivalent to Cubism.
struct TextureDomain {
    std::vector<uint8_t> rgba;
    int width=0,height=0;
    float logicalW=0,logicalH=0;
    bool premultiplied=false,conservativeFallback=false;
    double toUv[6]={1,0,0,0,1,0},fromUv[6]={1,0,0,0,1,0};
    std::pair<float,float> uv(float x,float y) const {
        return {float(toUv[0]*x+toUv[1]*y+toUv[2]),float(toUv[3]*x+toUv[4]*y+toUv[5])};
    }
    std::pair<float,float> pixel(float u,float v) const {
        return {float(fromUv[0]*u+fromUv[1]*v+fromUv[2]),float(fromUv[3]*u+fromUv[4]*v+fromUv[5])};
    }
    std::vector<std::vector<float>> contours; // source-image pixel corners
    std::vector<std::vector<float>> holes, holeSupports;
    int mergedHoleVertices=0;
    std::vector<float> stripCenterline;
    std::vector<float> strokeCenterline;
    int sharpTips=0;
};

bool simpleOutline(const std::vector<float>& points,const VertexLoop& loop,double tolerance);
bool nearBoundary(const std::vector<float>& pts,const VertexLoop& loop,float x,float y,float margin);
bool buildTextureContours(TextureDomain& domain,const MeshAutoGenParams& params,
                          std::string& error,double inset=0);
bool prepareTextureHoles(TextureDomain& domain,const MeshAutoGenParams& params,std::string& error);

struct CurvePoint {double x=0,y=0;};
// Material changes deformation support, not the painted alpha silhouette.
// Cloth keeps intermediate stations even on straight spans. Rigid outlines
// may remove redundant stations, but only if the original polygon is enclosed.
void materialOutline(std::vector<float>& path,MeshAutoGenStyle style,double spacing,
                     const meshSymmetry::Axis* symmetry=nullptr) {
    if(style==MeshAutoGenStyle::Soft){
        std::vector<float> sampled;
        for(size_t i=0;i<path.size();i+=2){const size_t j=(i+2)%path.size();
            const double length=std::hypot(path[j]-path[i],path[j+1]-path[i+1]);
            const int pieces=std::max(1,int(std::ceil(length/std::max(1.,spacing*1.25)-1e-5)));
            for(int k=0;k<pieces;++k){const float f=float(k)/pieces;sampled.push_back(path[i]+f*(path[j]-path[i]));sampled.push_back(path[i+1]+f*(path[j+1]-path[i+1]));}}
        path=std::move(sampled);return;
    }
    if(style!=MeshAutoGenStyle::Rigid)return;
    const auto original=path;const double tolerance=std::clamp(spacing*.02,.05,.6);
    for(int pass=0;pass<1000&&path.size()>6;++pass){bool changed=false;
        auto removable=[&](size_t i){const size_t n=path.size()/2,a=(i+n-1)%n,b=(i+1)%n;
            const double ux=path[2*i]-path[2*a],uy=path[2*i+1]-path[2*a+1],vx=path[2*b]-path[2*i],vy=path[2*b+1]-path[2*i+1];
            return (ux*vx+uy*vy)>std::hypot(ux,uy)*std::hypot(vx,vy)*.984807753&&std::hypot(ux+vx,uy+vy)<=spacing*2.5;};
        for(size_t i=0;i<path.size()/2;++i){if(!removable(i))continue;size_t mate=i;
            if(symmetry&&symmetry->detected){const auto q=symmetry->reflect(path[2*i],path[2*i+1]);
                for(size_t j=0;j<path.size()/2;++j)if(std::hypot(path[2*j]-q.x,path[2*j+1]-q.y)<.005){mate=j;break;}
                if(mate==i||!removable(mate))continue;}
            auto proposal=path;for(size_t k:std::set<size_t,std::greater<size_t>>{i,mate})proposal.erase(proposal.begin()+2*k,proposal.begin()+2*k+2);
            if(proposal.size()<6)continue;
            VertexLoop loop(proposal.size()/2);std::iota(loop.begin(),loop.end(),0);
            if(!simpleOutline(proposal,loop,1e-6))continue;
            bool valid=true;
            for(size_t k=0;k<original.size()&&valid;k+=2){double nearest=1e30;
                for(size_t a=0;a<proposal.size();a+=2){const size_t b=(a+2)%proposal.size();const double dx=proposal[b]-proposal[a],dy=proposal[b+1]-proposal[a+1],sq=dx*dx+dy*dy;
                    const double t=sq?std::clamp(((original[k]-proposal[a])*dx+(original[k+1]-proposal[a+1])*dy)/sq,0.,1.):0.;
                    nearest=std::min(nearest,std::hypot(original[k]-proposal[a]-t*dx,original[k+1]-proposal[a+1]-t*dy));}
                valid=nearest<=tolerance&&(nearest<1e-4||pointInLoop(proposal,loop,original[k],original[k+1]));}
            if(valid){path=std::move(proposal);changed=true;break;}
        }
        if(!changed)break;
    }
}
std::vector<float> pairedPolygon(const std::vector<float>& points,const meshSymmetry::Axis& axis,bool outside,double spacing) {
    if(points.size()<6)return {};
    const double pi=3.14159265358979323846;
    std::vector<meshSymmetry::Point> local;std::vector<double> angles{0,pi};
    for(size_t i=0;i<points.size();i+=2){const auto q=axis.local(points[i],points[i+1]);local.push_back(q);angles.push_back(std::atan2(std::abs(q.x),q.y));}
    std::sort(angles.begin(),angles.end());angles.erase(std::unique(angles.begin(),angles.end(),[](double a,double b){return std::abs(a-b)<1e-8;}),angles.end());
    auto radius=[&](double x,double y) {
        std::vector<double> hits;
        for(size_t i=0;i<local.size();++i) {
            const auto a=local[i],b=local[(i+1)%local.size()];const double dx=b.x-a.x,dy=b.y-a.y,den=x*dy-y*dx;
            if(std::abs(den)<1e-10)continue;
            const double t=(a.x*dy-a.y*dx)/den,u=(a.x*y-a.y*x)/den;
            if(t>0&&u>=-1e-7&&u<=1+1e-7)hits.push_back(t);
        }
        std::sort(hits.begin(),hits.end());hits.erase(std::unique(hits.begin(),hits.end(),[](double a,double b){return std::abs(a-b)<.001;}),hits.end());
        return hits.size()==1?hits[0]:-1.;
    };
    std::vector<meshSymmetry::Point> half;
    for(double angle:angles) {
        const double x=std::sin(angle),y=std::cos(angle),a=radius(x,y),b=radius(-x,y);
        if(a<0||b<0)return {};
        const double r=outside?std::max(a,b):std::min(a,b);half.push_back({x*r,y*r});
    }
    // Simplify the half curve once, then mirror exactly. End cap notches are
    // retained as turns rather than flattened by a horizontal slice envelope.
    std::vector<uint8_t> keep(half.size());keep.front()=keep.back()=1;
    std::vector<std::pair<size_t,size_t>> stack{{0,half.size()-1}};
    while(!stack.empty()) {
        const auto [a,b]=stack.back();stack.pop_back();const double dx=half[b].x-half[a].x,dy=half[b].y-half[a].y,squared=dx*dx+dy*dy;
        double worst=std::clamp(spacing*.005,.25,1.5);size_t at=a;
        for(size_t i=a+1;i<b;++i) {
            const double x=half[i].x-half[a].x,y=half[i].y-half[a].y,t=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
            const double d=std::hypot(x-t*dx,y-t*dy);if(d>worst){worst=d;at=i;}
        }
        if(at!=a){keep[at]=1;stack.emplace_back(a,at);stack.emplace_back(at,b);}
    }
    std::vector<meshSymmetry::Point> sampled;size_t previous=0;
    for(size_t i=1;i<half.size();++i)if(keep[i]) {
        const auto a=half[previous],b=half[i];const int n=std::max(1,int(std::lround(std::hypot(b.x-a.x,b.y-a.y)/std::max(1.,spacing))));
        for(int k=0;k<n;++k)sampled.push_back({a.x+(b.x-a.x)*k/n,a.y+(b.y-a.y)*k/n});previous=i;
    }
    sampled.push_back(half.back());std::vector<float> out;
    auto append=[&](double x,double y){const auto q=axis.world(x,y);out.push_back(float(q.x));out.push_back(float(q.y));};
    for(auto q:sampled)append(q.x,q.y);
    for(size_t i=sampled.size()-1;i>1;--i)append(-sampled[i-1].x,sampled[i-1].y);
    return out;
}
struct ClosedCurve {
    const std::vector<float>& points;
    std::vector<double> arc;
    explicit ClosedCurve(const std::vector<float>& p):points(p),arc(p.size()/2+1) {
        for(size_t i=0;i<p.size()/2;++i){const size_t j=(i+1)%(p.size()/2);arc[i+1]=arc[i]+std::hypot(p[2*j]-p[2*i],p[2*j+1]-p[2*i+1]);}
    }
    CurvePoint at(double distance) const {
        distance=std::fmod(distance,arc.back());if(distance<0)distance+=arc.back();
        const size_t i=std::min(arc.size()-2,size_t(std::upper_bound(arc.begin(),arc.end(),distance)-arc.begin()-1));
        const size_t j=(i+1)%(points.size()/2);const double f=(distance-arc[i])/std::max(1e-12,arc[i+1]-arc[i]);
        return {points[2*i]+f*(points[2*j]-points[2*i]),points[2*i+1]+f*(points[2*j+1]-points[2*i+1])};
    }
};
double curveCross(CurvePoint a,CurvePoint b,CurvePoint p) {
    return (b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x);
}

// Fit the actual contour with a periodic cubic B-spline at raster scale.
// Mesh density changes the sampling interval, never the amount of outward
// displacement. Refine individual chords until their curve error is small.
std::vector<float> smoothBoundary(const std::vector<float>& raw,int target,int cell) {
    if(raw.size()/2<=8)return raw; // intentional polygon corners / tiny islands
    const ClosedCurve original(raw);const double perimeter=original.arc.back();
    auto clean=[](std::vector<float>& polygon) {
        for(size_t i=0;i<polygon.size()/2&&polygon.size()>6;) {
            const size_t j=(i+1)%(polygon.size()/2);
            if(std::hypot(polygon[2*i]-polygon[2*j],polygon[2*i+1]-polygon[2*j+1])<1e-5)
                polygon.erase(polygon.begin()+2*i,polygon.begin()+2*i+2);
            else ++i;
        }
        for(size_t i=0;i<polygon.size()/2&&polygon.size()>6;) {
            const size_t n=polygon.size()/2,a=(i+n-1)%n,b=(i+1)%n;
            const double dx=polygon[2*b]-polygon[2*a],dy=polygon[2*b+1]-polygon[2*a+1];
            const double ex=polygon[2*i]-polygon[2*a],ey=polygon[2*i+1]-polygon[2*a+1];
            const double length=std::hypot(dx,dy);
            if(std::abs(ex*dy-ey*dx)<=1e-6*length&&ex*dx+ey*dy>=0&&ex*dx+ey*dy<=length*length) {
                polygon.erase(polygon.begin()+2*i,polygon.begin()+2*i+2);if(i)--i;
            } else ++i;
        }
    };
    auto encloses=[&](const std::vector<float>& polygon) {
        VertexLoop loop(polygon.size()/2);std::iota(loop.begin(),loop.end(),0);
        if(loop.size()<3||!simpleOutline(polygon,loop,1e-6)||loopArea(polygon,loop)<=0)return false;
        for(size_t i=0;i<raw.size()/2;++i) {
            const size_t j=(i+1)%(raw.size()/2);
            const CurvePoint a{raw[2*i],raw[2*i+1]},b{raw[2*j],raw[2*j+1]};
            if(!pointInLoop(polygon,loop,float(a.x),float(a.y)))return false;
            for(size_t k=0;k<polygon.size()/2;++k) {
                const size_t l=(k+1)%(polygon.size()/2);
                const CurvePoint c{polygon[2*k],polygon[2*k+1]},d{polygon[2*l],polygon[2*l+1]};
                if(curveCross(a,b,c)*curveCross(a,b,d)<-1e-10&&curveCross(c,d,a)*curveCross(c,d,b)<-1e-10)return false;
            }
        }
        return true;
    };
    const int controls=std::clamp(int(std::ceil(perimeter/(cell*3.))),12,4096);
    std::vector<CurvePoint> control;control.reserve(size_t(controls));
    for(int i=0;i<controls;++i)control.push_back(original.at(perimeter*i/controls));
    auto fittedAt=[&](double at) {
        at=std::fmod(at,double(controls));if(at<0)at+=controls;
        const int i=int(at);const double t=at-i,t2=t*t,t3=t2*t;
        const double basis[4]={std::pow(1-t,3)/6,(3*t3-6*t2+4)/6,(-3*t3+3*t2+3*t+1)/6,t3/6};
        CurvePoint point;
        for(int k=0;k<4;++k){const auto p=control[size_t((i+k-1+controls)%controls)];point.x+=basis[k]*p.x;point.y+=basis[k]*p.y;}
        return point;
    };
    // Outward offsets can form tiny loops at a sharp concavity. Cut those
    // loops at their actual intersection, retaining the larger closed region.
    // The full enclosing test below still decides whether the fit is usable.
    auto untangle=[&](std::vector<float>& polygon) {
        for(int pass=0;pass<64;++pass) {
            bool changed=false;const size_t count=polygon.size()/2;
            for(size_t i=0;i<count&&!changed;++i)for(size_t j=i+2;j<count&&!changed;++j) {
                const size_t k=(i+1)%count,l=(j+1)%count;if(l==i)continue;
                const CurvePoint a{polygon[2*i],polygon[2*i+1]},b{polygon[2*k],polygon[2*k+1]},c{polygon[2*j],polygon[2*j+1]},d{polygon[2*l],polygon[2*l+1]};
                const double rx=b.x-a.x,ry=b.y-a.y,sx=d.x-c.x,sy=d.y-c.y,den=rx*sy-ry*sx;
                if(std::abs(den)<1e-10)continue;
                const double t=((c.x-a.x)*sy-(c.y-a.y)*sx)/den,u=((c.x-a.x)*ry-(c.y-a.y)*rx)/den;
                if(t<=1e-6||t>=1-1e-6||u<=1e-6||u>=1-1e-6)continue;
                const float x=float(a.x+t*rx),y=float(a.y+t*ry);
                auto part=[&](size_t begin,size_t end) {
                    std::vector<float> path{x,y};
                    for(size_t v=begin;v!=end;v=(v+1)%count){path.push_back(polygon[2*v]);path.push_back(polygon[2*v+1]);}
                    return path;
                };
                auto one=part(k,l),two=part(l,k);
                VertexLoop p(one.size()/2),q(two.size()/2);std::iota(p.begin(),p.end(),0);std::iota(q.begin(),q.end(),0);
                polygon=std::abs(loopArea(one,p))>std::abs(loopArea(two,q))?std::move(one):std::move(two);changed=true;
            }
            if(!changed)break;
        }
    };
    const int dense=controls*4;
    // An angular fork may make the smooth normal offset self-intersect.
    // Prepare an error-bounded polygonal fit as an alternate curve, preserving
    // the fork instead of reverting the whole outline to pixel stair steps.
    auto angularFit=[&](double padding) {
        const size_t count=raw.size()/2;std::vector<uint8_t> keep(count);keep[0]=1;
        size_t opposite=1;double farthest=0;
        for(size_t i=1;i<count;++i){const double d=std::hypot(raw[2*i]-raw[0],raw[2*i+1]-raw[1]);if(d>farthest){farthest=d;opposite=i;}}
        keep[opposite]=1;std::vector<std::pair<size_t,size_t>> stack{{0,opposite},{opposite,count}};
        while(!stack.empty()) {
            const auto [a,b]=stack.back();stack.pop_back();const size_t end=b%count;
            const double dx=raw[2*end]-raw[2*a],dy=raw[2*end+1]-raw[2*a+1],squared=dx*dx+dy*dy;
            double error=.65*cell;size_t at=count;
            for(size_t i=a+1;i<b;++i) {
                const double x=raw[2*i]-raw[2*a],y=raw[2*i+1]-raw[2*a+1],t=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                const double d=std::hypot(x-t*dx,y-t*dy);if(d>error){error=d;at=i;}
            }
            if(at<count){keep[at]=1;stack.emplace_back(a,at);stack.emplace_back(at,b);}
        }
        std::vector<CurvePoint> polygon;
        for(size_t i=0;i<count;++i)if(keep[i])polygon.push_back({raw[2*i],raw[2*i+1]});
        std::vector<float> fit;if(polygon.size()<3)return fit;
        for(size_t i=0;i<polygon.size();++i) {
            const auto a=polygon[(i+polygon.size()-1)%polygon.size()],b=polygon[i],c=polygon[(i+1)%polygon.size()];
            const double l0=std::hypot(b.x-a.x,b.y-a.y),l1=std::hypot(c.x-b.x,c.y-b.y);
            if(l0<1e-8||l1<1e-8)return std::vector<float>{};
            const CurvePoint n0{(b.y-a.y)/l0,-(b.x-a.x)/l0},n1{(c.y-b.y)/l1,-(c.x-b.x)/l1};
            const double den=1+n0.x*n1.x+n0.y*n1.y;
            if(den<.25)return std::vector<float>{}; // no unbounded miters at a cusp
            const double x=(n0.x+n1.x)*padding*cell/den,y=(n0.y+n1.y)*padding*cell/den;
            if(std::hypot(x,y)>padding*cell*2)return std::vector<float>{};
            fit.push_back(float(b.x+x));fit.push_back(float(b.y+y));
        }
        return fit;
    };
    for(double padding:{1.25,1.75,2.25}) {
        std::vector<float> fitted;fitted.reserve(size_t(dense)*2);
        for(int i=0;i<dense;++i) {
            const double at=double(i)*controls/dense;
            auto point=fittedAt(at);const auto a=fittedAt(at-.04),b=fittedAt(at+.04);
            const double dx=b.x-a.x,dy=b.y-a.y,length=std::hypot(dx,dy);
            if(length>1e-8){point.x+=dy/length*padding*cell;point.y-=dx/length*padding*cell;}
            fitted.push_back(float(point.x));
            fitted.push_back(float(point.y));
        }
        clean(fitted);
        untangle(fitted);
        if(!encloses(fitted)) {
            auto alternate=angularFit(padding);
            untangle(alternate);
            if(alternate.empty()||!encloses(alternate))continue;
            fitted=std::move(alternate);
        }
        const ClosedCurve curve(fitted);const double length=curve.arc.back();
        auto weight=[&](double at) {
            const double window=std::max(cell*2.,std::min(cell*6.,length/target*.5));
            const auto a=curve.at(at-window),b=curve.at(at),c=curve.at(at+window);
            const double ax=b.x-a.x,ay=b.y-a.y,bx=c.x-b.x,by=c.y-b.y,product=std::hypot(ax,ay)*std::hypot(bx,by);
            const double angle=product>1e-10?std::acos(std::clamp((ax*bx+ay*by)/product,-1.,1.)):0;
            return 1+std::min(.8,angle*.7);
        };
        std::vector<double> samples;
        auto subdivide=[&](auto&& self,double start,double end,int depth)->void {
            const auto a=curve.at(start),b=curve.at(end);
            const double dx=b.x-a.x,dy=b.y-a.y,squared=dx*dx+dy*dy;
            double error=0;
            for(size_t i=size_t(std::upper_bound(curve.arc.begin(),curve.arc.end(),start)-curve.arc.begin());i+1<curve.arc.size()&&curve.arc[i]<end;++i) {
                const double x=fitted[2*i]-a.x,y=fitted[2*i+1]-a.y;
                const double f=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                error=std::max(error,std::hypot(x-f*dx,y-f*dy));
            }
            if(error>.6*cell&&end-start>cell*1.5&&depth<12) {
                const double middle=(start+end)*.5;self(self,start,middle,depth+1);self(self,middle,end,depth+1);
            } else samples.push_back(start);
        };
        const int measureCount=std::clamp(target*16,128,4096);
        std::vector<double> measure(size_t(measureCount)+1),starts(size_t(target)+1);
        for(int i=0;i<measureCount;++i)measure[size_t(i)+1]=measure[size_t(i)]+length/measureCount*weight(length*(i+.5)/measureCount);
        for(int i=0;i<target;++i) {
            const double at=measure.back()*i/target;
            const size_t k=std::min(size_t(measureCount-1),size_t(std::upper_bound(measure.begin(),measure.end(),at)-measure.begin()-1));
            starts[size_t(i)]=length/measureCount*(k+(at-measure[k])/(measure[k+1]-measure[k]));
        }
        starts.back()=length;
        for(int i=0;i<target;++i)subdivide(subdivide,starts[size_t(i)],starts[size_t(i)+1],0);
        for(int repair=0;repair<8;++repair) {
            std::sort(samples.begin(),samples.end());
            samples.erase(std::unique(samples.begin(),samples.end(),[](double a,double b){return std::abs(a-b)<1e-5;}),samples.end());
            std::vector<float> polygon;
            for(double at:samples){const auto p=curve.at(at);polygon.push_back(float(p.x));polygon.push_back(float(p.y));}
            clean(polygon);
            if(encloses(polygon)) {
                // Coverage repairs can leave a point almost on top of an
                // existing sample. Remove these short spans only while the
                // replacement chord stays close to the fitted curve and
                // still encloses the complete original raster contour.
                for(int pass=0;pass<512&&polygon.size()/2>size_t(target);++pass) {
                    std::vector<std::pair<double,size_t>> order;
                    const size_t n=polygon.size()/2;
                    for(size_t i=0;i<n;++i) {
                        const size_t a=(i+n-1)%n,b=(i+1)%n;
                        const double shortEdge=std::min(std::hypot(polygon[2*i]-polygon[2*a],polygon[2*i+1]-polygon[2*a+1]),
                            std::hypot(polygon[2*i]-polygon[2*b],polygon[2*i+1]-polygon[2*b+1]));
                        if(shortEdge<length/target*.5)order.emplace_back(shortEdge,i);
                    }
                    std::sort(order.begin(),order.end());bool removed=false;
                    for(auto [distance,i]:order) {
                        const size_t a=(i+n-1)%n,b=(i+1)%n;
                        const double dx=polygon[2*b]-polygon[2*a],dy=polygon[2*b+1]-polygon[2*a+1],squared=dx*dx+dy*dy;
                        const double x=polygon[2*i]-polygon[2*a],y=polygon[2*i+1]-polygon[2*a+1];
                        const double f=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                        if(std::hypot(x-f*dx,y-f*dy)>.75*cell)continue;
                        auto proposal=polygon;proposal.erase(proposal.begin()+2*i,proposal.begin()+2*i+2);
                        if(encloses(proposal)){polygon=std::move(proposal);removed=true;break;}
                    }
                    if(!removed)break;
                }
                // Balance arc intervals on the curve itself. Moving a point
                // is accepted only if the whole source outline stays covered.
                // Thus no averaged point can drift into the texture or create
                // an artificial bend between independently offset straight lines.
                std::vector<double> arcs(polygon.size()/2);
                for(size_t i=0;i<arcs.size();++i) {
                    double nearest=1e30;
                    for(size_t k=0;k<fitted.size()/2;++k) {
                        const size_t l=(k+1)%(fitted.size()/2);
                        const double dx=fitted[2*l]-fitted[2*k],dy=fitted[2*l+1]-fitted[2*k+1],squared=dx*dx+dy*dy;
                        const double f=squared?std::clamp(((polygon[2*i]-fitted[2*k])*dx+(polygon[2*i+1]-fitted[2*k+1])*dy)/squared,0.,1.):0;
                        const double distance=std::hypot(polygon[2*i]-fitted[2*k]-f*dx,polygon[2*i+1]-fitted[2*k+1]-f*dy);
                        if(distance<nearest){nearest=distance;arcs[i]=curve.arc[k]+f*std::sqrt(squared);}
                    }
                }
                for(int pass=0;pass<6;++pass)for(size_t i=0;i<arcs.size();++i) {
                    const size_t a=(i+arcs.size()-1)%arcs.size(),b=(i+1)%arcs.size();
                    double previous=arcs[a],next=arcs[b];while(previous>=arcs[i])previous-=length;while(next<=arcs[i])next+=length;
                    const double w0=weight((previous+arcs[i])*.5),w1=weight((next+arcs[i])*.5);
                    const double desired=previous+(next-previous)*w1/(w0+w1);
                    for(double fraction:{.5,.25,.125}) {
                        double at=arcs[i]+fraction*(desired-arcs[i]);at=std::fmod(at,length);if(at<0)at+=length;
                        const auto point=curve.at(at);auto proposal=polygon;proposal[2*i]=float(point.x);proposal[2*i+1]=float(point.y);
                        if(encloses(proposal)){polygon=std::move(proposal);arcs[i]=at;break;}
                    }
                }
                // Cleanup is for fitting only. Restore intermediate stations
                // on long straight spans after all simplification/relaxation;
                // otherwise a straight run collapses into one very long
                // edge and produces stretched adjacent triangles.
                const double interval=std::max(double(cell)*2,length/target);
                std::vector<float> spaced;
                for(size_t i=0;i<polygon.size()/2;++i) {
                    const size_t j=(i+1)%(polygon.size()/2);
                    const double dx=polygon[2*j]-polygon[2*i],dy=polygon[2*j+1]-polygon[2*i+1];
                    const int pieces=std::max(1,int(std::ceil(std::hypot(dx,dy)/(interval*1.35))));
                    for(int k=0;k<pieces;++k){spaced.push_back(float(polygon[2*i]+dx*k/pieces));spaced.push_back(float(polygon[2*i+1]+dy*k/pieces));}
                }
                if(std::getenv("L2DR_AUTOMESH_SMOOTH_PROFILE"))std::printf("[automesh-curve] raw=%zu target=%d fitted=%zu samples=%zu padding=%.2f\n",raw.size()/2,target,fitted.size()/2,spaced.size()/2,padding*cell);
                return spaced;
            }
            if(samples.size()>1600)break;
            const size_t previous=samples.size();VertexLoop loop(polygon.size()/2);std::iota(loop.begin(),loop.end(),0);
            // A coverage correction adds a point on the fitted curve near
            // only the offending source span. It never pushes a whole line
            // away from the contour or intersects extended neighboring lines.
            for(size_t i=0;i<raw.size()/2;++i)for(double f:{0.,.5}) {
                const size_t j=(i+1)%(raw.size()/2);
                const double x=raw[2*i]+f*(raw[2*j]-raw[2*i]),y=raw[2*i+1]+f*(raw[2*j+1]-raw[2*i+1]);
                if(pointInLoop(polygon,loop,float(x),float(y)))continue;
                size_t best=0;double nearest=std::numeric_limits<double>::max();
                for(size_t k=0;k<fitted.size()/2;++k){const double distance=std::hypot(fitted[2*k]-x,fitted[2*k+1]-y);if(distance<nearest){nearest=distance;best=k;}}
                samples.push_back(curve.arc[best]);
            }
            if(samples.size()==previous)break;
        }
    }
    // A sharp fork can reject all smooth outward fits. Keep its enclosing
    // raster contour, but remove staircase concavities and collinear stations
    // conservatively. Otherwise a single difficult corner would turn the
    // entire layer into hundreds of one-pixel edges, regardless of density.
    auto fallback=raw;
    std::vector<size_t> sourceStation(raw.size()/2);std::iota(sourceStation.begin(),sourceStation.end(),0);
    for(int pass=0;pass<64;++pass) {
        bool changed=false;
        for(size_t i=0;i<fallback.size()/2&&fallback.size()/2>size_t(target);) {
            const size_t count=fallback.size()/2,a=(i+count-1)%count,b=(i+1)%count;
            const CurvePoint x{fallback[2*a],fallback[2*a+1]},y{fallback[2*i],fallback[2*i+1]},z{fallback[2*b],fallback[2*b+1]};
            const double dx=z.x-x.x,dy=z.y-x.y,squared=dx*dx+dy*dy;
            const double t=squared?std::clamp(((y.x-x.x)*dx+(y.y-x.y)*dy)/squared,0.,1.):0;
            if(curveCross(x,y,z)>1e-8||std::hypot(y.x-x.x-t*dx,y.y-x.y-t*dy)>.75*cell){++i;continue;}
            bool tooFar=false;
            for(size_t k=(sourceStation[a]+1)%(raw.size()/2);k!=sourceStation[b];k=(k+1)%(raw.size()/2)) {
                const double rx=raw[2*k]-x.x,ry=raw[2*k+1]-x.y,f=squared?std::clamp((rx*dx+ry*dy)/squared,0.,1.):0;
                if(std::hypot(rx-f*dx,ry-f*dy)>.75*cell){tooFar=true;break;}
            }
            if(tooFar){++i;continue;}
            bool crossed=false;
            for(size_t k=0;k<count&&!crossed;++k) {
                const size_t l=(k+1)%count;if(k==a||k==i||k==b||l==a)continue;
                const CurvePoint u{fallback[2*k],fallback[2*k+1]},v{fallback[2*l],fallback[2*l+1]};
                if(std::max(std::min(x.x,z.x),std::min(u.x,v.x))<=std::min(std::max(x.x,z.x),std::max(u.x,v.x))&&
                   std::max(std::min(x.y,z.y),std::min(u.y,v.y))<=std::min(std::max(x.y,z.y),std::max(u.y,v.y))&&
                   curveCross(x,z,u)*curveCross(x,z,v)<=0&&curveCross(u,v,x)*curveCross(u,v,z)<=0)crossed=true;
            }
            if(crossed){++i;continue;}
            fallback.erase(fallback.begin()+2*i,fallback.begin()+2*i+2);sourceStation.erase(sourceStation.begin()+i);changed=true;
        }
        if(!changed)break;
    }
    if(!encloses(fallback))fallback=raw;
    std::vector<float> spaced;const ClosedCurve curve(fallback);const double interval=std::max(double(cell)*2,curve.arc.back()/target);
    for(size_t i=0;i<fallback.size()/2;++i) {
        const size_t j=(i+1)%(fallback.size()/2);const double dx=fallback[2*j]-fallback[2*i],dy=fallback[2*j+1]-fallback[2*i+1];
        const int pieces=std::max(1,int(std::ceil(std::hypot(dx,dy)/(interval*1.35))));
        for(int k=0;k<pieces;++k){spaced.push_back(float(fallback[2*i]+dx*k/pieces));spaced.push_back(float(fallback[2*i+1]+dy*k/pieces));}
    }
    if(std::getenv("L2DR_AUTOMESH_SMOOTH_PROFILE"))std::printf("[automesh-curve] raw=%zu target=%d fallback=%zu\n",raw.size()/2,target,spaced.size()/2);
    return spaced;
}

// A tapered stroke ends in one apex. Disk dilation otherwise rounds that apex
// into a cap and curvature sampling adds a needless fan of boundary vertices.
// Detect taper in the source alpha, not in the padded outline, and accept a
// replacement only if every affected native pixel corner remains enclosed.
void collapseSharpStripTips(TextureDomain& domain,double spacing,double margin) {
    const auto& guide=domain.strokeCenterline.empty()?domain.stripCenterline:domain.strokeCenterline;
    if(guide.size()<4||domain.contours.size()!=1)return;
    std::vector<double> arc(guide.size()/2);
    for(size_t i=1;i<arc.size();++i)arc[i]=arc[i-1]+std::hypot(guide[2*i]-guide[2*i-2],guide[2*i+1]-guide[2*i-1]);
    if(arc.back()<8)return;
    const double cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
    auto at=[&](double s) {
        s=std::clamp(s,0.,arc.back());
        const size_t i=std::clamp(size_t(std::upper_bound(arc.begin(),arc.end(),s)-arc.begin()),size_t(1),arc.size()-1);
        const double f=(s-arc[i-1])/std::max(1e-12,arc[i]-arc[i-1]);
        return CurvePoint{guide[2*i-2]+f*(guide[2*i]-guide[2*i-2]),guide[2*i-1]+f*(guide[2*i+1]-guide[2*i-1])};
    };
    for(int end=0;end<2;++end) {
        const double window=std::min(arc.back()*.2,std::max(spacing*.65,margin*2+cell*4));
        if(window<cell*6)continue;
        const auto p=at(end?arc.back():0),q=at(end?arc.back()-window*.25:window*.25);
        const double length=std::hypot(p.x-q.x,p.y-q.y);if(length<cell)continue;
        const CurvePoint d{(p.x-q.x)/length,(p.y-q.y)/length},n{-d.y,d.x};
        std::vector<CurvePoint> pixels;double maximum=-1e30;
        const int x0=std::max(0,int(std::floor(p.x-window*1.5))),x1=std::min(domain.width-1,int(std::ceil(p.x+window*1.5)));
        const int y0=std::max(0,int(std::floor(p.y-window*1.5))),y1=std::min(domain.height-1,int(std::ceil(p.y+window*1.5)));
        for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x)if(domain.rgba[(size_t(y)*domain.width+x)*4+3]) {
            const double dx=x+.5-p.x,dy=y+.5-p.y;
            if(std::hypot(dx,dy)>window*1.5)continue;
            pixels.push_back({dx*d.x+dy*d.y,dx*n.x+dy*n.y});maximum=std::max(maximum,pixels.back().x);
        }
        if(pixels.empty())continue;
        const double first=window*.22,last=window*.6,band=std::max(cell,window*.06);
        auto width=[&](double depth) {
            double lo=1e30,hi=-1e30;
            for(auto v:pixels)if(std::abs(maximum-v.x-depth)<band){lo=std::min(lo,v.y);hi=std::max(hi,v.y);}
            return hi>=lo?hi-lo+1:0.;
        };
        const double narrow=width(first),wide=width(last);
        if(narrow<=0||wide<=narrow*1.3+cell*.5||narrow>first*1.2||wide>last*1.2)continue;
        double lateral=0;int count=0;
        for(auto v:pixels)if(v.x>maximum-cell*1.5){lateral+=v.y;++count;}
        lateral/=std::max(1,count);
        maximum+=.5*(std::abs(d.x)+std::abs(d.y));
        auto& poly=domain.contours[0];const size_t size=poly.size()/2;
        auto projection=[&](size_t i){return (poly[2*i]-p.x)*d.x+(poly[2*i+1]-p.y)*d.y-maximum;};
        size_t peak=size;double highest=-1e30;
        for(size_t i=0;i<size;++i)if(std::hypot(poly[2*i]-p.x,poly[2*i+1]-p.y)<window*1.5+margin&&projection(i)>highest){peak=i;highest=projection(i);}
        if(peak==size)continue;
        bool collapsed=false;
        for(double fraction:{.55,.35,.2}) {
        const double cut=-std::min(window*fraction,std::max(cell*3,spacing*.45));
        size_t begin=peak,finish=peak;int removed=1;
        while(removed<int(size)-3&&projection((begin+size-1)%size)>cut){begin=(begin+size-1)%size;++removed;}
        while(removed<int(size)-3&&projection((finish+1)%size)>cut){finish=(finish+1)%size;++removed;}
        if(removed<2||removed>=int(size)-3)continue;
        const size_t before=(begin+size-1)%size,after=(finish+1)%size;
        const double affected=std::min(projection(before),projection(after))-2;
        for(double extension:{std::max(cell*1.5,std::min(margin*.35,spacing*.12)),std::max(cell*2,margin+cell)}) {
            const CurvePoint tip{p.x+d.x*(maximum+extension)+n.x*lateral,p.y+d.y*(maximum+extension)+n.y*lateral};
            std::vector<float> proposal{float(tip.x),float(tip.y)};
            for(size_t i=after;i!=begin;i=(i+1)%size){proposal.push_back(poly[2*i]);proposal.push_back(poly[2*i+1]);}
            VertexLoop loop(proposal.size()/2);std::iota(loop.begin(),loop.end(),0);
            if(!simpleOutline(proposal,loop,1e-6))continue;
            bool covered=true;
            for(int y=0;y<domain.height&&covered;++y)for(int x=0;x<domain.width&&covered;++x)if(domain.rgba[(size_t(y)*domain.width+x)*4+3]) {
                if((x+.5-p.x)*d.x+(y+.5-p.y)*d.y-maximum<affected-1)continue;
                for(auto offset:std::initializer_list<CurvePoint>{{0,0},{1,0},{0,1},{1,1}})
                    if(!pointInLoop(proposal,loop,float(x+offset.x),float(y+offset.y))){covered=false;break;}
            }
            if(covered){poly=std::move(proposal);++domain.sharpTips;collapsed=true;break;}
        }
        if(collapsed)break;
        }
    }
}

// Keep every occupied pixel cell, including alpha=1 and isolated pixels. OR
// reduction is conservative; point sampling/downscaling could erase thin tips.
bool textureOutline(const Cmo3Document& doc,const ArtMesh& mesh,
                    const MeshAutoGenParams& params,TextureDomain& domain,std::string& error) {
    const auto* texture=mesh.source?textureUvDefinition(doc.xml(),mesh.source->field("texture")):nullptr;
    const auto* resource=texture?textureUvDefinition(doc.xml(),texture->field("srcImageResource")):nullptr;
    int textureIndex=-1;
    if(resource)for(size_t i=0;i<doc.textures().size();++i)
        if(doc.textures()[i].imageResXsId==resource->id){textureIndex=int(i);break;}
    if(textureIndex<0)textureIndex=mesh.textureIndex;
    std::vector<uint8_t> png;PngDecodeResult decoded;
    if(textureIndex<0){error="网格缺少源纹理";return false;}
    // Resolve atlas ownership before decoding: the full page can be enormous,
    // while the authoritative per-layer source image is usually small.
    bool atlasResource=false;
    auto inspectAtlas=[&](auto&& self,const XsNode& node)->void {
        if(node.type=="CTextureAtlas") {
            const auto* cached=textureUvDefinition(doc.xml(),node.field("cachedAtlasImage"));
            if(cached&&resource&&cached->id==resource->id)atlasResource=true;
        }
        for(const auto& child:node.children)self(self,child);
    };
    inspectAtlas(inspectAtlas,doc.xml().root);

    TextureUvMapping mapping;
    const bool mapped=mesh.source&&sourceTextureUvMapping(doc.xml(),*mesh.source,mapping);
    if(!atlasResource) {
        if(!doc.extractTexture(textureIndex,png,error)||!pngDecodeToRgba(png.data(),png.size(),domain.rgba,&decoded,error)) {
            error="自动网格需要可读取的图层源纹理："+error;return false;
        }
        domain.width=decoded.width;domain.height=decoded.height;
    }
    const int imageWidth=atlasResource?(resource?resource->attrInt("width",0):0):domain.width;
    const int imageHeight=atlasResource?(resource?resource->attrInt("height",0):0):domain.height;
    domain.logicalW=imageWidth*(mapped?mapping.scaleU:1.f);
    domain.logicalH=imageHeight*(mapped?mapping.scaleV:1.f);
    if(!(domain.logicalW>0)||!(domain.logicalH>0)){error="源纹理尺寸无效";return false;}
    // pngDecodeToRgba returns unassociated RGBA. The GTexture2D flag describes
    // GPU upload; TextureCache premultiplies explicitly after PNG decoding.
    domain.premultiplied=false;
    domain.toUv[0]=1./domain.logicalW;domain.toUv[4]=1./domain.logicalH;
    if(atlasResource) {
        const XsNode* extension=nullptr;
        auto findExtension=[&](auto&& self,const XsNode& node)->void {
            if(node.type=="CTextureInputExtension")extension=textureUvDefinition(doc.xml(),&node);
            if(!extension)for(const auto& child:node.children)self(self,child);
        };
        findExtension(findExtension,*mesh.source);
        const XsNode* guid=nullptr;
        const auto* inputs=extension?textureUvDefinition(doc.xml(),extension->field("_textureInputs")):nullptr;
        if(inputs)for(const auto& item:inputs->children) {
            const auto* input=textureUvDefinition(doc.xml(),&item);
            if(input&&input->type=="CTextureInput_ModelImage")guid=textureUvDefinition(doc.xml(),input->field("_modelImageGuid"));
        }
        const XsNode *modelImage=nullptr,*placement=nullptr;int placements=0;
        auto findImage=[&](auto&& self,const XsNode& node)->void {
            if(node.ref<0&&(node.type=="CModelImage"||node.type=="ModelImageEntry")) {
                const auto* ownGuid=textureUvDefinition(doc.xml(),node.field(node.type=="CModelImage"?"guid":"modelImageGuid"));
                if(guid&&ownGuid&&!guid->uuid.empty()&&guid->uuid==ownGuid->uuid) {
                    if(node.type=="CModelImage")modelImage=&node;
                    else {
                        const auto* atlas=textureUvDefinition(doc.xml(),node.field("atlas"));
                        const auto* image=atlas?textureUvDefinition(doc.xml(),atlas->field("cachedAtlasImage")):nullptr;
                        if(image&&resource&&image->id==resource->id){placement=&node;++placements;}
                    }
                }
            }
            for(const auto& child:node.children)self(self,child);
        };
        findImage(findImage,doc.xml().root);
        const auto* image=modelImage?textureUvDefinition(doc.xml(),modelImage->field("_filteredImage")):nullptr;
        if(!image&&modelImage)image=textureUvDefinition(doc.xml(),modelImage->field("rawImage"));
        const auto* file=image?image->field("imageFileBuf"):nullptr;
        const auto* entry=file?doc.archive().findByPath(file->attr("path")):nullptr;
        const auto* transform=placement?placement->field("materialLocalToAtlasTransform"):nullptr;
        if(!entry||!transform||placements!=1||!CaffReader::extract(doc.archive(),*entry,png,error)||
           !pngDecodeToRgba(png.data(),png.size(),domain.rgba,&decoded,error)) {
            error="无法唯一对应图集中的完整图层源图："+error;return false;
        }
        domain.width=decoded.width;domain.height=decoded.height;
        auto fieldFloat=[](const XsNode* n,const char* field,double fallback) {
            const auto* f=n?n->field(field):nullptr;return f?std::strtod(f->text.c_str(),nullptr):fallback;
        };
        const auto* position=transform->field("position");const auto* scaling=transform->field("scale");
        const double tx=fieldFloat(position,"x",0),ty=fieldFloat(position,"y",0);
        const double sx=fieldFloat(scaling,"x",1),sy=fieldFloat(scaling,"y",1);
        const double angle=fieldFloat(transform,"eulerAngle",0)*3.14159265358979323846/180.;
        const double c=std::cos(angle),si=std::sin(angle);
        domain.toUv[0]=c*sx/domain.logicalW;domain.toUv[1]=-si*sy/domain.logicalW;domain.toUv[2]=tx/domain.logicalW;
        domain.toUv[3]=si*sx/domain.logicalH;domain.toUv[4]=c*sy/domain.logicalH;domain.toUv[5]=ty/domain.logicalH;
    }
    for(double value:domain.toUv)if(!std::isfinite(value)){error="图层纹理映射坐标无效";return false;}
    const auto* affine=domain.toUv;
    const double determinant=affine[0]*affine[4]-affine[1]*affine[3];
    if(!std::isfinite(determinant)||std::abs(determinant)<1e-15){error="图层纹理映射退化";return false;}
    domain.fromUv[0]=affine[4]/determinant;domain.fromUv[1]=-affine[1]/determinant;
    domain.fromUv[3]=-affine[3]/determinant;domain.fromUv[4]=affine[0]/determinant;
    domain.fromUv[2]=-(domain.fromUv[0]*affine[2]+domain.fromUv[1]*affine[5]);
    domain.fromUv[5]=-(domain.fromUv[3]*affine[2]+domain.fromUv[4]*affine[5]);
    if(!buildTextureContours(domain,params,error))return false;
    return !params.generateHoles||prepareTextureHoles(domain,params,error);
}

// Both sides of the border band come from the same alpha silhouette. Positive
// inset erodes that silhouette before fitting, rather than shrinking a polygon
// toward its center (which distorts concave shapes and narrow strips).
bool buildTextureContours(TextureDomain& domain,const MeshAutoGenParams& params,std::string& error,double inset) {
    const int cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
    const int sourceW=(domain.width+cell-1)/cell,sourceH=(domain.height+cell-1)/cell;
    const int radius=inset>0?int(std::min(double(std::max(sourceW,sourceH)),std::ceil(inset/cell))):
        int(std::min<int64_t>(std::max(sourceW,sourceH),(std::max<int64_t>(0,params.polygonMargin)+cell-1)/cell));
    // Extend the working raster before dilation. The source PNG bounds only
    // describe where pixels exist; they are not clipping planes for the mesh.
    // One empty cell beyond the requested disk keeps every offset contour closed.
    const int padding=radius+1,w=sourceW+2*padding,h=sourceH+2*padding;
    std::vector<uint8_t> mask(size_t(w)*h);
    for(int y=0;y<domain.height;++y)for(int x=0;x<domain.width;++x)
        if(domain.rgba[(size_t(y)*domain.width+x)*4+3])mask[size_t(y/cell+padding)*w+x/cell+padding]=1;
    if(std::find(mask.begin(),mask.end(),1)==mask.end()){error="图层源纹理完全透明";return false;}
    if(radius) {
        // Exact squared Euclidean distance via separable parabola envelopes.
        // A disk offset preserves round bends; a square max-filter distorts
        // the curvature differently along axes and diagonals.
        constexpr double infinity=1e12;
        std::vector<double> distance(mask.size()),temporary(mask.size());
        for(size_t i=0;i<mask.size();++i)distance[i]=(inset>0?!mask[i]:bool(mask[i]))?0:infinity;
        auto transform=[](const std::vector<double>& input,std::vector<double>& output,int n) {
            std::vector<int> sites(size_t(n),0);std::vector<double> cuts(size_t(n)+1);
            int last=0;cuts[0]=-1e30;cuts[1]=1e30;
            for(int q=1;q<n;++q) {
                double crossing;
                for(;;) {
                    const int p=sites[size_t(last)];
                    crossing=((input[size_t(q)]+double(q)*q)-(input[size_t(p)]+double(p)*p))/(2.*(q-p));
                    if(crossing>cuts[size_t(last)]||last==0)break;
                    --last;
                }
                sites[size_t(++last)]=q;cuts[size_t(last)]=crossing;cuts[size_t(last)+1]=1e30;
            }
            int k=0;
            for(int q=0;q<n;++q){while(cuts[size_t(k)+1]<q)++k;const int p=sites[size_t(k)];output[size_t(q)]=double(q-p)*(q-p)+input[size_t(p)];}
        };
        std::vector<double> input(size_t(std::max(w,h))),output(input.size());
        for(int y=0;y<h;++y) {
            for(int x=0;x<w;++x)input[size_t(x)]=distance[size_t(y)*w+x];
            transform(input,output,w);for(int x=0;x<w;++x)temporary[size_t(y)*w+x]=output[size_t(x)];
        }
        for(int x=0;x<w;++x) {
            for(int y=0;y<h;++y)input[size_t(y)]=temporary[size_t(y)*w+x];
            transform(input,output,h);for(int y=0;y<h;++y)mask[size_t(y)*w+x]=inset>0?
                output[size_t(y)]>double(radius)*radius:output[size_t(y)]<=double(radius)*radius;
        }
    }
    // Connect diagonal contacts with one extra empty cell. This avoids a
    // nonmanifold shared corner without removing either occupied pixel.
    for(int pass=0;pass<16;++pass) {
        bool changed=false;
        for(int y=0;y+1<h;++y)for(int x=0;x+1<w;++x) {
            const size_t a=size_t(y)*w+x,b=a+1,c=a+w,d=c+1;
            if(mask[a]&&mask[d]&&!mask[b]&&!mask[c]){mask[b]=1;changed=true;}
            if(mask[b]&&mask[c]&&!mask[a]&&!mask[d]){mask[a]=1;changed=true;}
        }
        if(!changed)break;
    }
    int minX=w,minY=h,maxX=0,maxY=0;
    const int stride=w+1;
    // Keep the last partial source cell exact while extending freely beyond
    // both image edges. This also preserves a zero-margin opaque rectangle.
    auto nativeCoordinate=[&](int grid,int cells,int extent) {
        const int local=grid-padding;
        return float(local<cells?local*cell:extent+(local-cells)*cell);
    };
    auto nativeX=[&](int x){return nativeCoordinate(x,sourceW,domain.width);};
    auto nativeY=[&](int y){return nativeCoordinate(y,sourceH,domain.height);};
    std::vector<uint8_t> edges(size_t(w+1)*(h+1));
    auto occupied=[&](int x,int y){return x>=0&&x<w&&y>=0&&y<h&&mask[size_t(y)*w+x];};
    for(int y=0;y<h;++y)for(int x=0;x<w;++x)if(occupied(x,y)) {
        minX=std::min(minX,x);minY=std::min(minY,y);maxX=std::max(maxX,x+1);maxY=std::max(maxY,y+1);
        if(!occupied(x,y-1))edges[size_t(y)*stride+x]|=1; // E
        if(!occupied(x+1,y))edges[size_t(y)*stride+x+1]|=2; // S
        if(!occupied(x,y+1))edges[size_t(y+1)*stride+x+1]|=4; // W
        if(!occupied(x-1,y))edges[size_t(y+1)*stride+x]|=8; // N
    }
    auto bboxFallback=[&] {
        domain.conservativeFallback=true;
        domain.holes.clear();
        domain.contours={{nativeX(minX),nativeY(minY),nativeX(maxX),nativeY(minY),
            nativeX(maxX),nativeY(maxY),nativeX(minX),nativeY(maxY)}};
    };
    const int walkDx[4]={1,0,-1,0},walkDy[4]={0,1,0,-1};
    size_t totalPoints=0;
    for(size_t start=0;start<edges.size();++start)while(edges[start]) {
        int dir=0;while(!(edges[start]&(1<<dir)))++dir;
        size_t at=start;std::vector<float> poly;
        bool closed=false;
        for(size_t walk=0;walk<edges.size()*4;++walk) {
            const int x=int(at%stride),y=int(at/stride);
            poly.push_back(nativeX(x));poly.push_back(nativeY(y));
            edges[at]&=uint8_t(~(1<<dir));
            const int nx=x+walkDx[dir],ny=y+walkDy[dir];at=size_t(ny)*stride+nx;
            if(at==start){closed=true;break;}
            int next=-1;for(int turn:{1,0,3,2})if(edges[at]&(1<<((dir+turn)%4))){next=(dir+turn)%4;break;}
            if(next<0)break;dir=next;
        }
        if(!closed){bboxFallback();return true;}
        VertexLoop ids(poly.size()/2);std::iota(ids.begin(),ids.end(),0);
        const bool hole=loopArea(poly,ids)<0;
        if(hole&&!params.generateHoles)continue;
        // Preserve curvature: discard only exact collinear raster corners.
        for(size_t i=0;i<poly.size()/2&&poly.size()>6;) {
            const size_t n=poly.size()/2,a=(i+n-1)%n,b=(i+1)%n;
            const double dx=poly[2*b]-poly[2*a],dy=poly[2*b+1]-poly[2*a+1];
            const double ex=poly[2*i]-poly[2*a],ey=poly[2*i+1]-poly[2*a+1];
            if(std::abs(ex*dy-ey*dx)<1e-8&&ex*dx+ey*dy>=0&&ex*dx+ey*dy<=dx*dx+dy*dy) {
                poly.erase(poly.begin()+2*i,poly.begin()+2*i+2);if(i)--i;
            } else ++i;
        }
        VertexLoop simplified(poly.size()/2);std::iota(simplified.begin(),simplified.end(),0);
        if(!simpleOutline(poly,simplified,1e-7)){bboxFallback();return true;}
        if(hole){domain.holes.push_back(std::move(poly));continue;}
        totalPoints+=poly.size()/2;domain.contours.push_back(std::move(poly));
        if(totalPoints>4000||domain.contours.size()>128){bboxFallback();return true;}
    }
    if(domain.contours.empty()){error="无法提取非透明纹理边界";return false;}
    for(size_t i=0;i<domain.contours.size();++i)for(size_t j=i+1;j<domain.contours.size();++j) {
        const auto& a=domain.contours[i];const auto& b=domain.contours[j];
        for(size_t u=0;u<a.size()/2;++u)for(size_t v=0;v<b.size()/2;++v) {
            const size_t un=(u+1)%(a.size()/2),vn=(v+1)%(b.size()/2);
            const double ax=a[2*u],ay=a[2*u+1],bx=a[2*un],by=a[2*un+1];
            const double cx=b[2*v],cy=b[2*v+1],dx=b[2*vn],dy=b[2*vn+1];
            auto cross=[](double x,double y,double xx,double yy,double px,double py){return (xx-x)*(py-y)-(yy-y)*(px-x);};
            const double ac=cross(ax,ay,bx,by,cx,cy),ad=cross(ax,ay,bx,by,dx,dy);
            const double ca=cross(cx,cy,dx,dy,ax,ay),cb=cross(cx,cy,dx,dy,bx,by);
            if(std::max(std::min(ax,bx),std::min(cx,dx))<=std::min(std::max(ax,bx),std::max(cx,dx))&&
               std::max(std::min(ay,by),std::min(cy,dy))<=std::min(std::max(ay,by),std::max(cy,dy))&&
               ac*ad<=0&&ca*cb<=0){bboxFallback();return true;}
        }
    }
    // Nested islands inside an already covered hole need no second mesh.
    std::vector<std::vector<float>> outer;
    for(size_t i=0;i<domain.contours.size();++i) {
        bool nested=false;
        for(size_t j=0;j<domain.contours.size()&&!nested;++j)if(i!=j) {
            VertexLoop ring(domain.contours[j].size()/2);std::iota(ring.begin(),ring.end(),0);
            nested=pointInLoop(domain.contours[j],ring,domain.contours[i][0],domain.contours[i][1]);
        }
        if(!nested||params.generateHoles)outer.push_back(domain.contours[i]);
    }
    domain.contours=std::move(outer);
    // Serialized "density" values are point intervals, not target counts.
    // A larger interval must produce a sparser mesh at the same resolution.
    const double outsideSpacing=meshSampling::pointInterval(domain.width,domain.height,params.polygonOuterDensity);
    const auto originalContours=domain.contours;
    for(auto& poly:domain.contours) {
        ClosedCurve original(poly);
        VertexLoop ring(poly.size()/2);std::iota(ring.begin(),ring.end(),0);
        // Long narrow regions need longitudinal stations on the scale of
        // their width, rather than a few extremely long triangles. Limit the
        // correction to 2x along the outline; never use it to densify the whole
        // interior. Tiny raster islands keep their mandatory corners below.
        const double meanWidth=2*std::abs(loopArea(poly,ring))/std::max(1.,original.arc.back());
        const double spacing=std::min(outsideSpacing,std::max(outsideSpacing*.5,meanWidth));
        const int minimum=poly.size()/2<=8&&original.arc.back()<outsideSpacing*1.5?3:std::clamp(params.polygonMinBoundsPt,3,1000);
        const int target=std::clamp(std::max(minimum,int(std::lround(original.arc.back()/spacing))),3,1000);
        poly=smoothBoundary(poly,target,cell);
        if(poly.size()/2>8)continue; // fitted closed curve is already uniformly sampled
        const ClosedCurve curve(poly);const double spanSpacing=curve.arc.back()/target;
        std::vector<float> sampled;
        // Mandatory polygon corners divide the ring into spans. Divide each
        // span evenly instead of adding corners to a separate global lattice.
        for(size_t i=0;i<poly.size()/2;++i) {
            const size_t j=(i+1)%(poly.size()/2);
            const double length=curve.arc[i+1]-curve.arc[i];
            const int pieces=std::max(1,int(std::lround(length/spanSpacing)));
            for(int k=0;k<pieces;++k){const float f=float(k)/pieces;sampled.push_back(poly[2*i]+f*(poly[2*j]-poly[2*i]));sampled.push_back(poly[2*i+1]+f*(poly[2*j+1]-poly[2*i+1]));}
        }
        poly=std::move(sampled);
    }
    // An outward fit must not create overlapping components. Closely packed
    // islands keep their exact enclosing contours if smoothing would touch.
    auto overlap=[](const std::vector<float>& a,const std::vector<float>& b) {
        VertexLoop ringA(a.size()/2),ringB(b.size()/2);std::iota(ringA.begin(),ringA.end(),0);std::iota(ringB.begin(),ringB.end(),0);
        if(pointInLoop(a,ringA,b[0],b[1])||pointInLoop(b,ringB,a[0],a[1]))return true;
        for(size_t i=0;i<a.size()/2;++i)for(size_t j=0;j<b.size()/2;++j) {
            const size_t k=(i+1)%(a.size()/2),l=(j+1)%(b.size()/2);
            const CurvePoint p{a[2*i],a[2*i+1]},q{a[2*k],a[2*k+1]},r{b[2*j],b[2*j+1]},t{b[2*l],b[2*l+1]};
            if(std::max(std::min(p.x,q.x),std::min(r.x,t.x))<=std::min(std::max(p.x,q.x),std::max(r.x,t.x))&&
               std::max(std::min(p.y,q.y),std::min(r.y,t.y))<=std::min(std::max(p.y,q.y),std::max(r.y,t.y))&&
               curveCross(p,q,r)*curveCross(p,q,t)<=0&&curveCross(r,t,p)*curveCross(r,t,q)<=0)return true;
        }
        return false;
    };
    for(size_t pass=0;pass<domain.contours.size();++pass) {
        bool changed=false;
        for(size_t i=0;i<domain.contours.size();++i)for(size_t j=i+1;j<domain.contours.size();++j)
            if(overlap(domain.contours[i],domain.contours[j])&&
               (domain.contours[i]!=originalContours[i]||domain.contours[j]!=originalContours[j])) {
                domain.contours[i]=originalContours[i];domain.contours[j]=originalContours[j];changed=true;
            }
        if(!changed)break;
    }
    if(inset==0&&domain.contours.size()==1) {
        domain.stripCenterline=meshSampling::thinCenterline(domain.rgba,domain.width,domain.height);
        domain.strokeCenterline=meshSampling::strokeCenterline(domain.rgba,domain.width,domain.height);
        collapseSharpStripTips(domain,outsideSpacing,std::max(0,params.polygonMargin));
    }
    if(inset==0)for(auto& contour:domain.contours)materialOutline(contour,params.style,outsideSpacing);
    return true;
}

// Sample an enclosed transparent boundary from inside the empty region. Any
// chord leaving that region is refined back to the exact raster stations.
std::vector<float> sampleHole(const std::vector<float>& raw,double spacing,const std::vector<float>* allowed=nullptr) {
    const auto& region=allowed?*allowed:raw;
    const ClosedCurve curve(raw);VertexLoop ids(region.size()/2);std::iota(ids.begin(),ids.end(),0);
    const int count=std::clamp(int(std::ceil(curve.arc.back()/std::max(1.,spacing))),3,1000);
    std::vector<double> stations;for(int i=0;i<count;++i)stations.push_back(curve.arc.back()*i/count);
    // Keep actual polygon corners (rectangular cutouts, screw slots, etc.).
    if(raw.size()/2<=8)stations.insert(stations.end(),curve.arc.begin(),curve.arc.end()-1);
    for(int pass=0;pass<12;++pass){std::sort(stations.begin(),stations.end());
        stations.erase(std::unique(stations.begin(),stations.end(),[](double a,double b){return std::abs(a-b)<1e-5;}),stations.end());
        std::vector<double> extra;
        for(size_t i=0;i<stations.size();++i){const double start=stations[i],end=i+1<stations.size()?stations[i+1]:curve.arc.back();
            const auto a=curve.at(start),b=curve.at(end);bool valid=true;
            const int probes=std::max(1,int(std::ceil(std::hypot(b.x-a.x,b.y-a.y)*2)));
            for(int k=0;k<=probes&&valid;++k){const float x=float(a.x+(b.x-a.x)*k/probes),y=float(a.y+(b.y-a.y)*k/probes);
                valid=pointInLoop(region,ids,x,y)||nearBoundary(region,ids,x,y,1e-4f);}
            // Pixel-spaced probes can miss a tiny excursion through a concave
            // corner. Test every interval cut by the exact source boundary.
            std::vector<double> cuts{0,1};const double dx=b.x-a.x,dy=b.y-a.y;
            for(size_t r=0;r<ids.size()&&valid;++r){const size_t next=(r+1)%ids.size();const double cx=region[2*r],cy=region[2*r+1],ex=region[2*next]-cx,ey=region[2*next+1]-cy;
                const double det=dx*ey-dy*ex;if(std::abs(det)<1e-12)continue;
                const double t=((cx-a.x)*ey-(cy-a.y)*ex)/det,u=((cx-a.x)*dy-(cy-a.y)*dx)/det;
                if(t>0&&t<1&&u>=-1e-8&&u<=1+1e-8)cuts.push_back(t);}
            std::sort(cuts.begin(),cuts.end());for(size_t k=1;k<cuts.size()&&valid;++k){const double t=(cuts[k-1]+cuts[k])*.5;const float x=float(a.x+t*dx),y=float(a.y+t*dy);
                valid=pointInLoop(region,ids,x,y)||nearBoundary(region,ids,x,y,1e-5f);}
            if(!valid)for(double arc:curve.arc)if(arc>start+1e-5&&arc<end-1e-5)extra.push_back(arc);
            if(valid){const double squared=dx*dx+dy*dy,tolerance=std::clamp(spacing*.04,.75,2.);double deviation=0;
                for(size_t k=0;k+1<curve.arc.size();++k)if(curve.arc[k]>start&&curve.arc[k]<end){
                    const double x=raw[2*k]-a.x,y=raw[2*k+1]-a.y,t=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                    deviation=std::max(deviation,std::hypot(x-t*dx,y-t*dy));}
                if(deviation>tolerance)extra.push_back((start+end)*.5);
            }
        }
        if(extra.empty())break;stations.insert(stations.end(),extra.begin(),extra.end());
    }
    std::vector<float> sampled;for(double arc:stations){const auto q=curve.at(arc);sampled.push_back(float(q.x));sampled.push_back(float(q.y));}
    VertexLoop ring(sampled.size()/2);std::iota(ring.begin(),ring.end(),0);
    return simpleOutline(sampled,ring,1e-6)?sampled:raw;
}

// Only adjacent stations of this one hole can collapse. Work before support
// rings and RGB constraints are built so preview and generation stay identical.
int mergeHoleNeighbors(std::vector<float>& path,double spacing,const TextureDomain& domain) {
    const auto original=path;const size_t count=original.size()/2;
    VertexLoop live(count);std::iota(live.begin(),live.end(),0);
    const double area=std::abs(loopArea(original,live));
    // Raster stair steps are not deformation corners. Judge the turn across
    // a short arc, rather than freezing every one-pixel right angle. Retain
    // the original stations (no averaging), with cumulative fit/area limits.
    const double limit=spacing*.65,tolerance=std::min(.75,std::max(.35,spacing*.06));
    const ClosedCurve curve(original);
    const double cornerSpan=std::min(curve.arc.back()*.08,std::clamp(spacing*.25,2.,6.));
    std::set<int> corners;
    for(size_t i=0;i<count;++i){const size_t a=(i+count-1)%count,b=(i+1)%count;
        const double ux=original[2*i]-original[2*a],uy=original[2*i+1]-original[2*a+1],vx=original[2*b]-original[2*i],vy=original[2*b+1]-original[2*i+1];
        const auto left=curve.at(curve.arc[i]-cornerSpan),right=curve.at(curve.arc[i]+cornerSpan);
        const double sx=original[2*i]-left.x,sy=original[2*i+1]-left.y,tx=right.x-original[2*i],ty=right.y-original[2*i+1];
        if(ux*vx+uy*vy<=.5*std::hypot(ux,uy)*std::hypot(vx,vy)&&
           sx*tx+sy*ty<=.5*std::hypot(sx,sy)*std::hypot(tx,ty))corners.insert(int(i));}
    auto transparentEar=[&](int a,int b,int c){
        const CurvePoint tri[3]={{original[2*a],original[2*a+1]},{original[2*b],original[2*b+1]},{original[2*c],original[2*c+1]}};
        if(std::abs(curveCross(tri[0],tri[1],tri[2]))<1e-10)return true;
        const int x0=std::max(0,int(std::floor(std::min({tri[0].x,tri[1].x,tri[2].x})))),x1=std::min(domain.width,int(std::ceil(std::max({tri[0].x,tri[1].x,tri[2].x}))));
        const int y0=std::max(0,int(std::floor(std::min({tri[0].y,tri[1].y,tri[2].y})))),y1=std::min(domain.height,int(std::ceil(std::max({tri[0].y,tri[1].y,tri[2].y}))));
        for(int y=y0;y<y1;++y)for(int x=x0;x<x1;++x)if(domain.rgba[(size_t(y)*domain.width+x)*4+3]){
            std::vector<CurvePoint> polygon(std::begin(tri),std::end(tri));
            for(int side=0;side<4&&!polygon.empty();++side){const double bound=side==0?x:side==1?x+1:side==2?y:y+1;
                auto distance=[&](CurvePoint q){return (side%2? -1.:1.)*((side<2?q.x:q.y)-bound);};std::vector<CurvePoint> clipped;
                for(size_t i=0;i<polygon.size();++i){const auto p=polygon[i],q=polygon[(i+1)%polygon.size()];const double dp=distance(p),dq=distance(q);
                    if(dp>=0)clipped.push_back(p);if((dp>=0)!=(dq>=0)){const double t=dp/(dp-dq);clipped.push_back({p.x+(q.x-p.x)*t,p.y+(q.y-p.y)*t});}}
                polygon=std::move(clipped);}
            double overlap=0;for(size_t i=0;i<polygon.size();++i){const auto p=polygon[i],q=polygon[(i+1)%polygon.size()];overlap+=(p.x-x)*(q.y-y)-(q.x-x)*(p.y-y);}
            if(std::abs(overlap)>1e-10)return false;
        }return true;
    };
    while(live.size()>3){std::vector<std::pair<double,std::pair<int,int>>> candidates;
        for(size_t i=0;i<live.size();++i){int a=live[i],b=live[(i+1)%live.size()];const double distance=std::hypot(original[2*a]-original[2*b],original[2*a+1]-original[2*b+1]);
            if(distance<limit)candidates.push_back({distance,{a,b}});}
        std::sort(candidates.begin(),candidates.end());bool changed=false;
        for(const auto& candidate:candidates){for(int remove:{candidate.second.first,candidate.second.second}){
            if(corners.count(remove))continue;
            const size_t at=size_t(std::find(live.begin(),live.end(),remove)-live.begin());const int a=live[(at+live.size()-1)%live.size()],b=live[(at+1)%live.size()];
            const double dx=original[2*b]-original[2*a],dy=original[2*b+1]-original[2*a+1],squared=dx*dx+dy*dy;if(squared<1e-12)continue;
            if(squared>spacing*spacing*1.35*1.35)continue;
            bool valid=true;
            // Check every original station on this arc, including previously
            // removed points: repeated collapses cannot accumulate drift.
            for(size_t v=(size_t(a)+1)%count;v!=size_t(b);v=(v+1)%count){const double x=original[2*v]-original[2*a],y=original[2*v+1]-original[2*a+1],t=std::clamp((x*dx+y*dy)/squared,0.,1.);
                if(std::hypot(x-t*dx,y-t*dy)>tolerance){valid=false;break;}}
            if(!valid||!transparentEar(a,remove,b))continue;
            auto proposal=live;proposal.erase(proposal.begin()+at);
            const double newArea=std::abs(loopArea(original,proposal));
            if(newArea<area*.98||newArea>area*1.02||!simpleOutline(original,proposal,1e-6))continue;
            live=std::move(proposal);changed=true;break;
        }if(changed)break;}
        if(!changed)break;
    }
    path.clear();for(int v:live){path.push_back(original[2*v]);path.push_back(original[2*v+1]);}
    return int(count-live.size());
}

bool prepareTextureHoles(TextureDomain& domain,const MeshAutoGenParams& params,std::string& error) {
    MeshAutoGenParams rawParams=params;rawParams.polygonMargin=0;rawParams.polygonMinBoundsPt=3;
    TextureDomain raw;raw.rgba=domain.rgba;raw.width=domain.width;raw.height=domain.height;
    if(!buildTextureContours(raw,rawParams,error)||raw.conservativeFallback){error="镂空轮廓无法安全提取，拒绝填平镂空";return false;}
    domain.holes.clear();domain.holeSupports.clear();domain.mergedHoleVertices=0;
    TextureDomain retreat;
    if(params.holeBoundaryMargin>0){retreat.rgba=domain.rgba;retreat.width=domain.width;retreat.height=domain.height;
        auto insetParams=rawParams;insetParams.polygonMargin=params.holeBoundaryMargin;std::string why;
        if(!buildTextureContours(retreat,insetParams,why)||retreat.conservativeFallback)retreat.holes.clear();}
    const double spacing=meshSampling::pointInterval(domain.width,domain.height,params.holePointSpacing);
    for(const auto& hole:raw.holes){VertexLoop ids(hole.size()/2);std::iota(ids.begin(),ids.end(),0);
        if(std::abs(loopArea(hole,ids))<std::max(1,params.holeMinArea))continue;
        const std::vector<float>* candidate=&hole;int matches=0;
        for(const auto& inset:retreat.holes)if(pointInLoop(hole,ids,inset[0],inset[1])){candidate=&inset;++matches;}
        if(matches!=1)candidate=&hole; // preserve narrow necks rather than splitting a hole
        auto sampled=sampleHole(*candidate,spacing,&hole);
        // Atlas UVs are stored as floats. Keep a subpixel safety inset so
        // roundoff cannot move a hole edge across an opaque pixel corner.
        VertexLoop vertices(sampled.size()/2);std::iota(vertices.begin(),vertices.end(),0);
        const double sign=loopArea(sampled,vertices)>0?1.:-1.;const auto original=sampled;
        for(size_t i=0;i<vertices.size();++i){const size_t a=(i+vertices.size()-1)%vertices.size(),b=(i+1)%vertices.size();
            const double dx0=original[2*i]-original[2*a],dy0=original[2*i+1]-original[2*a+1],dx1=original[2*b]-original[2*i],dy1=original[2*b+1]-original[2*i+1];
            const double l0=std::max(1e-8,std::hypot(dx0,dy0)),l1=std::max(1e-8,std::hypot(dx1,dy1));
            const double nx0=-dy0/l0*sign,ny0=dx0/l0*sign,nx1=-dy1/l1*sign,ny1=dx1/l1*sign;
            const double scale=.005/std::max(.1,1+nx0*nx1+ny0*ny1);
            sampled[2*i]+=float((nx0+nx1)*scale);sampled[2*i+1]+=float((ny0+ny1)*scale);}
        if(!simpleOutline(sampled,vertices,1e-6)){error="镂空过窄，无法保留安全的浮点边界";return false;}
        domain.mergedHoleVertices+=mergeHoleNeighbors(sampled,spacing,domain);
        domain.holes.push_back(std::move(sampled));}
    // Start from top-level opaque components. Reintroduce only islands inside
    // retained holes below; an area-filtered hole is filled, so its islands
    // must not create overlapping duplicate meshes.
    std::vector<std::vector<float>> roots;
    for(size_t i=0;i<domain.contours.size();++i){bool nested=false;
        for(size_t j=0;j<domain.contours.size()&&!nested;++j)if(i!=j){const auto& other=domain.contours[j];VertexLoop ids(other.size()/2);std::iota(ids.begin(),ids.end(),0);
            nested=pointInLoop(other,ids,domain.contours[i][0],domain.contours[i][1]);}
        if(!nested)roots.push_back(domain.contours[i]);}
    domain.contours=std::move(roots);
    if(domain.holes.empty())return true;
    // Outer dilation can absorb a tiny opaque island into the surrounding
    // component. Restore islands from the zero-margin alpha topology before
    // subtracting holes, otherwise those pixels would be cut out with the void.
    auto insideHole=[&](const std::vector<float>& contour){for(const auto& hole:domain.holes){VertexLoop ids(hole.size()/2);std::iota(ids.begin(),ids.end(),0);
        if(pointInLoop(hole,ids,contour[0],contour[1]))return true;}return false;};
    domain.contours.erase(std::remove_if(domain.contours.begin(),domain.contours.end(),insideHole),domain.contours.end());
    for(const auto& contour:raw.contours)if(insideHole(contour))domain.contours.push_back(contour);
    domain.stripCenterline.clear();domain.strokeCenterline.clear();
    if(domain.holes.size()>128){error="镂空数量超过 128，请提高最小镂空面积";return false;}
    if(!params.generateHoleBorder)return true;
    const double ringSpacing=meshSampling::pointInterval(domain.width,domain.height,params.holeBorderSpacing>0?params.holeBorderSpacing:params.holePointSpacing);
    for(size_t index=0;index<domain.holes.size();++index){bool found=false;
        for(double fraction:{1.,.5,.25}){
            TextureDomain inset;inset.rgba=domain.rgba;inset.width=domain.width;inset.height=domain.height;
            std::string why;if(!buildTextureContours(inset,rawParams,why,std::max(1.,params.holeBorderDistance*fraction))||inset.conservativeFallback)continue;
            for(const auto& expanded:inset.holes){VertexLoop rawIds(expanded.size()/2);std::iota(rawIds.begin(),rawIds.end(),0);
                const auto& hole=domain.holes[index];if(!pointInLoop(expanded,rawIds,hole[0],hole[1]))continue;
                int contained=0;for(const auto& other:domain.holes)contained+=pointInLoop(expanded,rawIds,other[0],other[1]);
                if(contained!=1)continue;
                for(double divisor:{-1.,1.,2.,4.,8.,0.}) {
                std::vector<float> path;
                if(divisor<0){std::vector<float> positive;for(size_t i=expanded.size();i>0;i-=2){positive.push_back(expanded[i-2]);positive.push_back(expanded[i-1]);}
                    path=smoothBoundary(positive,std::clamp(int(std::ceil(ClosedCurve(positive).arc.back()/ringSpacing)),3,1000),1);}
                else path=divisor?sampleHole(expanded,ringSpacing/divisor):expanded;
                VertexLoop ids(path.size()/2);std::iota(ids.begin(),ids.end(),0);bool valid=true;
                for(size_t i=0;i<hole.size()&&valid;i+=2)valid=pointInLoop(path,ids,hole[i],hole[i+1]);
                for(size_t i=0;i<ids.size()&&valid;++i){const size_t j=(i+1)%ids.size();
                    const double dx=path[2*j]-path[2*i],dy=path[2*j+1]-path[2*i+1];const int probes=std::max(1,int(std::ceil(std::hypot(dx,dy)*2)));
                    for(int k=0;k<=probes;++k){const int x=int(std::floor(path[2*i]+dx*k/probes)),y=int(std::floor(path[2*i+1]+dy*k/probes));
                        if(x<0||y<0||x>=domain.width||y>=domain.height||!domain.rgba[(size_t(y)*domain.width+x)*4+3]){valid=false;break;}}
                }
                if(valid){domain.holeSupports.push_back(std::move(path));found=true;break;}
                }
                if(found)break;
            }
            if(found)break;
        }
        if(!found){error="镂空周围空间不足以生成闭合锁边，请减小镂空锁边点距或关闭镂空锁边";return false;}
    }
    return true;
}

const XsNode* deep(const XsNode& n, const std::string& field) {
    if (n.fieldName == field) return &n;
    for (const auto& c : n.children)
        if (const XsNode* f = deep(c, field)) return f;
    return nullptr;
}
XsNode* deepMut(XsNode& n, const std::string& field) {
    if (n.fieldName == field) return &n;
    for (auto& c : n.children)
        if (XsNode* f = deepMut(c, field)) return f;
    return nullptr;
}
const XsNode* findType(const XsNode& n, const std::string& type) {
    if (n.type == type) return &n;
    for (const auto& c : n.children)
        if (const XsNode* f = findType(c, type)) return f;
    return nullptr;
}
XsNode* findTypeMut(XsNode& n, const std::string& type) {
    if (n.type == type) return &n;
    for (auto& c : n.children)
        if (XsNode* f = findTypeMut(c, type)) return f;
    return nullptr;
}
const XsNode* findId(const XsNode& n, int id) {
    if (n.id == id) return &n;
    for (const auto& c : n.children)
        if (const XsNode* f = findId(c, id)) return f;
    return nullptr;
}
std::string arrayText(const XsNode& root, const XsNode* n) {
    if (!n) return {};
    if (!n->text.empty()) return n->text;
    if (n->ref >= 0)
        if (const XsNode* t = findId(root, n->ref)) return t->text;
    return {};
}
std::vector<float> floats(const std::string& s) {
    std::vector<float> out;
    const char* p = s.c_str();
    while (*p) {
        char* e = nullptr;
        const float v = std::strtof(p, &e);
        if (e == p) break;
        out.push_back(v);
        p = e;
    }
    return out;
}
std::vector<int> ints(const std::string& s) {
    std::vector<int> out;
    const char* p = s.c_str();
    while (*p) {
        char* e = nullptr;
        const long v = std::strtol(p, &e, 10);
        if (e == p) break;
        out.push_back(static_cast<int>(v));
        p = e;
    }
    return out;
}
bool findPath(const XsNode& n, const XsNode* target, std::vector<int>& path) {
    if (&n == target) return true;
    for (size_t i = 0; i < n.children.size(); ++i) {
        path.push_back(static_cast<int>(i));
        if (findPath(n.children[i], target, path)) return true;
        path.pop_back();
    }
    return false;
}
XsNode* atPath(XsNode& root, const std::vector<int>& path) {
    XsNode* n = &root;
    for (int i : path) {
        if (i < 0 || i >= static_cast<int>(n->children.size())) return nullptr;
        n = &n->children[static_cast<size_t>(i)];
    }
    return n;
}
void setArray(XsNode& n, const std::string& text) {
    // A field may reference an array owned by another form or editable mesh.
    // Detach that field instead of leaving stale xs.ref alongside new contents.
    if(n.ref>=0) {
        n.ref=-1;
        n.attrs.erase(std::remove_if(n.attrs.begin(),n.attrs.end(),[](const auto& attr){return attr.first=="xs.ref";}),n.attrs.end());
    }
    n.selfClosing = false;
    n.text = text;
    const int i = n.attrIndex("count");
    if (i >= 0) n.attrs[static_cast<size_t>(i)].second = std::to_string(
        text.empty() ? 0 : static_cast<int>(std::count(text.begin(), text.end(), ' ') + 1));
}
std::string joinFloats(const std::vector<float>& v) {
    std::string s;
    char b[48];
    for (float f : v) {
        if (!s.empty()) s.push_back(' ');
        std::snprintf(b, sizeof(b), "%.9g", static_cast<double>(f));
        s += b;
    }
    return s;
}
std::string joinInts(const std::vector<int>& v) {
    std::string s;
    for (int x : v) {
        if (!s.empty()) s.push_back(' ');
        s += std::to_string(x);
    }
    return s;
}
bool barycentric(float px, float py, float ax, float ay, float bx, float by,
                 float cx, float cy, float& u, float& v, float& w, bool signedWeights=false) {
    const float d = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
    if (std::fabs(d) < 1e-9f) return false;
    u = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) / d;
    v = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) / d;
    w = 1.0f - u - v;
    return signedWeights || (u >= -1e-4f && v >= -1e-4f && w >= -1e-4f);
}

// The alpha contour can extend past the old triangulated polygon while staying
// in this layer's UV footprint. Use the nearest nondegenerate triangle's affine
// field there, preserving the same signed weights in every keyform.
int sampleTriangle(const std::vector<float>& points,const std::vector<int>& tris,
                   float x,float y,float& a,float& b,float& c,bool* extrapolated=nullptr) {
    double closest=std::numeric_limits<double>::max();int best=-1;
    for(size_t t=0;t+2<tris.size();t+=3) {
        const int i=tris[t],j=tris[t+1],k=tris[t+2];float u,v,w;
        if(!barycentric(x,y,points[2*i],points[2*i+1],points[2*j],points[2*j+1],points[2*k],points[2*k+1],u,v,w,true))continue;
        if(u>=-1e-4f&&v>=-1e-4f&&w>=-1e-4f) {
            a=u;b=v;c=w;if(extrapolated)*extrapolated=false;return int(t/3);
        }
        for(int edge=0;edge<3;++edge) {
            const int p=tris[t+edge],q=tris[t+(edge+1)%3];
            const double dx=points[2*q]-points[2*p],dy=points[2*q+1]-points[2*p+1];
            const double len=dx*dx+dy*dy;
            const double f=len?std::clamp(((x-points[2*p])*dx+(y-points[2*p+1])*dy)/len,0.,1.):0.;
            const double distance=std::pow(x-points[2*p]-f*dx,2)+std::pow(y-points[2*p+1]-f*dy,2);
            if(distance<closest){closest=distance;best=int(t/3);a=u;b=v;c=w;}
        }
    }
    if(extrapolated)*extrapolated=best>=0;
    return best;
}
bool nearBoundary(const std::vector<float>& pts, const VertexLoop& loop, float x, float y,
                  float margin) {
    const float m2 = margin * margin;
    for (size_t i = 0; i < loop.size(); ++i) {
        const int a = loop[i], b = loop[(i + 1) % loop.size()];
        const float ax = pts[static_cast<size_t>(2 * a)], ay = pts[static_cast<size_t>(2 * a + 1)];
        const float bx = pts[static_cast<size_t>(2 * b)], by = pts[static_cast<size_t>(2 * b + 1)];
        const float dx = bx - ax, dy = by - ay;
        const float len2 = dx * dx + dy * dy;
        float t = len2 > 1e-9f ? ((x - ax) * dx + (y - ay) * dy) / len2 : 0.0f;
        t = std::clamp(t, 0.0f, 1.0f);
        const float ex = ax + t * dx - x, ey = ay + t * dy - y;
        if (ex * ex + ey * ey < m2) return true;
    }
    return false;
}

bool simpleOutline(const std::vector<float>& points, const VertexLoop& loop, double tolerance) {
    if (std::set<int>(loop.begin(), loop.end()).size() != loop.size()) return false;
    auto cross = [&](int a, int b, int c) {
        return (double(points[2 * b]) - points[2 * a]) * (double(points[2 * c + 1]) - points[2 * a + 1]) -
               (double(points[2 * b + 1]) - points[2 * a + 1]) * (double(points[2 * c]) - points[2 * a]);
    };
    auto onSegment = [&](int a, int b, int c) {
        const double length = std::hypot(double(points[2 * b]) - points[2 * a], double(points[2 * b + 1]) - points[2 * a + 1]);
        return std::abs(cross(a, b, c)) <= tolerance * length &&
            points[2 * c] >= std::min(points[2 * a], points[2 * b]) - tolerance &&
            points[2 * c] <= std::max(points[2 * a], points[2 * b]) + tolerance &&
            points[2 * c + 1] >= std::min(points[2 * a + 1], points[2 * b + 1]) - tolerance &&
            points[2 * c + 1] <= std::max(points[2 * a + 1], points[2 * b + 1]) + tolerance;
    };
    for (size_t i = 0; i < loop.size(); ++i) {
        const int a = loop[i], b = loop[(i + 1) % loop.size()];
        for (size_t j = i + 2; j < loop.size(); ++j) {
            if (i == 0 && j + 1 == loop.size()) continue;
            const int c = loop[j], d = loop[(j + 1) % loop.size()];
            if (onSegment(a, b, c) || onSegment(a, b, d) || onSegment(c, d, a) || onSegment(c, d, b)) return false;
            const double abC = cross(a, b, c), abD = cross(a, b, d);
            const double cdA = cross(c, d, a), cdB = cross(c, d, b);
            if ((abC < 0) != (abD < 0) && (cdA < 0) != (cdB < 0)) return false;
        }
    }
    return true;
}

// Sum triangle/pixel-square intersection areas. Testing just pixel centers can
// miss a half-pixel cut at a silhouette edge. All nonzero alpha is mandatory.
bool verifyTextureCoverage(const TextureDomain& domain,const std::vector<float>& pts,
                           const std::vector<int>& tris,size_t& opaque,size_t& covered,std::string& error,float areaTolerance=1e-4f) {
    std::vector<float> coverage(size_t(domain.width)*domain.height);
    struct Point{double x,y;};
    for(size_t t=0;t<tris.size();t+=3) {
        Point triangle[3];double minX=domain.width,minY=domain.height,maxX=0,maxY=0;
        for(int k=0;k<3;++k) {
            triangle[k]={pts[2*tris[t+k]],pts[2*tris[t+k]+1]};
            minX=std::min(minX,triangle[k].x);maxX=std::max(maxX,triangle[k].x);
            minY=std::min(minY,triangle[k].y);maxY=std::max(maxY,triangle[k].y);
        }
        const int x0=std::clamp(int(std::floor(minX)),0,domain.width),x1=std::clamp(int(std::ceil(maxX)),0,domain.width);
        const int y0=std::clamp(int(std::floor(minY)),0,domain.height),y1=std::clamp(int(std::ceil(maxY)),0,domain.height);
        for(int y=y0;y<y1;++y)for(int x=x0;x<x1;++x) {
            const size_t pixel=size_t(y)*domain.width+x;
            if(!domain.rgba[pixel*4+3]||coverage[pixel]>=.99999f)continue;
            Point polygon[12],scratch[12];int n=3;std::copy_n(triangle,3,polygon);
            for(int side=0;side<4&&n;++side) {
                const bool horizontal=side>=2;const double bound=side==0?x:side==1?x+1:side==2?y:y+1;
                auto distance=[&](Point p){return (side%2==0?1.:-1.)*((horizontal?p.y:p.x)-bound);};
                int count=0;
                for(int i=0;i<n;++i) {
                    const Point a=polygon[i],b=polygon[(i+1)%n];const double da=distance(a),db=distance(b);
                    if(da>=0)scratch[count++]=a;
                    if((da>=0)!=(db>=0)) {const double f=da/(da-db);scratch[count++]={a.x+f*(b.x-a.x),a.y+f*(b.y-a.y)};}
                }
                n=count;std::copy_n(scratch,n,polygon);
            }
            double area=0;
            for(int i=0;i<n;++i){const auto a=polygon[i],b=polygon[(i+1)%n];area+=(a.x-x)*(b.y-y)-(b.x-x)*(a.y-y);}
            coverage[pixel]+=float(std::abs(area)*.5);
        }
    }
    opaque=covered=0;int firstX=-1,firstY=-1;
    for(int y=0;y<domain.height;++y)for(int x=0;x<domain.width;++x) {
        const size_t p=size_t(y)*domain.width+x;
        if(!domain.rgba[p*4+3])continue;
        ++opaque;
        if(coverage[p]>=1.f-areaTolerance)++covered;
        else if(firstX<0){firstX=x;firstY=y;}
    }
    if(opaque!=covered) {
        error="三角网格漏覆盖 "+std::to_string(opaque-covered)+" 个非透明像素，首处 ("+
            std::to_string(firstX)+", "+std::to_string(firstY)+")，拒绝写回";return false;
    }
    return true;
}

}  // namespace

bool MeshAutoGen::textureFeatures(const Cmo3Document& doc, int meshIdx,
    const MeshFeatureSettings& settings, MeshFeaturePreview& preview, std::string& error,
    const MeshAutoGenParams* params) {
    preview = {}; error.clear();
    if (meshIdx < 0 || meshIdx >= int(doc.artMeshes().size())) {error="网格不存在";return false;}
    const auto& mesh = doc.artMeshes()[size_t(meshIdx)];
    TextureDomain domain;
    if(!textureOutline(doc,mesh,params?*params:readParams(doc,meshIdx),domain,error))return false;
    preview.conservativeFallback=domain.conservativeFallback;
    for(const auto& contour:domain.contours) {
        auto uv=contour;for(size_t i=0;i<uv.size();i+=2){const auto p=domain.uv(uv[i],uv[i+1]);uv[i]=p.first;uv[i+1]=p.second;}
        preview.contours.push_back(std::move(uv));
    }
    for(const auto& hole:domain.holes){std::vector<float> uv;for(size_t i=0;i<hole.size();i+=2){const auto q=domain.uv(hole[i],hole[i+1]);uv.push_back(q.first);uv.push_back(q.second);}preview.holes.push_back(std::move(uv));}
    for(const auto& ring:domain.holeSupports){std::vector<float> uv;for(size_t i=0;i<ring.size();i+=2){const auto q=domain.uv(ring[i],ring[i+1]);uv.push_back(q.first);uv.push_back(q.second);}preview.holeBorders.push_back(std::move(uv));}
    // Include the complete mesh outline, including transparent space outside
    // the source PNG. The preview uses the same affine UV mapping as generation.
    float minX=float(domain.width),minY=float(domain.height),maxX=0,maxY=0;
    for(const auto& contour:domain.contours)for(size_t i=0;i<contour.size();i+=2) {
        minX=std::min(minX,contour[i]);maxX=std::max(maxX,contour[i]);minY=std::min(minY,contour[i+1]);maxY=std::max(maxY,contour[i+1]);
    }
    const int x0=int(std::floor(minX))-3,y0=int(std::floor(minY))-3;
    const int cw=int(std::ceil(maxX))+3-x0,ch=int(std::ceil(maxY))+3-y0;
    const float scale=std::min(1.f,1024.f/std::max(cw,ch));
    preview.width=std::max(7,int(cw*scale));preview.height=std::max(7,int(ch*scale));
    const int w=preview.width,h=preview.height;
    preview.u0=x0/domain.logicalW;preview.v0=y0/domain.logicalH;
    preview.uvToPixel={float(domain.fromUv[0]*w/cw),float(domain.fromUv[1]*w/cw),float((domain.fromUv[2]-x0)*w/cw),
        float(domain.fromUv[3]*h/ch),float(domain.fromUv[4]*h/ch),float((domain.fromUv[5]-y0)*h/ch)};
    preview.du=float(cw)/w/domain.logicalW;preview.dv=float(ch)/h/domain.logicalH;
    preview.premultiplied=domain.premultiplied;preview.rgba.resize(size_t(w)*h*4);
    std::vector<uint8_t> support(size_t(w)*h);
    for(int y=0;y<h;++y)for(int x=0;x<w;++x) {
        const int xx=int(std::floor(x0+(x+.5f)*cw/w)),yy=int(std::floor(y0+(y+.5f)*ch/h));
        if(xx<0||xx>=domain.width||yy<0||yy>=domain.height)continue; // transparent, never extrude edge texels
        const auto* pixel=domain.rgba.data()+(size_t(yy)*domain.width+xx)*4;
        std::copy_n(pixel,4,preview.rgba.data()+(size_t(y)*w+x)*4);support[size_t(y)*w+x]=pixel[3]>0;
    }
    auto detectorSettings=settings;
    const auto sampling=params?*params:readParams(doc,meshIdx);
    const double featureSpacing=meshSampling::pointInterval(domain.width,domain.height,sampling.polygonFeatureDensity);
    // At coarse mesh settings, short shading fragments must not force a dense
    // cluster of constraints. This only changes automatic candidate detection;
    // generate() still retains every enabled line supplied by the caller.
    detectorSettings.minLength=std::max(2,int(std::max(double(settings.minLength),featureSpacing*.5)*scale));
    detectorSettings.cornerTolerance=std::max(1,int(std::max(double(settings.cornerTolerance),featureSpacing*.05)*scale));
    detectorSettings.maxFoldWidth=std::max(0,int(settings.maxFoldWidth*scale));
    detectorSettings.nearbyMergeDistance=std::max(0,int(settings.nearbyMergeDistance*scale));
    auto featureStyle=sampling.style;
    if(sampling.material==MeshAutoGenMaterial::ByName)featureStyle=modeFromName(mesh.displayName()).style;
    if(sampling.material==MeshAutoGenMaterial::Soft)featureStyle=MeshAutoGenStyle::Soft;
    if(sampling.material==MeshAutoGenMaterial::Rigid)featureStyle=MeshAutoGenStyle::Rigid;
    const float featureScale=std::min(float(w)/cw,float(h)/ch);
    detectorSettings.curveSmoothing=(featureStyle==MeshAutoGenStyle::Rigid||featureStyle==MeshAutoGenStyle::FineLine)?0.f:settings.curveSmoothing*featureScale;
    detectorSettings.curveSpacing=float(featureSpacing)*featureScale;
    detectorSettings.curveBoundaryMargin=sampling.generateBorder?
        (std::max(0,sampling.polygonInnerMargin)+settings.curveSmoothing+2)*featureScale:0.f;
    detectorSettings.premultiplied=preview.premultiplied;
    auto detected=detectMeshTextureFeatures(preview.rgba.data(),w,h,detectorSettings,support.data());
    if(!detected.error.empty()){error=detected.error;preview={};return false;}
    preview.limited=detected.limited;
    if(sampling.generateBorder) {
        std::vector<double> lengths;
        for(const auto& contour:domain.contours)for(size_t i=0;i<contour.size()/2;++i) {
            const size_t j=(i+1)%(contour.size()/2);lengths.push_back(std::hypot(contour[2*j]-contour[2*i],contour[2*j+1]-contour[2*i+1]));
        }
        std::sort(lengths.begin(),lengths.end());
        const double cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
        const double outsideInterval=meshSampling::pointInterval(domain.width,domain.height,sampling.polygonOuterDensity);
        const double borderInterval=meshSampling::pointInterval(domain.width,domain.height,
            sampling.polygonBorderDensity>0?sampling.polygonBorderDensity:sampling.polygonOuterDensity);
        const double spacing=std::max(4.,std::min(borderInterval,
            std::max(lengths[lengths.size()/2]*borderInterval/outsideInterval,borderInterval*.5)));
        const double depth=std::max(1.,std::min(sampling.polygonInnerMargin>0?double(sampling.polygonInnerMargin):spacing*.55,spacing*.75));
        const double band=std::max(0,sampling.polygonMargin)+depth+3*cell;
        const double tolerance=std::clamp(featureSpacing*.05,cell,4*cell);
        // The inner ring already follows these silhouette-parallel color
        // transitions. Do not offer a duplicate row of automatic constraints
        // beside it. Explicit enabled paths passed to generate() are unchanged.
        detected.lines.erase(std::remove_if(detected.lines.begin(),detected.lines.end(),[&](const MeshFeatureLine& line) {
            int parallel=0,probes=0;
            for(size_t k=2;k<line.points.size();k+=2) {
                const double ax=x0+line.points[k-2]*cw/w,ay=y0+line.points[k-1]*ch/h;
                const double bx=x0+line.points[k]*cw/w,by=y0+line.points[k+1]*ch/h;
                const double dx=bx-ax,dy=by-ay,length=std::hypot(dx,dy);
                if(length<1e-5)continue;
                const int steps=std::max(1,int(std::ceil(length/std::max(cell*2,featureSpacing/16))));
                for(int s=0;s<steps;++s) {
                    const double x=ax+dx*(s+.5)/steps,y=ay+dy*(s+.5)/steps;
                    double nearest=1e30,cosine=0;
                    for(const auto& contour:domain.contours)for(size_t i=0;i<contour.size()/2;++i) {
                        const size_t j=(i+1)%(contour.size()/2);
                        const double vx=contour[2*j]-contour[2*i],vy=contour[2*j+1]-contour[2*i+1],squared=vx*vx+vy*vy;
                        const double t=squared?std::clamp(((x-contour[2*i])*vx+(y-contour[2*i+1])*vy)/squared,0.,1.):0;
                        const double distance=std::hypot(x-contour[2*i]-t*vx,y-contour[2*i+1]-t*vy);
                        if(distance<nearest){nearest=distance;cosine=squared?std::abs((dx*vx+dy*vy)/(length*std::sqrt(squared))):0;}
                    }
                    ++probes;if(std::abs(nearest-band)<=tolerance&&cosine>.94)++parallel;
                }
            }
            return probes>=4&&parallel>=probes*.85;
        }),detected.lines.end());
    }
    for(auto& line:detected.lines) {
        for(size_t i=0;i+1<line.points.size();i+=2) {
            const auto point=domain.uv(x0+line.points[i]*cw/w,y0+line.points[i+1]*ch/h);
            line.points[i]=point.first;line.points[i+1]=point.second;
        }
        line.length/=scale;
    }
    if(!domain.holes.empty()){
        // RGB tracing may connect a crease across a transparent opening. Split
        // automatic candidates at the actual mesh boundaries before presenting
        // them; manually supplied paths still fail validation if they cross.
        auto containingOuter=[&](float x,float y){int owner=-1;double smallest=1e30;
            for(size_t i=0;i<domain.contours.size();++i){const auto& path=domain.contours[i];VertexLoop ids(path.size()/2);std::iota(ids.begin(),ids.end(),0);
                if(pointInLoop(path,ids,x,y)){double area=std::abs(loopArea(path,ids));if(area<smallest){smallest=area;owner=int(i);}}}return owner;};
        std::vector<int> owners;for(const auto& hole:domain.holes)owners.push_back(containingOuter(hole[0],hole[1]));
        auto allowed=[&](double x,double y){const int owner=containingOuter(float(x),float(y));if(owner<0)return false;
            for(size_t i=0;i<domain.holes.size();++i)if(owners[i]==owner){const auto& path=domain.holes[i];VertexLoop ids(path.size()/2);std::iota(ids.begin(),ids.end(),0);
                if(pointInLoop(path,ids,float(x),float(y)))return false;}return true;};
        auto boundaries=domain.contours;boundaries.insert(boundaries.end(),domain.holes.begin(),domain.holes.end());
        std::vector<MeshFeatureLine> clipped;
        for(const auto& line:detected.lines){MeshFeatureLine piece=line;piece.points.clear();piece.length=0;
            auto flush=[&]{if(piece.points.size()>=4&&piece.length>=1)clipped.push_back(piece);piece.points.clear();piece.length=0;};
            for(size_t k=2;k<line.points.size();k+=2){const auto a=domain.pixel(line.points[k-2],line.points[k-1]),b=domain.pixel(line.points[k],line.points[k+1]);
                const double dx=b.first-a.first,dy=b.second-a.second,length=std::hypot(dx,dy);if(length<1e-8)continue;
                std::vector<double> cuts{0,1};
                for(const auto& path:boundaries)for(size_t i=0;i<path.size()/2;++i){const size_t j=(i+1)%(path.size()/2);
                    const double ex=path[2*j]-path[2*i],ey=path[2*j+1]-path[2*i+1],det=dx*ey-dy*ex;if(std::abs(det)<1e-12)continue;
                    const double cx=path[2*i]-a.first,cy=path[2*i+1]-a.second,t=(cx*ey-cy*ex)/det,u=(cx*dy-cy*dx)/det;
                    if(t>0&&t<1&&u>=0&&u<=1)cuts.push_back(t);}
                std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-8;}),cuts.end());
                for(size_t i=1;i<cuts.size();++i){const double mid=(cuts[i-1]+cuts[i])*.5;
                    if(!allowed(a.first+dx*mid,a.second+dy*mid)){flush();continue;}
                    const double epsilon=std::min(.01/length,(cuts[i]-cuts[i-1])*.1);
                    const double t0=cuts[i-1]+(cuts[i-1]>0?epsilon:0),t1=cuts[i]-(cuts[i]<1?epsilon:0);
                    const auto start=domain.uv(float(a.first+dx*t0),float(a.second+dy*t0)),end=domain.uv(float(a.first+dx*t1),float(a.second+dy*t1));
                    if(piece.points.empty()){piece.points.push_back(start.first);piece.points.push_back(start.second);}
                    piece.points.push_back(end.first);piece.points.push_back(end.second);piece.length+=float(length*(t1-t0));
                    if(cuts[i]<1)flush();}
            }flush();
        }
        detected.lines=std::move(clipped);
        if(settings.weakBoundaryDuplicates){
            // A painted edge following an already-supported alpha outline is
            // not a second crease. Require proximity AND parallel direction
            // along each retained segment. Split mixed candidates at existing
            // bends so an interior fold is not weakened with its rim-following
            // tail. Work after clipping, in source pixels, not UV units.
            struct Guide {const std::vector<float>* path;double band;};std::vector<Guide> guides;
            if(sampling.generateHoleBorder)for(const auto& hole:domain.holes)
                guides.push_back({&hole,double(std::max(0,sampling.holeBoundaryMargin)+std::max(1,sampling.holeBorderDistance)+3)});
            if(sampling.generateBorder)for(const auto& outer:domain.contours)
                guides.push_back({&outer,double(std::max(0,sampling.polygonMargin)+std::max(1,sampling.polygonInnerMargin)+5)});
            std::vector<MeshFeatureLine> classified;
            for(const auto& line:detected.lines){
                if(!line.automatic||line.weakReference||guides.empty()){classified.push_back(line);continue;}
                MeshFeatureLine piece=line;piece.points.clear();piece.length=0;
                auto flush=[&]{if(piece.points.size()>=4)classified.push_back(piece);piece.points.clear();piece.length=0;};
                for(size_t k=2;k<line.points.size();k+=2){const auto a=domain.pixel(line.points[k-2],line.points[k-1]),b=domain.pixel(line.points[k],line.points[k+1]);
                    const double dx=b.first-a.first,dy=b.second-a.second,length=std::hypot(dx,dy);if(length<1e-8)continue;
                    double total=0,matched=0,longestGap=0,gap=0;
                    const int steps=std::max(2,int(std::ceil(length/2)));const double weight=length/steps;
                    for(int s=0;s<steps;++s){const double x=a.first+dx*(s+.5)/steps,y=a.second+dy*(s+.5)/steps;bool follows=false;
                        for(const auto& guide:guides){const auto& path=*guide.path;double nearest=1e30,cosine=0;
                            for(size_t i=0;i<path.size()/2;++i){const size_t j=(i+1)%(path.size()/2);
                                const double vx=path[2*j]-path[2*i],vy=path[2*j+1]-path[2*i+1],squared=vx*vx+vy*vy;
                                const double t=squared?std::clamp(((x-path[2*i])*vx+(y-path[2*i+1])*vy)/squared,0.,1.):0;
                                const double distance=std::hypot(x-path[2*i]-t*vx,y-path[2*i+1]-t*vy);
                                if(distance<nearest){nearest=distance;cosine=squared?std::abs(dx*vx+dy*vy)/(length*std::sqrt(squared)):0;}
                            }
                            if(nearest<=guide.band&&cosine>=.9){follows=true;break;}
                        }
                        total+=weight;if(follows){matched+=weight;gap=0;}else{gap+=weight;longestGap=std::max(longestGap,gap);}
                    }
                    const bool weak=matched>=total*.9&&longestGap<=std::max(2.,total*.1);
                    if(!piece.points.empty()&&piece.boundaryReference!=weak)flush();
                    piece.boundaryReference=weak;piece.weakReference=weak;
                    if(piece.points.empty()){piece.points.push_back(line.points[k-2]);piece.points.push_back(line.points[k-1]);}
                    piece.points.push_back(line.points[k]);piece.points.push_back(line.points[k+1]);piece.length+=float(length);
                }
                flush();
            }
            detected.lines=std::move(classified);
        }
    }
    preview.mergedNearbyLines=detected.mergedNearbyLines;
    preview.lines=std::move(detected.lines);return true;
}

MeshAutoGenParams MeshAutoGen::readParams(const Cmo3Document& doc, int meshIdx) {
    MeshAutoGenParams p;
    if (meshIdx < 0 || meshIdx >= static_cast<int>(doc.artMeshes().size())) return p;
    const ArtMesh& m = doc.artMeshes()[static_cast<size_t>(meshIdx)];
    if (!m.source) return p;
    const XsNode* ext = findType(*m.source, "CMeshGeneratorExtension");
    const XsNode* s = ext ? deep(*ext, "meshGenerateSetting") : nullptr;
    if (!s) return p;
    auto get = [&](const char* name, int fallback) {
        const XsNode* n = deep(*s, name);
        return (n && !n->text.empty()) ? std::atoi(n->text.c_str()) : fallback;
    };
    p.polygonOuterDensity = get("polygonOuterDensity", p.polygonOuterDensity);
    p.polygonInnerDensity = get("polygonInnerDensity", p.polygonInnerDensity);
    p.polygonMargin = get("polygonMargin", p.polygonMargin);
    p.polygonInnerMargin = get("polygonInnerMargin", p.polygonInnerMargin);
    p.polygonMinMargin = get("polygonMinMargin", p.polygonMinMargin);
    p.polygonMinBoundsPt = get("polygonMinBoundsPt", p.polygonMinBoundsPt);
    p.thresholdAlpha = get("thresholdAlpha", p.thresholdAlpha);
    return p;
}

bool MeshAutoGen::writeParams(XsDocument& work, const Cmo3Document& doc, int meshIdx,
                              const MeshAutoGenParams& p, std::string& error) {
    error.clear();
    if (meshIdx < 0 || meshIdx >= static_cast<int>(doc.artMeshes().size())) {
        error = "mesh 下标越界";
        return false;
    }
    const int sourceId = doc.artMeshes()[static_cast<size_t>(meshIdx)].sourceXsId;
    XsNode* source = nullptr;
    for (auto& c : work.root.children) (void)c;
    {
        std::vector<XsNode*> stack{&work.root};
        while (!stack.empty()) {
            XsNode* cur = stack.back(); stack.pop_back();
            if (cur->id == sourceId && cur->type == "CArtMeshSource") { source = cur; break; }
            for (auto& c : cur->children) stack.push_back(&c);
        }
    }
    if (!source) { error = "找不到 ArtMesh source"; return false; }
    const std::vector<std::string> imports = {
        "<?import com.live2d.cubism.doc.model.extension.meshGenerator.CMeshGeneratorExtension?>",
        "<?import com.live2d.cubism.doc.model.extension.meshGenerator.MeshGenerateSetting?>"};
    const std::vector<std::string> versions = {
        "<?version CMeshGeneratorExtension:1?>", "<?version MeshGenerateSetting:1?>"};
    std::string declErr;
    xsp::ensureDecls(work, imports, versions, declErr);   // 已有声明时无副作用

    XsNode* list = nullptr;
    for (auto& c : source->children) if (c.fieldName == "_extensions") list = &c;
    if (!list) {
        XsNode l; l.type = "carray_list"; l.fieldName = "_extensions";
        l.attrs.emplace_back("xs.n", "_extensions"); l.attrs.emplace_back("count", "0");
        source->children.push_back(std::move(l));
        list = &source->children.back();
    }
    XsNode* ext = findTypeMut(*source, "CMeshGeneratorExtension");
    if (!ext) {
        XsNode e; e.type = "CMeshGeneratorExtension";
        XsNode setting; setting.type = "MeshGenerateSetting";
        setting.fieldName = "meshGenerateSetting"; setting.attrs.emplace_back("xs.n", "meshGenerateSetting");
        e.children.push_back(std::move(setting));
        list->children.push_back(std::move(e));
        ext = &list->children.back();
    }
    list->selfClosing = false;
    XsNode* setting = deepMut(*ext, "meshGenerateSetting");
    if (!setting) {
        XsNode s; s.type = "MeshGenerateSetting"; s.fieldName = "meshGenerateSetting";
        s.attrs.emplace_back("xs.n", "meshGenerateSetting");
        ext->children.push_back(std::move(s));
        setting = &ext->children.back();
    }
    auto put = [&](const char* name, int value) {
        XsNode* n = deepMut(*setting, name);
        if (!n) {
            XsNode v; v.type = "i"; v.fieldName = name; v.attrs.emplace_back("xs.n", name);
            setting->children.push_back(std::move(v));
            n = &setting->children.back();
        }
        n->text = std::to_string(value);
        n->selfClosing = false;
    };
    put("polygonOuterDensity", p.polygonOuterDensity);
    put("polygonInnerDensity", p.polygonInnerDensity);
    put("polygonMargin", p.polygonMargin);
    put("polygonInnerMargin", p.polygonInnerMargin);
    put("polygonMinMargin", p.polygonMinMargin);
    put("polygonMinBoundsPt", p.polygonMinBoundsPt);
    put("thresholdAlpha", p.thresholdAlpha);
    return true;
}

namespace {
bool generateInto(XsDocument& work, const Cmo3Document& doc, int meshIdx,
                           const MeshAutoGenParams& p, MeshAutoGenResult& result,
                           std::string& error, const std::vector<MeshFeatureLine>* features) {
    result = MeshAutoGenResult{};
    error.clear();
    if (meshIdx < 0 || meshIdx >= static_cast<int>(doc.artMeshes().size())) {
        error = "mesh 下标越界";
        return false;
    }
    const ArtMesh& m = doc.artMeshes()[static_cast<size_t>(meshIdx)];
    if (!m.source || !m.positionsNode || !m.uvsNode || !m.indicesNode) {
        error = "ArtMesh 缺 positions/uvs/indices";
        return false;
    }
    if (m.positionsNode->ref >= 0 || m.uvsNode->ref >= 0 || m.indicesNode->ref >= 0) {
        error = "positions/uvs/indices 是共享数组，拒绝自动生成（不做 COW）";
        return false;
    }
    std::vector<float> oldPts = floats(arrayText(work.root, m.positionsNode));
    std::vector<float> oldUv = floats(arrayText(work.root, m.uvsNode));
    std::vector<int> oldTris = ints(arrayText(work.root, m.indicesNode));
    const int nOld = m.vertexCount;
    if (nOld <= 0 || static_cast<int>(oldPts.size()) != nOld * 2 ||
        static_cast<int>(oldUv.size()) != nOld * 2 ||
        oldTris.empty() || oldTris.size() % 3 != 0) {
        error = "旧网格 positions/uvs/indices 不完整";
        return false;
    }
    for (int i : oldTris) {
        if (i < 0 || i >= nOld) { error = "旧网格三角形下标越界"; return false; }
    }
    for(float value:oldPts)if(!std::isfinite(value)){error="旧顶点坐标无效";return false;}
    for(float value:oldUv)if(!std::isfinite(value)){error="旧 UV 坐标无效";return false;}
    TextureDomain domain;
    if(!textureOutline(doc,m,p,domain,error))return false;
    meshSymmetry::Axis symmetry;
    if(p.autoSymmetry&&!domain.holes.empty())result.symmetryFallback="镂空轮廓使用完整约束三角化，暂不强制镜像采样";
    if(p.autoSymmetry&&domain.holes.empty()&&domain.contours.size()==1&&!domain.conservativeFallback) {
        symmetry=meshSymmetry::detect(domain.rgba,domain.width,domain.height);
        result.symmetryDetected=symmetry.detected;result.symmetryScore=symmetry.score;
        result.symmetryAxis={symmetry.x,symmetry.y,symmetry.dx,symmetry.dy};
        if(symmetry.detected) {
            const bool eligible=!meshSymmetry::pairedContour(domain.contours[0],symmetry,true,.5).empty();
            auto paired=eligible?pairedPolygon(domain.contours[0],symmetry,true,
                meshSampling::pointInterval(domain.width,domain.height,p.polygonOuterDensity)):std::vector<float>{};
            if(paired.empty()){symmetry.detected=false;result.symmetryFallback="轮廓不适合单轴截面";}
            else {domain.contours[0]=std::move(paired);materialOutline(domain.contours[0],p.style,meshSampling::pointInterval(domain.width,domain.height,p.polygonOuterDensity),&symmetry);result.symmetrySampling=true;}
        }
    }
    // Triangulate the source image in pixel UV space. A deformed source mesh
    // must not change which texture pixels the triangles cover.
    std::vector<float> oldPixel=oldUv;
    for(size_t i=0;i<oldPixel.size();i+=2){const auto pixelPoint=domain.pixel(oldUv[i],oldUv[i+1]);oldPixel[i]=pixelPoint.first;oldPixel[i+1]=pixelPoint.second;}
    std::vector<float> pts=oldPixel;
    std::vector<VertexLoop> loops;
    auto addPoint=[&](float x,float y,bool reuseOld=true) {
        for(size_t j=reuseOld?0:size_t(nOld);j<pts.size()/2;++j)if(std::hypot(pts[2*j]-x,pts[2*j+1]-y)<1e-5f)return int(j);
        const int index=int(pts.size()/2);pts.push_back(x);pts.push_back(y);return index;
    };
    for(const auto& contour:domain.contours) {
        VertexLoop loop;
        for(size_t i=0;i<contour.size();i+=2)loop.push_back(addPoint(contour[i],contour[i+1]));
        if(loop.size()<3||!simpleOutline(pts,loop,1e-6)){error="纹理边界自交或自接触，拒绝写回";return false;}
        loops.push_back(std::move(loop));
    }
    const size_t outerCount=loops.size();
    std::vector<int> holeOwners;
    for(const auto& hole:domain.holes){VertexLoop loop;for(size_t i=0;i<hole.size();i+=2)loop.push_back(addPoint(hole[i],hole[i+1]));
        int owner=-1;double area=std::numeric_limits<double>::max();
        for(size_t i=0;i<outerCount;++i)if(pointInLoop(pts,loops[i],hole[0],hole[1])){const double candidate=std::abs(loopArea(pts,loops[i]));if(candidate<area){area=candidate;owner=int(i);}}
        if(owner<0||!simpleOutline(pts,loop,1e-6)){error="镂空轮廓无法对应外轮廓";return false;}
        holeOwners.push_back(owner);loops.push_back(std::move(loop));}
    std::set<int> protectedOld;
    for(const auto& glue:doc.glues()) {
        if(glue.artMeshA==meshIdx)protectedOld.insert(glue.pairA.begin(),glue.pairA.end());
        if(glue.artMeshB==meshIdx)protectedOld.insert(glue.pairB.begin(),glue.pairB.end());
    }
    for(int v:protectedOld)if(v<0||v>=nOld){error="胶水绑定顶点下标无效";return false;}
    float minX=float(domain.width),minY=float(domain.height),maxX=0,maxY=0;
    for(const auto& loop:loops)for(int i:loop) {
        minX=std::min(minX,pts[2*i]);maxX=std::max(maxX,pts[2*i]);
        minY=std::min(minY,pts[2*i+1]);maxY=std::max(maxY,pts[2*i+1]);
    }
    const float w=maxX-minX,h=maxY-minY;
    const float boundaryTolerance=std::max(1e-6f,std::max(w,h)*1e-7f);
    auto component=[&](float x,float y) {
        int owner=-1;double area=std::numeric_limits<double>::max();
        for(size_t i=0;i<outerCount;++i)if(pointInLoop(pts,loops[i],x,y)){const double candidate=std::abs(loopArea(pts,loops[i]));if(candidate<area){area=candidate;owner=int(i);}}
        for(size_t i=0;i<holeOwners.size();++i)if(holeOwners[i]==owner&&pointInLoop(pts,loops[outerCount+i],x,y)&&!nearBoundary(pts,loops[outerCount+i],x,y,boundaryTolerance))return -1;
        return owner;
    };
    auto onBoundary=[&](float x,float y,float margin) {
        for(const auto& loop:loops)if(nearBoundary(pts,loop,x,y,margin))return true;
        return false;
    };
    double totalArea=0;std::vector<float> boundaryLengths;
    for(size_t li=0;li<loops.size();++li) {const auto& loop=loops[li];
        totalArea+=(li<outerCount?1:-1)*std::abs(loopArea(pts,loop));
        for(size_t i=0;i<loop.size();++i){const int a=loop[i],b=loop[(i+1)%loop.size()];boundaryLengths.push_back(std::hypot(pts[2*b]-pts[2*a],pts[2*b+1]-pts[2*a+1]));}
    }
    std::sort(boundaryLengths.begin(),boundaryLengths.end());
    const float medianBoundary=boundaryLengths[boundaryLengths.size()/2];
    const double nominalSpacing=p.interiorTriangleTarget>0?
        std::sqrt(totalArea/(.433012701892*std::clamp(p.interiorTriangleTarget,1,4000))):
        meshSampling::pointInterval(domain.width,domain.height,p.polygonInnerDensity);
    const int target=p.polygonInnerDensity>0?std::clamp(int(std::ceil(totalArea/(.866025403784*nominalSpacing*nominalSpacing))),1,2000):0;
    auto localSpacing=[&](float x,float y) {
        double nearest=std::numeric_limits<double>::max(),edgeLength=medianBoundary;
        for(const auto& loop:loops)for(size_t i=0;i<loop.size();++i) {
            const int a=loop[i],b=loop[(i+1)%loop.size()];
            const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],length=dx*dx+dy*dy;
            const double t=length?std::clamp(((x-pts[2*a])*dx+(y-pts[2*a+1])*dy)/length,0.,1.):0;
            const double distance=std::hypot(x-pts[2*a]-t*dx,y-pts[2*a+1]-t*dy);
            if(distance<nearest){nearest=distance;edgeLength=std::sqrt(length);}
        }
        const double factor=std::clamp(edgeLength/std::max(.01f,medianBoundary),.55,1.);
        return nominalSpacing*(factor+(1-factor)*std::min(1.,nearest/std::max(.01,nominalSpacing*2)));
    };
    std::vector<std::pair<int,int>> featureEdges,flowEdges;
    std::vector<CurvePoint> weakHints;
    if(features) {
        std::vector<std::vector<float>> paths;
        for(const auto& line:*features)if(line.enabled) {
            if(line.points.size()<4||line.points.size()%2){error="结构线坐标不完整";return false;}
            if(line.automatic&&line.weakReference){
                for(size_t i=0;i<line.points.size();i+=2){const auto q=domain.pixel(line.points[i],line.points[i+1]);
                    if(component(q.first,q.second)>=0&&!onBoundary(q.first,q.second,boundaryTolerance))weakHints.push_back({q.first,q.second});}
                continue;
            }
            auto path=line.points;int owner=-1;
            for(size_t i=0;i<path.size();i+=2) {
                const auto pixelPoint=domain.pixel(path[i],path[i+1]);path[i]=pixelPoint.first;path[i+1]=pixelPoint.second;
                if(!std::isfinite(path[i])||!std::isfinite(path[i+1])){error="结构线坐标无效";return false;}
                const int at=component(path[i],path[i+1]);
                if(at<0||(owner>=0&&at!=owner)||onBoundary(path[i],path[i+1],boundaryTolerance)) {
                    error="结构线跨出纹理轮廓，请忽略该线后重试（位置 "+std::to_string(path[i])+", "+std::to_string(path[i+1])+"）";return false;
                }
                owner=at;
            }
            if(line.automatic&&line.weakReference){for(size_t i=0;i<path.size();i+=2)weakHints.push_back({path[i],path[i+1]});}
            else paths.push_back(std::move(path));
        }
        std::set<int> retained=protectedOld;
        for(const auto& loop:loops)retained.insert(loop.begin(),loop.end());
        std::vector<float> constraintPoints;std::vector<int> constraintToOriginal;
        for(int v:retained){constraintPoints.push_back(pts[2*v]);constraintPoints.push_back(pts[2*v+1]);constraintToOriginal.push_back(v);}
        const size_t existing=constraintToOriginal.size();std::vector<std::pair<int,int>> constraints;
        if(!appendMeshFeatureConstraints(constraintPoints,paths,constraints,error))return false;
        for(size_t i=existing;i<constraintPoints.size()/2;++i)
            constraintToOriginal.push_back(addPoint(constraintPoints[2*i],constraintPoints[2*i+1]));
        for(auto [a,b]:constraints)featureEdges.push_back(std::minmax(constraintToOriginal[a],constraintToOriginal[b]));
        std::sort(featureEdges.begin(),featureEdges.end());featureEdges.erase(std::unique(featureEdges.begin(),featureEdges.end()),featureEdges.end());
        for(auto [a,b]:featureEdges)for(const auto& loop:loops)for(size_t k=0;k<loop.size();++k) {
            const int c=loop[k],d=loop[(k+1)%loop.size()];
            auto cross=[&](int i,int j,int v) {
                return (double(pts[2*j])-pts[2*i])*(double(pts[2*v+1])-pts[2*i+1])-
                       (double(pts[2*j+1])-pts[2*i+1])*(double(pts[2*v])-pts[2*i]);
            };
            const double ac=cross(a,b,c),ad=cross(a,b,d),ca=cross(c,d,a),cb=cross(c,d,b);
            const double length=std::hypot(double(pts[2*b])-pts[2*a],double(pts[2*b+1])-pts[2*a+1]);
            auto touches=[&](int v,double area) {
                return std::abs(area)<=boundaryTolerance*length&&
                    pts[2*v]>=std::min(pts[2*a],pts[2*b])-boundaryTolerance&&pts[2*v]<=std::max(pts[2*a],pts[2*b])+boundaryTolerance&&
                    pts[2*v+1]>=std::min(pts[2*a+1],pts[2*b+1])-boundaryTolerance&&pts[2*v+1]<=std::max(pts[2*a+1],pts[2*b+1])+boundaryTolerance;
            };
            if(((ac<0)!=(ad<0)&&(ca<0)!=(cb<0))||touches(c,ac)||touches(d,ad)) {
                error="结构线跨出或触及凹形轮廓，请忽略该线";return false;
            }
        }
    }
    // Long crease edges can have intermediate mesh vertices while their
    // detected semantic anchors remain endpoints/turns. With interior fill disabled,
    // retain the exact minimal endpoint-only topology.
    if(target&&p.style!=MeshAutoGenStyle::Rigid&&p.style!=MeshAutoGenStyle::FineLine) {
        std::vector<std::pair<int,int>> split;
        for(auto [a,b]:featureEdges) {
            const float ax=pts[2*a],ay=pts[2*a+1],bx=pts[2*b],by=pts[2*b+1];
            const double spacing=std::min(localSpacing((ax+bx)*.5f,(ay+by)*.5f)*1.2,double(medianBoundary)*2);
            const int pieces=std::clamp(int(std::ceil(std::hypot(bx-ax,by-ay)/std::max(.1,spacing))),1,256);
            int previous=a;
            for(int k=1;k<=pieces;++k) {
                const int vertex=k==pieces?b:addPoint(ax+(bx-ax)*k/pieces,ay+(by-ay)*k/pieces);
                split.push_back(std::minmax(previous,vertex));previous=vertex;
            }
        }
        featureEdges=std::move(split);std::sort(featureEdges.begin(),featureEdges.end());
        featureEdges.erase(std::unique(featureEdges.begin(),featureEdges.end()),featureEdges.end());
    }
    auto distanceToFeatures=[&](float x,float y) {
        double nearest=std::numeric_limits<double>::max();
        auto visit=[&](const std::vector<std::pair<int,int>>& edges) {for(auto [a,b]:edges) {
            const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],length=dx*dx+dy*dy;
            const double t=length?std::clamp(((x-pts[2*a])*dx+(y-pts[2*a+1])*dy)/length,0.,1.):0;
            nearest=std::min(nearest,std::hypot(x-pts[2*a]-t*dx,y-pts[2*a+1]-t*dy));
        }};visit(featureEdges);visit(flowEdges);
        return nearest;
    };
    const float preferredMargin=float(std::max(0,p.polygonInnerMargin)),minimumMargin=float(std::max(0,p.polygonMinMargin));
    std::set<int> anchors=protectedOld;
    for(const auto& loop:loops)anchors.insert(loop.begin(),loop.end());
    for(auto [a,b]:featureEdges){anchors.insert(a);anchors.insert(b);}
    std::vector<VertexLoop> borderLoops;
    std::vector<bool> borderIsHole;
    std::vector<std::pair<int,int>> borderEdges;
    const int beforeBorder=int(pts.size()/2);
    if(p.generateBorder||!domain.holeSupports.empty()) {
        const int cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
        const double outsideInterval=meshSampling::pointInterval(domain.width,domain.height,p.polygonOuterDensity);
        const double borderInterval=meshSampling::pointInterval(domain.width,domain.height,
            p.polygonBorderDensity>0?p.polygonBorderDensity:p.polygonOuterDensity);
        const double borderRatio=borderInterval/outsideInterval;
        const double spacing=std::max(4.,std::min(borderInterval,std::max(double(medianBoundary)*borderRatio,borderInterval*.5)));
        const double minimumRingInterval=std::max(cell*1.5,double(medianBoundary)*borderRatio*.45);
        const double depth=std::max(1.,std::min(preferredMargin>0?double(preferredMargin):spacing*.55,spacing*.75));
        MeshAutoGenParams ringParams=p;ringParams.polygonMargin=0;ringParams.polygonMinBoundsPt=3;
        ringParams.polygonOuterDensity=std::max(1,int(std::lround(spacing*2*meshSampling::workScale(domain.width,domain.height))));
        std::vector<std::vector<float>> ringPaths;
        auto segmentOnTexture=[&](CurvePoint a,CurvePoint b) {
            const int probes=std::max(1,int(std::ceil(std::hypot(b.x-a.x,b.y-a.y)/.75)));
            int owner=-1;
            for(int k=0;k<=probes;++k) {
                const float x=float(a.x+(b.x-a.x)*k/probes),y=float(a.y+(b.y-a.y)*k/probes);
                const int xx=int(std::floor(x)),yy=int(std::floor(y)),part=component(x,y);
                if(xx<0||xx>=domain.width||yy<0||yy>=domain.height||!domain.rgba[(size_t(yy)*domain.width+xx)*4+3]||
                   part<0||(owner>=0&&owner!=part)||onBoundary(x,y,boundaryTolerance))return false;
                owner=part;
            }
            return true;
        };
        auto validRing=[&](const std::vector<float>& path) {
            VertexLoop ring(path.size()/2);std::iota(ring.begin(),ring.end(),0);
            if(ring.size()<3||!simpleOutline(path,ring,1e-5))return false;
            for(size_t i=0;i<ring.size();++i) {
                const size_t j=(i+1)%ring.size();
                if(!segmentOnTexture({path[2*i],path[2*i+1]},{path[2*j],path[2*j+1]}))return false;
            }
            return true;
        };
        for(auto [fraction,clearance]:std::initializer_list<std::pair<double,double>>{{1.,2.25},{.5,2.25},{.25,2.25},{1.,3.5},{1.,5.}}) {
            if(!p.generateBorder)break;
            TextureDomain inset;inset.rgba=domain.rgba;inset.width=domain.width;inset.height=domain.height;
            std::string insetError;
            // The raster-scale fit can move outward by up to 2.25 cells.
            // Reserve that clearance so the resulting ring stays on texture.
            if(!buildTextureContours(inset,ringParams,insetError,depth*fraction+clearance*cell)||inset.conservativeFallback)continue;
            auto usableRing=[&](const std::vector<float>& contour) {
                const int owner=component(contour[0],contour[1]);if(owner<0)return false;
                VertexLoop ids(contour.size()/2);std::iota(ids.begin(),ids.end(),0);
                const double area=std::abs(loopArea(contour,ids));
                // Eroding a very narrow texture may leave isolated one-pixel
                // specks. They are not additional support loops. Their source
                // texture is still covered by the untouched outer contour.
                return std::abs(loopArea(pts,loops[size_t(owner)]))<cell*cell*50||
                    (area>=std::max(cell*cell*6.,spacing*spacing*.015)&&ClosedCurve(contour).arc.back()>=spacing*.75);
            };
            for(const auto& contour:inset.contours) {
                if(!usableRing(contour))continue;
                const int owner=component(contour[0],contour[1]);if(owner<0)continue;
                std::vector<float> outside;
                for(int v:loops[size_t(owner)]){outside.push_back(pts[2*v]);outside.push_back(pts[2*v+1]);}
                const ClosedCurve outer(outside);
                const ClosedCurve curve(contour);std::vector<float> path;
                // Carry the outer ring's arc stations to its inset curve.
                // Independently sampling by the shorter inset perimeter loses
                // both point correspondence and curvature detail at the ends.
                double phase=0,phaseDistance=std::numeric_limits<double>::max();
                for(size_t i=0;i<contour.size()/2;++i) {
                    const size_t j=(i+1)%(contour.size()/2);
                    const double dx=contour[2*j]-contour[2*i],dy=contour[2*j+1]-contour[2*i+1],squared=dx*dx+dy*dy;
                    const double t=squared?std::clamp(((outside[0]-contour[2*i])*dx+(outside[1]-contour[2*i+1])*dy)/squared,0.,1.):0;
                    const double distance=std::hypot(outside[0]-contour[2*i]-t*dx,outside[1]-contour[2*i+1]-t*dy);
                    if(distance<phaseDistance){phaseDistance=distance;phase=curve.arc[i]+t*std::sqrt(squared);}
                }
                double siblingLength=0;
                for(const auto& sibling:inset.contours)if(usableRing(sibling)&&component(sibling[0],sibling[1])==owner)siblingLength+=ClosedCurve(sibling).arc.back();
                const bool splitRegion=siblingLength>curve.arc.back()*1.001;
                const int wanted=std::clamp(int(std::lround(loops[size_t(owner)].size()/borderRatio*curve.arc.back()/siblingLength)),3,1000);
                const int stations=std::min(wanted,std::max(3,int(std::floor(curve.arc.back()/minimumRingInterval))));
                std::vector<double> starts;
                if(!splitRegion&&stations==int(loops[size_t(owner)].size()))for(size_t i=0;i+1<outer.arc.size();++i)starts.push_back(phase+outer.arc[i]/outer.arc.back()*curve.arc.back());
                else for(int k=0;k<stations;++k)starts.push_back(phase+curve.arc.back()*k/stations);
                starts.push_back(phase+curve.arc.back());
                auto sample=[&](auto&& self,double start,double end,int level)->void {
                    const auto a=curve.at(start),b=curve.at(end);const double dx=b.x-a.x,dy=b.y-a.y,squared=dx*dx+dy*dy;
                    double error=0;
                    for(size_t i=0;i<contour.size()/2;++i) {
                        double at=curve.arc[i];while(at<start)at+=curve.arc.back();
                        if(at<end) {
                            const double x=contour[2*i]-a.x,y=contour[2*i+1]-a.y,t=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                            error=std::max(error,std::hypot(x-t*dx,y-t*dy));
                        }
                    }
                    if((error>std::max(double(cell),std::min(2.,spacing*.1))||!segmentOnTexture(a,b))&&end-start>cell*1.5&&level<10) {
                        self(self,start,(start+end)*.5,level+1);self(self,(start+end)*.5,end,level+1);
                    } else {path.push_back(float(a.x));path.push_back(float(a.y));}
                };
                for(size_t k=1;k<starts.size();++k)sample(sample,starts[k-1],starts[k],0);
                if(symmetry.detected)path=pairedPolygon(path,symmetry,false,spacing);
                // Small asymmetric notches may lie inside the fitted ring.
                // Retry with paired inward clearance, bounded to two percent,
                // rather than keeping an off-texture or one-sided support row.
                for(int repair=0;symmetry.detected&&!path.empty()&&!validRing(path)&&repair<8;++repair)
                    for(size_t i=0;i<path.size();i+=2){path[i]=float(symmetry.x+(path[i]-symmetry.x)*.9975);path[i+1]=float(symmetry.y+(path[i+1]-symmetry.y)*.9975);}
                if(!validRing(path))continue;
                // Curvature/alpha repairs must not leave almost coincident
                // stations on a short inset. Remove a station only when the
                // replacement chord stays on texture and near the curve.
                for(int pass=0;!symmetry.detected&&pass<512&&path.size()>6;++pass) {
                    const size_t count=path.size()/2;bool removed=false;
                    for(size_t i=0;i<count;++i) {
                        const size_t a=(i+count-1)%count,b=(i+1)%count;
                        const double dx=path[2*b]-path[2*a],dy=path[2*b+1]-path[2*a+1],squared=dx*dx+dy*dy;
                        const double x=path[2*i]-path[2*a],y=path[2*i+1]-path[2*a+1];
                        if(std::min(std::hypot(x,y),double(std::hypot(path[2*i]-path[2*b],path[2*i+1]-path[2*b+1])))>=minimumRingInterval)continue;
                        const double t=squared?std::clamp((x*dx+y*dy)/squared,0.,1.):0;
                        if(std::hypot(x-t*dx,y-t*dy)>std::max(cell*.75,std::min(2.,minimumRingInterval*.12)))continue;
                        auto proposal=path;proposal.erase(proposal.begin()+2*i,proposal.begin()+2*i+2);
                        if(validRing(proposal)){path=std::move(proposal);removed=true;break;}
                    }
                    if(!removed)break;
                }
                // Move the sampling budget out of crowded bends onto long
                // spans. A normal inset keeps a population close to the
                // outside row; a genuinely short/narrow inset is capacity
                // limited by minimumRingInterval instead.
                while(!symmetry.detected&&path.size()/2<size_t(stations)) {
                    size_t best=0;double longest=minimumRingInterval*1.5;
                    for(size_t i=0;i<path.size()/2;++i) {
                        const size_t j=(i+1)%(path.size()/2);const double length=std::hypot(path[2*i]-path[2*j],path[2*i+1]-path[2*j+1]);
                        if(length>longest){longest=length;best=i;}
                    }
                    if(longest<=minimumRingInterval*1.5)break;
                    const size_t next=(best+1)%(path.size()/2);
                    const float x=(path[2*best]+path[2*next])*.5f,y=(path[2*best+1]+path[2*next+1])*.5f;
                    path.insert(path.begin()+2*best+2,{x,y});
                }
                // Rigid support favors straight spans and retains sharp turns.
                // Remove only near-collinear stations with bounded deviation
                // from the ORIGINAL ring. The alpha-contained closed ring and
                // requested spacing budget remain authoritative. Mirrored
                // stations are removed as pairs to preserve detected symmetry.
                if(p.style==MeshAutoGenStyle::Rigid) {
                    const auto original=path;
                    const double tolerance=std::clamp(spacing*.025,.35,1.25);
                    auto removable=[&](size_t i) {
                        const size_t count=path.size()/2,a=(i+count-1)%count,b=(i+1)%count;
                        const double ux=path[2*i]-path[2*a],uy=path[2*i+1]-path[2*a+1];
                        const double vx=path[2*b]-path[2*i],vy=path[2*b+1]-path[2*i+1];
                        const double l=std::hypot(ux,uy)*std::hypot(vx,vy);
                        return l>1e-8&&(ux*vx+uy*vy)/l>.9781476&&std::hypot(ux+vx,uy+vy)<=spacing*3;
                    };
                    auto nearOriginal=[&](const std::vector<float>& candidate) {
                        for(size_t k=0;k<original.size();k+=2) {
                            double nearest=1e30;
                            for(size_t i=0;i<candidate.size();i+=2) {
                                const size_t j=(i+2)%candidate.size();const double dx=candidate[j]-candidate[i],dy=candidate[j+1]-candidate[i+1],sq=dx*dx+dy*dy;
                                const double t=sq?std::clamp(((original[k]-candidate[i])*dx+(original[k+1]-candidate[i+1])*dy)/sq,0.,1.):0;
                                nearest=std::min(nearest,std::hypot(original[k]-candidate[i]-t*dx,original[k+1]-candidate[i+1]-t*dy));
                            }
                            if(nearest>tolerance)return false;
                        }
                        return true;
                    };
                    for(int pass=0;pass<1000&&path.size()>6;++pass) {
                        bool changed=false;
                        for(size_t i=0;i<path.size()/2;++i)if(removable(i)) {
                            size_t mate=i;
                            if(symmetry.detected) {
                                const auto reflected=symmetry.reflect(path[2*i],path[2*i+1]);double closest=.005;
                                for(size_t j=0;j<path.size()/2;++j){const double d=std::hypot(path[2*j]-reflected.x,path[2*j+1]-reflected.y);if(d<closest){closest=d;mate=j;}}
                                if(mate==i||!removable(mate))continue; // preserve axis tips
                            }
                            if(path.size()/2<=(mate==i?3:4))continue;
                            auto proposal=path;
                            for(size_t index:std::set<size_t,std::greater<size_t>>{i,mate})proposal.erase(proposal.begin()+index*2,proposal.begin()+index*2+2);
                            if(!nearOriginal(proposal)||!validRing(proposal))continue;
                            result.rigidBorderRemoved+=int((path.size()-proposal.size())/2);path=std::move(proposal);changed=true;break;
                        }
                        if(!changed)break;
                    }
                }
                // Align each regular support station with an outer-edge
                // midpoint at a common inward offset. Start from the alpha-safe
                // ring, retain corners/narrow exceptions, and accept only moves
                // that improve the local isosceles triangle without crowding
                // its neighboring support stations. Mirror moves are coupled.
                auto nearestEdge=[&](double x,double y) {
                    size_t edge=0;double nearest=1e30;
                    for(size_t k=0;k<outside.size();k+=2){const size_t j=(k+2)%outside.size();const double dx=outside[j]-outside[k],dy=outside[j+1]-outside[k+1],sq=dx*dx+dy*dy;
                        const double t=sq?std::clamp(((x-outside[k])*dx+(y-outside[k+1])*dy)/sq,0.,1.):0;
                        const double d=std::hypot(x-outside[k]-t*dx,y-outside[k+1]-t*dy);if(d<nearest){nearest=d;edge=k;}}
                    return std::make_pair(edge,nearest);
                };
                std::vector<double> offsets;for(size_t k=0;k<path.size();k+=2)offsets.push_back(nearestEdge(path[k],path[k+1]).second);
                std::sort(offsets.begin(),offsets.end());
                const double commonOffset=offsets[offsets.size()/2];
                auto triangleQuality=[&](size_t edge,double x,double y) {
                    const size_t j=(edge+2)%outside.size();const double dx=outside[j]-outside[edge],dy=outside[j+1]-outside[edge+1];
                    const double ax=x-outside[edge],ay=y-outside[edge+1],bx=x-outside[j],by=y-outside[j+1];
                    return 3.464101615137754*std::abs(dx*ay-dy*ax)/std::max(1e-12,dx*dx+dy*dy+ax*ax+ay*ay+bx*bx+by*by);
                };
                // Also protect the two oppositely oriented triangles based on
                // neighboring inner-ring edges. Improving only the outer-edge
                // triangle can squeeze these triangles into narrow slivers.
                auto adjacentQuality=[&](const std::vector<float>& ring,size_t a,size_t b,size_t outerIndex) {
                    const double dx=ring[2*b]-ring[2*a],dy=ring[2*b+1]-ring[2*a+1];
                    const double ax=outside[outerIndex]-ring[2*a],ay=outside[outerIndex+1]-ring[2*a+1];
                    const double bx=ax-dx,by=ay-dy;
                    return 3.464101615137754*std::abs(dx*ay-dy*ax)/std::max(1e-12,dx*dx+dy*dy+ax*ax+ay*ay+bx*bx+by*by);
                };
                auto corner=[&](size_t i) {const size_t n=path.size()/2,a=(i+n-1)%n,b=(i+1)%n;
                    const double ax=path[2*i]-path[2*a],ay=path[2*i+1]-path[2*a+1],bx=path[2*b]-path[2*i],by=path[2*b+1]-path[2*i+1];
                    return (ax*bx+ay*by)/std::max(1e-12,std::hypot(ax,ay)*std::hypot(bx,by))<.9781476;};
                std::set<size_t> aligned;
                for(size_t i=0;p.staggerBorder&&i<path.size()/2;++i) {
                    if(corner(i)||aligned.count(i))continue;
                    const auto [edge,distance]=nearestEdge(path[2*i],path[2*i+1]);const size_t end=(edge+2)%outside.size();
                    const double dx=outside[end]-outside[edge],dy=outside[end+1]-outside[edge+1],length=std::hypot(dx,dy);
                    if(length<1e-5||std::abs(distance-commonOffset)>std::max(double(cell)*2,commonOffset*.35))continue;
                    const double mx=(outside[edge]+outside[end])*.5,my=(outside[edge+1]+outside[end+1])*.5;
                    const double side=dx*(path[2*i+1]-my)-dy*(path[2*i]-mx)>0?1.:-1.;
                    // Cloth follows each local bend's offset rather than
                    // flattening the row onto one global median distance.
                    const double offset=p.style==MeshAutoGenStyle::Soft?distance:commonOffset;
                    const double tx=mx-side*dy/length*offset,ty=my+side*dx/length*offset;
                    const double before=triangleQuality(edge,path[2*i],path[2*i+1]);if(before>=.85)continue;
                    size_t mate=i;
                    if(symmetry.detected){const auto mirror=symmetry.reflect(path[2*i],path[2*i+1]);double nearest=.005;
                        for(size_t k=0;k<path.size()/2;++k){const double d=std::hypot(path[2*k]-mirror.x,path[2*k+1]-mirror.y);if(d<nearest){nearest=d;mate=k;}}
                        if(mate==i||corner(mate))continue;}
                    for(double amount:{1.,.5,.25}) {
                        auto proposal=path;proposal[2*i]=float(path[2*i]+(tx-path[2*i])*amount);proposal[2*i+1]=float(path[2*i+1]+(ty-path[2*i+1])*amount);
                        if(triangleQuality(edge,proposal[2*i],proposal[2*i+1])<before+.01)continue;
                        if(mate!=i){const auto mirrored=symmetry.reflect(proposal[2*i],proposal[2*i+1]);proposal[2*mate]=float(mirrored.x);proposal[2*mate+1]=float(mirrored.y);}
                        bool spaced=true;
                        for(size_t v:std::set<size_t>{i,mate}) {
                            const size_t n=path.size()/2,a=(v+n-1)%n,b=(v+1)%n;
                            const double oldA=std::hypot(path[2*v]-path[2*a],path[2*v+1]-path[2*a+1]),oldB=std::hypot(path[2*v]-path[2*b],path[2*v+1]-path[2*b+1]);
                            const double newA=std::hypot(proposal[2*v]-proposal[2*a],proposal[2*v+1]-proposal[2*a+1]),newB=std::hypot(proposal[2*v]-proposal[2*b],proposal[2*v+1]-proposal[2*b+1]);
                            const auto oldEdge=nearestEdge(path[2*v],path[2*v+1]).first;
                            spaced&=std::min(newA,newB)>=std::min(std::min(oldA,oldB),minimumRingInterval*.75)&&std::max(newA,newB)<=std::max(oldA,oldB)*1.15&&
                                triangleQuality(oldEdge,proposal[2*v],proposal[2*v+1])>=triangleQuality(oldEdge,path[2*v],path[2*v+1])+.009;
                            for(size_t neighbor:{a,b}) {
                                const double midX=(path[2*v]+path[2*neighbor])*.5,midY=(path[2*v+1]+path[2*neighbor+1])*.5;
                                size_t outerIndex=0;double closest=1e30;
                                for(size_t k=0;k<outside.size();k+=2){const double d=std::hypot(outside[k]-midX,outside[k+1]-midY);if(d<closest){closest=d;outerIndex=k;}}
                                spaced&=adjacentQuality(proposal,v,neighbor,outerIndex)+.005>=adjacentQuality(path,v,neighbor,outerIndex);
                            }
                        }
                        if(!spaced||!validRing(proposal))continue;
                        path=std::move(proposal);aligned.insert(i);aligned.insert(mate);break;
                    }
                }
                result.staggeredBorderVertices+=int(aligned.size());
                // Internal feature settings cannot move this fitted row.
                // Crossings are shared intersections on its existing segments.
                path.push_back(path[0]);path.push_back(path[1]);ringPaths.push_back(std::move(path));
            }
            if(!ringPaths.empty())break; // narrow regions retry with a shallower inset
        }
        for(auto path:domain.holeSupports){path.push_back(path[0]);path.push_back(path[1]);ringPaths.push_back(std::move(path));}
        if(!ringPaths.empty()) {
            std::vector<std::vector<float>> creasePaths;
            for(auto [a,b]:featureEdges)creasePaths.push_back({pts[2*a],pts[2*a+1],pts[2*b],pts[2*b+1]});
            auto allPaths=creasePaths;allPaths.insert(allPaths.end(),ringPaths.begin(),ringPaths.end());
            std::vector<float> constraintPoints;std::vector<int> constraintToOriginal;
            for(int v:anchors){constraintPoints.push_back(pts[2*v]);constraintPoints.push_back(pts[2*v+1]);constraintToOriginal.push_back(v);}
            const size_t existing=constraintToOriginal.size();std::vector<std::pair<int,int>> allEdges,creaseCuts,ringCuts;
            // Split crease/ring crossings into shared vertices before recovering
            // either set of constraints. A closed ring is an internal boundary,
            // never an outline to clip away the central texture.
            if(!appendMeshFeatureConstraints(constraintPoints,allPaths,allEdges,error)||
               !appendMeshFeatureConstraints(constraintPoints,creasePaths,creaseCuts,error)||
               !appendMeshFeatureConstraints(constraintPoints,ringPaths,ringCuts,error))return false;
            for(size_t i=existing;i<constraintPoints.size()/2;++i)
                constraintToOriginal.push_back(addPoint(constraintPoints[2*i],constraintPoints[2*i+1],false));
            featureEdges.clear();
            for(auto [a,b]:creaseCuts)featureEdges.push_back(std::minmax(constraintToOriginal[a],constraintToOriginal[b]));
            for(auto [a,b]:ringCuts)borderEdges.push_back(std::minmax(constraintToOriginal[a],constraintToOriginal[b]));
            for(auto [a,b]:allEdges){anchors.insert(constraintToOriginal[a]);anchors.insert(constraintToOriginal[b]);}
            for(size_t pathIndex=0;pathIndex<ringPaths.size();++pathIndex) {const auto& path=ringPaths[pathIndex];
                VertexLoop loop;
                for(size_t i=0;i+2<path.size();i+=2) {
                    for(size_t v=0;v<constraintToOriginal.size();++v)
                        if(std::hypot(path[i]-constraintPoints[2*v],path[i+1]-constraintPoints[2*v+1])<boundaryTolerance*2) {
                            loop.push_back(constraintToOriginal[v]);break;
                        }
                }
                if(loop.size()<3){error="锁边环顶点映射失败";return false;}
                borderLoops.push_back(std::move(loop));
                borderIsHole.push_back(pathIndex>=ringPaths.size()-domain.holeSupports.size());
            }
        }
    }
    if(symmetry.detected&&p.generateBorder&&borderLoops.empty()){error="对称锁边无法完整落在纹理内部";return false;}
    auto insideBorder=[&](float x,float y) {
        const int part=component(x,y);if(part<0)return false;bool hasRing=false,inside=false;
        for(size_t i=0;i<borderLoops.size();++i){const auto& loop=borderLoops[i];if(component(pts[2*loop[0]],pts[2*loop[0]+1])!=part)continue;
            const bool contains=pointInLoop(pts,loop,x,y);
            if(borderIsHole[i]){if(contains)return false;}
            else {hasRing=true;inside|=contains;}
        }
        return !hasRing||inside;
    };
    auto distanceToBorder=[&](float x,float y) {
        double nearest=std::numeric_limits<double>::max();
        for(auto [a,b]:borderEdges) {
            const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],squared=dx*dx+dy*dy;
            const double t=squared?std::clamp(((x-pts[2*a])*dx+(y-pts[2*a+1])*dy)/squared,0.,1.):0;
            nearest=std::min(nearest,std::hypot(x-pts[2*a]-t*dx,y-pts[2*a+1]-t*dy));
        }
        return nearest;
    };
    auto borderClearance=[&](double spacing) {
        return p.interiorBorderMargin>0?double(p.interiorBorderMargin):spacing*.4;
    };
    // Organize a single-stroke layer along its actual medial curve, including
    // wide/bent strips. Cross sections follow local width, not image X/Y axes.
    // The corresponding border rows already support a narrow strip: only add
    // interior columns when there is room for another ordinary point interval.
    const auto& guide=domain.strokeCenterline.empty()?domain.stripCenterline:domain.strokeCenterline;
    const bool strip=target&&!symmetry.detected&&!guide.empty()&&p.style!=MeshAutoGenStyle::Rigid;
    const int beforeFlow=int(pts.size()/2);
    if(strip&&p.style!=MeshAutoGenStyle::FineLine) {
        std::vector<double> arc(guide.size()/2);
        for(size_t i=1;i<arc.size();++i)arc[i]=arc[i-1]+std::hypot(guide[2*i]-guide[2*i-2],guide[2*i+1]-guide[2*i-1]);
        const double cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
        const double flowBudget=std::max(16.,(1900.-featureEdges.size()-borderEdges.size())/3);
        const double spacing=std::max({cell*3,nominalSpacing,double(medianBoundary),std::sqrt(totalArea/flowBudget)});
        auto at=[&](double distance) {
            distance=std::clamp(distance,0.,arc.back());
            const size_t i=std::clamp(size_t(std::upper_bound(arc.begin(),arc.end(),distance)-arc.begin()),size_t(1),arc.size()-1);
            const double f=(distance-arc[i-1])/std::max(1e-12,arc[i]-arc[i-1]);
            return CurvePoint{guide[2*i-2]+f*(guide[2*i]-guide[2*i-2]),guide[2*i-1]+f*(guide[2*i+1]-guide[2*i-1])};
        };
        auto onTexture=[&](CurvePoint c) {
            const int x=int(std::floor(c.x)),y=int(std::floor(c.y));
            return x>=0&&x<domain.width&&y>=0&&y<domain.height&&domain.rgba[(size_t(y)*domain.width+x)*4+3]&&insideBorder(float(c.x),float(c.y));
        };
        auto rayLength=[&](CurvePoint c,CurvePoint n,int owner) {
            double distance=1e30;
            auto visit=[&](const VertexLoop& loop) {
                for(size_t k=0;k<loop.size();++k) {
                    const int a=loop[k],b=loop[(k+1)%loop.size()];
                    const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],den=n.x*dy-n.y*dx;
                    if(std::abs(den)<1e-10)continue;
                    const double x=pts[2*a]-c.x,y=pts[2*a+1]-c.y;
                    const double t=(x*dy-y*dx)/den,u=(x*n.y-y*n.x)/den;
                    if(t>1e-5&&u>=-1e-7&&u<=1+1e-7)distance=std::min(distance,t);
                }
            };
            bool ring=false;
            for(const auto& loop:borderLoops)if(component(pts[2*loop[0]],pts[2*loop[0]+1])==owner) {visit(loop);ring=true;}
            if(!ring)visit(loops[size_t(owner)]);
            return distance;
        };
        std::vector<double> stations;
        const int sections=std::clamp(int(std::ceil(arc.back()/spacing)),2,512);
        // Refine turns by chord error while leaving straight runs sparse.
        auto sample=[&](auto&& self,double a,double b,int level)->void {
            const auto x=at(a),y=at(b);const double dx=y.x-x.x,dy=y.y-x.y,squared=dx*dx+dy*dy;
            double deviation=0;
            for(double f:{.25,.5,.75}) {
                const auto z=at(a+(b-a)*f);const double t=squared?std::clamp(((z.x-x.x)*dx+(z.y-x.y)*dy)/squared,0.,1.):0;
                deviation=std::max(deviation,std::hypot(z.x-x.x-t*dx,z.y-x.y-t*dy));
            }
            if(level<4&&b-a>spacing*.5&&deviation>std::max(cell*1.5,spacing*.08)) {
                self(self,a,(a+b)*.5,level+1);self(self,(a+b)*.5,b,level+1);
            } else stations.push_back((a+b)*.5);
        };
        for(int k=0;k<sections;++k)sample(sample,arc.back()*k/sections,arc.back()*(k+1)/sections,0);
        struct Station {int vertex;double fraction;};
        std::vector<Station> previous;
        std::vector<std::vector<float>> flowPaths;
        std::set<std::pair<int,int>> connected;
        auto connect=[&](int a,int b) {
            if(a==b||connected.count(std::minmax(a,b)))return;
            const CurvePoint x{pts[2*a],pts[2*a+1]},y{pts[2*b],pts[2*b+1]};
            const int probes=std::max(1,int(std::ceil(std::hypot(y.x-x.x,y.y-x.y)/std::max(1.,cell*.5))));
            for(int k=0;k<=probes;++k)if(!onTexture({x.x+(y.x-x.x)*k/probes,x.y+(y.y-x.y)*k/probes}))return;
            connected.insert(std::minmax(a,b));
            flowPaths.push_back({float(x.x),float(x.y),float(y.x),float(y.y)});
        };
        for(double s:stations) {
            const auto c=at(s),a=at(s-std::min(spacing*.25,arc.back()*.05)),b=at(s+std::min(spacing*.25,arc.back()*.05));
            const double length=std::hypot(b.x-a.x,b.y-a.y);const int owner=component(float(c.x),float(c.y));
            if(length<1e-5||owner<0||!onTexture(c)){previous.clear();continue;}
            const CurvePoint n{-(b.y-a.y)/length,(b.x-a.x)/length};
            const double left=rayLength(c,{-n.x,-n.y},owner),right=rayLength(c,n,owner),width=left+right;
            if(width>std::max(domain.width,domain.height)*2){previous.clear();continue;}
            const int columns=std::clamp(int(std::lround(width/spacing))-1,0,32);
            std::vector<Station> current;
            for(int k=1;k<=columns;++k) {
                const double fraction=double(k)/(columns+1),offset=-left+width*fraction;
                const float x=float(c.x+n.x*offset),y=float(c.y+n.y*offset);
                if(!onTexture({x,y})||distanceToFeatures(x,y)<spacing*.35||distanceToBorder(x,y)<borderClearance(spacing))continue;
                bool separated=true;
                for(int v:anchors)if(std::hypot(x-pts[2*v],y-pts[2*v+1])<spacing*.5){separated=false;break;}
                if(!separated)continue;
                const int v=addPoint(x,y,false);anchors.insert(v);current.push_back({v,fraction});
            }
            for(size_t i=1;i<current.size();++i)connect(current[i-1].vertex,current[i].vertex);
            // Ordered nearest-column matches avoid crisscrossing when width
            // changes and one interior column appears or disappears.
            auto join=[&](const std::vector<Station>& from,const std::vector<Station>& to) {
                if(to.empty())return;
                for(const auto& v:from) {
                    const auto q=std::min_element(to.begin(),to.end(),[&](const auto& x,const auto& y){return std::abs(x.fraction-v.fraction)<std::abs(y.fraction-v.fraction);});
                    connect(v.vertex,q->vertex);
                }
            };
            join(previous,current);join(current,previous);
            if(!current.empty())++result.crossSections;
            previous=std::move(current);
        }
        if(!flowPaths.empty()) {
            std::vector<std::vector<float>> creasePaths,ringPaths;
            for(auto [a,b]:featureEdges)creasePaths.push_back({pts[2*a],pts[2*a+1],pts[2*b],pts[2*b+1]});
            for(auto [a,b]:borderEdges)ringPaths.push_back({pts[2*a],pts[2*a+1],pts[2*b],pts[2*b+1]});
            auto allPaths=creasePaths;allPaths.insert(allPaths.end(),ringPaths.begin(),ringPaths.end());allPaths.insert(allPaths.end(),flowPaths.begin(),flowPaths.end());
            std::vector<float> graph;std::vector<int> mapping;
            for(int v:anchors){graph.push_back(pts[2*v]);graph.push_back(pts[2*v+1]);mapping.push_back(v);}
            const size_t existing=mapping.size();std::vector<std::pair<int,int>> all,cuts,ringCuts,flowCuts;
            if(!appendMeshFeatureConstraints(graph,allPaths,all,error)||!appendMeshFeatureConstraints(graph,creasePaths,cuts,error)||
               !appendMeshFeatureConstraints(graph,ringPaths,ringCuts,error)||!appendMeshFeatureConstraints(graph,flowPaths,flowCuts,error))return false;
            for(size_t i=existing;i<graph.size()/2;++i)mapping.push_back(addPoint(graph[2*i],graph[2*i+1],false));
            auto mapEdges=[&](const std::vector<std::pair<int,int>>& from,std::vector<std::pair<int,int>>& to) {
                to.clear();for(auto [a,b]:from)to.push_back(std::minmax(mapping[a],mapping[b]));
            };
            mapEdges(cuts,featureEdges);mapEdges(ringCuts,borderEdges);mapEdges(flowCuts,flowEdges);
            for(auto [a,b]:all){anchors.insert(mapping[a]);anchors.insert(mapping[b]);}
        }
    }
    std::set<int> discarded;
    const bool explicitFeatures=features&&std::any_of(features->begin(),features->end(),[](const auto& line){return line.enabled&&!line.automatic;});
    if(!symmetry.detected&&!explicitFeatures&&(!featureEdges.empty()||!flowEdges.empty())) {
        const double cell=std::max(1,(std::max(domain.width,domain.height)+511)/512);
        const double ratio=p.polygonBorderDensity>0?double(p.polygonBorderDensity)/std::max(1,p.polygonOuterDensity):1.;
        const double separation=std::max(cell*1.5,std::min(nominalSpacing*.15,double(medianBoundary)*.45*ratio));
        const double tolerance=std::max(cell*.75,std::min(cell*2.,separation*.25));
        std::set<int> fixed=protectedOld;for(const auto& loop:loops)fixed.insert(loop.begin(),loop.end());
        for(auto [a,b]:borderEdges){fixed.insert(a);fixed.insert(b);}
        for(auto [a,b]:featureEdges){fixed.insert(a);fixed.insert(b);}
        // Simplify only a degree-two station along the same constraint chain.
        // Junctions, texture-line anchors, both border rows and glue identities
        // stay fixed. Only generated flow chains may simplify here.
        for(int pass=0;pass<1000;++pass) {
            std::set<std::pair<int,int>> all(featureEdges.begin(),featureEdges.end());
            all.insert(borderEdges.begin(),borderEdges.end());all.insert(flowEdges.begin(),flowEdges.end());
            std::map<int,std::set<int>> adjacent;std::vector<std::pair<double,std::pair<int,int>>> shortEdges;
            for(auto [a,b]:all) {
                adjacent[a].insert(b);adjacent[b].insert(a);
                const double length=std::hypot(pts[2*a]-pts[2*b],pts[2*a+1]-pts[2*b+1]);
                if(length<separation)shortEdges.push_back({length,{a,b}});
            }
            std::sort(shortEdges.begin(),shortEdges.end());bool changed=false;
            for(const auto& candidate:shortEdges) {
                for(int direction=0;direction<2&&!changed;++direction) {
                    const int keep=direction?candidate.second.second:candidate.second.first;
                    const int remove=direction?candidate.second.first:candidate.second.second;
                    if(fixed.count(remove)||adjacent[remove].size()!=2)continue;
                    int next=-1;for(int v:adjacent[remove])if(v!=keep)next=v;
                    if(next<0||next==keep)continue;
                    const CurvePoint a{pts[2*keep],pts[2*keep+1]},b{pts[2*next],pts[2*next+1]},q{pts[2*remove],pts[2*remove+1]};
                    const double dx=b.x-a.x,dy=b.y-a.y,squared=dx*dx+dy*dy;
                    if(squared<1e-10)continue;
                    const double t=std::clamp(((q.x-a.x)*dx+(q.y-a.y)*dy)/squared,0.,1.);
                    if(std::hypot(q.x-a.x-t*dx,q.y-a.y-t*dy)>tolerance)continue;
                    bool valid=true;
                    for(auto [c,d]:all) {
                        if(c==remove||d==remove||c==keep||d==keep||c==next||d==next)continue;
                        const CurvePoint x{pts[2*c],pts[2*c+1]},y{pts[2*d],pts[2*d+1]};
                        if(curveCross(a,b,x)*curveCross(a,b,y)<-1e-10&&curveCross(x,y,a)*curveCross(x,y,b)<-1e-10){valid=false;break;}
                    }
                    if(!valid)continue;
                    for(const auto& [v,neighbors]:adjacent)if(v!=keep&&v!=next&&v!=remove) {
                        const double x=pts[2*v]-a.x,y=pts[2*v+1]-a.y,f=(x*dx+y*dy)/squared;
                        if(f>1e-7&&f<1-1e-7&&std::abs(x*dy-y*dx)<boundaryTolerance*std::sqrt(squared)){valid=false;break;}
                    }
                    if(!valid)continue;
                    const int probes=std::max(1,int(std::ceil(std::sqrt(squared)/.75)));
                    const int owner=component(float(a.x),float(a.y));
                    for(int k=0;k<=probes&&valid;++k) {
                        const float x=float(a.x+dx*k/probes),y=float(a.y+dy*k/probes);
                        valid=component(x,y)==owner&&!onBoundary(x,y,boundaryTolerance);
                    }
                    const bool onRing=std::any_of(borderEdges.begin(),borderEdges.end(),[&](auto e){return e.first==remove||e.second==remove;});
                    if(onRing) {
                        for(int k=0;k<=probes&&valid;++k) {
                            const int x=int(std::floor(a.x+dx*k/probes)),y=int(std::floor(a.y+dy*k/probes));
                            valid=x>=0&&x<domain.width&&y>=0&&y<domain.height&&domain.rgba[(size_t(y)*domain.width+x)*4+3];
                        }
                    }
                    auto proposal=borderLoops;
                    for(auto& loop:proposal) {
                        for(int& v:loop)if(v==remove)v=keep;
                        loop.erase(std::unique(loop.begin(),loop.end()),loop.end());
                        if(loop.size()>1&&loop.front()==loop.back())loop.pop_back();
                        valid&=loop.size()>=3&&simpleOutline(pts,loop,1e-6);
                    }
                    if(!valid)continue;
                    auto replace=[&](std::vector<std::pair<int,int>>& edges) {
                        std::set<std::pair<int,int>> unique;
                        for(auto [u,v]:edges){if(u==remove)u=keep;if(v==remove)v=keep;if(u!=v)unique.insert(std::minmax(u,v));}
                        edges.assign(unique.begin(),unique.end());
                    };
                    replace(featureEdges);replace(borderEdges);replace(flowEdges);borderLoops=std::move(proposal);
                    discarded.insert(remove);anchors.erase(remove);changed=true;++result.mergedConstraintVertices;
                }
                if(changed)break;
            }
            if(!changed)break;
        }
    }
    result.flowVertices=0;for(int v=beforeFlow;v<int(pts.size()/2);++v)if(!discarded.count(v))++result.flowVertices;
    result.flowEdges=int(flowEdges.size());
    struct Candidate{float x,y;double spacing,distance;bool reference=false;};std::vector<Candidate> candidates;
    const float grid=float(std::max(nominalSpacing*.45,std::sqrt(double(w)*h/16000)));
    int added=0;for(int v=beforeBorder;v<int(pts.size()/2);++v)if(!discarded.count(v))++added;
    std::vector<int> supportPoints;
    int row=0;
    if(target&&!strip&&p.style!=MeshAutoGenStyle::Rigid&&p.style!=MeshAutoGenStyle::FineLine)for(float y=minY+grid*.43f;y<maxY;y+=grid*.8660254f,++row)
        for(float x=minX+grid*(row%2?.75f:.25f);x<maxX;x+=grid) {
            if(component(x,y)<0||!insideBorder(x,y))continue;
            const double spacing=localSpacing(x,y);
            if(symmetry.detected&&symmetry.local(x,y).x<spacing*.4)continue;
            if(onBoundary(x,y,float(std::max(double(minimumMargin),std::min(double(preferredMargin),spacing*.55))))||
               distanceToFeatures(x,y)<spacing*.4||distanceToBorder(x,y)<borderClearance(spacing))continue;
            double nearest=std::numeric_limits<double>::max();
            for(int v:anchors) {
                nearest=std::min(nearest,std::pow(x-pts[2*v],2)+std::pow(y-pts[2*v+1],2));
                if(symmetry.detected){const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);nearest=std::min(nearest,std::pow(x-q.x,2)+std::pow(y-q.y,2));}
            }
            candidates.push_back({x,y,spacing,nearest});
        }
    if(symmetry.detected&&target&&p.style!=MeshAutoGenStyle::Rigid&&p.style!=MeshAutoGenStyle::FineLine) {
        double first=1e30,last=-1e30;
        for(const auto& loop:loops)for(int v:loop){const auto q=symmetry.local(pts[2*v],pts[2*v+1]);first=std::min(first,q.y);last=std::max(last,q.y);}
        const int sections=std::max(1,int(std::ceil((last-first)/nominalSpacing)));
        for(int k=0;k<sections;++k) {
            const auto q=symmetry.world(0,first+(last-first)*(k+.5)/sections);const float x=float(q.x),y=float(q.y);
            if(component(x,y)<0||!insideBorder(x,y)||onBoundary(x,y,minimumMargin)||distanceToBorder(x,y)<borderClearance(nominalSpacing)||distanceToFeatures(x,y)<nominalSpacing*.4)continue;
            double nearest=1e30;for(int v:anchors)nearest=std::min(nearest,std::pow(x-pts[2*v],2)+std::pow(y-pts[2*v+1],2));
            candidates.push_back({x,y,nominalSpacing,nearest});
        }
    }
    if(target)for(const auto& hint:weakHints){const float x=float(hint.x),y=float(hint.y);const double spacing=localSpacing(x,y);
        if(component(x,y)<0||!insideBorder(x,y)||onBoundary(x,y,minimumMargin)||distanceToBorder(x,y)<borderClearance(spacing)||distanceToFeatures(x,y)<spacing*.4)continue;
        if(symmetry.detected&&symmetry.local(x,y).x<spacing*.4)continue;
        double nearest=1e30;for(int v:anchors)nearest=std::min(nearest,std::pow(x-pts[2*v],2)+std::pow(y-pts[2*v+1],2));
        for(int v:supportPoints)nearest=std::min(nearest,std::pow(x-pts[2*v],2)+std::pow(y-pts[2*v+1],2));
        candidates.push_back({x,y,spacing,nearest,true});
    }
    int fillBudget=target;
    if(p.interiorTriangleTarget>0){std::set<int> outerSet,ringSet;
        for(const auto& loop:loops)outerSet.insert(loop.begin(),loop.end());
        for(auto [a,b]:borderEdges){ringSet.insert(a);ringSet.insert(b);}
        int fixedInside=0;for(int v:anchors)if(!outerSet.count(v)&&insideBorder(pts[2*v],pts[2*v+1]))++fixedInside;
        const int fixedFaces=std::max(0,2*fixedInside-int(ringSet.size())-2);
        fillBudget=std::max(0,(p.interiorTriangleTarget-fixedFaces)/2);
        if(symmetry.detected)fillBudget/=2;
    }
    const int bulkLimit=added+fillBudget;
    while(added<bulkLimit&&!candidates.empty()) {
        auto best=std::max_element(candidates.begin(),candidates.end(),[](const auto& a,const auto& b){auto score=[](const auto& c){return c.distance<c.spacing*c.spacing*.64?-1.:c.distance/(c.spacing*c.spacing)*(c.reference?1.15:1.);};return score(a)<score(b);});
        if(best->distance<best->spacing*best->spacing*.64)break;
        // Free vertices must be mapped again after relaxation. Reusing an old
        // UID here would leave its serialized UV/form coordinates unmoved.
        const float x=best->x,y=best->y;supportPoints.push_back(addPoint(x,y,false));++added;
        for(auto& c:candidates)c.distance=std::min(c.distance,std::pow(c.x-x,2)+std::pow(c.y-y,2));
    }
    if(symmetry.detected) {
        // Complete pairs for genuine crease anchors as well as automatic fill.
        // The reflected point is support geometry, not a fabricated RGB crease.
        // Keep explicit constraints and glue identities at their original positions.
        const size_t count=pts.size()/2;
        for(size_t v=0;v<count;++v) {
            if(discarded.count(int(v))||(v<size_t(nOld)&&!anchors.count(int(v))))continue;
            const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);
            if(component(float(q.x),float(q.y))<0&&!onBoundary(float(q.x),float(q.y),.002f)) {error="对称配对点超出轮廓";return false;}
            bool found=false;
            for(size_t k=0;k<pts.size()/2;++k)if(!discarded.count(int(k))&&(k>=size_t(nOld)||anchors.count(int(k)))&&std::hypot(pts[2*k]-q.x,pts[2*k+1]-q.y)<.002){found=true;break;}
            if(!found){addPoint(float(q.x),float(q.y),false);++added;}
        }
        // Split existing constraints at any newly introduced collinear mate.
        auto split=[&](std::vector<std::pair<int,int>>& edges) {
            std::vector<std::pair<int,int>> out;
            for(auto [a,b]:edges) {
                const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],squared=dx*dx+dy*dy;
                std::vector<std::pair<double,int>> stations{{0,a},{1,b}};
                for(size_t k=0;squared>0&&k<pts.size()/2;++k)if(int(k)!=a&&int(k)!=b&&!discarded.count(int(k))&&(k>=size_t(nOld)||anchors.count(int(k)))) {
                    const double x=pts[2*k]-pts[2*a],y=pts[2*k+1]-pts[2*a+1],t=(x*dx+y*dy)/squared;
                    if(t>1e-6&&t<1-1e-6&&std::abs(x*dy-y*dx)<1e-5*std::sqrt(squared))stations.emplace_back(t,int(k));
                }
                std::sort(stations.begin(),stations.end());for(size_t i=1;i<stations.size();++i)out.push_back(std::minmax(stations[i-1].second,stations[i].second));
            }
            edges=std::move(out);
        };
        split(featureEdges);split(borderEdges);split(flowEdges);
    }
    // Triangulate every outer component separately with polygon-first, so no
    // centroid-only clipping can silently drop a small island or narrow tip.
    MeshTopologyResult tri;
    auto triangulate=[&]() {
    tri={};
    for(size_t part=0;part<outerCount;++part) {
        MeshTopologyInput in;std::vector<int> activeToOriginal,originalToActive(pts.size()/2,-1);
        auto activate=[&](int original) {
            if(originalToActive[size_t(original)]>=0)return originalToActive[size_t(original)];
            const int index=int(activeToOriginal.size());originalToActive[size_t(original)]=index;activeToOriginal.push_back(original);
            in.pts.push_back(pts[2*original]);in.pts.push_back(pts[2*original+1]);return index;
        };
        VertexLoop activeLoop;
        for(int v:loops[part])activeLoop.push_back(activate(v));
        std::vector<VertexLoop> activeHoles;
        for(size_t i=0;i<holeOwners.size();++i)if(holeOwners[i]==int(part)){VertexLoop hole;
            for(int v:loops[outerCount+i])hole.push_back(activate(v));activeHoles.push_back(std::move(hole));}
        for(auto [a,b]:featureEdges)if(component(pts[2*a],pts[2*a+1])==int(part)){activate(a);activate(b);}
        for(int v=0;v<int(pts.size()/2);++v) {
            if(discarded.count(v)||(v<nOld&&!protectedOld.count(v))||originalToActive[size_t(v)]>=0)continue;
            if(component(pts[2*v],pts[2*v+1])==int(part)&&!onBoundary(pts[2*v],pts[2*v+1],boundaryTolerance))activate(v);
        }
        for(size_t i=0;i<activeLoop.size();++i)in.fixedEdges.push_back(std::minmax(activeLoop[i],activeLoop[(i+1)%activeLoop.size()]));
        for(const auto& hole:activeHoles)for(size_t i=0;i<hole.size();++i)in.fixedEdges.push_back(std::minmax(hole[i],hole[(i+1)%hole.size()]));
        const std::set<std::pair<int,int>> expectedBoundary(in.fixedEdges.begin(),in.fixedEdges.end());
        for(auto [a,b]:featureEdges)if(component(pts[2*a],pts[2*a+1])==int(part)) {
            if(originalToActive[size_t(a)]<0||originalToActive[size_t(b)]<0){error="结构线点未纳入三角化";return false;}
            in.fixedEdges.emplace_back(originalToActive[size_t(a)],originalToActive[size_t(b)]);
        }
        for(auto [a,b]:borderEdges)if(component(pts[2*a],pts[2*a+1])==int(part)) {
            if(originalToActive[size_t(a)]<0||originalToActive[size_t(b)]<0){error="锁边环点未纳入三角化";return false;}
            in.fixedEdges.emplace_back(originalToActive[size_t(a)],originalToActive[size_t(b)]);
        }
        for(auto [a,b]:flowEdges)if(component(pts[2*a],pts[2*a+1])==int(part)) {
            if(originalToActive[size_t(a)]<0||originalToActive[size_t(b)]<0){error="走向结构点未纳入三角化";return false;}
            in.fixedEdges.emplace_back(originalToActive[size_t(a)],originalToActive[size_t(b)]);
        }
        in.outlines={activeLoop};auto generated=rebuildMeshTriangles(in);
        if(!generated.ok||generated.tris.empty()){error="约束三角化失败："+std::string(generated.note);return false;}
        // Hole edges were recovered as hard constraints BEFORE classification,
        // so no retained face can straddle a hole. Boundary/area checks below
        // reject incomplete recovery instead of silently filling the opening.
        if(!activeHoles.empty()){std::vector<int> kept;
            for(size_t t=0;t<generated.tris.size();t+=3){const int a=generated.tris[t],b=generated.tris[t+1],c=generated.tris[t+2];
                const float x=(in.pts[2*a]+in.pts[2*b]+in.pts[2*c])/3,y=(in.pts[2*a+1]+in.pts[2*b+1]+in.pts[2*c+1])/3;
                bool empty=false;for(const auto& hole:activeHoles)empty|=pointInLoop(in.pts,hole,x,y);
                if(!empty)kept.insert(kept.end(),{a,b,c});}
            generated.tris=std::move(kept);}
        int nonManifold=0;const auto boundaries=boundaryLoops(generated.tris,int(activeToOriginal.size()),&nonManifold);
        std::set<std::pair<int,int>> actualBoundary;
        for(const auto& boundary:boundaries)for(size_t i=0;i<boundary.size();++i)
            actualBoundary.insert(std::minmax(boundary[i],boundary[(i+1)%boundary.size()]));
        double area=std::abs(loopArea(in.pts,activeLoop));for(const auto& hole:activeHoles)area-=std::abs(loopArea(in.pts,hole));
        if(nonManifold||boundaries.size()!=1+activeHoles.size()||actualBoundary!=expectedBoundary||std::abs(std::abs(trianglesArea(in.pts,generated.tris))-area)>area*1e-5) {
            error="生成三角形未完整覆盖纹理轮廓，拒绝写回";return false;
        }
        for(int index:generated.tris)tri.tris.push_back(activeToOriginal[size_t(index)]);
    }
    return true;
    };
    const bool trace=std::getenv("L2DR_AUTOMESH_TRACE")!=nullptr;
    if(trace)std::fprintf(stderr,"[hole-stage] triangulate points=%zu holes=%zu\n",pts.size()/2,holeOwners.size());
    if(!triangulate())return false;
    if(trace)std::fprintf(stderr,"[hole-stage] first triangulation done\n");
    // Collapse short mesh edges by role, never by averaging positions:
    // outer/glue/texture-line identities > support ring > flow > fill.
    // Ring stations may collapse, but the closed inset chain must survive.
    // Rebuild constraints and triangulation before
    // accepting a collapse; coupled mirror collapses preserve paired sampling.
    {
        std::set<int> immutable=protectedOld,ringVertices,structureVertices,outerVertices;
        for(const auto& loop:loops)immutable.insert(loop.begin(),loop.end());
        for(const auto& loop:loops)outerVertices.insert(loop.begin(),loop.end());
        for(auto [a,b]:borderEdges){ringVertices.insert(a);ringVertices.insert(b);}
        for(auto [a,b]:featureEdges){structureVertices.insert(a);structureVertices.insert(b);
            // Original crease stations remain immutable. In a hole mesh the
            // ring splitter also inserts derived collinear crease stations;
            // those may collapse along their unchanged crease, checked below.
            if(holeOwners.empty()||a<beforeBorder)immutable.insert(a);
            if(holeOwners.empty()||b<beforeBorder)immutable.insert(b);}
        for(auto [a,b]:flowEdges){structureVertices.insert(a);structureVertices.insert(b);}
        // Feature/flow intersections split borderEdges after borderLoops was
        // fitted. Include those shared stations before editing the ring chain.
        for(auto& loop:borderLoops){VertexLoop expanded;
            for(size_t i=0;i<loop.size();++i){const int a=loop[i],b=loop[(i+1)%loop.size()];
                const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1],squared=dx*dx+dy*dy;
                std::vector<std::pair<double,int>> stations{{0,a}};
                for(int v:ringVertices)if(v!=a&&v!=b&&squared>0){const double x=pts[2*v]-pts[2*a],y=pts[2*v+1]-pts[2*a+1],t=(x*dx+y*dy)/squared;
                    if(t>1e-7&&t<1-1e-7&&std::abs(x*dy-y*dx)<1e-5*std::sqrt(squared))stations.emplace_back(t,v);}
                std::sort(stations.begin(),stations.end());for(auto [t,v]:stations)expanded.push_back(v);
            }loop=std::move(expanded);}
        auto priority=[&](int v){return immutable.count(v)?3:ringVertices.count(v)?2:structureVertices.count(v)?1:0;};
        // A sparse interior must not coarsen an independently tuned lock row.
        // Compare its short edges to that row's own ordinary edge length.
        std::map<int,double> ringSpacing;
        for(const auto& loop:borderLoops){std::vector<double> lengths;
            for(size_t i=0;i<loop.size();++i){const int a=loop[i],b=loop[(i+1)%loop.size()];lengths.push_back(std::hypot(pts[2*b]-pts[2*a],pts[2*b+1]-pts[2*a+1]));}
            std::sort(lengths.begin(),lengths.end());for(int v:loop)ringSpacing[v]=lengths[lengths.size()/2];}
        auto mergeSpacing=[&](int v){return ringSpacing.count(v)?ringSpacing.at(v):localSpacing(pts[2*v],pts[2*v+1]);};
        // Evaluate the thin lock band separately: it must not hide an interior
        // sliver regression. Equal quality is enough when a redundant point goes.
        struct MergeScore {double minimum=1,loss=0;int severe=0;};
        auto mergeScore=[&](){std::array<MergeScore,2> scores;
            for(size_t t=0;t<tri.tris.size();t+=3){const int a=tri.tris[t],b=tri.tris[t+1],c=tri.tris[t+2];
                auto& score=scores[outerVertices.count(a)||outerVertices.count(b)||outerVertices.count(c)?0:1];
                const double ux=pts[2*b]-pts[2*a],uy=pts[2*b+1]-pts[2*a+1],vx=pts[2*c]-pts[2*a],vy=pts[2*c+1]-pts[2*a+1];
                const double q=3.464101615137754*std::abs(ux*vy-uy*vx)/std::max(1e-12,ux*ux+uy*uy+vx*vx+vy*vy+(ux-vx)*(ux-vx)+(uy-vy)*(uy-vy));
                const double deficit=std::max(0.,.35-q);score.minimum=std::min(score.minimum,q);score.loss+=deficit*deficit;score.severe+=q<.2;
            }return scores;};
        auto mirrorOf=[&](int v){const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);
            for(int k:std::set<int>(tri.tris.begin(),tri.tris.end()))if(std::hypot(pts[2*k]-q.x,pts[2*k+1]-q.y)<.005)return k;
            return -1;};
        // Many-hole meshes have far more hard constraints. Bound rejected
        // retriangulation attempts as well as successful collapses.
        int mergeTrials=0;const int mergeLimit=holeOwners.empty()?std::numeric_limits<int>::max():std::clamp(16000/int(pts.size()/2),12,64);
        for(int pass=0;pass<128&&mergeTrials<mergeLimit;++pass) {
            std::set<std::pair<int,int>> edges;
            for(size_t t=0;t<tri.tris.size();t+=3)for(int k=0;k<3;++k)edges.insert(std::minmax(tri.tris[t+k],tri.tris[t+(k+1)%3]));
            std::vector<std::pair<double,std::pair<int,int>>> close;
            for(auto [a,b]:edges){const double d=std::hypot(pts[2*a]-pts[2*b],pts[2*a+1]-pts[2*b+1]);
                const double limit=std::max(holeOwners.empty()?0.:1.,std::min(mergeSpacing(a),mergeSpacing(b))*.15);
                if(d<limit)close.push_back({d,{a,b}});}
            std::sort(close.begin(),close.end());bool changed=false;
            for(const auto& candidate:close) {
                auto [keep,remove]=candidate.second;
                if(priority(remove)>priority(keep))std::swap(keep,remove);
                if(immutable.count(remove))continue;
                if(mergeTrials>=mergeLimit)break;
                std::map<int,int> replacement{{remove,keep}};
                if(symmetry.detected){const int mr=mirrorOf(remove),mk=mirrorOf(keep);
                    if(mr<0||mk<0||mr==keep||mk==remove||immutable.count(mr)||priority(mk)<priority(mr))continue;
                    if(mr==remove&&mk!=keep)continue;
                    replacement[mr]=mk;}
                // Original crease endpoints/bends stay fixed; a generated
                // intersection may only collapse to its collinear neighbor.
                // This preserves the geometric union of all crease segments.
                bool valid=true;
                for(auto [from,to]:replacement){std::set<int> neighbors;
                    for(auto [a,b]:featureEdges){if(a==from)neighbors.insert(b);if(b==from)neighbors.insert(a);}
                    if(neighbors.empty())continue;
                    if(neighbors.size()!=2||!neighbors.count(to)){valid=false;break;}
                    const int other=*neighbors.begin()==to?*neighbors.rbegin():*neighbors.begin();
                    const double dx=pts[2*other]-pts[2*to],dy=pts[2*other+1]-pts[2*to+1],squared=dx*dx+dy*dy;
                    const double x=pts[2*from]-pts[2*to],y=pts[2*from+1]-pts[2*to+1],t=squared?(x*dx+y*dy)/squared:-1;
                    if(t<=0||t>=1||std::abs(x*dy-y*dx)>1e-5*std::sqrt(squared)){valid=false;break;}
                }
                // Hole and outer support rows may share a nearby station in a
                // narrow bridge, but neither may absorb a true boundary vertex.
                for(auto [from,to]:replacement)if(ringVertices.count(from)){
                    if(outerVertices.count(to)){valid=false;break;}
                    if(ringVertices.count(to)&&std::find(borderEdges.begin(),borderEdges.end(),std::minmax(from,to))==borderEdges.end()&&
                       (holeOwners.empty()||candidate.first>1.)){valid=false;break;}
                }
                if(!valid)continue;
                auto proposal=borderLoops;
                for(auto& loop:proposal){
                    for(int& v:loop)if(replacement.count(v))v=replacement.at(v);
                    loop.erase(std::unique(loop.begin(),loop.end()),loop.end());
                    if(loop.size()>1&&loop.front()==loop.back())loop.pop_back();
                    if(loop.size()<3||!simpleOutline(pts,loop,1e-6)){valid=false;break;}
                }
                size_t holeIndex=0;
                for(size_t i=0;i<proposal.size()&&valid;++i)if(borderIsHole[i]){
                    const auto& hole=domain.holes[holeIndex++];
                    for(size_t j=0;j<hole.size()&&valid;j+=2)valid=pointInLoop(pts,proposal[i],hole[j],hole[j+1]);
                }
                if(!valid)continue;
                const auto oldFeatures=featureEdges,oldFlow=flowEdges,oldBorder=borderEdges;
                const auto oldTri=tri;
                const auto beforeScore=mergeScore();
                auto replace=[&](std::vector<std::pair<int,int>>& source){std::set<std::pair<int,int>> unique;
                    for(auto [a,b]:source){if(replacement.count(a))a=replacement.at(a);if(replacement.count(b))b=replacement.at(b);if(a!=b)unique.insert(std::minmax(a,b));}
                    source.assign(unique.begin(),unique.end());};
                replace(featureEdges);replace(flowEdges);replace(borderEdges);
                // Check only changed ring chords against actual alpha. This
                // prevents a shortcut cutting across a transparent concavity.
                for(auto [a,b]:borderEdges)if(std::find(oldBorder.begin(),oldBorder.end(),std::pair<int,int>{a,b})==oldBorder.end()){
                    const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1];
                    const int probes=std::max(1,int(std::ceil(std::hypot(dx,dy)*2)));
                    for(int k=0;k<=probes;++k){const int x=int(std::floor(pts[2*a]+dx*k/probes)),y=int(std::floor(pts[2*a+1]+dy*k/probes));
                        if(x<0||y<0||x>=domain.width||y>=domain.height||!domain.rgba[(size_t(y)*domain.width+x)*4+3]){valid=false;break;}}
                    if(!valid)break;
                }
                for(auto [from,to]:replacement)discarded.insert(from);
                // Reject crossings rather than silently changing structural paths.
                auto constraints=featureEdges;constraints.insert(constraints.end(),flowEdges.begin(),flowEdges.end());constraints.insert(constraints.end(),borderEdges.begin(),borderEdges.end());
                for(size_t a=0;a<constraints.size()&&valid;++a)for(size_t b=a+1;b<constraints.size();++b){
                    const auto [u,v]=constraints[a];const auto [x,y]=constraints[b];if(u==x||u==y||v==x||v==y)continue;
                    const CurvePoint p0{pts[2*u],pts[2*u+1]},p1{pts[2*v],pts[2*v+1]},q0{pts[2*x],pts[2*x+1]},q1{pts[2*y],pts[2*y+1]};
                    if(curveCross(p0,p1,q0)*curveCross(p0,p1,q1)<-1e-10&&curveCross(q0,q1,p0)*curveCross(q0,q1,p1)<-1e-10){valid=false;break;}}
                // The budget bounds expensive retriangulations, not cheap
                // rejections of protected/crossing candidates ahead of them.
                if(valid){++mergeTrials;valid=triangulate();}
                if(valid){const auto afterScore=mergeScore();for(size_t i=0;i<beforeScore.size();++i)
                    valid&=afterScore[i].minimum+1e-6>=beforeScore[i].minimum&&afterScore[i].severe<=beforeScore[i].severe&&afterScore[i].loss<=beforeScore[i].loss+1e-5;}
                if(!valid){
                    featureEdges=oldFeatures;flowEdges=oldFlow;borderEdges=oldBorder;tri=oldTri;
                    for(auto [from,to]:replacement)discarded.erase(from);
                    error.clear();continue;
                }
                borderLoops=std::move(proposal);
                for(auto [from,to]:replacement)if(std::any_of(oldFeatures.begin(),oldFeatures.end(),[&](auto edge){return edge.first==from||edge.second==from;}))
                    ++result.mergedDerivedFeatureVertices;
                for(auto [from,to]:replacement){if(anchors.erase(from))anchors.insert(to);structureVertices.erase(from);if(from>=beforeBorder)--added;
                    if(ringVertices.erase(from)){ringVertices.insert(to);if(!ringSpacing.count(to))ringSpacing[to]=mergeSpacing(from);ringSpacing.erase(from);++result.mergedBorderVertices;}}
                result.priorityMergedVertices+=int(replacement.size());changed=true;break;
            }
            if(!changed)break;
        }
        supportPoints.erase(std::remove_if(supportPoints.begin(),supportPoints.end(),[&](int v){return discarded.count(v)||anchors.count(v);}),supportPoints.end());
        result.flowEdges=int(flowEdges.size());
        result.flowVertices=0;for(int v=beforeFlow;v<int(pts.size()/2);++v)if(!discarded.count(v)&&structureVertices.count(v))++result.flowVertices;
    }
    if(trace)std::fprintf(stderr,"[hole-stage] merges done\n");
    // Balance adjacent edge lengths using a common local rest length. Only
    // free support vertices move; boundary, border ring, crease and glue anchors stay put.
    // Reject moves that crowd another vertex or cross a structural boundary.
    for(int pass=0;!symmetry.detected&&pass<8&&!supportPoints.empty();++pass) {
        const auto before=pts;const auto previousTri=tri;
        std::vector<std::vector<int>> adjacent(pts.size()/2);
        for(size_t t=0;t<tri.tris.size();t+=3)for(int k=0;k<3;++k) {
            const int a=tri.tris[t+k],b=tri.tris[t+(k+1)%3];adjacent[size_t(a)].push_back(b);adjacent[size_t(b)].push_back(a);
        }
        for(int v:supportPoints) {
            auto& neighbors=adjacent[size_t(v)];std::sort(neighbors.begin(),neighbors.end());neighbors.erase(std::unique(neighbors.begin(),neighbors.end()),neighbors.end());
            if(neighbors.empty())continue;
            double fx=0,fy=0;
            for(int n:neighbors) {
                const double dx=before[2*n]-before[2*v],dy=before[2*n+1]-before[2*v+1],length=std::hypot(dx,dy);
                if(length<1e-8)continue;
                const double rest=std::min(localSpacing((before[2*v]+before[2*n])*.5f,(before[2*v+1]+before[2*n+1])*.5f),double(medianBoundary)*1.25);
                const double force=(length-rest)/length;
                fx+=force*dx;fy+=force*dy;
            }
            const float x=float(before[2*v]+fx*.6/neighbors.size()),y=float(before[2*v+1]+fy*.6/neighbors.size());
            const double spacing=std::min(localSpacing(x,y),double(medianBoundary)*1.25);
            bool separated=true;
            for(size_t n=0;n<pts.size()/2&&separated;++n)if(int(n)!=v&&!discarded.count(int(n))&&(int(n)>=nOld||anchors.count(int(n))))
                separated=std::hypot(x-pts[2*n],y-pts[2*n+1])>=spacing*.65;
            if(separated&&insideBorder(x,y)&&component(x,y)==component(before[2*v],before[2*v+1])&&
               !onBoundary(x,y,float(std::max(double(minimumMargin),spacing*.35)))&&distanceToFeatures(x,y)>spacing*.3&&
               distanceToBorder(x,y)>(p.interiorBorderMargin>0?double(p.interiorBorderMargin):std::min(spacing*.3,double(medianBoundary)*.4))) {
                pts[2*v]=x;pts[2*v+1]=y;
            }
        }
        if(!triangulate()){pts=before;tri=previousTri;error.clear();break;}
    }
    // Interior quality is optimized AFTER crease and ring constraints are fixed.
    // Fill oversized cells first, then move only unconstrained vertices toward
    // equilateral opposite-edge targets. Narrow constrained cells may stay skewed.
    if(p.regularInterior&&target&&p.style!=MeshAutoGenStyle::FineLine) {
        int qualityAttempts=0;
        const int qualityLimit=holeOwners.empty()?std::numeric_limits<int>::max():std::clamp(10000/int(pts.size()/2),8,24);
        auto rebuildQuality=[&]{if(qualityAttempts>=qualityLimit)return false;++qualityAttempts;return triangulate();};
        auto spacingAt=[&](float x,float y){return p.interiorTriangleTarget>0?nominalSpacing:std::min(localSpacing(x,y),double(medianBoundary)*1.75);};
        auto quality=[&](int a,int b,int c){const double ux=pts[2*b]-pts[2*a],uy=pts[2*b+1]-pts[2*a+1],vx=pts[2*c]-pts[2*a],vy=pts[2*c+1]-pts[2*a+1];
            return 3.464101615137754*std::abs(ux*vy-uy*vx)/std::max(1e-12,ux*ux+uy*uy+vx*vx+vy*vy+(ux-vx)*(ux-vx)+(uy-vy)*(uy-vy));};
        const int budget=std::min(holeOwners.empty()?256:12,target*8);
        for(int pass=0;pass<budget;++pass) {
            if(p.interiorTriangleTarget>0){std::set<int> outerSet;for(const auto& loop:loops)outerSet.insert(loop.begin(),loop.end());
                int faces=0;for(size_t t=0;t<tri.tris.size();t+=3)if(!outerSet.count(tri.tris[t])&&!outerSet.count(tri.tris[t+1])&&!outerSet.count(tri.tris[t+2]))++faces;
                if(faces+(symmetry.detected?4:2)>p.interiorTriangleTarget)break;}
            double bestScore=0;float bestX=0,bestY=0;
            const std::set<int> active(tri.tris.begin(),tri.tris.end());
            for(size_t t=0;t<tri.tris.size();t+=3){const int a=tri.tris[t],b=tri.tris[t+1],c=tri.tris[t+2];
                const float x=(pts[2*a]+pts[2*b]+pts[2*c])/3,y=(pts[2*a+1]+pts[2*b+1]+pts[2*c+1])/3;
                const double spacing=spacingAt(x,y);
                if(symmetry.detected&&symmetry.local(x,y).x<spacing*.3)continue;
                if(!insideBorder(x,y)||onBoundary(x,y,minimumMargin)||distanceToBorder(x,y)<std::max(double(p.interiorBorderMargin),spacing*.35)||distanceToFeatures(x,y)<spacing*.3)continue;
                const double twice=std::abs((pts[2*b]-pts[2*a])*(pts[2*c+1]-pts[2*a+1])-(pts[2*b+1]-pts[2*a+1])*(pts[2*c]-pts[2*a]));
                const double score=twice/(spacing*spacing);
                if(score<1.3||score<=bestScore)continue;
                bool separated=true;for(int v:active)if(std::hypot(x-pts[2*v],y-pts[2*v+1])<spacing*.6){separated=false;break;}
                if(separated){bestScore=score;bestX=x;bestY=y;}
            }
            if(!bestScore)break;
            const auto previousTri=tri;const size_t previousSize=pts.size();
            addPoint(bestX,bestY,false);
            if(symmetry.detected){const auto q=symmetry.reflect(bestX,bestY);if(!insideBorder(float(q.x),float(q.y))){pts.resize(previousSize);break;}addPoint(float(q.x),float(q.y),false);}
            if(!rebuildQuality()){pts.resize(previousSize);tri=previousTri;error.clear();break;}
            const int count=int((pts.size()-previousSize)/2);added+=count;result.refinedInteriorVertices+=count;
        }
        std::set<int> fixed=protectedOld;
        for(const auto& loop:loops)fixed.insert(loop.begin(),loop.end());
        for(const auto& list:{featureEdges,borderEdges,flowEdges})for(auto [a,b]:list){fixed.insert(a);fixed.insert(b);}
        const std::set<int> active(tri.tris.begin(),tri.tris.end());
        std::map<int,int> mates;
        if(symmetry.detected){for(int v:active){const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);for(int k:active)if(std::hypot(pts[2*k]-q.x,pts[2*k+1]-q.y)<.005){mates[v]=k;break;}}
            const auto originalFixed=fixed;for(int v:originalFixed)if(mates.count(v))fixed.insert(mates[v]);}
        for(int pass=0;pass<16;++pass){bool changed=false;const auto previousPoints=pts;const auto previousTri=tri;
            std::vector<std::vector<size_t>> incident(pts.size()/2);
            for(size_t t=0;t<tri.tris.size();t+=3)for(int k=0;k<3;++k)incident[tri.tris[t+k]].push_back(t);
            for(int v:active){if(fixed.count(v))continue;int mate=v;
                if(symmetry.detected){if(!mates.count(v))continue;mate=mates[v];if(mate<v||fixed.count(mate))continue;}
                double tx=0,ty=0;int count=0;
                for(size_t t:incident[v]){int a=-1,b=-1;for(int k=0;k<3;++k){const int n=tri.tris[t+k];if(n!=v){if(a<0)a=n;else b=n;}}
                    const double dx=pts[2*b]-pts[2*a],dy=pts[2*b+1]-pts[2*a+1];
                    const double side=dx*(pts[2*v+1]-pts[2*a+1])-dy*(pts[2*v]-pts[2*a])>0?1.:-1.;
                    tx+=(pts[2*a]+pts[2*b])*.5-side*dy*.866025403784;ty+=(pts[2*a+1]+pts[2*b+1])*.5+side*dx*.866025403784;++count;}
                if(!count)continue;tx/=count;ty/=count;
                if(symmetry.detected&&mate==v){const auto q=symmetry.local(tx,ty),onAxis=symmetry.world(0,q.y);tx=onAxis.x;ty=onAxis.y;}
                std::set<size_t> affected(incident[v].begin(),incident[v].end());affected.insert(incident[mate].begin(),incident[mate].end());
                double oldMin=1,oldSum=0;std::vector<double> signs;
                for(size_t t:affected){const int a=tri.tris[t],b=tri.tris[t+1],c=tri.tris[t+2];const double q=quality(a,b,c);oldMin=std::min(oldMin,q);oldSum+=q;signs.push_back((pts[2*b]-pts[2*a])*(pts[2*c+1]-pts[2*a+1])-(pts[2*b+1]-pts[2*a+1])*(pts[2*c]-pts[2*a]));}
                const float vx=pts[2*v],vy=pts[2*v+1],mx=pts[2*mate],my=pts[2*mate+1];
                for(double amount:{.5,.25,.125}){
                    pts[2*v]=float(vx+(tx-vx)*amount);pts[2*v+1]=float(vy+(ty-vy)*amount);
                    if(mate!=v){const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);pts[2*mate]=float(q.x);pts[2*mate+1]=float(q.y);}
                    bool valid=true;double newMin=1,newSum=0;size_t index=0;
                    for(int k:std::set<int>{v,mate}){valid&=insideBorder(pts[2*k],pts[2*k+1])&&!onBoundary(pts[2*k],pts[2*k+1],minimumMargin)&&distanceToBorder(pts[2*k],pts[2*k+1])>=p.interiorBorderMargin;
                        for(int n:active)if(n!=k&&std::hypot(pts[2*k]-pts[2*n],pts[2*k+1]-pts[2*n+1])<spacingAt(pts[2*k],pts[2*k+1])*.25){valid=false;break;}}
                    for(size_t t:affected){const int a=tri.tris[t],b=tri.tris[t+1],c=tri.tris[t+2];const double sign=(pts[2*b]-pts[2*a])*(pts[2*c+1]-pts[2*a+1])-(pts[2*b+1]-pts[2*a+1])*(pts[2*c]-pts[2*a]);
                        valid&=sign*signs[index++]>1e-12;const double q=quality(a,b,c);newMin=std::min(newMin,q);newSum+=q;}
                    if(valid&&newMin+1e-6>=oldMin&&newSum>oldSum+1e-4){changed=true;break;}
                    pts[2*v]=vx;pts[2*v+1]=vy;pts[2*mate]=mx;pts[2*mate+1]=my;
                }
            }
            if(!changed)break;
            if(!rebuildQuality()){pts=previousPoints;tri=previousTri;error.clear();break;}
        }
        // Penalize shape, not just area/edge length. Keep the deliberately thin
        // outer lock band out of this objective. Fixed corners may impose acute
        // angles; never erase their constraints just to improve a score.
        std::set<int> outerSet;for(const auto& loop:loops)outerSet.insert(loop.begin(),loop.end());
        struct ShapeScore { double minimum=1,loss=0,maxEdge=0;int severe=0; };
        auto shapeScore=[&](){ShapeScore score;
            for(size_t t=0;t<tri.tris.size();t+=3){const int a=tri.tris[t],b=tri.tris[t+1],c=tri.tris[t+2];
                if(outerSet.count(a)||outerSet.count(b)||outerSet.count(c))continue;
                const double q=quality(a,b,c),deficit=std::max(0.,.35-q);
                score.minimum=std::min(score.minimum,q);score.loss+=deficit*deficit;score.severe+=q<.2;
                for(auto [u,v]:{std::pair<int,int>{a,b},{b,c},{c,a}})score.maxEdge=std::max(score.maxEdge,double(std::hypot(pts[2*u]-pts[2*v],pts[2*u+1]-pts[2*v+1])));
            }return score;};
        auto improvesShape=[](const ShapeScore& before,const ShapeScore& after){
            return after.minimum+1e-6>=before.minimum&&after.severe<=before.severe&&
                after.loss+1e-5<before.loss&&after.maxEdge<=before.maxEdge*1.02+1e-4;};
        // Only free mirrored fill points may weld onto the symmetry axis.
        // Detected texture lines stay exact, even when close to that axis.
        if(symmetry.detected&&!explicitFeatures){
            std::set<int> immutable=protectedOld;immutable.insert(outerSet.begin(),outerSet.end());
            for(const auto& edges:{featureEdges,borderEdges,flowEdges})for(auto [a,b]:edges){immutable.insert(a);immutable.insert(b);}
            for(int v:active){if(discarded.count(v)||immutable.count(v)||!mates.count(v))continue;
                const int mate=mates[v];if(mate<=v||discarded.count(mate)||immutable.count(mate))continue;
                const double separation=std::hypot(pts[2*v]-pts[2*mate],pts[2*v+1]-pts[2*mate+1]);
                if(separation>std::min(2.,spacingAt(pts[2*v],pts[2*v+1])*.15))continue;
                const auto before=shapeScore();const auto oldTri=tri;const auto oldFeatures=featureEdges;
                const float x=pts[2*v],y=pts[2*v+1];const auto q=symmetry.local(x,y),center=symmetry.world(0,q.y);
                pts[2*v]=float(center.x);pts[2*v+1]=float(center.y);discarded.insert(mate);
                std::set<std::pair<int,int>> replacement;
                for(auto [a,b]:featureEdges){if(a==mate)a=v;if(b==mate)b=v;if(a!=b)replacement.insert(std::minmax(a,b));}
                featureEdges.assign(replacement.begin(),replacement.end());
                if(rebuildQuality()&&improvesShape(before,shapeScore())){
                    mates[v]=v;mates.erase(mate);if(mate>=beforeBorder)--added;++result.qualityRemovedVertices;continue;}
                pts[2*v]=x;pts[2*v+1]=y;discarded.erase(mate);featureEdges=oldFeatures;tri=oldTri;error.clear();
            }
        }
        // Try removing redundant free points first. Ring/crease/manual anchors
        // are immutable, and a symmetric pair is accepted or rolled back together.
        int removalTrials=0;
        for(int pass=0;pass<16&&removalTrials<64;++pass){
            std::map<int,double> candidates;
            for(size_t t=0;t<tri.tris.size();t+=3){const double q=quality(tri.tris[t],tri.tris[t+1],tri.tris[t+2]);if(q>=.3)continue;
                for(int k=0;k<3;++k){int v=tri.tris[t+k];if(fixed.count(v))continue;
                    auto it=candidates.find(v);if(it==candidates.end())candidates[v]=q;else it->second=std::min(it->second,q);}}
            std::vector<std::pair<double,int>> ordered;for(auto [v,q]:candidates)ordered.push_back({q,v});std::sort(ordered.begin(),ordered.end());
            bool accepted=false;const auto before=shapeScore();
            for(auto [q,v]:ordered){if(++removalTrials>64)break;std::set<int> remove{v};
                if(symmetry.detected){if(!mates.count(v)||fixed.count(mates[v]))continue;remove.insert(mates[v]);}
                const auto oldTri=tri;for(int k:remove)discarded.insert(k);
                if(rebuildQuality()&&improvesShape(before,shapeScore())){
                    for(int k:remove)if(k>=beforeBorder)--added;result.qualityRemovedVertices+=int(remove.size());accepted=true;break;}
                for(int k:remove)discarded.erase(k);tri=oldTri;error.clear();
            }
            if(!accepted)break;
        }
        // Area-only refinement misses long thin triangles and cannot shorten
        // fixed crease edges. Split these edges without moving their line.
        for(int pass=0;pass<std::min(128,target*4);++pass){
            if(p.interiorTriangleTarget>0){int faces=0;for(size_t t=0;t<tri.tris.size();t+=3)if(!outerSet.count(tri.tris[t])&&!outerSet.count(tri.tris[t+1])&&!outerSet.count(tri.tris[t+2]))++faces;
                if(faces+(symmetry.detected?4:2)>p.interiorTriangleTarget)break;}
            std::map<std::pair<int,int>,double> edgeQuality;
            std::set<std::pair<int,int>> sliverLongest;
            for(size_t t=0;t<tri.tris.size();t+=3){const double q=quality(tri.tris[t],tri.tris[t+1],tri.tris[t+2]);
                if(q<.3&&!outerSet.count(tri.tris[t])&&!outerSet.count(tri.tris[t+1])&&!outerSet.count(tri.tris[t+2])){
                    double longest=0;std::pair<int,int> edge;
                    for(int k=0;k<3;++k){const int a=tri.tris[t+k],b=tri.tris[t+(k+1)%3];const double length=std::hypot(pts[2*a]-pts[2*b],pts[2*a+1]-pts[2*b+1]);
                        if(length>longest){longest=length;edge=std::minmax(a,b);}}
                    sliverLongest.insert(edge);}
                for(int k=0;k<3;++k){const auto edge=std::minmax(tri.tris[t+k],tri.tris[t+(k+1)%3]);
                    auto it=edgeQuality.find(edge);if(it==edgeQuality.end())edgeQuality[edge]=q;else it->second=std::min(it->second,q);}}
            const std::set<std::pair<int,int>> ringEdges(borderEdges.begin(),borderEdges.end());
            std::set<std::pair<int,int>> structure(featureEdges.begin(),featureEdges.end());structure.insert(flowEdges.begin(),flowEdges.end());
            const std::set<int> active(tri.tris.begin(),tri.tris.end());
            const int qualityBudget=std::min(24,std::max(4,int(active.size()/5)));
            std::vector<std::pair<double,std::pair<int,int>>> candidates;
            for(const auto& [edge,q]:edgeQuality){const auto [a,b]=edge;if(ringEdges.count(edge)||(outerSet.count(a)&&outerSet.count(b)))continue;
                const float x=(pts[2*a]+pts[2*b])*.5f,y=(pts[2*a+1]+pts[2*b+1])*.5f;
                const double ratio=std::hypot(pts[2*a]-pts[2*b],pts[2*a+1]-pts[2*b+1])/spacingAt(x,y);
                const bool repair=sliverLongest.count(edge)&&result.qualityRefinedVertices<qualityBudget;
                if(repair||(ratio>1.65&&(structure.count(edge)||q<.55||ratio>2.)))candidates.push_back({-(repair?3.+(.3-q):ratio),edge});}
            std::sort(candidates.begin(),candidates.end());bool accepted=false;
            for(const auto& candidate:candidates){const auto [a,b]=candidate.second;
                const bool repair=sliverLongest.count(candidate.second)&&result.qualityRefinedVertices<qualityBudget;
                const auto beforeShape=shapeScore();
                const float x=(pts[2*a]+pts[2*b])*.5f,y=(pts[2*a+1]+pts[2*b+1])*.5f;
                std::vector<CurvePoint> insertions{{x,y}};
                if(symmetry.detected){const auto q=symmetry.reflect(x,y);if(std::hypot(q.x-x,q.y-y)>.005)insertions.push_back({q.x,q.y});}
                if(repair&&result.qualityRefinedVertices+int(insertions.size())>qualityBudget)continue;
                bool valid=true;
                for(const auto& q:insertions){const double spacing=spacingAt(float(q.x),float(q.y));
                    valid&=insideBorder(float(q.x),float(q.y))&&!onBoundary(float(q.x),float(q.y),minimumMargin)&&distanceToBorder(float(q.x),float(q.y))>std::max(double(p.interiorBorderMargin),spacing*.1);
                    for(int v:active)if(std::hypot(q.x-pts[2*v],q.y-pts[2*v+1])<std::max(1.,spacing*(repair?.15:.3))){valid=false;break;}}
                if(insertions.size()==2&&std::hypot(insertions[0].x-insertions[1].x,insertions[0].y-insertions[1].y)<spacingAt(x,y)*.3)valid=false;
                if(!valid)continue;
                const size_t oldSize=pts.size();const auto oldTri=tri;const auto oldFeature=featureEdges,oldFlow=flowEdges;
                // A previously collapsed vertex may occupy this location but
                // is still in discarded. Allocate fresh IDs for the split.
                std::vector<int> addedIds;for(const auto& q:insertions){addedIds.push_back(int(pts.size()/2));pts.push_back(float(q.x));pts.push_back(float(q.y));}
                auto split=[&](std::vector<std::pair<int,int>>& edges){std::vector<std::pair<int,int>> resultEdges;
                    for(auto [u,v]:edges){const double dx=pts[2*v]-pts[2*u],dy=pts[2*v+1]-pts[2*u+1],sq=dx*dx+dy*dy;
                        std::vector<std::pair<double,int>> stations{{0,u},{1,v}};
                        for(int k:addedIds){const double qx=pts[2*k]-pts[2*u],qy=pts[2*k+1]-pts[2*u+1],t=sq?(qx*dx+qy*dy)/sq:0;
                            if(t>1e-6&&t<1-1e-6&&std::abs(qx*dy-qy*dx)<.001*std::sqrt(sq))stations.push_back({t,k});}
                        std::sort(stations.begin(),stations.end());for(size_t k=1;k<stations.size();++k)resultEdges.push_back(std::minmax(stations[k-1].second,stations[k].second));}
                    edges=std::move(resultEdges);};
                split(featureEdges);split(flowEdges);
                if(!rebuildQuality()||(repair&&!improvesShape(beforeShape,shapeScore()))){pts.resize(oldSize);tri=oldTri;featureEdges=oldFeature;flowEdges=oldFlow;error.clear();continue;}
                added+=int(addedIds.size());
                if(repair)result.qualityRefinedVertices+=int(addedIds.size());else result.longEdgeVertices+=int(addedIds.size());
                accepted=true;break;
            }
            if(!accepted)break;
        }
    }
    if(trace)std::fprintf(stderr,"[hole-stage] quality done\n");
    const auto generatedEdges=edgesOfTriangles(tri.tris);
    {std::set<int> outerSet;for(const auto& loop:loops)outerSet.insert(loop.begin(),loop.end());
     for(size_t t=0;t<tri.tris.size();t+=3)if(!outerSet.count(tri.tris[t])&&!outerSet.count(tri.tris[t+1])&&!outerSet.count(tri.tris[t+2]))++result.interiorTriangles;}
    for(auto edge:featureEdges)if(!std::binary_search(generatedEdges.begin(),generatedEdges.end(),edge)) {
        error="三角化未保留结构线约束，拒绝写回";return false;
    }
    for(auto edge:borderEdges)if(!std::binary_search(generatedEdges.begin(),generatedEdges.end(),edge)) {
        error="三角化未保留锁边环约束，拒绝写回";return false;
    }
    for(auto edge:flowEdges)if(!std::binary_search(generatedEdges.begin(),generatedEdges.end(),edge)) {
        error="三角化未保留长条走向约束，拒绝写回";return false;
    }
    if(!verifyTextureCoverage(domain,pts,tri.tris,result.opaquePixels,result.coveredPixels,error))return false;
    if(symmetry.detected) {
        std::set<int> active;for(int v:tri.tris)active.insert(v);
        for(int v:active) {
            const auto q=symmetry.reflect(pts[2*v],pts[2*v+1]);bool paired=false;
            for(int other:active)if(std::hypot(pts[2*other]-q.x,pts[2*other+1]-q.y)<.003){paired=true;break;}
            if(!paired){error="对称点未完整纳入三角网格";return false;}
        }
    }
    std::vector<int> oldTriOf(pts.size()/2,-1);
    std::vector<float> bu(pts.size()/2),bv(pts.size()/2),bw(pts.size()/2),sourcePoints=oldPts,uv=oldUv;
    for(size_t i=size_t(nOld);i<pts.size()/2;++i) {
        bool extrapolated=false;
        oldTriOf[i]=sampleTriangle(oldPixel,oldTris,pts[2*i],pts[2*i+1],bu[i],bv[i],bw[i],&extrapolated);
        if(oldTriOf[i]<0){error="旧 UV 三角形全部退化，无法映射关键形";return false;}
        if(extrapolated)++result.extrapolatedVertices;
        const int t=3*oldTriOf[i],a=oldTris[t],b=oldTris[t+1],c=oldTris[t+2];
        sourcePoints.push_back(bu[i]*oldPts[2*a]+bv[i]*oldPts[2*b]+bw[i]*oldPts[2*c]);
        sourcePoints.push_back(bu[i]*oldPts[2*a+1]+bv[i]*oldPts[2*b+1]+bw[i]*oldPts[2*c+1]);
        const auto mapped=domain.uv(pts[2*i],pts[2*i+1]);uv.push_back(mapped.first);uv.push_back(mapped.second);
    }
    // Remove unreferenced legacy points instead of carrying a dense old mesh
    // into every new result. Bound identities stay available to existing glue.
    std::vector<uint8_t> used(pts.size()/2);
    for(int i:tri.tris)used[size_t(i)]=1;
    for(int i:protectedOld)used[size_t(i)]=1;
    std::vector<int> kept,remap(used.size(),-1);
    for(size_t i=0;i<used.size();++i)if(used[i]){remap[i]=int(kept.size());kept.push_back(int(i));}
    auto compact=[&](const std::vector<float>& values) {
        std::vector<float> result;result.reserve(kept.size()*2);
        for(int i:kept){result.push_back(values[2*i]);result.push_back(values[2*i+1]);}return result;
    };
    const auto generatedPoints=compact(sourcePoints),generatedUv=compact(uv);
    for(int& i:tri.tris)i=remap[size_t(i)];
    for(auto& loop:loops)for(int& i:loop)i=remap[size_t(i)];
    for(auto& edge:featureEdges){edge.first=remap[size_t(edge.first)];edge.second=remap[size_t(edge.second)];}
    for(auto& edge:borderEdges){edge.first=remap[size_t(edge.first)];edge.second=remap[size_t(edge.second)];}
    for(auto& edge:flowEdges){edge.first=remap[size_t(edge.first)];edge.second=remap[size_t(edge.second)];}
    // Validate the actual float UVs that will be serialized, too. The 0.001
    // pixel-area tolerance only accommodates UV roundoff on large atlas pages.
    auto serializedPixels=generatedUv;
    for(size_t i=0;i<serializedPixels.size();i+=2) {
        const auto point=domain.pixel(serializedPixels[i],serializedPixels[i+1]);
        serializedPixels[i]=point.first;serializedPixels[i+1]=point.second;
    }
    if(!verifyTextureCoverage(domain,serializedPixels,tri.tris,result.opaquePixels,result.coveredPixels,error,.001f))return false;
    for(float value:generatedPoints)if(!std::isfinite(value)){error="源图外推产生无效位置，拒绝写回";return false;}
    // 写 source positions / uvs / indices。
    XsNode* source = nullptr;
    {
        std::vector<XsNode*> stack{&work.root};
        while (!stack.empty()) {
            XsNode* cur = stack.back(); stack.pop_back();
            if (cur->id == m.sourceXsId && cur->type == "CArtMeshSource") { source = cur; break; }
            for (auto& c : cur->children) stack.push_back(&c);
        }
    }
    if (!source) { error = "work 里找不到 ArtMesh source"; return false; }
    std::vector<int> path;
    const XsNode& docRoot = doc.xml().root;
    path.clear();
    if (!findPath(docRoot, m.positionsNode, path)) { error = "positions 路径失败"; return false; }
    XsNode* positions = atPath(work.root, path);
    path.clear();
    if (!findPath(docRoot, m.uvsNode, path)) { error = "uvs 路径失败"; return false; }
    XsNode* uvs = atPath(work.root, path);
    path.clear();
    if (!findPath(docRoot, m.indicesNode, path)) { error = "indices 路径失败"; return false; }
    XsNode* indices = atPath(work.root, path);
    if (!positions || !uvs || !indices) {
        error = "work 路径解析为空";
        return false;
    }
    setArray(*positions, joinFloats(generatedPoints));
    setArray(*uvs, joinFloats(generatedUv));
    setArray(*indices, joinInts(tri.tris));

    // GEditableMesh2：point / pointUid / edge / edgePriority。
    XsNode* editable = findTypeMut(*source, "GEditableMesh2");
    if (editable) {
        if (XsNode* point = deepMut(*editable, "point")) {
            setArray(*point, joinFloats(generatedPoints));
        }
        if (XsNode* puid = deepMut(*editable, "pointUid")) {
            std::vector<int> uids = ints(arrayText(work.root, puid));
            if(uids.size()!=size_t(nOld)){error="编辑顶点 UID 数量与旧网格不一致";return false;}
            int next = 1;
            for (int u : uids) next = std::max(next, u + 1);
            const int ai = editable->attrIndex("nextPointUid");
            if (ai >= 0) next = std::max(next, std::atoi(editable->attrs[static_cast<size_t>(ai)].second.c_str()));
            for (int i = nOld; i < static_cast<int>(pts.size() / 2); ++i) uids.push_back(next++);
            std::vector<int> compactUids;compactUids.reserve(kept.size());
            for(int i:kept)compactUids.push_back(uids[size_t(i)]);
            setArray(*puid, joinInts(compactUids));
            if (ai >= 0) editable->attrs[static_cast<size_t>(ai)].second = std::to_string(next);
            else editable->attrs.emplace_back("nextPointUid", std::to_string(next));
        }
        if (XsNode* edge = deepMut(*editable, "edge")) {
            std::set<std::pair<int, int>> seen;
            std::vector<int> ev, ep;
            std::set<std::pair<int, int>> boundary;
            boundary.insert(featureEdges.begin(),featureEdges.end());
            boundary.insert(borderEdges.begin(),borderEdges.end());
            boundary.insert(flowEdges.begin(),flowEdges.end());
            for (const auto& loop:loops)for (size_t i = 0; i < loop.size(); ++i) {
                int a = loop[i], b = loop[(i + 1) % loop.size()];
                if (a > b) std::swap(a, b);
                boundary.emplace(a, b);
            }
            for (size_t t = 0; t + 2 < tri.tris.size(); t += 3) {
                const int ii[3] = {tri.tris[t], tri.tris[t + 1], tri.tris[t + 2]};
                for (int k = 0; k < 3; ++k) {
                    int a = ii[k], b = ii[(k + 1) % 3];
                    if (a > b) std::swap(a, b);
                    if (!seen.insert({a, b}).second) continue;
                    ev.push_back(a); ev.push_back(b);
                    ep.push_back(boundary.count({a, b}) ? 30 : 10);
                }
            }
            setArray(*edge, joinInts(ev));
            if (XsNode* prio = deepMut(*editable, "edgePriority")) setArray(*prio, joinInts(ep));
        }
    }

    // Include ordinary, extended-grid, and BlendShape forms; shared form
    // references are written once, while shared position arrays use COW below.
    std::vector<std::pair<const XsNode*,std::vector<float>>> forms;
    std::set<const XsNode*> seenForms;
    auto addForm=[&](const XsNode* ref,const std::vector<float>* norm=nullptr) {
        const XsNode* form=ref?doc.resolveRef(ref):nullptr;
        if(form&&seenForms.insert(form).second)
            forms.emplace_back(form,norm?*norm:floats(arrayText(doc.xml().root,deep(*form,"positions"))));
    };
    for(const auto& cell:m.grid.cells)addForm(cell.form,&cell.norm);
    for(const auto* form:m.keyformForms)addForm(form);
    for(const auto& morph:m.morphTargets)addForm(morph.form,&morph.norm);
    for (const auto& formData : forms) {
        if (formData.second.size() != static_cast<size_t>(nOld) * 2) {
            error = "某个 keyform 的 norm 与旧顶点数不一致，拒绝只写一部分";
            return false;
        }
        std::vector<float> np(formData.second);
        for (int pi = nOld; pi < static_cast<int>(pts.size() / 2); ++pi) {
            const int t = oldTriOf[static_cast<size_t>(pi)] * 3;
            const int a = oldTris[static_cast<size_t>(t)], b = oldTris[static_cast<size_t>(t + 1)],
                      cc = oldTris[static_cast<size_t>(t + 2)];
            for (int k = 0; k < 2; ++k)
                np.push_back(bu[static_cast<size_t>(pi)] * formData.second[2 * a + k] +
                             bv[static_cast<size_t>(pi)] * formData.second[2 * b + k] +
                             bw[static_cast<size_t>(pi)] * formData.second[2 * cc + k]);
        }
        std::vector<int> fpath;
        if (!findPath(doc.xml().root, formData.first, fpath)) { error = "keyform 路径失败"; return false; }
        XsNode* form = atPath(work.root, fpath);
        XsNode* fpos = form ? deepMut(*form, "positions") : nullptr;
        if (!fpos && form) {
            XsNode position; position.type = "float-array"; position.fieldName = "positions";
            position.attrs.emplace_back("xs.n", "positions");
            form->children.push_back(std::move(position));
            fpos = &form->children.back();
        }
        if (!fpos) { error = "keyform positions 缺失，拒绝"; return false; }
        if (fpos->ref >= 0) {
            // 共享数组 COW：把该 form 的引用改成自己的 inline 数组，不污染其它 form。
            fpos->ref = -1;
            for (size_t ai = 0; ai < fpos->attrs.size();) {
                if (fpos->attrs[ai].first == "xs.ref") fpos->attrs.erase(fpos->attrs.begin() + static_cast<std::ptrdiff_t>(ai));
                else ++ai;
            }
        }
        fpos->selfClosing = false;
        const auto formPositions=compact(np);
        for(float value:formPositions)if(!std::isfinite(value)){error="关键形外推产生无效位置，拒绝写回";return false;}
        setArray(*fpos, joinFloats(formPositions));
    }
    result.oldVertices = nOld;
    result.newVertices = static_cast<int>(kept.size());
    result.oldTriangles = static_cast<int>(oldTris.size() / 3);
    result.newTriangles = static_cast<int>(tri.tris.size() / 3);
    result.changed = generatedPoints != oldPts || generatedUv != oldUv || tri.tris != oldTris;
    result.featureEdges = int(featureEdges.size());
    result.flowEdges = int(flowEdges.size());
    std::set<int> borderVertices;for(auto [a,b]:borderEdges){borderVertices.insert(a);borderVertices.insert(b);}
    result.borderVertices=int(borderVertices.size());result.borderEdges=int(borderEdges.size());result.borderRings=int(borderLoops.size());
    for(const auto& loop:loops)result.boundaryVertices+=int(loop.size());
    result.holeContours=int(domain.holes.size());result.holeBorderRings=int(domain.holeSupports.size());
    result.mergedHoleVertices=domain.mergedHoleVertices;
    for(size_t i=outerCount;i<loops.size();++i)result.holeVertices+=int(loops[i].size());
    result.supportVertices=added;result.stripSampling=strip;result.components=int(outerCount);result.conservativeFallback=domain.conservativeFallback;
    result.sharpTips=domain.sharpTips;
    return true;
}
} // namespace

MeshNameMode MeshAutoGen::modeFromName(const std::string& name) {
    // ASCII word boundaries prevent "hair" matching "chair". Split CamelCase
    // before folding case; Chinese terms use substring matching. More specific
    // compounds consume their span before generic terms (sword tassel / 剑穗).
    auto lower=[](unsigned char c){return c>='a'&&c<='z';};
    auto upper=[](unsigned char c){return c>='A'&&c<='Z';};
    auto letter=[&](unsigned char c){return lower(c)||upper(c);};
    std::string text;
    for(size_t i=0;i<name.size();++i) {
        const unsigned char c=name[i];
        if(i&&upper(c)&&(lower(static_cast<unsigned char>(name[i-1]))||
            (upper(static_cast<unsigned char>(name[i-1]))&&i+1<name.size()&&lower(static_cast<unsigned char>(name[i+1])))))text+=' ';
        text+=upper(c)?char(c-'A'+'a'):char(c);
    }
    struct Rule {std::string term;MeshAutoGenStyle style;};
    std::vector<Rule> rules;
    auto add=[&](MeshAutoGenStyle style,std::initializer_list<const char*> terms){for(auto term:terms)rules.push_back({term,style});};
    add(MeshAutoGenStyle::Soft,{"软质","裙摆","裙子","裙","布料","布条","披风","头发","刘海","侧发","后发","发丝","袖子","飘带","丝带","剑穗",
        "soft","skirt","cloth","fabric","hair","bangs","cape","sleeve","ribbon","tassel","sword tassel"});
    add(MeshAutoGenStyle::Rigid,{"硬质","金属","发夹","发卡","剑","刀","盔甲","纽扣","拉链","镜框",
        "rigid","sword","blade","metal","armor","armour","button","zipper","buckle","glasses","hairclip","hairpin","hair clip","hair pin"});
    add(MeshAutoGenStyle::FineLine,{"眼线","睫毛","eyeliner","eye liner","eyelash","eyelashes","lash","lashes"});
    // These describe facial features, not weapons; keep material unspecified.
    add(MeshAutoGenStyle::Automatic,{"剑眉","刀疤"});
    std::stable_sort(rules.begin(),rules.end(),[](const Rule& a,const Rule& b){return a.term.size()>b.term.size();});
    // Treat common separators identically in multiword English compounds.
    std::string normalized;
    for(char c:text) {
        const bool separator=c=='_'||c=='-'||c==' '||c=='\t'||c=='\r'||c=='\n';
        if(separator){if(!normalized.empty()&&normalized.back()!=' ')normalized+=' ';}
        else normalized+=c;
    }
    text=std::move(normalized);
    std::vector<bool> consumed(text.size(),false);std::set<MeshAutoGenStyle> modes;std::string hits;
    for(const auto& rule:rules)for(size_t pos=text.find(rule.term);pos!=std::string::npos;pos=text.find(rule.term,pos+1)) {
        const size_t end=pos+rule.term.size();const bool ascii=static_cast<unsigned char>(rule.term[0])<128;
        if(ascii&&((pos&&letter(text[pos-1]))||(end<text.size()&&letter(text[end]))))continue;
        if(std::any_of(consumed.begin()+pos,consumed.begin()+end,[](bool v){return v;}))continue;
        std::fill(consumed.begin()+pos,consumed.begin()+end,true);modes.insert(rule.style);
        if(!hits.empty())hits+="、";hits+=rule.term;
    }
    if(modes.empty())return {MeshAutoGenStyle::Automatic,"未命中名称关键词，采用常规布点"};
    if(modes.size()>1)return {MeshAutoGenStyle::Automatic,"名称含不同模式关键词（"+hits+"），采用常规布点；可手动指定"};
    return {*modes.begin(),"名称命中："+hits};
}

bool MeshAutoGen::generate(XsDocument& work, const Cmo3Document& doc, int meshIdx,
                           const MeshAutoGenParams& params, MeshAutoGenResult& result,
                           std::string& error, const std::vector<MeshFeatureLine>* features) {
    result = {};
    XsDocument staged = work;
    MeshAutoGenResult generated;
    // A painted straight fold says nothing about rigidity. Resolve the user's
    // material intent before any fine-line, strip, crease or symmetry branch.
    auto effective=params;
    MeshNameMode nameMode;
    if(params.material==MeshAutoGenMaterial::ByName&&meshIdx>=0&&meshIdx<int(doc.artMeshes().size())) {
        nameMode=modeFromName(doc.artMeshes()[size_t(meshIdx)].displayName());effective.style=nameMode.style;
    }
    if(params.material==MeshAutoGenMaterial::Soft)effective.style=MeshAutoGenStyle::Soft;
    if(params.material==MeshAutoGenMaterial::Rigid)effective.style=MeshAutoGenStyle::Rigid;
    if (!generateInto(staged, doc, meshIdx, effective, generated, error, features)) {
        if(!generated.symmetrySampling)return false;
        const auto attempt=generated;const auto reason=error;auto fallback=effective;fallback.autoSymmetry=false;staged=work;
        if(!generateInto(staged,doc,meshIdx,fallback,generated,error,features))return false;
        generated.symmetryDetected=attempt.symmetryDetected;generated.symmetryScore=attempt.symmetryScore;
        generated.symmetryAxis=attempt.symmetryAxis;generated.symmetryFallback=reason;
    }
    work = std::move(staged);
    work.idIndex.clear();
    auto index = [&](auto&& self, XsNode& node) -> void {
        if (node.id >= 0) work.idIndex[node.id] = &node;
        for (auto& child : node.children) self(self, child);
    };
    index(index, work.root);
    result = generated;
    result.effectiveStyle = effective.style;
    result.nameModeReason = nameMode.reason;
    return true;
}

}  // namespace l2d
