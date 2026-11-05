#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace l2d::meshSymmetry {
struct Point { double x=0,y=0; };
struct Axis {
    bool detected=false;
    double x=0,y=0,dx=0,dy=1,score=0;
    Point local(double px,double py) const {return {(px-x)*dy-(py-y)*dx,(px-x)*dx+(py-y)*dy};}
    Point world(double u,double v) const {return {x+u*dy+v*dx,y-u*dx+v*dy};}
    Point reflect(double px,double py) const {const auto q=local(px,py);return world(-q.x,q.y);}
};
// Geometry only: RGB shading is deliberately excluded from shape detection.
inline Axis detect(const std::vector<uint8_t>& rgba,int w,int h) {
    Axis best;if(w<3||h<3||rgba.size()!=size_t(w)*h*4)return best;
    std::vector<Point> samples;
    const int step=std::max(1,int(std::ceil(std::max(w,h)/160.)));
    double cx=0,cy=0;
    for(int y=0;y<h;y+=step)for(int x=0;x<w;x+=step)if(rgba[(size_t(y)*w+x)*4+3]) {
        samples.push_back({x+.5,y+.5});cx+=x+.5;cy+=y+.5;
    }
    if(samples.size()<20)return best;
    cx/=samples.size();cy/=samples.size();
    auto evaluate=[&](double angle,double offset) {
        Axis a;a.dx=std::sin(angle);a.dy=std::cos(angle);a.x=cx+offset*a.dy;a.y=cy-offset*a.dx;
        int matched=0;
        for(auto p:samples) {
            const auto r=a.reflect(p.x,p.y);const int x=int(std::floor(r.x)),y=int(std::floor(r.y));
            matched+=x>=0&&x<w&&y>=0&&y<h&&rgba[(size_t(y)*w+x)*4+3]>0;
        }
        a.score=double(matched)/(2*samples.size()-matched);
        if(a.score>best.score+1e-9)best=a;
    };
    const double pi=3.14159265358979323846,range=std::min(w,h)*.025;
    for(int degree=0;degree<180;degree+=10)for(int offset=-2;offset<=2;++offset)evaluate(degree*pi/180,offset*range/2);
    for(double angular:{1.,.1}) {
        const Axis previous=best;const double angle=std::atan2(previous.dx,previous.dy);
        const double offset=(previous.x-cx)*previous.dy-(previous.y-cy)*previous.dx;
        for(int d=-10;d<=10;++d)for(int k=-2;k<=2;++k)evaluate(angle+d*angular*pi/180,offset+k*range*(angular==1?.1:.02));
    }
    best.detected=best.score>=.90;
    return best;
}

// A single interval per cross section excludes branches/holes in this branch.
// Outside uses the union with the reflection; inside uses their intersection.
// This retains every original texture pixel without forcing the image itself
// to become symmetric. The paired profile is reduced only conservatively.
inline std::vector<float> pairedContour(const std::vector<float>& polygon,const Axis& axis,bool outside,double tolerance,double spacing=0) {
    std::vector<Point> p;std::vector<double> rows;
    for(size_t i=0;i<polygon.size();i+=2){auto q=axis.local(polygon[i],polygon[i+1]);p.push_back(q);rows.push_back(q.y);}
    std::sort(rows.begin(),rows.end());rows.erase(std::unique(rows.begin(),rows.end(),[](double a,double b){return std::abs(a-b)<1e-5;}),rows.end());
    if(rows.size()<2)return {};
    auto crossings=[&](double y) {
        std::vector<double> x;
        for(size_t i=0;i<p.size();++i) {
            const auto a=p[i],b=p[(i+1)%p.size()];
            if(y<std::min(a.y,b.y)-1e-7||y>std::max(a.y,b.y)+1e-7)continue;
            if(std::abs(a.y-b.y)<1e-8){if(std::abs(y-a.y)<1e-7){x.push_back(a.x);x.push_back(b.x);}continue;}
            const double t=(y-a.y)/(b.y-a.y);if(t>=-1e-7&&t<=1+1e-7)x.push_back(a.x+t*(b.x-a.x));
        }
        std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end(),[](double a,double b){return std::abs(a-b)<1e-4;}),x.end());return x;
    };
    for(size_t i=1;i<rows.size();++i) {
        const double y=(rows[i-1]+rows[i])*.5;
        const auto x=crossings(y);
        const bool cap=std::min(y-rows.front(),rows.back()-y)<std::max(2.,(rows.back()-rows.front())*.03);
        if(x.empty()||(!cap&&(x.size()!=2||x.front()>1e-4||x.back()<-1e-4)))return {};
    }
    std::vector<Point> profile;
    for(double y:rows) {
        const auto x=crossings(y);if(x.empty())return {};
        double radius=outside?std::max(std::abs(x.front()),std::abs(x.back())):std::min(-x.front(),x.back());
        if(radius<-.001){if(!outside)continue;return {};}
        profile.push_back({std::max(0.,radius),y});
    }
    if(profile.size()<2)return {};
    // Remove only stations whose replacement stays on the conservative side.
    for(bool changed=true;changed&&profile.size()>2;) {
        changed=false;
        for(size_t i=1;i+1<profile.size();++i) {
            const auto a=profile[i-1],b=profile[i+1],q=profile[i];
            const double r=a.x+(b.x-a.x)*(q.y-a.y)/(b.y-a.y),delta=r-q.x;
            if(std::abs(delta)<=tolerance&&(outside||delta<=1e-7)) {
                profile.erase(profile.begin()+i);changed=true;break;
            }
        }
    }
    if(spacing>0) {
        std::vector<Point> sampled;
        for(size_t i=1;i<profile.size();++i) {
            const auto a=profile[i-1],b=profile[i];
            const int n=std::max(1,int(std::lround(std::hypot(a.x-b.x,a.y-b.y)/spacing)));
            for(int k=0;k<n;++k)sampled.push_back({a.x+(b.x-a.x)*k/n,a.y+(b.y-a.y)*k/n});
        }
        sampled.push_back(profile.back());profile=std::move(sampled);
    }
    std::vector<float> out;
    auto append=[&](double x,double y){const auto q=axis.world(x,y);if(out.empty()||std::hypot(out[out.size()-2]-q.x,out.back()-q.y)>1e-5){out.push_back(float(q.x));out.push_back(float(q.y));}};
    for(auto q:profile)append(q.x,q.y);
    for(auto i=profile.rbegin();i!=profile.rend();++i)append(-i->x,i->y);
    if(out.size()>4&&std::hypot(out[0]-out[out.size()-2],out[1]-out.back())<1e-5)out.resize(out.size()-2);
    return out.size()>=6?out:std::vector<float>{};
}
} // namespace l2d::meshSymmetry
