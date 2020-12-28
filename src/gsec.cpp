#include<tuple>
#include<vector>
#include<cmath>
#include "gsec.hpp"

// template<typename... args>
float* gss(
    float (&func)(std::vector<std::tuple<float, int, int>>&, float), std::vector<std::tuple<float, int, int>> &arg1, 
    float a, float b,
    float tol)
{
    // std::vector<std::tuple<float, int, int>>
    a = std::min(a,b);
    b = std::max(a,b);
    float h = b - a;
    return gss(func, arg1, a, b, tol, h, true, 0.f, 0.f, true, 0.f, 0.f);
}

// template<class T>
float* gss(
    // float (&func)(T&, float), T& arg1, 
    float (&func)(std::vector<std::tuple<float, int, int>>&, float), std::vector<std::tuple<float, int, int>> &arg1,
    float a, float b,
    float tol,
    float h,
    bool noC, float c, float fc,
    bool noD, float d, float fd)
{
    if(std::abs(h) <= tol) 
    {
        if(noD) return new float[2]{c,c};
        if(noC) return new float[2]{d,d};
        // return new float[2]{a,b};
    }
    
    if(noC){
        c = a + INVPHI2 * h;
        fc = func(arg1, c);
    }
    if(noD){
        d = a + INVPHI * h;
        fd = func(arg1, d);
    }
    if(fc < fd){
        return gss(func, arg1, a, d, tol, h*INVPHI, true, 0.f,0.f, false,c,fc);
    }
    else{
        return gss(func, arg1, c, b, tol, h*INVPHI, false, d,fd, true ,0.f,0.f);
    }
    // rerun float[] {1,2}
}