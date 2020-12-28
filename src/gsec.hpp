#pragma once

const float GOLDEN_RATIO = (std::sqrt(5.0) + 1) / 2.0; // phi
const float INVPHI = (std::sqrt(5.0) - 1) / 2.0;  //# 1 / phi
const float INVPHI2 = (3 - std::sqrt(5.0)) / 2.0;  //# 1 / phi^2

// template<typename... args>
float* gss(
    float (&func)(std::vector<std::tuple<float, int, int>>&, float), std::vector<std::tuple<float, int, int>>& arg1, 
    float a, float b,
    float tol = 1e-5);

// template<typename... args>
float* gss(
    float (&func)(std::vector<std::tuple<float, int, int>>&, float), std::vector<std::tuple<float, int, int>>& arg1, 
    float a, float b,
    float tol,
    float h,
    bool noC, float c, float fc,
    bool noD, float d, float fd
);