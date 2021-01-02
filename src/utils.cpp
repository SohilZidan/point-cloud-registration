#include <vector>
#include <tuple>
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include "../happly.h"
#include "../nanoflann.hpp"
#include "utils.hpp"



float icp::mse_cost(std::vector<std::tuple<float, int, int>>& distances, float xi)
{
    size_t Np = distances.size();
    size_t Npo = xi*Np;
    
    float Sts = std::accumulate(distances.begin(), distances.begin() + Npo, 0, tupleAccumulateOP);
    xi = (float)std::pow(xi, 3);
    float e = Sts/Npo;
    
    return e/xi;
} 

bool icp::comparePairs(const std::tuple<float, int, int>& lhs, const std::tuple<float, int, int>& rhs)
{
    return std::get<0>(lhs) < std::get<0>(rhs);
}

float icp::tupleAccumulateOP(const float &a, const std::tuple<float, int, int> &b)
{
    return a+std::get<0>(b);
}

icp::PointCloud icp::readPointCloud(const std::string filename)
{
    //"../data/fountain_a.ply"
    happly::PLYData plyData(filename);
    // auto propsNames = plyData.getElement("vertex").getPropertyNames();
    std::vector<std::array<double, 3>> vPos = plyData.getVertexPositions();
    size_t dim = 3;
    size_t N = vPos.size();
    icp::PointCloud mat(N, dim);
  
    mat.resize(N, dim);
    for (size_t i = 0; i < N; i++)
        for (size_t d = 0; d < dim; d++)
            mat(i, d) = vPos[i][d];
    // for(const auto& prop: propsNames)
    //     std::cout<<prop<<std::endl;
    return mat;
}

icp::PointCloud icp::addNoise(icp::PointCloud M)
{
    // get dimentions
    size_t rows = M.rows();
    size_t cols = M.cols();
    
    return icp::PointCloud::Random(rows, cols);
}


// rotation
Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> icp::rotate(icp::PointCloud& M, float deg, Eigen::Vector3f axis)
{
    // float deg = 10.f;
    float rad = deg * M_PI / 180;
    // Eigen::Vector3f axis(0,1,1);
    Eigen::Transform<float, 3, Eigen::Affine> t;
    Eigen::AngleAxis<float> rot(rad, axis);
    t = rot;

    // std::cout<<t.linear()<<std::endl;
    // std::cout<<t.linear().inverse()<<std::endl;
    // std::cout<<t.linear().inverse().matrix()<<std::endl;
    M = (t.linear() * M.transpose()).transpose();
    return t.linear().matrix();
}