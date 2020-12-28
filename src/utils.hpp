#pragma once


namespace icp {
    // datatypes
    typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> PointCloud;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<PointCloud> EigenMatrixKDTree;    

    // MSE cost function
    float mse_cost(std::vector<std::tuple<float, int, int>>& distances, float xi);

    // pair rejection
    // reject pair p,m if (m is closest point to p) and (p is not closest to p)
    inline
    bool rejectWrongPair(PointCloud &M, size_t j, EigenMatrixKDTree& P, size_t i)
    {
        float out_dist_sqr;
        size_t ret_idx;
        nanoflann::KNNResultSet<float> resultSet(1);
        resultSet.init(&ret_idx, &out_dist_sqr);

        std::vector<float> query_pt = {M(j,0), M(j,1), M(j,2)};

        P.index->findNeighbors(
            resultSet, 
            &query_pt[0],
            nanoflann::SearchParams(10));
        // if(j!=i) std::cout<<i<<", "<<j<<std::endl;
        // P.index[]
        return ret_idx!=i;
    }

    // reject pair p,m if (m is closest point to p) and (p is not closest to p)
    inline
    bool rejectWrongPair(
        // point m
        std::vector<float>& m_query_pt, 
        PointCloud& P, EigenMatrixKDTree& p_index, size_t i, 
        float eps = 0.01)
    {
        float out_dist_sqr;
        size_t ret_idx;
        nanoflann::KNNResultSet<float> resultSet(1);
        resultSet.init(&ret_idx, &out_dist_sqr);

        // std::vector<float> query_pt = {M(j,0), M(j,1), M(j,2)};

        p_index.index->findNeighbors(
            resultSet, 
            &m_query_pt[0],
            nanoflann::SearchParams(10));
        // if(j!=i) std::cout<<i<<", "<<j<<std::endl;
        // P.index[]
        Eigen::Vector3f p_query_pt = P.row(i);
        Eigen::Vector3f p_query_pt_prime = P.row(ret_idx);
        float norm = (p_query_pt - p_query_pt_prime).squaredNorm();

        return norm > eps;
    }

    bool comparePairs(const std::tuple<float, int, int>& lhs, const std::tuple<float, int, int>& rhs);

    float tupleAccumulateOP(const float &a, const std::tuple<float, int, int> &b);

    PointCloud readPointCloud(const std::string filename);

    // noise
    PointCloud addNoise(PointCloud M);

    // rotation
    PointCloud rotate(PointCloud &M, float deg = 20.f, Eigen::Vector3f axis = Eigen::Vector3f(0,1,0));
    
}
