#pragma once

namespace icp 
{
    // data types
    typedef std::vector<std::tuple<float, int, int>> PointsMapper;

    class Icp
    {
    protected:
        icp::PointCloud* P;
        icp::PointCloud* M;
        const size_t dim = 3;
        PointsMapper distances;

        EigenMatrixKDTree m_index;// = EigenMatrixKDTree();
        nanoflann::KNNResultSet<float> resultSet;
        size_t m_ret_index;
        float m_out_dist_sqr;
        

        EigenMatrixKDTree p_index;//(this->dim, this->M, 10);
        std::vector<float> p_query_pt;

        size_t iter = 0;
        float Sts;
        float e, e_prev= FLT_MAX; // Trimmed MSE
        float Stsprime = FLT_MAX;

        size_t Np, Npo;
        // Stopping Conditions
        size_t Niter = 100;
        float mse_tol = 0.01;
        float tol = 0.01, stsTolerance = 20.0;

        // motion affine transformation matrix
        icp::PointCloud R;
        icp::PointCloud Rtotal = icp::PointCloud::Identity(3,3);
        Eigen::Vector3f t;
        Eigen::Vector3f ttotal = Eigen::Vector3f(0,0,0);
        Eigen::Matrix<float, Eigen::Dynamic , Eigen::Dynamic> T = Eigen::Matrix<float, Eigen::Dynamic , Eigen::Dynamic>::Identity(4,4);

        // wrongPair rejection algo
        bool reciprocal = false;
        void (Icp::*closestPtsFinder)();
        float eps = 0.01f;
    public:
        // icp(/* args */);
        Icp(icp::PointCloud& P, icp::PointCloud& M);
        Icp(icp::PointCloud& , icp::PointCloud& , bool);
        ~Icp();
        //
        void init();
        void constructKDTree();
        void initdists();
        virtual void closestpoints();
        void closestpointsReciprocal();
        virtual void updateNpo();
        void calculateErrors();
        bool convergencetest();
        void calculateMotion();
        void transform();
        virtual void updateParameters();
        void run();
        // getters
        size_t getIter();
        size_t getNpo();
        float getSts();
        float getMse();
        icp::PointCloud getR();
        float rotationError(Eigen::Matrix3f R);
        Eigen::Vector3f getT();
        
    };
    

    class TrIcp: public Icp
    {
    private:
        bool isUpdatedNpo = false;
        float xi=0.5f;
    public:
        TrIcp(icp::PointCloud& P, icp::PointCloud& M);
        TrIcp(icp::PointCloud& P, icp::PointCloud& M, float xi);
        ~TrIcp();
        //
        void closestpoints();
        void updateNpo();
        void updateParameters();

    };
    
    
    
}
