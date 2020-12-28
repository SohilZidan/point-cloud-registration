#include <iostream>
#include <string>
#include <numeric>
#include <float.h>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "../happly.h"
#include "../nanoflann.hpp"
#include "gsec.hpp"
#include "utils.hpp"
#include "icp.hpp"

using namespace nanoflann;
using namespace icp;


// Refs:
// ICP: http://www-evasion.inrialpes.fr/people/Franck.Hetroy/Teaching/ProjetsImage/2007/Bib/besl_mckay-pami1992.pdf
// TrICP: https://www.researchgate.net/publication/3974183_The_Trimmed_Iterative_Closest_Point_algorithm
// http://people.csail.mit.edu/bkph/articles/Quaternions.pdf
// https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
// Golden Ratio: https://en.wikipedia.org/wiki/Golden_ratio
// Golden Section Search: https://en.wikipedia.org/wiki/Golden-section_search
// https://github.com/nmwsharp/happly
// https://github.com/ClayFlannigan/icp/blob/master/icp.py
//
// ToLookUp
// https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
// https://ieeexplore.ieee.org/document/466913

int main1(int argc, char* argv[])
{
    
    // 
    // Reading Two Data clouds
    // 
    std::cout << "Hello nanoFlann!" <<std::endl;
    PointCloud M, P;
    M = readPointCloud("../data/fountain_a.ply");
    P = readPointCloud("../data/fountain_b.ply");
    std::cout << P.row(0)<<std::endl;
    P += addNoise(P);
    // get dimention
    const size_t dim = 3;
    std::cout << P.row(0)<<std::endl;
    // rotation
    P = icp::rotate(P);
    std::cout << P.row(0)<<std::endl;


    std::cout << "The matrix mat_a is of size "
            << M.rows() << "x" << M.cols() << std::endl;
    std::cout << "It has " << M.size() << " coefficients" << std::endl;

    std::cout << "The matrix mat_b is of size "
                << P.rows() << "x" << P.cols() << std::endl;
    std::cout << "It has " << P.size() << " coefficients" << std::endl;

    // Initialize icp object
    // Icp tricp(P,M);

    // std::cout << Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>::Random(M.rows(), M.cols())<<std::endl;
    

    // 
    // Construct KDTree
    // 
    EigenMatrixKDTree mat_index(dim, std::cref(M), 10 /* max leaf */);
    mat_index.index->buildIndex();
    // for Geometric criteria: rejecting wrong point pairs
    // EigenMatrixKDTree p_index(dim, P, 10 /* max leaf */);
    // p_index.index->buildIndex();
    // tricp.constructKDTree();
    // 
    // Do a KNN search
    // 
    // query point
    // std::vector<float> query_pt(dim);
    // query_pt[0] = -13.0945;
    // query_pt[1] = -12.5863;
    // query_pt[2] = -3.24795;
    // std::vector<float> query_pt = {-15.7116, -10.9587, 0.650531};
    // {-13.5072, -12.4584, -3.282};//{-13.0945, -12.5863, -3.24795};
    // do a knn search
    const size_t num_results = 1;
    size_t ret_index;
    float out_dist_sqr;
    std::vector<float> query_pt;
    // distances: vector of pair<dist2, Pindex, Mindex>
    std::vector<std::tuple<float, int, int>> distances;

    // Step 0: parameters
    nanoflann::KNNResultSet<float> resultSet(num_results);
    //
    float xi = 0.5f; // 1 -> ICP
    bool unknownXi = true;
    size_t Np = P.rows();
    size_t Npo = xi*Np;
    float Sts;
    float e, e_prev= FLT_MAX; // Trimmed MSE
    float Stsprime = FLT_MAX;
    // Stopping Conditions
    size_t Niter = 100;
    float tol = 0.01, stsTolerance = 20;
    size_t iter = 0;

    do
    {
        std::cout << "------------------ Iteration "<< iter <<" -------------------" << std::endl;
        
        // Npo = xi*Np;
        
        // std::cout << "-------------------------------------------------" << std::endl;
        iter++;
        // TODO: set xi automatically: done
        // points rejection
        // 
    
    
        // initialize heap
        distances.clear();
        // std::cout<< std::distance(distances.begin(), distances.begin()+Npo) << std::endl;
        // std::cout<< mse_cost(distances, xi) << std::endl;
        std::make_heap(distances.begin(), distances.end(), comparePairs);
        // tricp.initdists();

        

        // // Step 1: Closest Point
        // tricp.closestpoints();
        for (size_t i = 0; i < Np; i++)
        {
            
            resultSet.init(&ret_index, &out_dist_sqr);
            //auto query_ptt = P.row(i);
            query_pt = {P(i,0), P(i,1), P(i,2)};//{-15.7116, -10.9587, 0.650531};
            // 
            mat_index.index->findNeighbors(
                resultSet, 
                &query_pt[0],
                nanoflann::SearchParams(10));
            // store distance, Pindex, Mindex into heap
            // std::cout << "knnSearch(nn="<<num_results<<"): \n";
            // std::cout << "ret_index=" << ret_index << " out_dist_sqr=" << out_dist_sqr << std::endl;
            // std::cout << "-------------------------------------------"<< std::endl;
            
            // wrong pair rejection
            // if(icp::rejectWrongPair(M, ret_index, p_index, i)) continue;
            Eigen::Vector3f t1 = P.row(i);
            Eigen::Vector3f t2 = M.row(ret_index);
            float norm = (t1-t2).squaredNorm();
            std::tuple<float, int, int> tmpDist = std::make_tuple(norm, i, ret_index);
            distances.push_back(tmpDist);
            // if(out_dist_sqr < 10)
            //     std::cout<< i << ", " << ret_index << ", " << norm << ", " << out_dist_sqr <<std::endl;
            std::push_heap(distances.begin(), distances.end(), comparePairs);
            
        }
        std::sort_heap(distances.begin(), distances.end(),  comparePairs);
        // std::pop_heap(v.begin(), v.end(), comparePairs);
        // int largest = v.back(); 
        // v.pop_back();
        //  for (size_t i = 0; i < 10; i++)
        // {
        //     auto [ dist, pidx, midx ] = distances[i];
        //     std::cout<<"dist = "<<dist<<std::endl;
        //     // A.row(i) = P.row(pidx);
        //     // B.row(i) = M.row(midx);
        // }
        if(unknownXi)
        {
            float* rng = gss(mse_cost, distances, 0.4, 1.0);
            xi = rng[1]; //(rng[0] + rng[1])/ 2.f; //0.5; // 1 -> ICP
            // std::cout << "xi = " << xi << std::endl;
            // std::cout << "xil = " << rng[0] << ", xir = " << rng[1] << std::endl;
            Npo = xi*Np;
            unknownXi = false;
        }
        
        // Step 2: Trimmed Squares
        // distances.resize(Npo);
        // Npo = distances.size();//< Npo? distances.size():Npo;
        std::cout << "Npo = " << Npo << std::endl;
        // std::cout << "dists size = " << distances.size() << std::endl;
        Sts =  std::accumulate(distances.begin(), distances.begin()+Npo, 0.0, tupleAccumulateOP);
        // for (auto it = distances.begin(); it != distances.begin()+Npo; ++it)
        // {
        //     Sts += std::get<0>(*it);
        // }
        
        e = Sts / Npo;
        std::cout << "Trimmed Squares: " << Sts << std::endl;
        std::cout << "Trimmed MSE: " << e << std::endl;
        // tricp.calculateErrors();
            
        // Step 3: Convergence test
        if (
            e <= 0.1 || 
            // std::abs(e-e_prev)/e <= tol || 
            iter == Niter)// ||
            //std::abs(Sts-Stsprime) <= stsTolerance) 
            break;
        // if(tricp.convergencetest(iter)) break;

        // Step 4: Motion Calculation
        // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
        // select Npo points
        PointCloud A(Npo, dim), B(Npo, dim);
        for (size_t i = 0; i < Npo; i++)
        {
            auto [ dist, pidx, midx ] = distances[i];
            // std::cout<<"dist = "<<dist<<std::endl;
            A.row(i) = P.row(pidx);
            B.row(i) = M.row(midx);
        }

        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> umeyama = Eigen::umeyama(A.transpose(), B.transpose(), false);
        // std::cout << "The matrix umeyama is of size "
        //         << umeyama.rows() << "x" << umeyama.cols() << std::endl;
        // std::cout << umeyama.block(0,0,3,3) << std::endl;
        auto R = umeyama.block(0,0,3,3);
        auto t  = umeyama.block(0,3,3,1);
        // std::cout << "The matrix tt is of size "
        //         << tt.rows() << "x" << tt.cols() << std::endl;
        // std::cout << tt << std::endl<<std::endl;

        // auto centerA = A.colwise().mean();
        // auto centerB = B.colwise().mean();
        // // std::cout << "mean before centralizing: " <<  << std::endl;
        // A.rowwise() -= centerA;
        // B.rowwise() -= centerB;
        // std::cout << "mean after centralizing: " << A.colwise().mean() << std::endl;
        // PointCloud H = A.transpose()*B;
        // const auto AA = A;
        // const auto BB = B;
        
    // std::cout << "It has " << P.size() << " coefficients" << std::endl;
        // Eigen::BDCSVD<Eigen::Matrix<float, -1, -1> > svd(H,Eigen::ComputeFullU | Eigen::ComputeFullV);
        // // Eigen::JacobiSVD<Eigen::Matrix<float, -1, -1>> svd(H,Eigen::ComputeFullU | Eigen::ComputeFullV);
        // Eigen::MatrixXf Vt;
        // Vt = svd.matrixV().transpose();
        // Rotation
        // PointCloud R = Vt.transpose() * svd.matrixU().transpose();

        // if (R.determinant() < 0 ){
        //     Vt.block<1,3>(2,0) *= -1;
        //     R = Vt.transpose() * svd.matrixU().transpose();
        // }

        // std::cout << R << std::endl;
        
        // Translation: giving different coefs from 
        // auto t = centerB.transpose() - R * centerA.transpose();
        // std::cout << t << std::endl;
        // tricp.calculateMotion();

        // Step 5: Transformation
        // affine transformation
        /* for (size_t i = 0; i < Npo; i++)
        {
            P.row(i) = (R * P.row(i).transpose() + t).transpose();
        }*/
        // umeyama(0,0);
        // umeyama(2,3);
        // R << umeyama(0,0), umeyama(0,1), umeyama(0,2),
        //     umeyama(1,0), umeyama(1,1), umeyama(1,2),
        //     umeyama(2,0), umeyama(2,1), umeyama(2,2);
        // t[0] = umeyama(0,3);
        // t[1] = umeyama(1,3);
        // t[2] = umeyama(2,3);
        // t=tt;
        // tt.transpose();
        P = (R * (P.transpose())).transpose();
        P.rowwise() += t.col(0).transpose();
        // tricp.transform();
        
        // Step 6: update parameters
        // float* rng = gss(mse_cost, distances, 0.4, 1.0);
        // xi = (rng[0] + rng[1])/ 2.f; //0.5; // 1 -> ICP
        // Npo = xi*Np;
        // p_index.index->buildIndex();
        e_prev = e;
        Stsprime = Sts;
        // tricp.updateParameters();

    } while (true);


    return 0;
}

int main(int argc, char* argv[])
{
    
    // 
    // Reading Two Data clouds
    // 
    PointCloud M, P;
    M = readPointCloud("../data/fountain_a.ply");
    P = readPointCloud("../data/fountain_b.ply");
    // add noise
    P += addNoise(P);
    // rotation
    float deg = 20.f;
    Eigen::Vector3f axis(0,1,0);
    P = icp::rotate(P, deg, axis);


    //
    // Run Robust ICP
    //
    // Initialize icp object
    // std::cout<<"&P = "<<&P<<std::endl;
    Icp robusticp(P,M, true);
    TrIcp tricp(P,M);
    // tricp.setStoppingConditions();


    auto start = std::chrono::high_resolution_clock::now();
    robusticp.run();    
    auto stop = std::chrono::high_resolution_clock::now();


    // millis
    std::chrono::duration<float, std::milli> fp_ms = stop - start;
    auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    // secs
    std::chrono::duration<float, std::ratio<1,1>> fp_s = stop - start;
    auto int_s = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
    //
    std::cout << "robust Icp took " << int_s.count() << "(" << fp_s.count() <<")" << "s" << " / " 
                                    << int_ms.count() << "(" << fp_ms.count() <<")" << "ms." << std::endl;
    std::cout<<std::endl;
    std::cout << "---initial configs---"<<std::endl;
    std::cout << "rotation: "<< deg << "ْ around axis " << axis.transpose() <<  std::endl;
    std::cout<<std::endl;
    std::cout << "---icp results---"<<std::endl;
    std::cout << "#iterations "<< robusticp.getIter() << std::endl;
    std::cout << "valid Npo = " << robusticp.getNpo() << std::endl;
    std::cout << "TriSquares: " << robusticp.getSts() << std::endl;
    std::cout << "TriMSE: " << robusticp.getMse() << std::endl;
    std::cout << std::endl;


    // TrIcp
    start = std::chrono::high_resolution_clock::now();
    tricp.run();    
    stop = std::chrono::high_resolution_clock::now();
    // millis
    fp_ms = stop - start;
    int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    // secs
    fp_s = stop - start;
    int_s = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
    //
    std::cout << "-----------------------" << std::endl;
    std::cout << "robust TrIcp took " << int_s.count() << "(" << fp_s.count() <<")" << "s" << " / " 
                                    << int_ms.count() << "(" << fp_ms.count() <<")" << "ms." << std::endl;
    std::cout<<std::endl;
    std::cout << "---initial configs---"<<std::endl;
    std::cout << "rotation: "<< deg << "ْ around axis " << axis.transpose() <<  std::endl;
    std::cout<<std::endl;
    std::cout << "---tricp results---"<<std::endl;
    std::cout << "#iterations "<< tricp.getIter() << std::endl;
    std::cout << "valid Npo = " << tricp.getNpo() << std::endl;
    std::cout << "TriSquares: " << tricp.getSts() << std::endl;
    std::cout << "TriMSE: " << tricp.getMse() << std::endl;

    return 0;
}