#include<iostream>
#include <vector>
#include <tuple>
#include <float.h>
#include <numeric>
#include <Eigen/Dense>
#include "../nanoflann.hpp"
#include "gsec.hpp"
#include "utils.hpp"
#include "icp.hpp"

using namespace nanoflann;

//
// Robust ICP
//

icp::Icp::Icp(icp::PointCloud &P, icp::PointCloud &M): 
    resultSet(1), 
    m_index(dim, std::cref(M), 10),
    p_index(dim, P, 10)
{
    // std::cout<<"&P = "<<&P<<std::endl;
    this->P = &P;
    // std::cout<<"&this->P = "<<(this->P)<<std::endl;
    this->Np = P.rows();
    this->M = &M;
    // 
    this->constructKDTree();
    //
    this->resultSet.init(&(this->m_ret_index), &(this->m_out_dist_sqr));
    this->closestPtsFinder = &icp::Icp::closestpoints;
}

icp::Icp::Icp(icp::PointCloud &P, icp::PointCloud &M, bool reciprocal): icp::Icp(P, M)
{
    this->reciprocal = reciprocal;
    if(this->reciprocal) this->closestPtsFinder = &icp::Icp::closestpointsReciprocal;
}

void icp::Icp::init()
{
    this->iter++;
    this->initdists();
}

void icp::Icp::initdists()
{
    // initialize heap
    this->distances.clear();
    // std::cout<< std::distance(distances.begin(), distances.begin()+Npo) << std::endl;
    // std::cout<< mse_cost(distances, xi) << std::endl;
    std::make_heap(this->distances.begin(), this->distances.end(), icp::comparePairs);
}

void icp::Icp::closestpoints()
{
    // Step 1: Closest Point
    for (size_t i = 0; i < this->Np; i++)
    {
        this->resultSet.init(&(this->m_ret_index), &(this->m_out_dist_sqr));
        this->p_query_pt = {(*this->P)(i,0), (*this->P)(i,1), (*this->P)(i,2)};//{-15.7116, -10.9587, 0.650531};
        // 
        this->m_index.index->findNeighbors(
            this->resultSet, 
            &this->p_query_pt[0],
            nanoflann::SearchParams(10));

        // wrong pair rejection
        if(icp::rejectWrongPair(*this->M, this->m_ret_index, this->p_index, i)) continue;
    
        // update heap(distance)
        std::tuple<float, int, int> tmpDist = std::make_tuple(this->m_out_dist_sqr, i, this->m_ret_index);
        this->distances.push_back(tmpDist);
        std::push_heap(this->distances.begin(), this->distances.end(), icp::comparePairs);
        
    }
}

void icp::Icp::closestpointsReciprocal()
{
    // Step 1: Closest Point
    for (size_t i = 0; i < this->Np; i++)
    {
        
        this->resultSet.init(&(this->m_ret_index), &(this->m_out_dist_sqr));
        //auto query_ptt = P.row(i);
        this->p_query_pt = {(*this->P)(i,0), (*this->P)(i,1), (*this->P)(i,2)};
        // 
        this->m_index.index->findNeighbors(
            this->resultSet, 
            &this->p_query_pt[0],
            nanoflann::SearchParams(10));

        std::vector<float> m_query_pt = {(*this->M)(this->m_ret_index,0), (*this->M)(this->m_ret_index,1), (*this->M)(this->m_ret_index,2)};
        // wrong pair rejection
        if(icp::rejectWrongPair(
            m_query_pt, *this->P, this->p_index, i, 
            this->eps)) 
            continue;
    
        std::tuple<float, int, int> tmpDist = std::make_tuple(this->m_out_dist_sqr, i, this->m_ret_index);
        this->distances.push_back(tmpDist);
        std::push_heap(this->distances.begin(), this->distances.end(), icp::comparePairs);
        
    }
}

void icp::Icp::updateNpo()
{
    this->Npo = this->distances.size();//< Npo? distances.size():Npo;
    // std::cout << "Npo = " << Npo << std::endl;
}

void icp::Icp::calculateErrors()
{
    this->updateNpo();
    // this->Npo = this->distances.size();
    // Step 2: Squares
    this->Sts = std::accumulate(this->distances.begin(), this->distances.begin()+this->Npo, 0.0, icp::tupleAccumulateOP);
    // MSE
    this->e = this->Sts / this->Npo;

}

bool icp::Icp::convergencetest()
{
    // Step 3: Convergence test
    return (
        this->e <= this->mse_tol || 
        std::abs(this->e - this->e_prev)/this->e <= this->tol || 
        this->iter > this->Niter);// ||
        //std::abs(Sts-Stsprime) <= stsTolerance) 
    
}

void icp::Icp::constructKDTree()
{
    // // this->m_index(this->dim, std::cref(this->M), 10);
    // this->m_index.index->buildIndex();
    // // for Geometric criteria: rejecting wrong point pairs
    // this->p_index.index->buildIndex();
}

void icp::Icp::calculateMotion()
{
    // Step 4: Motion Calculation
    // https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    // select Npo points
    icp::PointCloud A(this->Npo, this->dim), B(this->Npo, this->dim);
    for (size_t i = 0; i < (this->Npo); i++)
    {
        auto [ dist, pidx, midx ] = this->distances[i];
        // std::cout<<"dist = "<<dist<<std::endl;
        A.row(i) = this->P->row(pidx);
        B.row(i) = this->M->row(midx);
    }

    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> umeyama = Eigen::umeyama(A.transpose(), B.transpose(), false);
    // std::cout << "The matrix umeyama is of size "
    //         << umeyama.rows() << "x" << umeyama.cols() << std::endl;
    // std::cout << umeyama.block(0,0,3,3) << std::endl;
    this->R = umeyama.block(0,0,3,3);
    this->t  = umeyama.block(0,3,3,1);

    // auto centerA = A.colwise().mean();
    // auto centerB = B.colwise().mean();
    
    // A.rowwise() -= centerA;
    // B.rowwise() -= centerB;
    
    // icp::PointCloud H = A.transpose()*B;
    // Eigen::BDCSVD<Eigen::Matrix<float, -1, -1> > svd(H,Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    // // Rotation
    // this->R = svd.matrixV() * svd.matrixU().transpose();
    // // Translation
    // this->t = centerB.transpose() - this->R * centerA.transpose();
}

void icp::Icp::transform()
{
    // affine transformation
    (*this->P) = (R * (this->P->transpose())).transpose();
    (*this->P).rowwise() += this->t.transpose();
    // t.col(0)
}

void icp::Icp::updateParameters()
{
    // update parameters
    // float* rng = gss(mse_cost, distances, 0.4, 1.0);
    // xi = (rng[0] + rng[1])/ 2.f; //0.5; // 1 -> ICP
    // Npo = xi*Np;
    
    (this->p_index).index->buildIndex();
    (this->e_prev) = this->e;
    (this->Stsprime) = this->Sts;
}


void icp::Icp::run()
{
    do
    {
        // initialize heap
        this->init();

        // Step 1: Closest Point
        // if(this->reciprocal)
        //     this->closestpoints(this->eps);
        // else this->closestpoints();
        (this->*closestPtsFinder)();
        
        // Step 2: Trimmed Squares
        this->calculateErrors();
            
        // Step 3: Convergence test
        if(this->convergencetest()) break;

        // Step 4: Motion Calculation
        this->calculateMotion();

        // Step 5: Transformation
        this->transform();
        
        // Step 6: update parameters
        this->updateParameters();

    } while (true);
}

// getters
size_t icp::Icp::getIter()
{
    return this->iter;
}
size_t icp::Icp::getNpo()
{
    return this->Npo;
}
float icp::Icp::getSts()
{
    return this->Sts;
}
float icp::Icp::getMse()
{
    return this->e;
}

//
icp::Icp::~Icp()
{
    this->P = NULL;
    delete this->P;
    this->M = NULL;
    delete this->M;
    // this->closestPtsFinder = NULL;
    // delete this->closestPtsFinder;
}


//
// Trimmed Icp
//
// Constructors
icp::TrIcp::TrIcp(icp::PointCloud& P, icp::PointCloud& M): icp::Icp(P, M)
{
    this->Npo = this->xi * this->Np;
}

icp::TrIcp::TrIcp(icp::PointCloud& P, icp::PointCloud& M, float xi): icp::Icp(P, M)
{
    this->xi = xi;
    this->Npo = xi * this->Np;
}

// public methods
void icp::TrIcp::closestpoints()
{
    // Step 1: Closest Point
    for (size_t i = 0; i < this->Np; i++)
    {
        this->resultSet.init(&(this->m_ret_index), &(this->m_out_dist_sqr));
        this->p_query_pt = {(*this->P)(i,0), (*this->P)(i,1), (*this->P)(i,2)};//{-15.7116, -10.9587, 0.650531};
        // 
        this->m_index.index->findNeighbors(
            this->resultSet, 
            &this->p_query_pt[0],
            nanoflann::SearchParams(10));

    
        // update heap(distance)
        std::tuple<float, int, int> tmpDist = std::make_tuple(this->m_out_dist_sqr, i, this->m_ret_index);
        this->distances.push_back(tmpDist);
        std::push_heap(this->distances.begin(), this->distances.end(), icp::comparePairs);
        
    }
    std::sort_heap(this->distances.begin(), this->distances.end(),  icp::comparePairs);
}

void icp::TrIcp::updateNpo()
{
    if(!this->isUpdatedNpo)
    {
        float* rng = gss(mse_cost, this->distances, 0.4, 1.0);
        this->xi = (rng[0] + rng[1])/ 2.f; //0.5; // 1 -> ICP
        this->Npo = this->xi * this->Np;
        //
        // this->isUpdatedNpo = true;
    }
}

void icp::TrIcp::updateParameters()
{
    (this->e_prev) = this->e;
    (this->Stsprime) = this->Sts;
}


// destructor   
icp::TrIcp::~TrIcp()
{
}