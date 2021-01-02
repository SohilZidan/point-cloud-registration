# point-cloud-registration
A C++ implementation for Iterative Closest Point algorithm:
* Robust ICP
* Trimmed ICP

## reading point cloud:
It is only working with .ply extension using the function `readPointCloud(const std::string filename)` which returns an Eigen::Matrix<float, nSamples, 3> object.

## ICP:
1. Create an object of `icp::Icp` or `icp::TrIcp` with the required arguments:  
```
icp::PointCloud P, M;  
icp::Icp robusticp(P,M);  
```
2. run the algorithm
```
robusticp.run();
```
3. performance measurement:  
`icp::Icp.getIter()`: #iteration taken for until convergence  
`icp::Icp.getNpo()`: last #trimmedpoints considered in computations  
`icp::Icp.getSts()`: sum of squared distances  
`icp::Icp.getMse()`: mean squared error or `Sts/Npo`  
`icp::Icp.getR()`: an accumulated 3X3 rotation matrix  
`icp::Icp.rotationError(Eigen::Matrix3f R)`: computes the rotation error given a ground truth rotation  
`icp::Icp.getT()`: an accumulated `Eigen::Vector3f` translation vector  
