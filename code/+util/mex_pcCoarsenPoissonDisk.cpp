/*
 * mex_pcCoarsenPoissonDisk Coarsens a 3D point cloud using Poisson Disk Sampling as 
 *    implemented in PoissonPruning function of the VCGLIB package.
 *    xc = mex_pcCoarsenPoissonDisk(x,M) coarsens the the point cloud x (stored as a
 *    N-by-3 array of doubles) to approximately M points. Note that M is a very
 *    rough a approximation for the size of the coarsened set.  Typically, the
 *    coarsened point cloud has many more than M points.
 *    
 * To create the mex file for the function:
 * 1. Download the vcglib: https://github.com/cnr-isti-vclab/vcglib
 * 2. Change directories to cfpu/code/+util and type
 *       mex mex_pcCoarsenPoissonDisk.cpp -Ipath_to_vcglib -Ipath_to_vcglib_eigenlib
 *    where path_to_vcglib is the path to where vcglib is installed, e.g., 
 *    ~/packages/vcglib and path_to_vcglib_eigenlib is the path to where the 
 *    Eigen libray is installed, e.g., ~/packages/vcglib/eigenlib
 */

#include <iostream>
#include <vector>
#include <random>

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/create/resampler.h>
#include <vcg/complex/algorithms/clustering.h>
#include <vcg/simplex/face/distance.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/voronoi_processing.h>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};
class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj,  vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
		typedef vcg::tri::TrivialSampler<MyMesh> BaseSampler;
		MyMesh m;

        checkArguments(outputs, inputs);
        uint Nc = (uint) inputs[1][0];

        TypedArray<double> doubleArray = std::move(inputs[0]);
        uint N = doubleArray.getNumberOfElements()/3;

		Allocator<MyMesh>::AddVertices(m,N);
		int idx = 0; 
		for ( idx = 0; idx < N; idx++ ) {
			m.vert[idx].P() = MyMesh::CoordType(doubleArray[idx][0], doubleArray[idx][1], doubleArray[idx][2]);
		}

		vcg::tri::UpdateBounding<MyMesh>::Box(m); // updates bounding box

		MyMesh::ScalarType radius = vcg::tri::SurfaceSampling<MyMesh,BaseSampler>::ComputePoissonDiskRadius(m,Nc);
		std::vector<MyMesh::VertexPointer> prunedPts;
		vcg::tri::PoissonPruning(m,prunedPts,radius,71205);
		Nc = prunedPts.size();

      	ArrayFactory f;
      	TypedArray<double> xc = f.createArray<double>({Nc, 3});
		for ( idx = 0; idx < Nc; idx++ ) {
			// xc[idx][0] = -m.vert[idx].N().dot(MyMesh::CoordType(1.0, 0.0, 0.0));
			// xc[idx][1] = -m.vert[idx].N().dot(MyMesh::CoordType(0.0, 1.0, 0.0));
			// xc[idx][2] = -m.vert[idx].N().dot(MyMesh::CoordType(0.0, 0.0, 1.0));
			xc[idx][0] = prunedPts[idx]->P()[0];
			xc[idx][1] = prunedPts[idx]->P()[1];
			xc[idx][2] = prunedPts[idx]->P()[2];
		}

      	outputs[0] = xc;
	}

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check first input argument
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE ||
            inputs[0].getNumberOfElements() == 1)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("First input must be an array of doubles") }));
        }

        // Check second input argument
        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be scalar integer") }));
        }
        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
};
