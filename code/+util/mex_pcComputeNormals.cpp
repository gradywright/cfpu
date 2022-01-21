/*
 * mex_pcComputeNormals Computes the normals for a 3D point cloud using 
 *    PointCloudNormal function of the VCGLIB package.
 *
 *    xc = mex_pcComputeNormals(x,nn,siter) computes normals for the point cloud
 *    x (stored as a N-by-3 array of doubles) using nn neighbors for each point
 *    and siter smoothing operations. 
 *    
 * To create the mex file for the function:
 * 1. Download the vcglib: https://github.com/cnr-isti-vclab/vcglib
 * 2. Change directories to cfpu/code/+util and type
 *       mex mex_pcComputeNormals.cpp -Ipath_to_vcglib -Ipath_to_vcglib_eigenlib
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
		MyMesh m;

        checkArguments(outputs, inputs);
        uint fittingAdjNum = (uint) inputs[1][0];
        uint smoothingIterNum = (uint) inputs[2][0];

        TypedArray<double> doubleArray = std::move(inputs[0]);
        uint N = doubleArray.getNumberOfElements()/3;

		Allocator<MyMesh>::AddVertices(m,N);
		int idx = 0; 
		for ( idx = 0; idx < N; idx++ ) {
			m.vert[idx].P() = MyMesh::CoordType(doubleArray[idx][0], doubleArray[idx][1], doubleArray[idx][2]);
		}

		vcg::tri::UpdateBounding<MyMesh>::Box(m); // updates bounding box

		vcg::tri::PointCloudNormal<MyMesh>::Param p;
		p.fittingAdjNum = fittingAdjNum;
		p.smoothingIterNum = smoothingIterNum;
		p.useViewPoint = 0;
		vcg::tri::PointCloudNormal<MyMesh>::Compute(m, p, 0);

      	ArrayFactory f;
      	TypedArray<double> xc = f.createArray<double>({N, 3});
		for ( idx = 0; idx < N; idx++ ) {
			xc[idx][0] = -m.vert[idx].N()[0];
			xc[idx][1] = -m.vert[idx].N()[1];
			xc[idx][2] = -m.vert[idx].N()[2];
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

        // Check second input argument
        if (inputs[2].getType() != ArrayType::DOUBLE ||
            inputs[2].getNumberOfElements() != 1)
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
