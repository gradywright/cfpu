/* 
 * mex_pcCoarsenWse Coarsens a 3D point cloud using Weighted Sample Elimination (WSE)
 *    as implemented in Eliminate function of the cyCodeBase package.
 *    xc = mex_pcCoarsenWse(x,M,area) coarsens the the point cloud x (stored as a
 *    N-by-3 array of doubles) to exactly M points. area is an estimate for the
 *    surface area of the object.  If set to zero then this value is estimated
 *    from the bounding box.
 *    
 * To create the mex file for the function:
 * 1. Download the cyCodeBase: https://github.com/cemyuksel/cyCodeBase
 * 2. Change directories to cfpu/code/+util and type
 *       mex mex_pcCoarsenWse.cpp -Ipath_to_cyCodeBase 
 *    where path_to_cyCodeBase is the path to where cyCodeBase is installed, 
 *    e.g., ~/packages/cyCodeBase
*/
#include <iostream>
#include <vector>
#include <random>
#include "cySampleElim.h"
#include "cyPointCloud.h"
#include "cyPoint.h"
#include "cyHeap.h"
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace cy;

// # define M_PI	3.14159265358979323846

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // checkArguments(outputs, inputs);
        uint Nc = (uint) inputs[1][0];

        TypedArray<double> doubleArray = std::move(inputs[0]);
        int N = doubleArray.getNumberOfElements()/3;
        std::vector< cy::Point3d > inputPoints(N);
		int idx = 0; 
		for ( idx = 0; idx < N; idx++ ) {
		    inputPoints[idx].x = doubleArray[idx][0];
		    inputPoints[idx].y = doubleArray[idx][1];
		    inputPoints[idx].z = doubleArray[idx][2];
		}        

		// create weighted elimination object
		cy::WeightedSampleElimination< Point3d, double, 3, int > wse;

		// execute weighted elimination
		std::vector< cy::Point3d > outputPoints(Nc);
        float area = inputs[2][0];
		float d_max = 2 * wse.GetMaxPoissonDiskRadius( 
									2, outputPoints.size(), area );
		
		bool isProgressive = true;
		wse.Eliminate( inputPoints.data(), inputPoints.size(), 
		               outputPoints.data(), outputPoints.size(), 
		               isProgressive, d_max, 2 );

		Nc = outputPoints.size();

		ArrayFactory f;
		TypedArray<double> xc = f.createArray<double>({Nc, 3});
		for ( idx = 0; idx < Nc; idx++ ) {
		    xc[idx][0] = outputPoints[idx].x;
		    xc[idx][1] = outputPoints[idx].y;
		    xc[idx][2] = outputPoints[idx].z;		    
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
                std::vector<Array>({ factory.createScalar("Input must be scalar double") }));
        }

        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
};
