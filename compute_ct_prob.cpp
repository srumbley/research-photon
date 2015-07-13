#include <random>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include "mex.h"
#include "matrix.h"

using namespace std;

struct Detection {
	size_t idx;
	int t;

	Detection() {

	}

	Detection(size_t pixelIndex, int detectionTime) {
		idx = pixelIndex;
		t = detectionTime;
	}

	friend bool operator<(Detection d1, Detection d2) {
		//mexPrintf("comparing %d > %d: %d \n", d1.t, d2.t, (d1.t > d2.t));
		return d1.t < d2.t;
	}
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Load data
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("simulate_crosstalk:invalidNumOutputs", "One output argument required.");
	}
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("simulate_crosstalk:invalidNumInputs", "Two input arguments required.");
	}
	if (!mxIsSparse(prhs[0])) {
		mexErrMsgIdAndTxt("simulate_crosstalk:detectionsNotSparse", "Input detections must be sparse.");
	}

	const int numPixels = (int)mxGetM(prhs[0]);
	const int numPulses = (int)mxGetN(prhs[0]);
	const mwIndex* detectionsJc = mxGetJc(prhs[0]);
	if (detectionsJc[numPulses] == 0) {
		plhs[0] = mxCreateSparseLogicalMatrix(numPixels, numPulses, 0);
		return;
	}
	double* detectionsPr = mxGetPr(prhs[0]);
	mwIndex* detectionsIr = mxGetIr(prhs[0]);

	if (mxGetNumberOfDimensions(prhs[1]) != 3) {
		mexErrMsgIdAndTxt("simulate_crosstalk:pctNot3D", "Crosstalk probability distribution should be 3-D.");
	}
	const mwSize* pctDims = mxGetDimensions(prhs[1]);
	const int ctMaxDelay = (int)pctDims[0];
	const mwSize nrows = (pctDims[2] + 1) / 2;
	const mwSize ncols = (pctDims[1] + 1) / 2;
	if (nrows*ncols != numPixels) {
		mexErrMsgIdAndTxt("simulate_crosstalk:pctWrongDims", "Crosstalk probability matrix should be size [ctMaxDelay, ncols, nrows].");
	}
	const double* pct = mxGetPr(prhs[1]);

	plhs[0] = mxCreateSparse(numPixels, numPulses, detectionsJc[numPulses], mxREAL);
	double* lctPr = mxGetPr(plhs[0]);
	mwIndex* lctIr = mxGetIr(plhs[0]);
	mwIndex* lctJc = mxGetJc(plhs[0]);

	int frameCounts;
	size_t cIdx = 0;
	size_t prIdx = 0;
	for (int f = 0; f < numPulses; f++)
	{
		frameCounts = detectionsJc[f + 1] - detectionsJc[f];

		Detection* frameDetections = new Detection[frameCounts]();
		for (int c = 0; c < frameCounts; c++)
		{
			frameDetections[c] = { detectionsIr[cIdx], (int)detectionsPr[cIdx] };
			//mexPrintf("%d, %d - %d\n", cIdx, detectionsIr[cIdx], (int)detectionsPr[cIdx]);
			cIdx++;
		}
		sort(frameDetections, frameDetections + frameCounts);

		size_t pr1 = 0;
		lctJc[f + 1] = lctJc[f];
		int di, dj, dt;
		Detection d1, d2;
		size_t pr;
		double dLct = 0;
        double pct_xyt = 0;
		for (size_t pr2 = 0; pr2 < frameCounts; pr2++)
		{
			d1 = frameDetections[pr1];
			d2 = frameDetections[pr2];
			dt = d2.t - d1.t;
			//mexPrintf("d1 = (%d, %d), d2 = (%d, %d)\n", d1.idx, d1.t, d2.idx, d2.t);
			//mexEvalString("drawnow;");
			while (dt > ctMaxDelay && pr1 < pr2)
			{
				pr1++;
				d1 = frameDetections[pr1];
				dt = d2.t - d1.t;
				//mexPrintf("d1 = (%d, %d)\n", d1.idx, d1.t);
				//mexEvalString("drawnow;");
			}
			pr = pr1;
			while (dt >= 0)
			{
				di = (d2.idx % nrows) - (d1.idx % nrows);
				dj = (d2.idx / nrows) - (d1.idx / nrows);
                pct_xyt = pct[ctMaxDelay*(pctDims[1] * (nrows + di - 1) + ncols + dj - 1) + dt];
                //if (dt == 0)
                //{
                //    pct_xyt = pct_xyt * 0.5;
                //}
				dLct += pct_xyt;
				//mexPrintf("pct[%d][%d][%d] = %f. dLct = %f\n", di, dj, dt, pct[ctMaxDelay*(pctDims[1] * (nrows + di - 1) + ncols + dj - 1) + dt], dLct);
				pr++;
				//mexPrintf("%d\n", pr);
				//mexEvalString("drawnow;");
				if (pr >= frameCounts)
				{
					break;
				}
				d1 = frameDetections[pr];
				dt = d2.t - d1.t;
				//mexPrintf("d1 = (%d, %d)\n", d1.idx, d1.t);
				//mexEvalString("drawnow;");
			}
			if (dLct != 0)
			{
				(lctJc[f + 1])++;
				lctIr[prIdx] = d2.idx;
				lctPr[prIdx] = dLct;

				prIdx++;
				dLct = 0;
			}
		}

		delete[] frameDetections;
	}

}

