#include <random>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include "mex.h"
#include "matrix.h"

using namespace std;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dis(0.0, 1.0);
double randnum = dis(gen);

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
		return d1.t > d2.t;
	}
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Load data
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("simulate_crosstalk:invalidNumOutputs", "One output argument required.");
	}
	if (nrhs != 5) {
		mexErrMsgIdAndTxt("simulate_crosstalk:invalidNumInputs", "Five input arguments required.");
	}
	if (!mxIsSparse(prhs[0])) {
		mexErrMsgIdAndTxt("simulate_crosstalk:detectionsNotSparse", "Input detections must be sparse.");
	}

	const int numPixels = (int)mxGetM(prhs[0]);
	const int numPulses = (int)mxGetN(prhs[0]);
	const mwIndex* detectionsJc = mxGetJc(prhs[0]);
	if (detectionsJc[numPulses] == 0) {
		plhs[0] = mxCreateSparse(numPixels, numPulses, 0, mxREAL);
	}
	double* detectionsPr = mxGetPr(prhs[0]);
	mwIndex* detectionsIr = mxGetIr(prhs[0]);

	//mexPrintf("numPixels = %d \n", numPixels);
	//mexPrintf("numPulses = %d \n", numPulses);


	if (mxGetNumberOfDimensions(prhs[1]) != 3) {
		mexErrMsgIdAndTxt("simulate_crosstalk:pctNot3D", "Crosstalk probability distribution should be 3-D.");
	}

	const mwSize* pctDims = mxGetDimensions(prhs[1]);
	//const int ctMaxDelay = (int)pctDims[2];
	//const mwSize nrows = (pctDims[0] + 1) / 2;
	//const mwSize ncols = (pctDims[1] + 1) / 2;
	//mexPrintf("nrows = %d \n ncols = %d, ctMaxDelay = %d \n", nrows, ncols, ctMaxDelay);
	//double*** pct = new double**[ctMaxDelay];
	//double* pctPr = mxGetPr(prhs[1]);
	size_t counter = 0;
	//for (size_t t = 0; t < ctMaxDelay; t++)
	//{
	//	pct[t] = new double*[pctDims[0]];

	//	for (size_t x = 0; x < pctDims[0]; x++)
	//	{
	//		pct[t][x] = new double[pctDims[1]];
	//		for (size_t y = 0; y < pctDims[1]; y++)
	//		{
	//			pct[t][x][y] = pctPr[counter++];
	//			//mexPrintf("%d, %f\n",counter, pct[t][x][y]);
	//		}
	//	}
	//}

	const int ctMaxDelay = (int)pctDims[0];
	const mwSize nrows = (pctDims[2] + 1) / 2;
	const mwSize ncols = (pctDims[1] + 1) / 2;
	const double* pct = mxGetPr(prhs[1]);
	if (nrows*ncols != numPixels) {
		mexErrMsgIdAndTxt("simulate_crosstalk:pctWrongDims", "Crosstalk probability matrix should be size [ctMaxDelay, ncols, nrows].");
	}

	//mexPrintf("nrows = %d \n ncols = %d, ctMaxDelay = %d \n", nrows, ncols, ctMaxDelay);


	if (mxGetNumberOfElements(prhs[2]) != 1) {
		mexErrMsgIdAndTxt("simulate_crosstalk:maxTimebinNotScalar", "Max timebin should be scalar integer.");
	}
	const int maxTimebin = (int)mxGetScalar(prhs[2]);

	if (mxGetNumberOfDimensions(prhs[3]) != 2) {
		mexErrMsgIdAndTxt("simulate_crosstalk:qNot2D", "Quantum efficiency should be 2D.");
	}
	const double* q = mxGetPr(prhs[3]);

	if (mxGetNumberOfElements(prhs[4]) != numPixels || !mxIsLogical(prhs[4])) {
		mexErrMsgIdAndTxt("simulate_crosstalk:deadpixNot2D", "Deadpix should be logical with numPixel elements.");
	}
	const mxLogical* deadpix = (mxLogical*)mxGetData(prhs[4]);

	// Crosstalk 
	size_t startIdx, frameCounts, finalFrameCounts;
	size_t cIdx, pixIdx;
	int detectionTime;
	int i, j;
	int t, tPixel;
	mwIndex* finalJc = (mwIndex*)mxCalloc(numPulses + 1, sizeof(mwIndex));
	int** finalIrArrays = new int*[numPulses];
	int** finalPrArrays = new int*[numPulses];
	priority_queue<Detection> frameDetections;
	for (int f = 0; f < numPulses; f++)
	{
		finalFrameCounts = 0;
		startIdx = detectionsJc[f];
		frameCounts = detectionsJc[f + 1] - startIdx;
		//mexPrintf("frame %d: %d detections\n", f, frameCounts);
		//mexEvalString("drawnow;");

		// Get full matrix (1-d array) of detections in current frame.
		// Also populate vector of detection indices - to be sorted.
		int* frameDetectionsFull = new int[numPixels]();
		for (size_t c = 0; c < frameCounts; c++) 
		{
			cIdx = startIdx + c;
			pixIdx = detectionsIr[cIdx];
			detectionTime = (int)detectionsPr[cIdx];
			//mexPrintf("enqueuing %d, %d \n", pixIdx, detectionTime);
			
			frameDetections.emplace(pixIdx, detectionTime);
		}

		// Generate crosstalk for each detection in frame, in order of increasing time.
		//for (auto &dPrimary : frameDetections)
		Detection dPrimary;
		//priority_queue<Detection> ctFrame;
		//int* ctFrameFull = new int[numPixels];
		while (!frameDetections.empty())
		{
			dPrimary = frameDetections.top();
			frameDetections.pop();

			if (frameDetectionsFull[dPrimary.idx] == 0 && deadpix[dPrimary.idx])
			{
				// at this point, it is safe to store dPrimary.t because it cannot get blocked. 
				// All previous detections have already been processed, because we go in 
				// sequential order and crosstalk can only happen later in time. (??)
				frameDetectionsFull[dPrimary.idx] = dPrimary.t;

				finalFrameCounts++;

				i = dPrimary.idx % nrows;
				j = dPrimary.idx / nrows;
				t = dPrimary.t;
				//mexPrintf("i = %d, j = %d, t = %d \n", i, j, t);
				//mexEvalString("drawnow;");

				//const double* pixPct = pct;
				//for (int ii = 0; ii < 25; ii++) {
				//	mexPrintf("%d, %d: %f, %f\n", &(pct[ii*3+1]), pixPct+1, pct[ii*3+1], pixPct[1]);
				//	pixPct = pixPct + ctMaxDelay;
				//}
				//mexPrintf("pct[41] = %f, Ptr = %d\n", pct[41], pixPct);
				
				double pctq;
				const double* pixPct = pct + ctMaxDelay*(pctDims[1] * (nrows - i - 1) + ncols - j - 1);
				for (int di = -i; di < (int)nrows - i; di++)
				{
					for (int dj = -j; dj < (int)ncols - j; dj++)
					{
						pixIdx = i + di + (j + dj)*nrows;
						tPixel = frameDetectionsFull[pixIdx];
						//mexPrintf("di = %d, dj = %d, pixIdx = %d, tPixel = %d \n", di, dj, pixIdx, tPixel);
						//mexEvalString("drawnow;");
						if (tPixel == 0) 
						{
							for (int dt = 0; dt < min(ctMaxDelay, maxTimebin - t + 1); dt++)
							{
								//pctq = pct[dt][dj + ncols - 1][di + nrows - 1];
								//pctq = pct[di + nrows - 1][dj + ncols - 1][dt] * q[pixIdx];
								pctq = pixPct[dt] * q[pixIdx];
								//pctq = pixPct[dt];

								/*mexPrintf("%d pct[%d][%d][%d] = %f\n", pixPct+dt, di, dj, dt, pixPct[dt]);
								mexEvalString("drawnow;");*/
								if (dis(gen) < pctq)
								{
									//mexprintf("crosstalk at pixidx %d, timebin %d \n", i + di + (j + dj)*nrows, t + dt);
									
									//framedetectionsfull[pixIdx] = t + dt;
									//finalframecounts++;
									frameDetections.emplace(pixIdx, t + dt);
									break;
								}
							}
						}

						pixPct = pixPct + ctMaxDelay;
					}
					pixPct = pixPct + ctMaxDelay*(ncols - 1);
				}
			}
			
		}

		// Populate sparse arrays for current frame.
		finalJc[f + 1] = finalJc[f] + finalFrameCounts;
		int* frameFinalIr = new int[finalFrameCounts]();
		int* frameFinalPr = new int[finalFrameCounts]();
		counter = 0;
		for (size_t pIdx = 0; pIdx < numPixels; pIdx++)
		{
			tPixel = frameDetectionsFull[pIdx];
			//mexPrintf("pIdx = %d, tPixel = %d, ", pIdx, tPixel);
			//mexEvalString("drawnow;");
			if (tPixel != 0)
			{
				//mexPrintf("!= 0");
				//mexEvalString("drawnow;");
				frameFinalIr[counter] = pIdx;
				frameFinalPr[counter] = tPixel;
				counter++;
			}
			//mexPrintf("\n");
			//mexEvalString("drawnow;");
		}
		finalIrArrays[f] = frameFinalIr;
		finalPrArrays[f] = frameFinalPr;


		delete[] frameDetectionsFull;
	}

	plhs[0] = mxCreateSparse(numPixels, numPulses, finalJc[numPulses], mxREAL);
	double* finalPr = mxGetPr(plhs[0]);
	mwIndex* finalIr = mxGetIr(plhs[0]);
	mxSetJc(plhs[0], finalJc);
	counter = 0;
	for (size_t f = 0; f < numPulses; f++)
	{
		frameCounts = finalJc[f + 1] - finalJc[f];
		for (size_t c = 0; c < frameCounts; c++)
		{
			finalIr[counter] = finalIrArrays[f][c];
			finalPr[counter] = (double)finalPrArrays[f][c];
			counter++;
			//mexPrintf("finalIrArrays[%d][%d] = %d\n", f, c, finalIrArrays[f][c]);
			//mexEvalString("drawnow;");
			//mexPrintf("finalPrArrays[%d][%d] = %f\n", f, c, (double)finalPrArrays[f][c]);
		}
		delete[] finalIrArrays[f];
		delete[] finalPrArrays[f];
	}

	delete[] finalIrArrays;
	delete[] finalPrArrays;
	//for (size_t f = 0; f < numPulses+1; f++) {
	//	mexPrintf("finalJc[%d] = %d\n", f, finalJc[f]);
	//	mexEvalString("drawnow;");
	//}

}
