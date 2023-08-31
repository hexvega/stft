
#include "stft.h"

// Standard lib
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/*
@brief:
		STFT所需的数据加载内存分配.
@oaram:
		ptrRawArray:	  批量原始数据指针;
		rawArraylength	  批量原始数据长度；
@return:
		ctxt     上下文
		
*/
void* stft_init(float* ptrRawArray, int rawArraylength) {

	STFT_CONTEXT* ctxt = NULL;
	ctxt = calloc(1, sizeof(kiss_fft_cpx**) + sizeof(float**) * 2 + sizeof(float*) * 2 + sizeof(int*) * 3);
	if (!ctxt) return NULL;

	ctxt->_ptrRawArrayData = ptrRawArray;
	ctxt->_rawArraylength = rawArraylength;
	ctxt->_proMatrixHeight = N_FFT / 2 + 1;
	ctxt->_proMatrixWidth = rawArraylength / N_HOP + 1;
	////createHanningWindow(ctxt->_hannwin);
	ctxt->_hannwin = (float*)calloc(N_FFT, sizeof(float*));
	createHanningWindow(ctxt->_hannwin);

	ctxt->_complexMatrix = (kiss_fft_cpx**)calloc(ctxt->_proMatrixHeight, sizeof(kiss_fft_cpx*));
	for (int i = 0;i < ctxt->_proMatrixHeight; i++) {
		ctxt->_complexMatrix[i] = (kiss_fft_cpx*)calloc(ctxt->_proMatrixWidth, sizeof(kiss_fft_cpx));
	}

	ctxt->_floatMatrix = (float**)calloc(ctxt->_proMatrixHeight, sizeof(float*));
	for (int i = 0;i < ctxt->_proMatrixHeight; i++) {
		ctxt->_floatMatrix[i] = (float*)calloc(ctxt->_proMatrixWidth, sizeof(float));
	}

	ctxt->_dbMatrix = (float**)calloc(ctxt->_proMatrixHeight, sizeof(float*));
	for (int i = 0;i < ctxt->_proMatrixHeight; i++) {
		ctxt->_dbMatrix[i] = (float*)calloc(ctxt->_proMatrixWidth, sizeof(float));
	}
	return ctxt;
}

//    释放复数矩阵堆内存
void freeComplexMatrix(kiss_fft_cpx** matrixInput, int height) {
	if (matrixInput != NULL) {
		for (int i = 0; i < height; i++)
			free(matrixInput[i]);
		free(matrixInput);
	}
}

//  释放复数模矩阵堆内存
void freeFloatMatrix(float** floatMatrix, int height) {
	if (floatMatrix != NULL) {
		for (int i = 0; i < height; i++)
			free(floatMatrix[i]);
		free(floatMatrix);
	}
}

void freeDbMatrix(float** dbMatrix, int height) {
	if (dbMatrix != NULL) {
		for (int i = 0; i < height; i++)
			free(dbMatrix[i]);
		free(dbMatrix);
	}
}

void freeHanningWindow(float* hannwin) {
	if (hannwin != NULL) {
		free(hannwin);
	}
}


/*
@brief:
		256点汉明窗
*/
static void createHanningWindow(float* win)
{
	for (int i = 0; i < N_FFT; i++) { win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / N_FFT); }
}

/*
@brief:
		STFT_PAD 使用reflect模式
@oaram:
		in:	  批量原始数据指针;
		out	  pad延伸后的数据指针；
		arrayLength   输入的原始数据长度
*/
static void pad(float* in, float* out, int arrayLength) {
	static int pad_len = N_FFT / 2;

	memcpy(&out[pad_len], in, arrayLength * sizeof(float));

	for (int i = 0; i < pad_len; ++i) {
		out[i] = in[pad_len - i];
	}

	for (int i = pad_len; i < N_FFT; ++i) {
		out[i + arrayLength] = in[arrayLength - 2 - i + pad_len];
	}
}

/*
@brief:
		复数取模
*/
static float INL_MOD_COMPLEX(float r, float i)
{
	float temp;
	temp = r * r + i * i;
	temp = sqrt(temp);
	return temp;
}

/*
@brief:
		STFT变换
@oaram:
		_ptrRawArrayData：   原始数据指针   
		_complexMatrix:	     STFT输出的复数矩阵指针;
*/
void stft(STFT_CONTEXT* ctxt) {
	kiss_fftr_cfg cfg = kiss_fftr_alloc(N_FFT, 0, NULL, NULL);
	kiss_fft_cpx freqdata[N_FFT / 2 + 1];

	float* arr = (float*)calloc(ctxt->_rawArraylength + N_FFT, sizeof(float));
	pad(ctxt->_ptrRawArrayData, arr, ctxt->_rawArraylength);

	for (int i = 0; i < ctxt->_proMatrixWidth; ++i) {
		float tempDataRaw[N_FFT];
		for (int j = 0; j < N_FFT; j++) {
			tempDataRaw[j] = arr[i * N_HOP + j] * ctxt->_hannwin[j];
		}
		kiss_fftr(cfg, tempDataRaw, freqdata);
		for (int k = 0; k < ctxt->_proMatrixHeight; k++) {
			ctxt->_complexMatrix[k][i] = freqdata[k];
		}
	}
	//及时释放汉明窗堆内存。
	freeHanningWindow(ctxt->_hannwin);
}


void complexAbs(kiss_fft_cpx** matrixInput, float** matrixOutput, int height, int width) {

	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			matrixOutput[i][j] = INL_MOD_COMPLEX(matrixInput[i][j].r, matrixInput[i][j].i);
		}
	}

	freeComplexMatrix(matrixInput, height);
}

/*
@brief:
		获取矩阵数据中最大值，作为功率谱转换的形参
@oaram:
		matrixInput:	  STFT变换后数据矩阵复指针;
		height            matrixInput矩阵高
		width             matrixInput矩阵宽
@return:
		maxValue     矩阵中最大值。
*/
float getMaxVaule(float** matrixInput, int height, int width) {
	float maxValue = 0;
	for (int i = 0;i < height;i++) {
		for (int j = 0;j < width;j++) {
			float sss = matrixInput[i][j];
			maxValue = (matrixInput[i][j] > maxValue) ? matrixInput[i][j] : maxValue;
		}
	}
	return maxValue;
}

/*
@brief:
		功率谱转换.
@oaram:
		matrixInput:	  STFT变换后数据矩阵复指针;
		matrixOutput	  转换成DB分量的矩阵复指针；
		height            matrixInput矩阵高
		width             matrixInput矩阵宽
		refValue          log计算所需的参考值。默认1.0，当前取matrixInput最大值。
*/
void power2db(float** matrixInput, float** matrixOutput, int height, int width, float refValue) {

	float maxMatrixValue = 0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			float magnitude = (matrixInput[i][j] < 0) ? matrixInput[i][j] * -1 : matrixInput[i][j];
			float logSpec = (magnitude > 1e-10) ? 10.0 * log10(magnitude) : 10.0 * (1e-10);
			float logSub = (refValue > (1e-10)) ? 10.0 * log10(refValue) : 10.0 * (1e-10);
			matrixOutput[i][j] = logSpec - logSub;
			maxMatrixValue = (matrixOutput[i][j] > maxMatrixValue) ? matrixOutput[i][j] : maxMatrixValue;
		}
	}
	// 设置80.0为 top_db
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (matrixOutput[i][j] < (maxMatrixValue - TOP_DB)) {
				matrixOutput[i][j] = (maxMatrixValue - TOP_DB);
			}
		}
	}

	freeFloatMatrix(matrixInput, height);
}

