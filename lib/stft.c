
#include "stft.h"

// Standard lib
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/*
@brief:
		STFT��������ݼ����ڴ����.
@oaram:
		ptrRawArray:	  ����ԭʼ����ָ��;
		rawArraylength	  ����ԭʼ���ݳ��ȣ�
@return:
		ctxt     ������
		
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

//    �ͷŸ���������ڴ�
void freeComplexMatrix(kiss_fft_cpx** matrixInput, int height) {
	if (matrixInput != NULL) {
		for (int i = 0; i < height; i++)
			free(matrixInput[i]);
		free(matrixInput);
	}
}

//  �ͷŸ���ģ������ڴ�
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
		256�㺺����
*/
static void createHanningWindow(float* win)
{
	for (int i = 0; i < N_FFT; i++) { win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / N_FFT); }
}

/*
@brief:
		STFT_PAD ʹ��reflectģʽ
@oaram:
		in:	  ����ԭʼ����ָ��;
		out	  pad����������ָ�룻
		arrayLength   �����ԭʼ���ݳ���
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
		����ȡģ
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
		STFT�任
@oaram:
		_ptrRawArrayData��   ԭʼ����ָ��   
		_complexMatrix:	     STFT����ĸ�������ָ��;
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
	//��ʱ�ͷź��������ڴ档
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
		��ȡ�������������ֵ����Ϊ������ת�����β�
@oaram:
		matrixInput:	  STFT�任�����ݾ���ָ��;
		height            matrixInput�����
		width             matrixInput�����
@return:
		maxValue     ���������ֵ��
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
		������ת��.
@oaram:
		matrixInput:	  STFT�任�����ݾ���ָ��;
		matrixOutput	  ת����DB�����ľ���ָ�룻
		height            matrixInput�����
		width             matrixInput�����
		refValue          log��������Ĳο�ֵ��Ĭ��1.0����ǰȡmatrixInput���ֵ��
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
	// ����80.0Ϊ top_db
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (matrixOutput[i][j] < (maxMatrixValue - TOP_DB)) {
				matrixOutput[i][j] = (maxMatrixValue - TOP_DB);
			}
		}
	}

	freeFloatMatrix(matrixInput, height);
}

