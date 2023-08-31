/*
	��ʾ��Ҫ˵��������ʵ��ʹ��csv����������Ϊ���룬���һ��bool�����

	�㷨ִ�й�����ѭPythonʾ���㷨�����gsensor_vad.py��

	��������������Pythonʾ���㷨���һ�¡�
*/
#include "testData/csvdata.h"
#include "lib/stft.h"

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

int main() {

	STFT_CONTEXT* context = NULL;

	static float* ptrDataSlice;

	ptrDataSlice = &accs_423_left[0][0];
	int arrayLength = sizeof(accs_423_left) / sizeof(float) / DATA_RAW_DIM;

	float* fusionArray = (float*)calloc(arrayLength, sizeof(float*));
	for (int i = 0;i < arrayLength; i++) {
		float sum = 0.0f;
		sum += *(ptrDataSlice++);
		sum += *(ptrDataSlice++);
		sum += *(ptrDataSlice++);
		fusionArray[i] = sum;
	}

	context = stft_init(fusionArray, arrayLength);
	
	stft(context);

	complexAbs(context->_complexMatrix, context->_floatMatrix, context->_proMatrixHeight, context->_proMatrixWidth);

	float maxValue = getMaxVaule(context->_floatMatrix, context->_proMatrixHeight, context->_proMatrixWidth);

	power2db(context->_floatMatrix, context->_dbMatrix, context->_proMatrixHeight, context->_proMatrixWidth, maxValue);

	for (int i = 0; i < context->_proMatrixWidth; i++) {
		float sum1 = 0.0f;
		float sum2 = 0.0f;
		float rate = 0.0f;
		for (int j = 15;j < 45;j++) {
			sum1 += context->_dbMatrix[j][i];
		}
		for (int j = 92;j < 122;j++) {
			sum2 += context->_dbMatrix[j][i];
		}
		rate = sum2 / sum1;
		if (rate > 1.05) {
			printf("step %d :   1\r\n", i);
		}
		else {
			printf("step %d :   0\r\n", i);
		}
	}

	freeDbMatrix(context->_dbMatrix, context->_proMatrixHeight);
	return 0;
}