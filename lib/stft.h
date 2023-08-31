#pragma once
#ifndef __STFT_H__
#define __STFT_H__
// Independent lib
#include "kiss_fft.h"
#include "kiss_fftr.h"
// Standard lib
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define M_PI 3.14159265358979323846
#define DATA_RAW_DIM			3
#define N_FFT 256
#define N_HOP (N_FFT/4)
#define TOP_DB 80.0

/*
@oaram:
		_ptrRawArrayData��   ԭʼ����ָ��
		_rawArraylength:	 ԭʼ���ݳ���
		_proMatrixHeight:	 STFT����������ݸ�
		_proMatrixWidth:	 STFT����������ݿ�
		_hannwin             ������
		_complexMatrix       ָ��STFT����ĸ�������
		_floatMatrix         ָ��_complexMatrix�ĸ���ģ
		_dbMatrix            ָ�����ױ任�����
*/
typedef struct {
	float* _ptrRawArrayData;
	int    _rawArraylength;
	int    _proMatrixHeight;
	int    _proMatrixWidth;
	float* _hannwin;
	kiss_fft_cpx** _complexMatrix;
	float** _floatMatrix;
	float** _dbMatrix;
} STFT_CONTEXT;

void* stft_init(float* ptrRawArray, int rawArraylength);
void power2db(float** matrixInput, float** matrixOutput, int height, int width, float refValue);
static void createHanningWindow(float* win);
static void pad(float* in, float* out, int arrayLength);
static float INL_MOD_COMPLEX(float r, float i);
void stft(STFT_CONTEXT* ctxt);
void complexAbs(kiss_fft_cpx** matrixInput, float** matrixOutput, int height, int width);
float getMaxVaule(float** matrixInput, int height, int width);

#endif
