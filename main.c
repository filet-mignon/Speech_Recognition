/*
 * main.c
 *
 *  Created on: Nov 7, 2016
 *      Author: Mignon Huang and Rosemary Lee
 *
 *      Framelength of 64 ms, 1024 samples at 16 kHz
 *
 */

#define FRAMELENGTH 1024
#define K 513
#define SAMPLELEN 16000
#define PI 3.141592654
#define N 26

#include <math.h>
#include "L138_LCDK_aic3106_init.h"
//#include "L138_LCDK_switch_led.h"
#include "fft.h"
#include "stdio.h"

int sum;
int counter;
int playback;
int reset;
int yesno;
int input;
int sel_index;
int record; 	//-1 no recording available
					//1 currently recording
					//0 recording processed
char* sound;
char* sound_0 = "ah";
char* sound_1 = "eh";
char* sound_2 = "ee";
char* sound_3 = "oo";
char* sound_4 = "xx";

int16_t left_sample;
int16_t sample[SAMPLELEN];
double hamming[FRAMELENGTH];
int PHI[N+2];
int MEL[N+2];
double MU[N+2];
double Y[N];
double X[N];
double MFCC[N/2];
double DATABASE[4][14];				//ah = 0
									//eh = 1
									//ee = 2
									//oo = 3
//double ARR[FRAMELENGTH];

COMPLEX twiddle[FRAMELENGTH];
COMPLEX frame[FRAMELENGTH];

double mag[FRAMELENGTH];

char* state_0 = "No active task";
char* state_1 = "Recording";
char* state_2 = "Analyzing";
char* state_3 = "Is this sound correct?";
char* state_4 = "Please select correct sound";
char* state_5 = "Database updated";

char* state_message;


interrupt void interrupt4(void)
{
	if (reset == 1){
		record = -1;
		playback = 0;
		counter = -1;
		input = -1;
		reset = 0;
		sound = sound_4;
		state_message = state_0;
	}
	left_sample = 0;
	sum = record + playback;
	if (sum > 1){
		record = -1;
		playback = 0;
		counter = -1;
		input = -1;
		reset = 0;
	}
	else if (sum == 1){
		if (record == 1){
			state_message = state_1;
			counter++;
			//Collects a 1 second sample
			if(counter > 100 && counter< (100+SAMPLELEN))
			{
//				LCDK_LED_on(4);
				left_sample = input_left_sample();
				sample[(counter-101)] = left_sample;
			}
			else if(counter >= (100+SAMPLELEN)){
				record = 0;
				counter = -1;
//				LCDK_LED_off(4);
			}
		}
		if (playback == 1){
			counter++;
			if(counter < SAMPLELEN){
				left_sample = sample[counter];
			}
			else if (counter >= SAMPLELEN){
				playback = 0;
				counter = -1;
			}
		}
	}

	output_left_sample(left_sample);
	return;
}

void main()
{
	sum = 0;
	counter = 0;
	playback = 0;
	record = -1;
	input = -1;
	yesno = -1;

//	LCDK_LED_init();

	int i, m, k;

	//Set up FFT twiddle factors
	for(i = 0; i < FRAMELENGTH; i++)
	{
		twiddle[i].real = cos(PI*i/FRAMELENGTH);
		twiddle[i].imag = -sin(PI*i/FRAMELENGTH);
		hamming[i] = 0.46 - 0.54*cos(2*PI*i/(FRAMELENGTH-1));
	}

	for(i = 0; i < (N+2); i++){
		MEL[i] = 345 + (2840-345)*i/27;
		MU[i] = (pow(10.0, (MEL[i]/2595.0)) - 1)*700.0;
		PHI[i] = (MU[i]/8000)*K;
		//printf("MEL %d, MU %f, PHI %d \n", MEL[i], MU[i], PHI[i]);
	}

	for(i = 0; i < N; i++)
		Y[i] = 0;
	for(i = 0; i < (N/2); i++)
		MFCC[i] = 0;

	for(m = 0; m < 4; m++){
		for(k = 0; k < 14; k++)
			DATABASE[m][k] = 0;
	}

	L138_initialise_intr(FS_16000_HZ,ADC_GAIN_24DB,DAC_ATTEN_0DB,LCDK_MIC_INPUT);

	while(1){
		if(record == 0){
			state_message = state_2;
			for(i = 0; i < FRAMELENGTH; i++){
				frame[i].real = sample[i+100] * hamming[i];
				frame[i].imag = 0.0;
			}

			fft(frame, FRAMELENGTH, twiddle);

			for(i = 0; i < FRAMELENGTH; i++){
				mag[i] = (frame[i].real*frame[i].real + frame[i].imag*frame[i].imag);
			}

			for(m = 1; m < N+1; m++){
				for(k = PHI[m-1]; k < PHI[m]; k++){
					Y[m] += (k-PHI[m-1])/(PHI[m]-PHI[m-1]);
				}
				for(k = PHI[m]; k < PHI[m+1]; k++){
					Y[m] += (PHI[m+1]-k)/(PHI[m+1]-PHI[m]);
				}
				X[m] = log10(Y[m]);
			}

			for(k = 1; k <= 13; k++){
				for(m = 1; m <= N; m++){
					MFCC[k] += X[m]*cos((m-0.5)*k*PI/26);
				}
			}

			double minDist = 30000;
			int min_index = -1;
			double distance;
			for(i = 0; i < 4; i++){
				//Compare to existing values
				//find ID with minimum distance
				for(k = 0; k <13; k++){
					distance += (MFCC[k] - DATABASE[i][k])*(MFCC[k] - DATABASE[i][k]);
				}
				distance = sqrt(distance);

				if(distance < minDist){
					minDist = distance;
					min_index = i;
				}
			}

			switch(min_index){
			//Using min_index, set the sound identifier
			case 0:
				sound = sound_0;
				break;
			case 1:
				sound = sound_1;
				break;
			case 2:
				sound = sound_2;
				break;
			case 3:
				sound = sound_3;
				break;
			default:
				sound = sound_4;
				break;
			}

			state_message = state_3;

			while(input == -1){}
			if(input == 1){
				for(i = 0; i < 13; i++){
					DATABASE[min_index][i] =
							(DATABASE[min_index][13]*DATABASE[min_index][i] + MFCC[i]);
					DATABASE[min_index][i] /= (DATABASE[min_index][13]+1);
				}
				DATABASE[min_index][13]++;
				state_message = state_5;
				input = -1;
				record = -1;
			}
			if(input == 0){
				state_message = state_4;
				//Ask for correct sound identifier and update data for that sound
				while(sel_index == -1){}
				switch(sel_index){
				//Using min_index, set the sound identifier
				case 0:
					sound = sound_0;
					break;
				case 1:
					sound = sound_1;
					break;
				case 2:
					sound = sound_2;
					break;
				case 3:
					sound = sound_3;
					break;
				default:
					sound = sound_4;
					break;
				}

				sound[2] = '\0';
				for(i = 0; i < 13; i++){
					DATABASE[sel_index][i] =
							(DATABASE[sel_index][13]*DATABASE[sel_index][i] + MFCC[i]);
					DATABASE[sel_index][i] /= (DATABASE[sel_index][13]+1);
				}
				DATABASE[sel_index][13]++;
				sel_index = -1;
				state_message = state_5;
				input = -1;
				record = -1;
			}
		}
		else
			state_message = state_0;
	}

}
