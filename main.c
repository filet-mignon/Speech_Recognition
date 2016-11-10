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

#include <math.h>
#include "L138_LCDK_aic3106_init.h"
//#include "L138_LCDK_switch_led.h"
#include "fft.h"

int sum;
int counter;
int playback;
int record; 	//-1 no recording available
					//1 currently recording
					//0 recording processed

int16_t left_sample;
int16_t sample[SAMPLELEN];
double hamming[FRAMELENGTH];
//double ARR[FRAMELENGTH];

COMPLEX twiddle[FRAMELENGTH];
COMPLEX frame[FRAMELENGTH];

double mag[FRAMELENGTH];

interrupt void interrupt4(void)
{
	left_sample = 0;
	sum = record + playback;
	if (sum > 1){
		record = 0;
		playback = 0;
		counter = -1;
	}
	else if (sum == 1){
		if (record == 1){
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

//	LCDK_LED_init();

	int i;

	//Set up FFT twiddle factors
	for(i = 0; i < FRAMELENGTH; i++)
	{
		twiddle[i].real = cos(PI*i/FRAMELENGTH);
		twiddle[i].imag = -sin(PI*i/FRAMELENGTH);
	}

	for(i=0; i < FRAMELENGTH; i++){
		hamming[i] = 0.46 - 0.54*cos(2*PI*i/(FRAMELENGTH-1));
	}

	L138_initialise_intr(FS_16000_HZ,ADC_GAIN_24DB,DAC_ATTEN_0DB,LCDK_MIC_INPUT);

	while(1){
		if(record == 0){
			for(i = 0; i < FRAMELENGTH; i++){
				frame[i].real = sample[i+100] * hamming[i];
				frame[i].imag = 0.0;
			}

			fft(frame, FRAMELENGTH, twiddle);

			for(i = 0; i < FRAMELENGTH; i++){
				mag[i] = sqrt(frame[i].real*frame[i].real + frame[i].imag*frame[i].imag);
			}

			record = -1;
		}
	}

}
