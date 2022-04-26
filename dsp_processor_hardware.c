/*
In this first experiment, we will add source code that blinks LED #0 at a rate of about 2.5 times per second using the LED module of the DSK6713 Board Support Library.  The example also reads the state of DIP switch #3 and lights LED #3 if the switch is depressed or turns it off if the switch is not depressed. The purpose of this experiment is to demonstrate basic Board Support Library usage as well as provide a project base for your own code.  The BSL is divided into several modules, each of which has its own include file.  The file dsk6713.h must be included in every program that uses the BSL.  This example also includes  dsk6713_led.h and dsk6713_dip.h because it uses the LED and DIP modules. 
*/
/*  ======== led.c ========*/
#include <dsk6713.h>
#include <dsk6713_led.h>
#include <dsk6713_dip.h>
void main()
{
/* Initialize the board support library, must be first BSL call */ DSK6713_init();

   		 /* Initialize the LED and DIP switch modules of the BSL */
    		DSK6713_LED_init();
    		DSK6713_DIP_init();

    		while(1)
    		{
        		/* Toggle LED #0 */
        		DSK6713_LED_toggle(0);

/* Check DIP switch #3 and light LED #3 accordingly, */
//0 = switch pressed 
        		if (DSK6713_DIP_get(3) == 0)
            		/* Switch pressed, turn LED #3 on */
            			DSK6713_LED_on(3);
        		else
            		/* Switch not pressed, turn LED #3 off */
            			DSK6713_LED_off(3);

        		/* Spin in a software delay loop for about 200ms */
        		DSK6713_waitusec(200000);
    		}
}









/*
In the code, right click on DSK6713_AIC23_DEFAULTCONFIG and click Open Declaration. This will open the header file which contains detailed information about how the codec is configured. (eg: how to decrease headphone volume?). The codec is started by calling the BSL function DSK6713_AIC23_openCodec(). 
Each iteration of the infinite loop generates a sample and writes it to the codec. Note that the float value is converted to int and the upper 16 bits are set to 0 before outputting. This is because the BSL function that configures the codec sets McBSP1 to send and receive 32-bit words, with the left sample in the upper 16 bits and right sample in the lower 16 bits. The 16-bit samples are in signed 2s complement form. Since the upper 16 bits of out_sample are set to 0,  the tone will be heard on the right channel only. If we want the output to be on the left channel, we can use the statement out_sample = (int)sample << 16; instead, which puts the 16-bit value in the top half, and sets the lower 16 bits to 0s. We can also pack two 16-bit samples in out_sample for output on both L and R channels. 
The function DSK6713_AIC23_write() is used to write a pair of samples to the DAC. The function uses polling to write samples and returns 0 if codec is not ready and returns 1 if write is successful. The while loop continues till write is successful. Build and Debug the program. Connect your headphones to the headphone and listen to the tone(beware of volume). You should hear the tone on R channel only. Using bitwise operators in C, try to output the tone on both channels. Change frequency F0 to 500 Hz and listen. Next, change F0 to 7500 Hz and listen. Explain what you observe. 
*/

#include <math.h>
#include <dsk6713.h>
#include <dsk6713_aic23.h>
int main()
{
    float Fs = 8000.;
    float F0 = 1000.;
    float pi = 3.141592653589;
    float theta = 0.;
    float delta = 2. * pi * F0 / Fs; // increment for theta
    float sample;
    unsigned out_sample;
    /* Initialize the board support library, must be called first */
    DSK6713_init();

    DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;
    DSK6713_AIC23_CodecHandle hCodec;
    /* Start the codec */
    hCodec = DSK6713_AIC23_openCodec(0, &config);

    /* Change the sampling rate to 16 kHz */
    DSK6713_AIC23_setFreq(hCodec, DSK6713_AIC23_FREQ_8KHZ);

    for (;;)
    { /* Infinite loop */
        sample = 15000.0 * sin(theta); /* Scale for DAC */
        out_sample = (int)sample & 0x0000ffff; // Put in lower half (R) 

        /* Poll XRDY bit until true, then write to DXR */
        while (!DSK6713_AIC23_write(hCodec, out_sample))
            ;
        theta += delta;
        if (theta > 2 * pi)
            theta -= 2 * pi;

    }
}














/*
A structure COMPLEX array is used to store the real and imaginary parts of x[n] and X[k].  The computations are performed in-place with the input array over-written by the output array. The program computes the 64-point DFT on the 64 samples of a 1KHz signal sampled at 8000Hz. 

Procedure: Insert a breakpoint in the code on the line calling the dft() function (To add a breakpoint, double click on the line number OR right-click anywhere on the line and click Breakpoint(CCS) >Breakpoint. You should see a small blue mark on the line where the breakpoint is inserted). 
Build the project and Debug. 
The program will first halt at entry to main(). Click Resume and program will halt at the  breakpoint. The array samples now contain the time-domain signal. We can visualize the signal using the Graph tool in CCS. Click on Tools>Graph>Single-Time. Set the Graph Properties as follows: 
Acquisition buffer size: 64 (since we have stored 64 samples of x[n]
DSP Data type: 32 bit floating point 
Index increment: 2 (The nature of the structure array samples is such that it comprises 2N float values ordered so that the first value is the real part of x[0],  the second is the imaginary part of x[0], the third is the real part of x[1], and so on. Since x[n] is purely real, we take alternate values only)
Start Address: samples
Data plot style: Bar
Display Data Size: 64
Click Ok. The graph of x[n] should appear. 
Click Step Over and the program and will now halt after the function dft() returns. At this point, the array samples contain the 64 DFT coefficients X[k]. The graph now displays the real part of the X[k]. If needed, click on Reset the Graph button and then Refresh button in the graph window. You should see two distinct peaks at k= 8 and k=56. 
*/

#include <math.h>
#define PI 3.1415926535897
#define M 64 //signal length
#define N 64 //DFT length
typedef struct
{
    float real;
    float imag;
} COMPLEX;
void dft(COMPLEX *);

void main()
{
    int n;
    float F = 1000.0, Fs = 8000.0;
    COMPLEX samples[N]={0.0};

//Generate time-domain signal
    for (n = 0; n < M; n++) //M samples of x[n]
    {
        samples[n].real = cos(2 * PI * F * n / Fs);
        samples[n].imag = 0.0;
    }

    
    dft(samples); //call DFT function
    
}

void dft(COMPLEX *x)
{
    COMPLEX result[N];
    int k, n;
    for (k = 0; k < N; k++) // N point DFT
    {
        result[k].real = 0.0;
        result[k].imag = 0.0;
        for (n = 0; n < N; n++)
        {
            result[k].real += x[n].real * cos(2 * PI * k * n / N) +
                              x[n].imag * sin(2 * PI * k * n / N);
            result[k].imag += x[n].imag * cos(2 * PI * k * n / N) -
                              x[n].real * sin(2 * PI * k * n / N);
        }
    }
    for (k = 0; k < N; k++)
    {
        x[k] = result[k];
    }
}








/*
 In the code below, we use the window method to design an FIR low pass filter for a cut-off frequency of 1 KHz and test it on real-time signals on the DSK. To convert the cut-off frequency in Hz to cut-off frequency in rad/sample, we use the relation omega = 2*pi*F/Fs.  Assuming a sampling frequency Fs=8kHz, the cut-off frequency in rad/sample is omega_c = pi/4. 

Test the filter with different frequencies in the pass-band and stop-band (connect headphone out of PC to line-in of DSK using aux cable and play tones from youtube).  Observe the effect of the filter on a music signal
*/
#include <dsk6713.h>
#include <dsk6713_aic23.h>
#include <math.h>
#define L 11 //length of filter
int main(void)
{
    Uint32 sample_pair;
    float pi = 3.141592653589;
    float hamming[L], h[L], x[L] = { 0 }, y;
    int i;

    //Generate Hamming window sequence
    for (i = 0; i < L; i++)
        hamming[i] = 0.54 - 0.46 * cos(2 * pi * i / (L - 1));

    //cut-off frequency of filter in rad/sample
    float wc = pi / 4.0;

    //compute filter coeffs
    for (i = 0; i < L; i++)
        //avoid division by 0 when i=(L-1)/2
        if (i == (L - 1) / 2)
            h[i] = wc / pi * hamming[i];
        else
            h[i] = sin(wc * (i - (L - 1) / 2.0)) / (pi*(i - (L - 1)/2.0))
                    * hamming[i];

    DSK6713_init();
    DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;
    DSK6713_AIC23_CodecHandle hCodec;
    hCodec = DSK6713_AIC23_openCodec(0, &config);

    /* Change the sampling rate to 8 kHz */
    DSK6713_AIC23_setFreq(hCodec, DSK6713_AIC23_FREQ_8KHZ);

    while (1)
    {
        while (!DSK6713_AIC23_read(hCodec, &sample_pair))
            ;
        //store top-half of sample from codec in x[0]
        x[0] = (int)sample_pair >>16;

        //process input sample:
        y = 0.0;

        for (i = 0; i < L; i++) //compute filter output
            y += h[i] * x[i];

        //shift delay line contents
        for (i = (L - 1); i > 0; i--)
            x[i] = x[i - 1];

        //output y to left channel
        sample_pair = (int)y <<16;
        while (!DSK6713_AIC23_write(hCodec, sample_pair))
            ;

    }

    /* Close the codec */
    DSK6713_AIC23_closeCodec(hCodec);

    return 0;
}







/*
C Library File IO: CCS also supports standard C library file I/O functions such as fopen( ), fclose( ), fread( ), fwrite( ), and so on. These functions not only provide the ability of operating on different file formats, but also allow users to directly use the functions on computers. Comparing with the memory load/save method introduced in the previous part, these file I/O functions are portable to other development environments. An example of C program that uses fopen( ), fclose( ), fread( ), and fwrite( ) functions is included below. Verify the working of the program in CCS.
*/

#include <stdio.h>
void writedatafile(void);
void conv(float *x, float *h, float *y, short l, short m);
void main()
{
    writedatafile();//put some data in a file
    FILE *fp;
    fp = fopen("../../myfilebin", "rb");//open file for binary read
    float x[5];
//read 5 floats from file and store in x
fread(x, sizeof(float), 5, fp);
    float h[3]={1, 1, 1};
    float y[7];
    conv(x, h, y, 5, 3);

}

void writedatafile()
{
    float a[] = { 1, 2, 3, 4, 5}; //some data
    FILE *fp;
    fp = fopen("../../myfilebin", "wb");//open for binary write

    //write 5 floats in array a to myfilebin
    fwrite(a, sizeof(float), 5, fp); 
    fclose(fp);

}

void conv(float *x, float *h, float *y, short l, short m)
{
    short k, kmin, kmax, n;
    for (n = 0; n < l + m - 1; n++)
    {
        y[n] = 0;
        kmin = (n > l - 1) ? n - l + 1 : 0;
        kmax = (n > m - 1) ? m - 1 : n;
        // printf("%d %d\n", kmin, kmax);
        for (k = kmin; k <= kmax; k++)
        {
            y[n] = y[n] + h[k] * x[n - k];
        }

    }

}

