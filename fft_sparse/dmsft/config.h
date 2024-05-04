// Credit to MSU, markiwen08 & Ruochuan Zhang
// Modified by team
// DMSFT

#include <complex>

//kappa is the number of points we use in the convolution for each side, the actual number of points in convolution is 2*kappa + 1
//usually choose kappa between 4 to 6 can achive very good balance between accuracy and running time.
int kappa = 4;

//determine the first main prime in sampleing scheme. 
//Theoretically, increase this number can increase the accuracy, but will also increase the running time
double first_main_prime = 1.2;

//control the number of sets of prime numbers in sampling scheme.
//Theoretically, increase this number can increase the accuracy, but will also increase the running time
double number_of_sample_sets = 1.3;

//the parameter c in periodic Gaussian function
//the value of c determines how repidly (slowly) the periodic Gaussian function decrease in time (Fourier) space
//increase the value of c means need less number in convolutio part to keep the same accuracy
double c_param = 3;

//number of repeated numerical test, need to be at least 1
int num_test = 1;

//number of different numerical test, need to be at least 1
int rounds = 10;

//N is the bandwidth, later will take the value of N from outside (N is the smallest power of two number that greater than sample size)
//e.g. if sample size is 1000, then N is 2^10 = 1024
const int N = pow(2, 16);

//Number of significant frequencies, later will also take this number from outside
const int k = 50;

//If run in the noiseless case, set noiseless to be true, otherwoise, set it to be false.
bool noiseless = true;

//noise level: 1000000 means SNR = 60; 
//100000 means SNR = 50;
//10000 means SNR = 40;
//1000 means SNR = 30;
//100 means SNR = 20;
//10 means SNR = 10;
//1 means SNR = 0;
//increase noiselevel will decrease the SNR
//e.g. SNR 0 means the power of noise is as much as the power of signal
double noiselevel = 1;