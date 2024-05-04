#include <complex>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "primes.h"
#include <ctime>
#include <cstdlib>
#include <fftw3.h>

#include "config.h"

using namespace std;

//use this index to track the correct location of sig_sample
int sig_sample_index = 0;

double SNR = log(noiselevel)/log(10) * 10;

double std_noise = noiseless ? 0 : sqrt(k/(2*noiselevel));

const double pi = acos(-1);
const complex<double> I(0, 1);


//calculate a^b
int power(int a, int b) {
	int temp1 = 1;
	for (int i = 0; i < b; ++i) {
		temp1 = temp1*a;
	}
	return temp1;
}


//calculate pow(a, 2)
double pow_double(double a){
	return a*a;
}



//In this function, we try to find the index of the smallest main prime number
//this number is determined by N, k and first_main_prime
int index_of_smallest_main_prime(int N, int k) {
	int prime_num = log(N) / log(2) + 1;  //we will need at most log(N) / log(2) + 1 prime numbers.
	vector<int> p(prime_num);
	//this vector record which prime number shows up previously
	//0 indicate does not show up; otherwise indicate shows up
	vector<int> show_up(prime_num);
	//alpha record the power of each prime numbers
	vector<int> alpha(prime_num);
	//We will do not want q shows up in our prime sequence again, since we already have it in our prime sequence
	for (int i = 0; i < prime_num; ++i) {
		p[i] = PRIMES[i];
		show_up[i] = 0;
		alpha[i] = 0;
	}
	
	int temp = 1;
	int currentmin;
	int possiblemin;
	int l;
	int l_min;
	
	int S = 1;
	//2 must shows up in the first round
	show_up[0] = 1;
	
	while (temp < N / p[S]) {
		temp = 1;
		currentmin = N;
		l = 0;
		l_min = 0;
		while (alpha[l] > 0) {
			possiblemin = power(p[l], alpha[l])*(p[l] - 1);
			if (possiblemin < currentmin) {
				currentmin = possiblemin;
				l_min = l;
			}
			++l;
		}
		if (p[l] < currentmin) {
			l_min = l;
			alpha[l_min] = 1;
			++show_up[l_min];
		}
		else {
			++alpha[l_min];
			++show_up[l_min];
		}
		for (int i = 0; i < prime_num; ++i) {
			temp = temp*power(p[i], alpha[i]);
		}
		//update S the be the smallest index in show_up that is equal to 0
		for (int i = 0; i < show_up.size(); ++i){
			if (show_up[i] == 0){
				S = i;
				break;
			}
		}
	}
	
	//we want the first main prime number to be at least first_main_prime*k
	if (PRIMES[S] >= first_main_prime*k){
		return S;
	}
	else{
		while (PRIMES[S] < first_main_prime*k){
			++S;
		}
		return S;
	}
	
}



//In this part we try to build a sequence of prime numbers that we will use in our sampling scheme.
//Here we follow the principle that we want product of the sequence p_l^(alpha_l) goes to N/S as slow as possible.
//More specifically, the scheme we have will be q_j^(alpha_j), q_j^(alpha_j)*p_j1^(beta_j1), q_j^(alpha_j)*p_j1^(beta_j1)*p_j2^(beta_j2)....
//This prime sequence should satisfy q_j^(alpha_j)*p_j1^(beta_j1)*p_j2^(beta_j2)...>=N
//Also we want sum(q_j^(alpha_j), p_j1^(beta_j1), p_j2^(beta_j2)...) to ba as small as possible
//And for each p_j1^(beta_j1), p_j2^(beta_j2)... they should be close to each other as much as possible.
//The last property will be very useful later when we make our algorithm parallel.
void PrimeSequence(vector<int>& prime, vector<int>& power_of_prime, int q) {
	
	prime.push_back(q);
	power_of_prime.push_back(1);
	
	int S = power(q, power_of_prime[0]);
	int prime_num = log(N / S) / log(2) + 1;  //we will need at most log(N / S) / log(2) + 1 prime numbers.
	vector<int> p(prime_num);


	//We will do not want q shows up in our prime sequence again, since we already have it in our prime sequence
	int q_occurs = 0;
	for (int i = 0; i < prime_num; ++i) {
		if (PRIMES[i] != q) {
			p[i] = PRIMES[i + q_occurs];
		}
		else {
			++q_occurs;
			p[i] = PRIMES[i + q_occurs];
		}
	}
	vector<int> alpha(prime_num);
	for (int i = 0; i < prime_num; ++i) {
		alpha[i] = 0;
	}
	int temp = 1;
	int currentmin;
	int possiblemin;
	int l;
	int l_min;
	while (temp < N / S) {
		temp = 1;
		currentmin = N;
		l = 0;
		l_min = 0;
		while (alpha[l] > 0) {
			possiblemin = power(p[l], alpha[l])*(p[l] - 1);
			if (possiblemin < currentmin) {
				currentmin = possiblemin;
				l_min = l;
			}
			++l;
		}
		if (p[l] < currentmin) {
			l_min = l;
			alpha[l_min] = 1;
		}
		else {
			++alpha[l_min];
		}
		for (int i = 0; i < prime_num; ++i) {
			temp = temp*power(p[i], alpha[i]);
		}
	}
	//for (int i = 0; i < prime_num; ++i) {
	//std::cout << alpha[i] << std::endl;
	//}
	for (int i = 0; i < prime_num; ++i) {
		if (alpha[i] != 0) {
			prime.push_back(p[i]);
			power_of_prime.push_back(alpha[i]);
		}
	}
	//return alpha;
}


//In this function we calculate the value for phi function in Fermat-Euler Theorem by using Euler's product formula first, 
//Then return the value of phi function - 1, we will then use this value to calculate the modular inverse.
//Here the input n = q^(alpha)*p_1^(beta_1)*p_2^(beta_2)...*p_j^(beta_j)
//prime is a vector contains q, p_1, p_2, ..., p_j
int Phi_n_1(int n, vector<int> prime) {
	int length = prime.size();
	int answer = n;
	for (int i = 0; i < length; ++i) {
		answer = answer / prime[i] * (prime[i] - 1);
	}
	return answer - 1;
}





//In this part we use Fermat�Euler theorem to find the inverse (modulo ) for 'a mod n'. 
//The inverse (modulo) for 'a mod n' is defined in this way: if 'b' is inverse (modulo) for 'a mod n', 
//then b*a = 1 mod n. 
//The value of phi_n_1 comes from Phi_n_1 function, here the phi function is the one in Fermat-Euler Theorem
//In this part, we use long long int to prevent overflow. 
long long int InverseModP(int a, int phi_n_1, int n)
{
	long long int y = a;
	long long int temp = 1;
	
	//Here we apply the property of Euler's totient function. 
	//When q is a prime nbumber, the Euler's totient function of q^alpha is equal to q^(alpha)-q^(alpha-1)
	long long int temp1 = phi_n_1;
	//Since it will overflow if we calculate a^(temp1), so we will express temp1 into base 2 and then use the
	//equation: ab mod n = (a mod n)(b mod n) mod n
	while (temp1 != 0)
	{
		if (temp1 % 2 == 1)   //check whether it's current digit in base 2 is 1
		{
			temp = (temp*y) % n;
		}
		y = (y*y) % n;
		temp1 = temp1 / 2;
	}
	return temp%n;
}


//In this part we create one row for sample scheme matrix, it is in the form:
//q_i^(alpha_i)  p_i1^(beta_i1)  p_i2^(beta_i2), ..., p_ik^(beta_ik)
//where q_i, p_ij are prime numbers
//beta_ij are the power of prime numbers
void row_sampling_scheme(vector<vector <int> >& sample_scheme, vector<int> prime, vector<int>alpha) {
	vector<int> temp;   //use temp to store value for each row and then push_back temp in to sample_scheme
	int length_prime = prime.size();
	int temp1 = power(prime[0], alpha[0]);
	temp.push_back(temp1);
	for (int i = 1; i < length_prime; ++i) {
		temp.push_back(power(prime[i], alpha[i]));
	}
	//Add a row in sample_scheme.
	sample_scheme.push_back(temp);
	// cout << temp.size() << endl;
}


//create struct, struct contains four elements: 
//one is p, it is equal to the value of p_ij^(beta_ij)
//one is qp, it is equal to the vaue of  q_i^(alpha_i)*p_i1^(beta_i1)*...*p_i(j-1)^(beta_i(j-1))
//one is pInverseqp, it is equal to p_ij^(beta_ij) inverse mod q_i^(alpha_i)*p_i1^(beta_i1)*...*p_i(j-1)^(beta_i(j-1))
struct recovery_info
{
	int p;
	int qp;
	long long int pInverseqp;
};



//Generate sample_scheme, which is in the form:
//q_1^(alpha_1)  p_11^(beta_11)  p_12^(beta_12), ..., p_1k^(beta_1k)
//.
//.
//.
//q_m^(alpha_m)  p_m1^(beta_m1)  p_m2^(beta_m2), ..., p_mk^(beta_mk')
//The matrix contains m rows, but for each row, the number of elements may be different.
//Here m is a number relative to sparsity, N is the bandwidth
void sampling_scheme_recovery_info(vector<vector <int> >& sample_scheme, vector<vector <recovery_info> >& Recovery_Info, int m) {
	int phi;
	fftw_complex *in, *out;
	
	int smallest_main_index = index_of_smallest_main_prime(N, k);
	
	for (int i = 0; i < m; ++i) {
		//we will use temp1 to store subvector of prime, for each loop, we initialize temp1 first in order to set it as an empty vector.
		vector<int> temp;
		//generate prime sequence and their power first
		vector<int> prime;
		vector<int> power_of_prime;
		//Here we use PrimeSequence function generate prime and power_of_prime, they are vectors with same length
		PrimeSequence(prime, power_of_prime, PRIMES[smallest_main_index + i]);
		//use row_sampling_scheme generate a row of sample scheme matrix and push it back into sample_scheme
		row_sampling_scheme(sample_scheme, prime, power_of_prime);

		Recovery_Info.push_back(vector<recovery_info>());

		//Build Recovery_Info
		//temp1 is the pointer to the i th row of Recovery_Info
		vector<recovery_info> &temp1 = Recovery_Info[i];
		temp1.push_back(recovery_info());
		//p_i1^(beta_i1)
		temp1[0].p = sample_scheme[i][1];
		//q_i^(alpha_i)
		temp1[0].qp = sample_scheme[i][0];
		temp.push_back(prime[0]);
		//here phi is the short version for phi(n) - 1
		phi = Phi_n_1(sample_scheme[i][0], temp);
		//p_i1^(beta_i1) inverse mod q_i^(alpha_i)
		temp1[0].pInverseqp = InverseModP(sample_scheme[i][1], phi, sample_scheme[i][0]);
		
		/*
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * temp1[0].qp);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * temp1[0].qp);
		temp1[0].plan = fftw_plan_dft_1d(temp1[0].qp, in, out, FFTW_FORWARD, FFTW_MEASURE);
		*/
		
		//So far, we have first struct in row i
		//In this for loop, we will push_back other structs in row i of Recoery_Info
		for (int j = 1; j < prime.size() - 1; ++j) {
			temp1.push_back(recovery_info());
			temp1[j].p = sample_scheme[i][j+1];
			temp1[j].qp = temp1[j - 1].qp*temp1[j - 1].p;
			temp.push_back(prime[j]);
			phi = Phi_n_1(temp1[j].qp, temp);
			temp1[j].pInverseqp = InverseModP(temp1[j].p, phi, temp1[j].qp);
			
			/*
			in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * temp1[j].qp);
			out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * temp1[j].qp);
			temp1[j].plan = fftw_plan_dft_1d(temp1[j].qp, in, out, FFTW_FORWARD,FFTW_MEASURE);
			*/
			
			//temp1[j].plan = p;
			//paque pointer
		}
	}
}




//generate fftw_plans as a dictionary: map <int, fftw_plan> fftw_plans
//the element of fftw_plans is a pair, for which contains two part:
//int: corresponding to the length of DFT
//fftw_plan: store the pointer of fftw_plan with corresponding length
//It is worth to mention that we can use length of sequence as the key is because
//we use pairwise prime number to construct our sampling scheme
void create_plans(map <int, fftw_plan>& fftw_plans, vector<vector <int> >& sample_scheme){
	fftw_complex *in, *out;
	int dft_len;
	for (int i = 0; i < sample_scheme.size(); ++i){
		in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sample_scheme[i][0]);
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sample_scheme[i][0]);
		fftw_plans.insert(pair<int, fftw_plan>(sample_scheme[i][0], fftw_plan_dft_1d(sample_scheme[i][0], in, out, FFTW_FORWARD, FFTW_MEASURE)));
		for (int j = 1; j < sample_scheme[i].size(); ++j){
			dft_len = sample_scheme[i][0]*sample_scheme[i][j];
			in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dft_len);
		    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dft_len);
		    fftw_plans.insert(pair<int, fftw_plan>(dft_len, fftw_plan_dft_1d(dft_len, in, out, FFTW_FORWARD, FFTW_MEASURE)));
		}
	}
}



/*
//exp function is very slow, since in create_sample part, the x in exp(x) always in the interval [-2*pi, 0]
//so we generate the value of exp function on equal spaced points (the number of points determined by precision we want) on interval [-2*pi, 0]
//then everytime we need to calculate exp(x), we just use linear interpolation to approximate the exp(x)
//the precison can be guaranteed because: 1. the exp function in [-2*pi, 0] is smooth; 
//2. we can always get the precision we need by increasing number of pre-calculating equal spaced points.
void exp_equal_pts(vector<double>& exp_pts, int N){
	int pts_num = exp_pts.size();
	//create samples on [-2*pi, 0], include end points on both side
	for (int i = 0; i < pts_num; ++i){
		exp_pts[i] = exp(-2.1 * pi * log(N) / log(2) + i * 2.1 * pi * log(N) / log(2) / (pts_num - 1));
	}
}
*/



//create value of G(x_j - x) that needed in convolution part (create_sample function)
//we precalculate all the values that need in the convolution part and store them in map<pair<int, int>, vector<complex <double>>> conv_G
//In map<pair<int, int>, vector<complex <double>>> conv_G, the pair<int, int> is the key, which contains two part:
//the first int is the sample size, the second int means the shift of filter function
void create_conv_G(map<pair<int, int>, fftw_complex*>& conv_G, double sample_size, double c, int kappa, int filter_shift) {

	fftw_complex *temp;
	temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (int(sample_size)*(2*kappa + 1)));
	
	//change all the exp function to interpolation form
	int current_index = 0;
	
	double c_2 = pow_double(c)*2.0;
	
	double log_2_N = log(N)/log(2);
	
	int j;
	for (double i = 0.0; i < sample_size; ++i) {
		j = ceil(i * N / sample_size);
		if (kappa <= j && j < N - kappa) {
			//cout << "j" << j << endl;
			for (int l = j - kappa; l <= j + kappa; ++l) {
			    temp[current_index][0] = real(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				temp[current_index][1] = imag(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
			    ++current_index;
			}
		}
		else if (j < kappa) {
			for (int l = N - kappa + j; l <= N - 1; ++l) {
				temp[current_index][0] = real(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size - 2 * pi) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				temp[current_index][1] = imag(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size - 2 * pi) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				++current_index;
			}
			for (int l = 0; l <= j + kappa; ++l) {
				temp[current_index][0] = real(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				temp[current_index][1] = imag(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				++current_index;
			}
		}
		else if (j >= N - kappa) {
			for (int l = j - kappa - 1; l <= N - 1; ++l) {
				temp[current_index][0] = real(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				temp[current_index][1] = imag(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				++current_index;
			}
			for (int l = 0; l <= j + kappa - N; ++l) {
				temp[current_index][0] = real(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size + 2 * pi) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				temp[current_index][1] = imag(exp(-pow_double(2 * pi*l / N - 2 * pi*i / sample_size + 2 * pi) / c_2) * exp(I * double(filter_shift) * (2 * pi*i / sample_size - 2 * pi*l / N))/c);
				++current_index;
			}
		}
	}
	conv_G.insert(make_pair(make_pair(sample_size, filter_shift), temp));
}



//generate samples
//input: sig_sample(equispaced sample from signal), samples_voctor(store a single needed sample sequence)
//samples: a vector with elements as vector, this will be uesed to store all the samples needed in SFT algorithm
//conv_G: pre-computed values for periodic Gaussian function G
//sample_size: size of current sequence
//filter_shift: determine how to shift the filter function
//result: add one vector to samples

//create different equispace samples from pre-setting simple testing function
void create_sample2(fftw_complex*& sig_sample, fftw_complex*& samples_vector, map<pair<int, int>, fftw_complex*>& conv_G, double sample_size, int filter_shift) {
	complex<double> temp1;
	
	double log_2_N = log(N)/log(2);
	
	int j;
	
	int temp_conv_G = 0;
	
	fftw_complex *crt_G = conv_G[make_pair(sample_size, filter_shift)];
	
	//cout << "FLAG!!!!" << sample_size << endl;
	for (int i = 0; i < sample_size; ++i) {
		
		//cout << "The value of i is: " << i << endl;
		//cout << "sig_sample_index is : " << sig_sample_index << endl;
		//cout << "temp_conv_G is: " << temp_conv_G << endl;
		
		samples_vector[i][0] = 0;
		samples_vector[i][1] = 0;
		
		for (int l = 0; l < 2*kappa + 1; ++l) {
			samples_vector[i][0] += sig_sample[sig_sample_index][0] * (*(crt_G + temp_conv_G))[0] - sig_sample[sig_sample_index][1] * (*(crt_G + temp_conv_G))[1];
			samples_vector[i][1] += sig_sample[sig_sample_index][1] * (*(crt_G + temp_conv_G))[0] + sig_sample[sig_sample_index][0] * (*(crt_G + temp_conv_G))[1];
			
			++temp_conv_G;
			++sig_sample_index;
				
		}
	}
}



//create signal sample at x. 
//It is worth to notice that x is always on the N uniform sample scheme, i.e. x belongs to {2*pi*j/N for j = 0,...,N-1}.
complex<double> single_sig_sample(vector<double>& frequency, vector<complex<double>>& freq_coeff, double x) {

	complex<double> temp1;
	complex<double> noise;
	double u1, v1;
	u1 = drand48();
    v1 = drand48();
	//create random noise
	noise = std_noise * sqrt(-2*log(u1)) * exp(2*M_PI * I * v1);
	
	temp1 = 0;
	for (int j = 0; j < frequency.size(); ++j) {
		temp1 = temp1 + freq_coeff[j] * exp(I * frequency[j] * x);
	}
	//add noise to signal sample
	temp1 += noise;
	return temp1;
}




//generate needed uniform samples only
//the samples generated from this funciton are a subset of {*pi*j/N for j = 0,...,N-1}
//input: 
//sig_sample: contains all the needed sig_samples
//sample_size: size of needed sequence (a number much less than N)
//frequency: vector contains significant frequencies
//freq_coeff: vector contains the coefficients of significant frequencies

//result: add one vector to sig_samples

//create different equispace samples from pre-setting simple testing function
//Later need to set sample_size as duoble in order to use complex computing
void create_sample(fftw_complex*& sig_sample, double sample_size, vector<double>& frequency, vector<complex<double>>& freq_coeff) {
	
	double log_2_N = log(N)/log(2);
	complex<double> temp1;
	int j;
	for (double i = 0.0; i < sample_size; ++i) {
		j = ceil(i * N / sample_size);
		if (kappa <= j && j < N - kappa) {
			for (int l = j - kappa; l <= j + kappa; ++l) {
				//temp1 += sig_sample[l] * (exp_pts[ceil(desire_pt_index)] + exp_pts[floor(desire_pt_index)]) / 2.0 * (cos(filter_shift * (2 * pi*i / sample_size - 2 * pi*l / NN)) + I * sin(filter_shift * (2 * pi*i / sample_size - 2 * pi*l / NN)));
				temp1 = single_sig_sample(frequency, freq_coeff, 2 * pi*l / N);
				sig_sample[sig_sample_index][0] = real(temp1);
				sig_sample[sig_sample_index][1] = imag(temp1);
				++sig_sample_index;
			}
		}
		else if (j < kappa) {
			for (int l = N - kappa + j; l <= N - 1; ++l) {
				temp1 = single_sig_sample(frequency, freq_coeff, 2 * pi*l / N);
				sig_sample[sig_sample_index][0] = real(temp1);
				sig_sample[sig_sample_index][1] = imag(temp1);
				++sig_sample_index;
			}
			for (int l = 0; l <= j + kappa; ++l) {
				temp1 = single_sig_sample(frequency, freq_coeff, 2 * pi*l / N);
				sig_sample[sig_sample_index][0] = real(temp1);
				sig_sample[sig_sample_index][1] = imag(temp1);
				++sig_sample_index;
			}
		}
		else if (j >= N - kappa) {
			for (int l = j - kappa - 1; l <= N - 1; ++l) {
				temp1 = single_sig_sample(frequency, freq_coeff, 2 * pi*l / N);
				sig_sample[sig_sample_index][0] = real(temp1);
				sig_sample[sig_sample_index][1] = imag(temp1);
				++sig_sample_index;
			}
			for (int l = 0; l <= j + kappa - N; ++l) {
				temp1 = single_sig_sample(frequency, freq_coeff, 2 * pi*l / N);
				sig_sample[sig_sample_index][0] = real(temp1);
				sig_sample[sig_sample_index][1] = imag(temp1);
				++sig_sample_index;
			}
		}
	}
}




//calculate Discrete Fourier Transform by using FFTW
//input: vector<fftw_complex*>& samples (vector that cotains all the samples)
//input: map <int, fftw_plan>& fftw_plans, this is the 
//inout: vector<int>& sample_scheme_flat: a vector stores the length of all needed sequence in MSFT
//output: map <int, fftw_complex*>& DFT_samples (vector that contains all the samples after Fourier TRansform)
void DFT(vector<fftw_complex*>& samples, map <int, fftw_plan>& fftw_plans, map <int, fftw_complex*>& DFT_samples, vector<int>& sample_scheme_flat) {
	fftw_complex *in, *out;
	//for each DFT define
	//1. define in and out
	//2. record the length of each vector that will take DFT
	//3. execute corresponding fftw_plan to get out
	//4. store out in DFT_samples
	int dft_length;
	for (int i = 0; i < samples.size(); ++i) {
		dft_length = sample_scheme_flat[i];
		//1. define in and out
		in = samples[i];
		out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * dft_length);
		
		//2. fill in with complex numbers
		for (int l = 0; l < dft_length; ++l){
			in[l][0] = in[l][0]/dft_length;
			in[l][1] = in[l][1]/dft_length;
		}
		fftw_execute_dft(fftw_plans[dft_length], in, out);
		DFT_samples.insert(make_pair(dft_length, out));
	}
}


//define own function to compare two pairs <double number, index(integer)>
bool mag_complex(pair <complex<double>, int> a, pair <complex<double>, int> b)
{
	return abs(a.first) > abs(b.first);
}

//convert fftw to complex<double>
complex<double> fftw_to_std(fftw_complex c){
	return {c[0], c[1]};
}


//find index of t_min in MSFT
//input: A_sort[b], DFT_samples[p[l]*p[q+j]], p[l], p[q+j], r[0]
int t_m(complex <double> sorted_element, fftw_complex* DFT_vec, double p_l, double p_qj, int r_0) {
	double temp_min = abs(sorted_element - fftw_to_std(DFT_vec[r_0]));
	int t_m = 0;
	for (int t = 1; t < p_l; ++t) {
		if (abs(sorted_element - fftw_to_std(DFT_vec[t*(int)p_qj + r_0])) < temp_min) {
			temp_min = abs(sorted_element - fftw_to_std(DFT_vec[t*(int)p_qj + r_0]));
			t_m = t;
		}
	}
	return t_m;
}


// Returns modulo inverse of a with respect to mm using extended
int inv(int a, int mm)
{
	int m0 = mm, t, q;
	int x0 = 0, x1 = 1;
	if (mm == 1)
		return 0;
	// Apply extended Euclid Algorithm
	while (a > 1)
	{
		// q is quotient
		q = a / mm;
		t = mm;
		// m is remainder now, process same as
		// euclid's algo
		mm = a % mm, a = t;
		t = x0;
		x0 = x1 - q * x0;
		x1 = t;
	}
	// Make x1 positive
	if (x1 < 0)
		x1 += m0;
	return x1;
}


// k is size of num[] and rem[].  Returns the smallest
// number x such that:
//  x % num[0] = rem[0],
//  x % num[1] = rem[1],
//  ..................
//  x % num[k-1] = rem[k-1]
// Assumption: Numbers in num[] are pairwise coprime
// (gcd for every pair is 1)
int findMinX(vector<int> num, vector<int> rem, int k)
{
	// Compute product of all numbers
	int prod = 1;
	for (int i = 0; i < k; i++)
		prod *= num[i];

	// Initialize result
	int result = 0;

	// Apply above formula
	for (int i = 0; i < k; i++)
	{
		int pp = prod / num[i];
		result += rem[i] * inv(pp, num[i]) * pp;
	}

	return result % prod;
}


//find median of vector according to absolute value of elements
double findMedian(vector<double> vec) {
	double median;
	int size_vec = vec.size();

	sort(vec.begin(), vec.end());

	if (size_vec % 2 == 0) {
		median = (vec[size_vec / 2 - 1] + vec[size_vec / 2]) / 2.0;
	}
	else {
		median = vec[size_vec / 2];
	}
	return median;
}


//Recover frequency using information in Recovery_Info
//try torecover frequency x satisfy: rem[i] = freq % num[i]
//line_num tells which line of information of Recovery_Info should be used here
int findMinFreq(vector<int>& rem, vector<vector <recovery_info> >& Recovery_Info, int line_num){
	int x = rem[0];
	for (int i = 1; i < rem.size(); ++i){
		x = (rem[i] + (x - rem[i]) * Recovery_Info[line_num][i-1].pInverseqp * Recovery_Info[line_num][i-1].p) % (Recovery_Info[line_num][i-1].qp * Recovery_Info[line_num][i-1].p);
		if (x < 0){
			x = (Recovery_Info[line_num][i-1].qp * Recovery_Info[line_num][i-1].p) + x;
		}
	}
	return x;
}


//Energetic Frequencies Identification
//Find energetic frequencies in effective support of current window function
//input: 
//vector<vector <recovery_info> >& Recovery_Info: stores the inverse modular infomation that cna be used to recover significant frequencies according to Chinese Reminder Theorem
//map <int, int>& C: stores the information of all recovered frequencies. The first numnber of the frequency; the second number is how many times it has been recovered.
//vector<vector <int> >& sample_scheme: smapling scheme of MSFT. 
//map <int, fftw_complex*>& DFT_samples: DSFT transform of needed sequences
//int filter_shift: current center of periodic Gaussian function
//int FilterShift: floor(N / double(filter_num))
void EFI(vector<vector <recovery_info> >& Recovery_Info, map <int, int>& C, vector<vector <int> >& sample_scheme, map <int, fftw_complex*>& DFT_samples, int filter_shift, int FilterShift){
	//the element of A_sort is a pair. In the pair the first element is a complex number, the second element is the index based on sort result
	vector <pair <complex<double>, int>> A_sort;
	vector <int> r;
	int t_min;
	int temp_C;
	
	for (int j = 0; j < sample_scheme.size(); ++j) {
		//initialize a new vector A_sort to store the sorted complex vector
		A_sort.clear();
		for (int i = 0; i < sample_scheme[j][0]; ++i) {
			const auto &c = DFT_samples[sample_scheme[j][0]][i];
			A_sort.push_back(make_pair(complex<double>{c[0], c[1]}, i));
		}
		sort(A_sort.begin(), A_sort.end(), mag_complex);
		for (int b = 0; b < k; ++b) {
			r.clear();
			//no need to add index by 1, since index start from zero
			r.push_back(A_sort[b].second);
			for (int l = 1; l < sample_scheme[j].size(); ++l) {
				t_min = t_m(A_sort[b].first, DFT_samples[sample_scheme[j][0] * sample_scheme[j][l]], sample_scheme[j][l], sample_scheme[j][0], r[0]);
				r.push_back(int(r[0] + t_min*sample_scheme[j][0]) % int(sample_scheme[j][l]));
			}
			
			//Use Chinese Reminder Theorem to reconstruct the frequency
			//temp_C = findMinX(sample_scheme[j], r, sample_scheme[j].size()) % N;
			temp_C = findMinFreq(r, Recovery_Info, j) % N;
			//only recover frequencies in the interval [filter_shift - filter_shift/2, filter_shift + filter_shift/2]
			if (temp_C <= filter_shift - FilterShift/2.0 or temp_C >= filter_shift + FilterShift/2.0){
				continue;
			}
			
			//count how many times this frequency has been found
			if (C.find(temp_C) == C.end()) {
				//not found yet
				C.insert(pair<int, int>(temp_C, 1));
			}
			else {
				//found before
				C[temp_C] += 1;
			}
		}
	}
}


//Coefficient Estimation
//CE(R, C, sample_scheme, DFT_samples, c, filter_shift)
//R is the place store the final result, i.e. the frequnecy and coefficient pairs
//C contains the <frequency, show up times> pair
//sample_scheme contains the information of sampling scheme, here we use it to find the correct position of DFT_smaples
//DFT_samples contains the sequences in Fourier space
//c is the parameter in periodic Gaussian function
//filter_shift contains the information of corrent center of periodic Gaussian funciton, we use this to scale the coefficient from G*f to f
void CE(vector <pair<int, complex<double>>>& R, map <int, int>& C, vector<vector <int> >& sample_scheme, map <int, fftw_complex*>& DFT_samples, double c, int filter_shift){
	
	//create a vector to store the selected elements from A[p[m-1]*p[q+h]]
	vector<double> real_part, imag_part;
	//create a vector with pair as elements to store the significant frequencies
	
	complex<double> temp_coeff;
	
	for (const auto &item : C) {
		if (item.second >= 2) {
			real_part.clear();
			imag_part.clear();
			
			for (int h = 0; h < sample_scheme.size(); ++h) {
				real_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-1]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-1])][0]);
				imag_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-1]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-1])][1]);
				real_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-2]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-2])][0]);
				imag_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-2]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-2])][1]);
				real_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-3]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-3])][0]);
				imag_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-3]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-3])][1]);
				real_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-4]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-4])][0]);
				imag_part.push_back(DFT_samples[sample_scheme[h][0] * sample_scheme[h].end()[-4]][item.first % int(sample_scheme[h][0] * sample_scheme[h].end()[-4])][1]);
			}
			
			//find the median of real_part and imag_part and then get the coefficient
			temp_coeff = (findMedian(real_part) + findMedian(imag_part)*I)/sqrt(2*pi)/exp(-pow_double(c*(item.first - filter_shift))/2)*2.0*pi/double(N);
						
			//set the coefficient to be the median of real_part and imag_part
			R.push_back(pair<int, complex<double>>(item.first, temp_coeff));
		}
	}
}


//define own function to compare two pairs in R
bool R_mag_complex(pair<int, complex<double>> a, pair<int, complex<double>> b)
{
	return abs(a.second) > abs(b.second);
}

//define own function to compare two pairs in R
bool R_freq(pair<int, complex<double>> a, pair<int, complex<double>> b)
{
	return abs(a.first) > abs(b.first);
}


int main()
{	
	//in rounds number of numerical experiments, err_num of them fails to find all the significant frequencies
	int err_num = 0;
	//later will be use to store the largest error (abs(true - extimate)/abs(true)) in current round of numerical experiments
	double err = 0;
	//average L1 error: defined as the summation of all errors per frequencies divide by the number of significant frequencies
	double avg_l1_err = 0;
	
	int temp = 1;
	//find the index of first main prime number according to bandwidth, sparsity and first_main_prime
	int m_start = index_of_smallest_main_prime(N, k);
	
	//m is the number of set of prime numbers that we will need in sampling_scheme	
	int m = 0;
	
	//assign m the correct value
	while (temp < N) {
		temp = temp*PRIMES[m_start + m];
		++m;
	}
	//adjust the smallest number of set of samples to increase stability
	//the adjustment is according to number_of_sample_sets
	m = (m < ceil(number_of_sample_sets*log(k))) ? ceil(number_of_sample_sets*log(k)) : m;
	
	//Sample scheme
	vector<vector <int> > sample_scheme;
	//Recovery information
	vector<vector <recovery_info> > Recovery_Info;
	//q_1^(alpha_1)  q_1^(alpha_1)*p_11^(beta_11)  q_1^(alpha_1)*p_12^(beta_12), ..., q_1^(alpha_1)*p_1k^(beta_1k)
	//is one set of sample scheme
	sampling_scheme_recovery_info(sample_scheme, Recovery_Info, m);
	
	//count how many samples we needed in the MSFT algorithm for one filter function
	int total_samples = 0;
	cout << "Below is the sampling scheme matrix:" << '\n';
	for (int i = 0;i < sample_scheme.size(); ++i){
		total_samples = total_samples + sample_scheme[i][0];
		for (int j = 0; j < sample_scheme[i].size(); ++j){
			cout << sample_scheme[i][j] << " ";
			if (j > 0){
				total_samples = total_samples + sample_scheme[i][0]*sample_scheme[i][j];
			}
		}
		cout << '\n';
	}
	cout << "Number of samples need in one filter function is: " << total_samples << endl;
	
	//flatten the sample scheme. Later will be used as the set of key to find the correct position in some data structures like DFT_smaples
	vector<int> sample_scheme_flat;
	for (int i = 0;i < sample_scheme.size(); ++i){
		sample_scheme_flat.push_back(sample_scheme[i][0]);
		for (int j = 1; j < sample_scheme[i].size(); ++j){
			sample_scheme_flat.push_back(sample_scheme[i][0]*sample_scheme[i][j]);
		}
	}
	
	//sample_scheme.size() is number of rows in sample_scheme
	//sample_scheme[i].size() is number of elements in i th row of sample_scheme
	
	//define fftw_plans, this will later be used to store all the fftw_plan in FFTW
	//the position for each element in fftw_plans is corresponding to the elements in sample_scheme
	//for example: fftw_plans[0][0] corresponding to FFT with size sample_scheme[0][0]
	//fftw_plans[0][1] corresponding to FFT with size sample_scheme[0][0]*sample_scheme[0][1]
	//for j not equal to zero, fftw_plans[i][j] corresponding to FFT with size sample_scheme[i][0]*sample_scheme[i][j]
	map <int, fftw_plan> fftw_plans;
    cout << "Start create fftw_plans..." << '\n';
	create_plans(fftw_plans, sample_scheme);
	cout << "Done!" << '\n';
	
	
	//samples contains all the samples that will be needed in MSFT
	//at this point we only initialize samples based on sample_scheme
	//the value of each element in samples will be filled in later
	//each element of samples is a fftw_complex pointer, fftw_complex is a high efficiency data structure of fftw3
	vector<fftw_complex*> samples;
	//initialize samples once fix sampling_scheme
	fftw_complex *samples_temp;
	for (int j = 0; j < sample_scheme.size(); ++j) {
		samples_temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * sample_scheme[j][0]);
		samples.push_back(samples_temp);
		for (int l = 1; l < sample_scheme[j].size(); ++l) {
			samples_temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (sample_scheme[j][0]*sample_scheme[j][l]));
		    samples.push_back(samples_temp);
		}
	}
	//cout << "Size of samples is: " << samples.size() << endl;
	int samples_index = 0;
	//mmm is the k in our paper. mmm determines how many filter function we use to cover the whole bandwidth
	int mmm = 1;
	
	double BB = sqrt(24 * log(2));
	
	//c is the key variable of filter function
	double c = BB * sqrt(log(N) / log(2)) / (c_param*N);
	
	//double BB_prime = BB / sqrt(pi);
	double BB_prime = BB * sqrt(24*log(2)) / pi;
	
	//number of filter functions use to cover the whole bandwidth in Fourier space
	//for mmm = 1
	//if bandwidth is smaller than 2^16, the number of filter function we need is 2
	//if bandwidth is between 2^16 to 2^36, the number of filter function we need is 3
	//if bandwidth is between 2^36 to 2^64, the number of filter function we need is 4
	int filter_num = ceil(mmm * sqrt(log(N) / log(2))/2);
	cout << "filter number is: " << filter_num << endl;
	
	//length of a single shift
	int FilterShift = floor(N / double(filter_num));
	cout << "FilterShift is: " << FilterShift << endl;
	
	
	cout << "Start pre-compute convolution function G..." << '\n';
	
	//Since MSFT is a deterministic algorithm, this means that once we fixed bandwidth N and sparsity k
	//the sample_scheme can be determined immediately
	//once the sample_scheme has been determined, we can precompute the values of periodic Gaussian functions that will be used in the sample creating function
	//for fixed bandwidth and sparsity, the conv_G can be used repeatedly for many times
	//conv_G is a map. Its key is a pair of two integers: the first one is the size of create smaples, the second one is the current center of periodic Gaussian function
	map<pair<int, int>, fftw_complex*> conv_G;
	int filter_shift_current;
	for (int i = 0; i < filter_num; ++i){
		filter_shift_current = floor(N / double(filter_num) / 2.0) + i*FilterShift;
	    for (int j = 0; j < sample_scheme.size(); ++j) {
		    create_conv_G(conv_G, double(sample_scheme[j][0]), c, kappa, filter_shift_current);
		    for (int l = 1; l < sample_scheme[j].size(); ++l) {
		        //cout << p[l]*p[q+j] << '\n';
				create_conv_G(conv_G, double(sample_scheme[j][0]*sample_scheme[j][l]), c, kappa, filter_shift_current);
		    }
	    }
	}
	cout << "Done!" << '\n';
	
	//For large bandwidth N, create the whole uniform sample vector can be very time consuming
	//only need those uniform samples that will be used in convolution
	//This will not impact the running time and at the same time dramatically decrease the needed storage.
	//we can repeatly use sig_sample for different filter function, since the shift of filter function will not impact the location of uniform samples.
	fftw_complex *sig_sample;
	sig_sample = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (total_samples * (2* kappa + 1)));
	
	//record the error for each round of numerical test.
	vector<double> err_record;
	//record how long time it takes for num_test repeated numerical experiments
	double duration;
	//random seed
	int seed = 123;
	
	//start different numerical test
	for (int round = 0; round < rounds; ++round){
	
	srand(seed);
	++seed;
	vector<double> freq;
	vector<complex<double>> coeff (k,0);
	//create random frequencies and coefficients
	//the way we create coefficients guarantee that they are being nomrliazed already
	vector<double>::iterator it;
	while (freq.size() < k)
	{
		for (int i = 0; i < k - freq.size(); i++)
			freq.push_back(rand() % N);
		sort(freq.begin(), freq.end());
		it = unique(freq.begin(), freq.end());
		freq.resize(it - freq.begin());
	}
	
	double arg;
	for (int i = 0; i < k; ++i){
		arg = 2 * pi * drand48();
		//it is worth to mention that below we do not need to normalize the coefficient with abs(complex<double>(cos(arg), sin(arg)))
		coeff[i] = /*rad*/complex<double>(cos(arg), sin(arg))/abs(complex<double>(cos(arg), sin(arg)));
	}
	
	//create needed signal samples by using function create_sample
	//during this part, we fill values into sig_sample
	//two things need to mention here:
	//first: sig_sample contains only needed uniform samples
	//second: usually, with large bandwith, i.e. 2^26, the size of sig_sample is less than 1% of bandwidth, which dramatically decrease the signal generation time
	for (int j = 0; j < sample_scheme.size(); ++j) {
		create_sample(sig_sample, double(sample_scheme[j][0]), freq, coeff);
	    for (int l = 1; l < sample_scheme[j].size(); ++l) {
			create_sample(sig_sample, double(sample_scheme[j][0]*sample_scheme[j][l]), freq, coeff);
	    }
	}
	//cout << "Size of sig_sample is:" << sig_sample_index << endl;
	
	sig_sample_index = 0;
	
	
	
	//R is where store the frequencies and coefficients recovered from one filter (window) function
	//pair<frequency, coefficient>
	vector <pair<int, complex<double>>> R;
	//pick k frequencies with largest magniture from R
	//Store it in R_final
	vector <pair<int, complex<double>>> R_final;
	//convert freq and coeff into the same form of R_final
	vector <pair<int, complex<double>>> R_true;
	
	
	clock_t start;
	
	start = std::clock();
	//start repeated numerical test
	for (int boo = 0; boo < num_test; ++boo){
	
	R.clear();
	R_final.clear();
	R_true.clear();
	
	//shift filter function for filter_num times to cover the whole bandwidth
	for (int i = 0; i < filter_num; ++i){
		
		filter_shift_current = floor(N / double(filter_num) / 2.0) + i*FilterShift;
				
		//create needed samples
		//now, these samples are created by using function create_sample2
		//with sig_sample and pre-computing conv_G
	    for (int j = 0; j < sample_scheme.size(); ++j) {
			create_sample2(sig_sample, samples[samples_index], conv_G, double(sample_scheme[j][0]), filter_shift_current);
			++samples_index;
		    for (int l = 1; l < sample_scheme[j].size(); ++l) {
			    create_sample2(sig_sample, samples[samples_index], conv_G, double(sample_scheme[j][0]*sample_scheme[j][l]), filter_shift_current);
				++samples_index;
		    }
	    }
		
		sig_sample_index = 0;
		samples_index = 0;
		
		//DFT_samples stores the DFT of sequence in samples
	    map <int, fftw_complex*> DFT_samples;
	    //Use FFTW to get DFT of sequence		
		
	    DFT(samples, fftw_plans, DFT_samples, sample_scheme_flat);	    
	    
	    //Set C as a map, for each item in C
	    //item.first is the reconstructed frequency, item.second is the counted number of this frequency
	    map <int, int> C;
	    
	    //Energetic frequency identification
	    EFI(Recovery_Info, C, sample_scheme, DFT_samples, filter_shift_current, FilterShift);
	    
	    //coefficient estimation
	    CE(R, C, sample_scheme, DFT_samples, c, filter_shift_current);
		
	}
	//end filter_num
	
	
	}
	//end repeat numerical test
	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	
	//choose k frequencies with largest magnitude
	sort(R.begin(), R.end(), R_mag_complex);
	for (int i = 0; i < k; ++i) {
		R_final.push_back(R[i]);
	}
	
	
	//convert true freq and coeff into the same form of R_final
	//varify results
	for (int i = 0; i < k; ++i){
		R_true.push_back(pair<int, complex<double>>(freq[i], coeff[i]));
	}
	
	
	//sort R_final and R_true
	//R_final contains the results from DMSFT algorithm
	//R_true contains the true frequency and coefficient pairs
	sort(R_final.begin(), R_final.end(), R_freq);
	sort(R_true.begin(), R_true.end(), R_freq);
	
	//check whether we recovered all the significant frequencies
	bool TRUE_SET = true;
	//err and err_num
	for (int i = 0; i < k; ++i){
		//if one of the frequencies is different with what we have in R_true, err_num += 1
		if (R_final[i].first != R_true[i].first){
			err_num += 1;
			//It is not a true set
			TRUE_SET = false;
			break;
		}
	}
	//only compute error and average L1 error for correct recovery
	if (TRUE_SET){
		for (int i = 0; i < k; ++i){
			avg_l1_err = avg_l1_err + abs(R_final[i].second - R_true[i].second)/k;
			if (abs(R_final[i].second - R_true[i].second)/abs(R_true[i].second) > err){
				err = abs(R_final[i].second - R_true[i].second)/abs(R_true[i].second);
			}
		}
		err_record.push_back(err);
		err = 0;
	}
	
	freq.clear();
	coeff.clear();
	//at the end of each round, clear the sig_sample
	sig_sample_index = 0;
	
	}
	//end rounds
	
	//compute the average L1 error for each correct recovery
	double avg_err = 0;
	for (int i = 0; i < err_record.size(); ++i){
		//cout << err_record[i] << endl;
		avg_err += err_record[i];
	}
	avg_err = avg_err/err_record.size();
	
	
	//output results
	cout << "Bandwidth is: " << N << '\n';
	cout << "Sparsity is: " << k << '\n';
	cout << "Number of terms in each dicrete convolution is: " << 2*kappa + 1 << '\n';
	if (!noiseless){
	    cout << "The noise level (SNR) is: " << SNR << endl;
	}
	else{
		cout << "There is no noise in the original signal." << endl;
	}
	cout << "Running time (w/o create equispaced samples) is: " << duration/num_test << "seconds" << '\n';
	cout << "Amoung " << rounds << " of numerical experiments, " << err_num << " of them are wrong." << '\n';
	cout << "The average infinite norm error is: " << avg_err << endl;
	cout << "The average L1 norm error is: " << avg_l1_err / rounds << endl;
}