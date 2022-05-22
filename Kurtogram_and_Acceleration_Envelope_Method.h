
                            ////////////////////////////////////////////
							////////////* Macro Definition *////////////
							////////////////////////////////////////////
//////* kurtogram *//////
#define N_LEVEL 7
#define PI 3.14159265
#define MAX_ORDER 25

//////* acc envelope method *//////
//////////////////////////////*band pass filter*//////////////////////////////
//bandpass
#define BANDPASS_HALF_ORDER 50

//////////////////////////////*envelope*//////////////////////////////
//downsampling
#define DOWN_SAMPLE_FACTOR 2

//lowpass
#define LOWPASS_HALF_ORDER 8
#define LOWPASS_CUTOFF 3500.0

//////////////////////////////*FFT*//////////////////////////////
//FFT
#define TRSIZ 8
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

//////////////////////////////*Bearing Fault Detector*//////////////////////////////
#define TOLERANCE 0.25
#define TREE_HEIGHT 6
#define CTR_THRESHOLD 5

							  ///////////////////////////////////////////////
							  ////////////* Variable Definition *////////////
							  ///////////////////////////////////////////////
//////* variables for file path *//////
const char* input_file_path;

//////* variables for kurtogram *//////
int time_length;
int sampling_rate;
int row = 2;
int raw_data_length;
int data_length;
int base_number; //for kurtosis table construction
double sub_kurtosis_table[5] = { 0 };
double* source_real_ptr = 0;
double* source_imag_ptr = 0;
double* target_real_ptr = 0;
double* target_imag_ptr = 0;

int sub_index;
int filter_length;
int down_factor;

double* filter_real_select = 0;
double* filter_imag_select = 0;

int down_cnt;
int idx_cnt;

int source_start_point;
int target_start_point;


double	h_real[17] = { -0.0019,0.0022,0,0,0.0409,0.0316,0,0.2067,0.4003,0.2067,0,0.0316,0.0409,0,0,0.0022,-0.0019 },
h_imag[17] = { 0,0.0022,0.0108,0,0,0.0316,-0.081,-0.2067,0,0.2067,0.081,-0.0316,0,0,-0.0108,-0.0022,0 },
g_real[16] = { -0.0022,0,0,0.0409,-0.0316,0,-0.2067,0.4003,-0.2067,0,-0.0316,0.0409,0,0,-0.0022,-0.0019 },
g_imag[16] = { 0.0022,-0.0108,0,0,0.0316,0.081,-0.2067,0,0.2067,-0.081,-0.0316,0,0,0.0108,-0.0022,0 },
h1_real[25] = { -0.0012,0.0005,0.002,0,-0.0025,0.0067,0.0272,0.0315,0.0064,0,0.0743,0.2016,0.2666,0.2016,0.0743,0,0.0064,0.0315,0.0272,0.0067,-0.0025,0,0.002,0.0005,-0.0012 },
h1_imag[25] = { 0,0.0003,0.0034,0.0072,0.0043,-0.004,0,0.0182,0.011,-0.054,-0.1286,-0.1164,0,0.1164,0.1286,0.054,-0.011,-0.0182,0,0.0039,-0.0043,-0.0072,-0.0034,-0.0003,0 },
h2_real[25] = { -0.0012,0,-0.0039,0,0.005,0,0.0272,0,-0.0127,0,-0.1485,0,0.2666,0,-0.1485,0,-0.0127,0,0.0272,0,0.005,0,-0.0039,0,-0.0012 },
h2_imag[25] = { 0,0.0006,0,-0.0072,0,-0.0078,0,0.0363,0,0.054,0,-0.2328,0,0.2328,0,-0.054,0,-0.0363,0,0.0078,0,0.0072,0,-0.0006,0 },
h3_real[25] = { -0.0012,-0.0005,0.002,0,-0.0025,-0.0067,0.0272,-0.0315,0.0064,0,0.0743,-0.2016,0.2666,-0.2016,0.0743,0,0.0064,-0.0315,0.0272,-0.0067,-0.0025,0,0.002,-0.0005,-0.0012 },
h3_imag[25] = { 0,0.0003,-0.0034,0.0072,-0.0043,-0.0039,0,0.0182,-0.011,-0.054,0.1286,-0.1164,0,0.1164,-0.1286,0.054,0.011,-0.0182,0,0.0039,0.0043,-0.0072,0.0034,-0.0003,0 };

double h_data[2][17] = { 0 }, g_data[2][16] = { 0 }, h1_data[2][25] = { 0 }, h2_data[2][25] = { 0 }, h3_data[2][25] = { 0 };
double* h_real_ptr = &h_real[0];
double* h_imag_ptr = &h_imag[0];
double* g_real_ptr = &g_real[0];
double* g_imag_ptr = &g_imag[0];
double* h1_real_ptr = &h1_real[0];
double* h1_imag_ptr = &h1_imag[0];
double* h2_real_ptr = &h2_real[0];
double* h2_imag_ptr = &h2_imag[0];
double* h3_real_ptr = &h3_real[0];
double* h3_imag_ptr = &h3_imag[0];


struct str_filter {
	int order_num;
	int idx;
	int filter_cnt;
	double data[2][MAX_ORDER];
	double filtered[2];
};

struct str_filter comp_filter;

////////////* variables for searching maximum kurtosis, center frequency and bandwidth *////////////
double max_value;
double bandwidth;
double center_frequency;

//////* variables for acc_envelope_method *//////
int down_length;
//////////////////////////////*band pass filter*//////////////////////////////
double transfreq_low;
double transfreq_high;
double ft_low;
double ft_high;

double bandpass_coefficient[BANDPASS_HALF_ORDER * 2 + 1];
double bandpass_output;

struct str_bandpass_filter {
	int bandpass_order_num;
	int bandpass_idx;
	int bandpass_filter_cnt;
	double bandpass_data[BANDPASS_HALF_ORDER * 2 + 1];
	double bandpass_filtered;
};

struct str_bandpass_filter comp_bandpass_filter;

//////////////////////////////*envelope*//////////////////////////////
double squaring_output;
double lowpass_output;

struct str_lowpass_filter {
	int lowpass_order_num;
	int lowpass_idx;
	int lowpass_filter_cnt;
	double lowpass_data[LOWPASS_HALF_ORDER * 2 + 1];
	double lowpass_filtered;
};
struct str_lowpass_filter comp_lowpass_filter;
double normalized_freq;
double sqrt_output;
int envelope_idx_cnt;
int envelope_down_cnt;
double lowpass_coefficient[LOWPASS_HALF_ORDER * 2 + 1];

//////////////////////////////*FFT*//////////////////////////////
int fft_max_size;
struct str_FFT_data {
	double fs;
	int size;
};

struct str_FFT_data FFT;

//////* variables for Bearing Fault Detector *//////
int diagnosis_result_flag;
double n, pd, bd, theta, omega;
const char* log_text;