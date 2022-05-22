
#include <iostream>
#include <fstream>
#include <math.h> 

#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "Kurtogram_and_Acceleration_Envelope_Method.h"

using namespace std;

                            
                            ///////////////////////////////////////////////
                            ////////////* Function Definition *////////////
                            ///////////////////////////////////////////////

//////////////////////////////////////////////
//////* function for parameter setting *//////
//////////////////////////////////////////////
void set_bearing_spec(double bn, double ball_diameter, double pitch_diameter, double contact_angle, double motor_speed) {
    n = bn; //滾珠數(N)
    bd = ball_diameter; //滾珠直徑(bd)(mm) 
    pd = pitch_diameter; //軸承節徑(pd)(mm)
    theta = contact_angle; //接觸角(theta)(degree)
    omega = motor_speed; //轉速(omega)(RPM)
}
void set_exp_condition(int time_len, int samplerate) {
    time_length = time_len; //test time length(sec.)
    sampling_rate = samplerate; //sampling rate(per sec.)
}
void set_input_filepath(const char* argv_in){
    input_file_path = argv_in;
}

//////////////////////////////////////////////////////////
//////* Kurtogram_and_Acceleration_Envelope_Method *//////
//////////////////////////////////////////////////////////


//////* function for kurtogram *//////
/////////////////////////* position transfer to acceleration */////////////////////////
void pos_to_acc(vector<vector<double>> &x)
{
    int ppr;
    int max_value_tmp = -20000000;

    for (int j = 0; j < raw_data_length; j++)
    {
        if (x[0][j] > max_value_tmp)
        {
            max_value_tmp = x[0][j];
        }
    }

    if (max_value_tmp > 16000000) {
        ppr = 16777216;
    }
    else if (max_value_tmp > 4200000) {
        ppr = 5120000;
    }
    else {
        ppr = 4194304;
    }

    for (int i = 0; i < raw_data_length - 1; i++) {
        x[0][i] = x[0][i + 1] - x[0][i];
        if (fabs(x[0][i]) > ppr * 0.5) {
            if (x[0][i] > 0) {
                x[0][i] = x[0][i] - ppr;
            }
            else {
                x[0][i] = x[0][i] + ppr;
            }
        }
        x[0][i] = x[0][i] * 60 * sampling_rate / ppr;
    }

    for (int i = 0; i < raw_data_length - 2; i++) {
        x[0][i] = (x[0][i + 1] - x[0][i]) * 60 * sampling_rate;
    }

}

/////////////////////////* moving acceleration signal to remove axis frequency */////////////////////////
void acc_moving(vector<vector<double>>& x, int movmean_length)
{
    double acc_movmean_sum;
    double acc_movmean_mean;

    for (int i = 0; i < ((raw_data_length - 2) - movmean_length + 1); i++) {
        acc_movmean_sum = 0;
        for (int j = i; j < (i + movmean_length - 1); j++)
        {
            acc_movmean_sum = acc_movmean_sum + x[0][j];
        }
        acc_movmean_mean = (acc_movmean_sum / movmean_length);
        x[0][i] = acc_movmean_mean;
    }

}

void filter_init(struct str_filter* p, int point)
{
    int n;

    p->idx = 0;
    p->filter_cnt = 1;
    p->filtered[0] = 0;
    p->filtered[1] = 0;
    p->order_num = point;

    for (n = 0; n < point; n++) {
        p->data[0][n] = 0;
        p->data[1][n] = 0;
    }
}

double* filter(struct str_filter* p, double x_real, double x_imag, double* filter_real_type, double* filter_imag_type)
{
    int idx = p->idx;

    p->data[0][idx] = x_real;
    p->data[1][idx] = x_imag;

    p->filtered[0] = 0;
    p->filtered[1] = 0;

    for (int i = 0; i < p->filter_cnt; i++)
    {
        int j;
        j = p->filter_cnt - i - 1;

        p->filtered[0] = p->filtered[0] +
            p->data[0][j] * *(filter_real_type + i) - p->data[1][j] * *(filter_imag_type + i);
        p->filtered[1] = p->filtered[1] +
            p->data[0][j] * *(filter_imag_type + i) + p->data[1][j] * *(filter_real_type + i);
    }


    filter_real_type = filter_real_type + p->filter_cnt;
    filter_imag_type = filter_imag_type + p->filter_cnt;

    for (int i = 0; i < (p->order_num - p->filter_cnt); i++)
    {
        int j;
        j = p->order_num - i - 1;

        p->filtered[0] = p->filtered[0] +
            p->data[0][j] * *(filter_real_type + i) - p->data[1][j] * *(filter_imag_type + i);
        p->filtered[1] = p->filtered[1] +
            p->data[0][j] * *(filter_imag_type + i) + p->data[1][j] * *(filter_real_type + i);
    }

    p->idx++; //update idx
    p->filter_cnt++;

    if (p->filter_cnt > p->order_num)
    {
        p->filter_cnt = 1;
    }

    if (p->idx >= p->order_num) {
        p->idx = 0;
    }

    return(p->filtered);
}

int power(int base, int exponent) {
    int result;
    result = 1;
    for (int i = 0; i < exponent; i++)
    {
        result *= base;
    }
    return result;
}

double kurt(double* real_ptr, double* imag_ptr, int layer, int down_factor, int filter_length) {
    double a = 0, b = 0;
    int count = 0, count2 = 0;
    int complex_amount = 0;
    int start_sum = 0, end_sum = 0;
    double K = 0, real_sum = 0, imag_sum = 0, real_mean = 0, imag_mean = 0, abs = 0, abs_square_sum = 0, abs_fourth_sum = 0, E = 0;
    if (down_factor == 1) {
        complex_amount = data_length;
        start_sum = 0;
        end_sum = complex_amount;
    }
    else if (down_factor == 2 && filter_length == 16)
    {
        complex_amount = (int)data_length / (power(2, layer) * down_factor) - 15;
        start_sum = 15;
        end_sum = start_sum + complex_amount;
    }
    else
    {
        complex_amount = (int)data_length / (power(2, layer) * down_factor) - 16;
        start_sum = 16;
        end_sum = start_sum + complex_amount;
    }

    for (int j = start_sum; j < end_sum; j++)
    {
        if (*(real_ptr + j) == 0.0 && *(imag_ptr + j) == 0.0)
        {
            count2++;
        }
        a = *(real_ptr + j);
        b = *(imag_ptr + j);
        real_sum = real_sum + a;
        imag_sum = imag_sum + b;
    }
    if (count2 == complex_amount)
    {
        K = 0;
    }
    real_mean = real_sum / complex_amount;
    imag_mean = imag_sum / complex_amount;
    for (int j = start_sum; j < end_sum; j++)
    {
        abs = sqrt((*(real_ptr + j) - real_mean) * (*(real_ptr + j) - real_mean) +
            (*(imag_ptr + j) - imag_mean) * (*(imag_ptr + j) - imag_mean));
        abs_square_sum += abs * abs;
        abs_fourth_sum += abs * abs * abs * abs;
    }
    E = abs_square_sum / complex_amount;
    if (E < 2.2204e-16)
    {
        K = 0;
    }
    K = (abs_fourth_sum / complex_amount) / (E * E);
    for (int j = start_sum; j < end_sum; j++)
    {
        if (*(imag_ptr + j) == 0.0)
            count++;
    }
    if (count == complex_amount)
        K -= 3;
    else
        K -= 2;

    return K;
}

double* build_kurtosis_table(int layer, int index, vector<vector<double>>& data_A, vector<vector<double>>& data_B)
{

    int data_length_now = (int)data_length / power(2, layer);

    if (index == 0) {
        source_start_point = 0;
    }
    else
    {
        source_start_point = (int)data_length / power(2, layer) * index;
    }
    for (sub_index = 0; sub_index < 5; sub_index++)
    {
        // determine start point of target
        if (sub_index == 0)
        {
            target_start_point = source_start_point;
            filter_length = 25;
            down_factor = 3;

            filter_real_select = h1_real_ptr;
            filter_imag_select = h1_imag_ptr;
        }
        else if (sub_index == 1)
        {
            target_start_point = (int)source_start_point + data_length_now / 3;
            filter_length = 25;
            down_factor = 3;

            filter_real_select = h2_real_ptr;
            filter_imag_select = h2_imag_ptr;
        }
        else if (sub_index == 2)
        {
            target_start_point = (int)source_start_point + data_length_now / 3 * 2;
            filter_length = 25;
            down_factor = 3;

            filter_real_select = h3_real_ptr;
            filter_imag_select = h3_imag_ptr;
        }
        else if (sub_index == 3)
        {
            target_start_point = source_start_point;
            filter_length = 17;
            down_factor = 2;

            filter_real_select = h_real_ptr;
            filter_imag_select = h_imag_ptr;
        }
        else if (sub_index == 4)
        {
            target_start_point = (int)source_start_point + data_length_now / 2;
            filter_length = 16;
            down_factor = 2;

            filter_real_select = g_real_ptr;
            filter_imag_select = g_imag_ptr;
        }

        // determine start point of source by layer and index
        if (layer % 2 == 0) {
           
            source_real_ptr = &data_A[0][source_start_point];
            source_imag_ptr = &data_A[1][source_start_point];
            target_real_ptr = &data_B[0][target_start_point];
            target_imag_ptr = &data_B[1][target_start_point];

        }
        else
        {
            source_real_ptr = &data_B[0][source_start_point];
            source_imag_ptr = &data_B[1][source_start_point];
            target_real_ptr = &data_A[0][target_start_point];
            target_imag_ptr = &data_A[1][target_start_point];

        }

        down_cnt = down_factor;
        idx_cnt = 0;

        filter_init(&comp_filter, filter_length);

        for (int i = 0; i < data_length_now; i++)
        {
            double* filtered_ptr = 0;
            double input_real;
            double input_imag;

            input_real = *(source_real_ptr + i);
            input_imag = *(source_imag_ptr + i);

            // FILTERING
            filtered_ptr = filter(&comp_filter, input_real, input_imag, filter_real_select, filter_imag_select);
            // DOWN-SAMPLING
            if (down_cnt == down_factor)
            {
                *(target_real_ptr + idx_cnt) = *(filtered_ptr + 0);
                *(target_imag_ptr + idx_cnt) = *(filtered_ptr + 1);
                down_cnt = 0;
                idx_cnt++;
            }
            down_cnt++;
        }
        // kurtosis calculation
        double sign = -1;
        if (filter_length == 16) {
            for (int i = 0; i < (int)data_length_now / 2; i++)
            {
                *(target_real_ptr + i) = *(target_real_ptr + i) * sign;
                *(target_imag_ptr + i) = *(target_imag_ptr + i) * sign;
                sign = sign * -1;
            }

        }
        sub_kurtosis_table[sub_index] = kurt(target_real_ptr, target_imag_ptr, layer, down_factor, filter_length);
    }

    return sub_kurtosis_table;

}

void maximum_kurtosis_search(vector<vector<double>>& table) {
    //searching index of max value 
    double normality_Freq = sampling_rate / 2;
    double bandwidth_all[2 * N_LEVEL + 1] = { 0 };
    double two_way_bandwidth[N_LEVEL] = { 0 };
    double three_way_bandwidth[N_LEVEL] = { 0 };
    double max_layer = 0, max_index = 0;

    for (int i = 1; i <= N_LEVEL; i++)
    {
        two_way_bandwidth[i - 1] = normality_Freq / power(2, i);
        three_way_bandwidth[i - 1] = normality_Freq / (power(2, i - 1) * 3.0);
    }

    bandwidth_all[0] = normality_Freq;

    int count_two = 1, count_three = 2;

    for (int i = 1; i <= N_LEVEL; i++)
    {
        bandwidth_all[count_two] = two_way_bandwidth[i - 1];
        bandwidth_all[count_three] = three_way_bandwidth[i - 1];
        count_two += 2;
        count_three += 2;
    }

    for (int i = 0; i < 2 * N_LEVEL; i++)
    {
        for (int j = 0; j < base_number - 1; j++)
        {
            if (table[i][j] > max_value)
            {
                max_value = table[i][j];
                max_layer = i;
                max_index = j;

                bandwidth = bandwidth_all[i];
                center_frequency = bandwidth * j + bandwidth / 2;
            }
        }
    }
}


//////* function for acc_envelope_method *//////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////*band pass filter*//////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void bandpass_filter_init(struct str_bandpass_filter* p, int point)
{
    int n;

    p->bandpass_idx = 0;
    p->bandpass_filter_cnt = 1;
    p->bandpass_filtered = 0;
    p->bandpass_order_num = point;

    for (n = 0; n < point; n++) {
        p->bandpass_data[n] = 0;
    }
}

double bandpass_filter(struct str_bandpass_filter* p, double x, double bandpass_coef[BANDPASS_HALF_ORDER * 2 + 1], double FT1, double FT2)
{
    double bandpass_m_2;
    double bandpass_val1;
    double bandpass_val2;
    bandpass_m_2 = 0.5 * (BANDPASS_HALF_ORDER * 2.0);
    bandpass_coef[BANDPASS_HALF_ORDER] = 2 * (FT2 - FT1);
    for (int n = 0; n < BANDPASS_HALF_ORDER; n++)
    {
        bandpass_val1 = (sin(2.0 * PI * FT1 * (n - bandpass_m_2)) /
            (PI * (n - bandpass_m_2))) * (0.54 - (0.46 * (cos((2 * PI * n) / (BANDPASS_HALF_ORDER * 2.0)))));
        bandpass_val2 = (sin(2.0 * PI * FT2 * (n - bandpass_m_2)) /
            (PI * (n - bandpass_m_2))) * (0.54 - (0.46 * (cos((2 * PI * n) / (BANDPASS_HALF_ORDER * 2.0)))));
        bandpass_coef[n] = bandpass_val2 - bandpass_val1;
        bandpass_coef[BANDPASS_HALF_ORDER * 2 + 1 - n - 1] = bandpass_val2 - bandpass_val1;
    }

    int idx = p->bandpass_idx;

    p->bandpass_data[idx] = x;

    p->bandpass_filtered = 0;

    double* bandpass_ptr = &bandpass_coef[0];

    for (int i = 0; i < p->bandpass_filter_cnt; i++)
    {
        int j;
        j = p->bandpass_filter_cnt - i - 1;

        p->bandpass_filtered = p->bandpass_filtered + p->bandpass_data[j] * *(bandpass_ptr + i);
    }

    bandpass_ptr = bandpass_ptr + p->bandpass_filter_cnt;

    for (int i = 0; i < (p->bandpass_order_num - p->bandpass_filter_cnt); i++)
    {
        int j;
        j = p->bandpass_order_num - i - 1;

        p->bandpass_filtered = p->bandpass_filtered + p->bandpass_data[j] * *(bandpass_ptr + i);
    }

    p->bandpass_idx++; //update idx
    p->bandpass_filter_cnt++;

    if (p->bandpass_filter_cnt > p->bandpass_order_num)
    {
        p->bandpass_filter_cnt = 1;
    }

    if (p->bandpass_idx >= p->bandpass_order_num)
    {
        p->bandpass_idx = 0;
    }

    return(p->bandpass_filtered);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////*envelope*//////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////*squaring*//////////////////////////////
double squaring_one(double x)
{
    double y;
    y = x * x * 2;
    return y;
}

//////////////////////////////*low pass filter*//////////////////////////////

void lowpass_filter_init(struct str_lowpass_filter* p, int point)
{
    int n;

    p->lowpass_idx = 0;
    p->lowpass_filter_cnt = 1;
    p->lowpass_filtered = 0;
    p->lowpass_order_num = point;

    for (n = 0; n < point; n++) {
        p->lowpass_data[n] = 0;
    }
}

double lowpass_filter(struct str_lowpass_filter* p, double x, double lowpass_coef[LOWPASS_HALF_ORDER * 2 + 1])
{
    double lowpass_val;
    double lowpass_m_2;
    lowpass_m_2 = 0.5 * (LOWPASS_HALF_ORDER * 2.0);
    lowpass_coef[LOWPASS_HALF_ORDER] = 2 * normalized_freq;
    for (int n = 0; n < LOWPASS_HALF_ORDER; n++)
    {
        lowpass_val = ((sin(2 * PI * normalized_freq * (n - lowpass_m_2)) / (PI * (n - lowpass_m_2))));
        lowpass_coef[n] = lowpass_val;
        lowpass_coef[LOWPASS_HALF_ORDER * 2 + 1 - n - 1] = lowpass_val;
    }

    int idx = p->lowpass_idx;

    p->lowpass_data[idx] = x;

    p->lowpass_filtered = 0;

    double* lowpass_ptr = &lowpass_coef[0];

    for (int i = 0; i < p->lowpass_filter_cnt; i++)
    {
        int j;
        j = p->lowpass_filter_cnt - i - 1;

        p->lowpass_filtered = p->lowpass_filtered + p->lowpass_data[j] * *(lowpass_ptr + i);
    }


    lowpass_ptr = lowpass_ptr + p->lowpass_filter_cnt;

    for (int i = 0; i < (p->lowpass_order_num - p->lowpass_filter_cnt); i++)
    {
        int j;
        j = p->lowpass_order_num - i - 1;

        p->lowpass_filtered = p->lowpass_filtered + p->lowpass_data[j] * *(lowpass_ptr + i);
    }

    p->lowpass_idx++; //update idx
    p->lowpass_filter_cnt++;

    if (p->lowpass_filter_cnt > p->lowpass_order_num)
    {
        p->lowpass_filter_cnt = 1;
    }

    if (p->lowpass_idx >= p->lowpass_order_num) {
        p->lowpass_idx = 0;
    }

    return(p->lowpass_filtered);
}

//////////////////////////////*SQRT*//////////////////////////////
double SQRT_one(double x)
{
    double y;
    y = sqrt(fabs(x));
    return y;
}


//////////////////////////////*envelope_integration*//////////////////////////////
void envelope(vector <double> &data) {
    double* data_ptr = &data[0];
    //squaring and downsampling
    for (int i = 0; i < data_length; i++)
    {
        squaring_output = squaring_one(*(data_ptr + i));
        if (envelope_down_cnt == DOWN_SAMPLE_FACTOR)
        {
            *(data_ptr + envelope_idx_cnt) = squaring_output;
            envelope_down_cnt = 0;
            envelope_idx_cnt++;
        }
        envelope_down_cnt++;
    }
    //lowpass filter and square root
    for (int i = 0; i < down_length; i++)
    {
        lowpass_output = lowpass_filter(&comp_lowpass_filter, *(data_ptr + i), lowpass_coefficient);
        sqrt_output = SQRT_one(lowpass_output);
        *(data_ptr + i) = sqrt_output;
    }

}

/////////////////////////////////////////////////////////////////
//////////////////////////////*FFT*//////////////////////////////
/////////////////////////////////////////////////////////////////

void diftt(double* data1, int N)
{
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    int i = 0, j = 0, n = 0, k = 0, m = 0, isign = -1, istep, mmax;
    double* data;
    data = &data1[0] - 1;
    n = N * 2;
    mmax = n / 2;
    // calculate the FFT
    while (mmax >= 2) {
        istep = mmax << 1;
        theta = isign * (2 * PI / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                int debug333 = 0;
                if (j >= n)
                {
                    debug333++;
                }
                tempr = data[i];
                tempi = data[i + 1];
                data[i] = data[i] + data[j];
                data[i + 1] = data[i + 1] + data[j + 1];
                data[j] = tempr - data[j];
                data[j + 1] = tempi - data[j + 1];
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = tempr;
                data[j + 1] = tempi;
            }
            wtemp = wr;
            wr += wtemp * wpr - wi * wpi;
            wi += wtemp * wpi + wi * wpr;
        }
        mmax = mmax / 2;
    }
    // do the bit-reversal
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

bool do_FFT(vector <double> &data, int size, int sampling_freq)
{
    int j_cnt = 0;
    for (int i = 0; i < size; i++)
    {
        int j;
        j = data_length - i - 1;
        data[j - j_cnt] = data[i];
        data[j - 1 - j_cnt] = 0;

        j_cnt++;
    }
    //Reverse data from start to end
    int start = 0;
    int end = data_length - 1;
    while (start < end)
    {
        double temp = data[start];
        data[start] = data[end];
        data[end] = temp;
        start++;
        end--;
    }

    FFT.size = size;
    FFT.fs = sampling_freq;

    diftt(&data[0], size);

    return true;
}

void FFT_post_processing(vector <double> &data, int size, int sampling_freq, int down_sample_factor) {
    for (int i = 0; i < (size + 1); i++)
    {
        if (i == 0 || i == size) {
            data[i] = fabs(data[i] / (2.0 * size));
        }
        else
        {
            data[i] = 2 * fabs(data[i] / (2.0 * size));
        }
    }

    for (int i = 0; i < (size + 1); i++)
    {
        data[i + size + 1] = ((double)sampling_freq / down_sample_factor) * i / (2.0 * size);
    }
}

//////* function for bearing fault detector *//////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////*Moving Average*//////////////////////////////
////////////////////////////////////////////////////////////////////////////
vector<double> moving_avg(vector<double> amp, int l_size) {

    double sum = 0;
    vector<double> mAvg;

    for (int i = 0; i <= (amp.size() - l_size); i++) {
        sum = 0;
        //cout << "Numbers ";

        for (int j = i; j < i + l_size; j++) {
            sum += amp[j];
            //cout << amplitude[j] << " ";
        }

        mAvg.push_back(sum / l_size);
        //cout << endl << "Moving Average: " << mAvg << endl;
    }
    return mAvg;
}

/////////////////////////////////////////////////////////////////////////////
//////////////////////////////*failure_confirm*//////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void failure_confirm(vector<double> amp, int Characteristic_f_points[], int window[], int size, int sampling_freq, int down_sample_factor)
{
    vector<double>mAvg;
    vector<double>mAvg1 = moving_avg(amp, (Characteristic_f_points[1] + Characteristic_f_points[2] + Characteristic_f_points[3]) / 6);
    string Characteristic_f_name[4] = { "BSF: ", "BPFO: ", "BPFI: ", "FTF: " };
    //four Characteristic f
    diagnosis_result_flag = 0;
    for (int k = 0; k < 4; k++) {

        cout << Characteristic_f_name[k] << "\n";

        if (k < 3) {
            mAvg = mAvg1;
        }
        else {
            mAvg = moving_avg(amp, (Characteristic_f_points[0] / 2));
        }


        //int Characteristic_f_points = round(Characteristic_f[k] / (hzindex[1] - hzindex[0])); // find f points
        int ctr = 0;          // counter to count the amount of peaks


        //jump
        //cout << "window: " << window;
        for (int i = Characteristic_f_points[k]; i < Characteristic_f_points[k] * 11;) {

            //start Interval point :  num.begin() + i - window
            //end Interval point :  num.begin() + i + window

            double maxelement = amp[i - window[k]];
            double sum = 0;
            vector<double> temp;
            int pos = 0;
            for (int j = i - window[k]; j < i + window[k]; j++) {
                int revolution_freq = round(omega / 60 /
                    (((double)sampling_freq / down_sample_factor) / (2.0 * size)));
           
                if (amp[j] > maxelement && (j % revolution_freq != 0) &&
                    (j % revolution_freq != 1) && (j % revolution_freq != 2) &&
                    (j % revolution_freq != 3) && (j % revolution_freq != 4) &&
                    (j % revolution_freq != (revolution_freq - 1)) &&
                    (j % revolution_freq != (revolution_freq - 2)) &&
                    (j % revolution_freq != (revolution_freq - 3)) &&
                    (j % revolution_freq != (revolution_freq - 4))){
                    maxelement = amp[j];
                }
                pos = j;
                sum += amp[j];
                temp.push_back(amp[j]);
            }

            double maxelement_mov = mAvg[i - window[k]];
            double sum_mov = 0;
            vector<double> temp_mov;
            int pos_mov = 0;
            for (int j = i - window[k]; j < i + window[k]; j++) {
                if (mAvg[j] > maxelement_mov)
                    maxelement_mov = mAvg[j];
                pos_mov = j;
                sum_mov += mAvg[j];
                temp_mov.push_back(mAvg[j]);
            }

            cout << "pos: " << i << "  window: " << window[k]
                << "  amp: " << maxelement << "  amp_mov: " << maxelement_mov
                //<< "  avg: " << avg << "  std: " << std
                //<< "  sum: " << sum 
                << "\n";
            if (maxelement - maxelement_mov > TREE_HEIGHT * maxelement_mov) {   //threshold
                //cout << "max element: " << maxelement << "\n" << "position: " << pos << "\n";
                ctr++;
            }
            //cout <<"pos: " << i << "  amp: " << maxelement << "\n";

            if (pos - i > 0.5 * window[k]) {
                i = i + Characteristic_f_points[k];
            }
            else {
                i = pos + Characteristic_f_points[k];
            }
        }


        if (ctr >= CTR_THRESHOLD) //adjustable

        {
            cout << "ctr=  " << ctr << "\n";
            printf("fault \n\n");
            diagnosis_result_flag = 1;
        }
        else
        {
            cout << "ctr=  " << ctr << "\n";
            printf("here is not fault \n\n");
        }
    }
    if (diagnosis_result_flag == 1) {
        printf("Bearing is fault!!! \n\n");
    }
    else {
        printf("Bearing is health!!! \n\n");
    }

}

void Kurtogram_and_Acceleraton_Envelope_Method() {
    //////*kurtogram*//////
    //////////////////////////////* fool proof *//////////////////////////////
    raw_data_length = time_length * sampling_rate;
    data_length = raw_data_length - 2;
    base_number = power(2, (N_LEVEL - 1)) * 3;
    vector<vector<double>> data_test_A(row, vector<double>(raw_data_length));
    vector<vector<double>> data_test_B(row, vector<double>(raw_data_length));
    vector<vector<double>> kurtosis_table((2 * N_LEVEL + 1), vector<double>(base_number));
    vector<double> data_test(raw_data_length);
    double* data_test_ptr = &data_test[0];
    int fool_proof;
    fool_proof = log2(data_length) - 7;  //setup maximum numbers of level 

    if (N_LEVEL > fool_proof) {
        cout << ("Please enter a smaller number of decomposition levels");
    }
    else {
        //////////////////////////////* load data(position data) *//////////////////////////////
        //ifstream fin("D://summer_intern//Bearing_Diagnosis_Kurtogram_AccEnvelopeMethod_C++code_ver3//25-1_pos.txt");
        ifstream fin(input_file_path);
        // ckeck whether the file reading is successful or not
        if (!fin) {
            cout << "false to load data" << endl;
        }

        double number;
        int cnt = 0;

        while (fin >> number) {
            data_test_A[0][cnt] = number;
            data_test_A[1][cnt] = 0;
            cnt++;
        }
        fin.close();   // close file
        /////////////////////////* position transfer to acceleration */////////////////////////
        pos_to_acc(data_test_A);
#if 0
        /////////////////////////* moving acceleration signal to remove axis frequency */////////////////////////
        int moving_length;
        moving_length = (sampling_rate * 60 / omega);
        acc_moving(data_test_A, moving_length);
#endif
        for (int i = 0; i < raw_data_length; i++) {
            data_test[i] = data_test_A[0][i];
        }

        //////////////////////////////* calculation *//////////////////////////////
        double* all_kurt_source_ptr = &data_test_A[0][0];
        double* all_kurt_target_ptr = &data_test_A[1][0];

        kurtosis_table[0][0] = kurt(all_kurt_source_ptr, all_kurt_target_ptr, 0, 1, 0);


        double* kurtosis_table_ptr = 0;
        int layer_cnt = -1;
        for (int Layer = 0; Layer < N_LEVEL; Layer++)
        {
            for (int Index = 0; Index < power(2, Layer); Index++)
            {
                kurtosis_table_ptr = build_kurtosis_table(Layer, Index, data_test_A, data_test_B); // Kurtogram construction

                for (int Sub_index = 0; Sub_index < 5; Sub_index++)
                {
                    if (Layer == 0)
                    {
                        if (Sub_index == 0 || Sub_index == 1 || Sub_index == 2)
                        {
                            kurtosis_table[Layer + 2][Sub_index] = *(kurtosis_table_ptr + Sub_index);
                        }
                        else if (Sub_index == 3 || Sub_index == 4)
                        {
                            kurtosis_table[Layer + 1][Sub_index - 3] = *(kurtosis_table_ptr + Sub_index);
                        }

                    }
                    else
                    {
                        if (Sub_index == 0 || Sub_index == 1 || Sub_index == 2)
                        {
                            kurtosis_table[Layer + 3 + layer_cnt][Index * 3 + Sub_index] = *(kurtosis_table_ptr + Sub_index);
                        }
                        else if (Sub_index == 3 || Sub_index == 4)
                        {
                            kurtosis_table[Layer + 2 + layer_cnt][Index * 2 + (Sub_index - 3)] = *(kurtosis_table_ptr + Sub_index);
                        }
                    }

                }
            }
            layer_cnt++;
        }

        maximum_kurtosis_search(kurtosis_table); //maximum kurtosis search

        //////////////////////////////* write data *//////////////////////////////
#if 0
        ofstream fout("D://summer_intern//Bearing_Diagnosis_Kurtogram_AccEnvelopeMethod_C++code_ver3//25-1_kurtogram_table.txt");
        // ckeck whether the file reading is successful or not
        if (!fout)
        {
            cout << "false to write data" << endl;
        }
        //write data to file
        for (int i = 0; i < (2 * N_LEVEL + 1); i++) {
            for (int j = 0; j < base_number; j++)
            {
                fout << kurtosis_table[i][j];
                if (j == base_number - 1) {
                    fout << "\n";
                }
                else {
                    fout << " ";
                }
            }
        }

        fout << "\n";
        fout << "Maximum Kurtosis : " << max_value << endl;
        fout << "\n";
        fout << "Bandwidth : " << bandwidth << endl;
        fout << "Center Frequency : " << center_frequency;

        fout.close();   // close file
#endif
    }

    //////*acc_envelope_method*//////
    down_length = (int)data_length / DOWN_SAMPLE_FACTOR;
    transfreq_high = (center_frequency + (bandwidth / 2));
    transfreq_low = (center_frequency - (bandwidth / 2));
    ft_high = transfreq_high / sampling_rate;
    ft_low = transfreq_low / sampling_rate;
    normalized_freq = LOWPASS_CUTOFF / sampling_rate;
    ///////////////////////////////////////////////////////////
    /////////////////////*bandpass filter*/////////////////////
    ///////////////////////////////////////////////////////////
    bandpass_filter_init(&comp_bandpass_filter, BANDPASS_HALF_ORDER * 2 + 1);
    for (int i = 0; i < data_length; i++)
    {
        bandpass_output = bandpass_filter(&comp_bandpass_filter, *(data_test_ptr + i), bandpass_coefficient, ft_high, ft_low);
        *(data_test_ptr + i) = bandpass_output;
    }

    ////////////////////////////////////////////////////
    /////////////////////*envelope*/////////////////////
    ////////////////////////////////////////////////////
    envelope_down_cnt = DOWN_SAMPLE_FACTOR;
    envelope_idx_cnt = 0;
    lowpass_filter_init(&comp_lowpass_filter, LOWPASS_HALF_ORDER * 2 + 1);
    envelope(data_test);

    ///////////////////////////////////////////////
    /////////////////////*FFT*/////////////////////
    ///////////////////////////////////////////////
    fft_max_size = down_length;//data size is a power if 2 which is near raw data size after downsampling(data_length/DOWN_SAMPLE_FACTOR)
    int size = 1 << ((int)log2(down_length));
    do_FFT(data_test, size, sampling_rate);

    // spectrum post processing(frequency and amplitude calculation)
    FFT_post_processing(data_test, FFT.size, sampling_rate, DOWN_SAMPLE_FACTOR);


    //////////////////////////////*write data for amplitude index*//////////////////////////////
    vector<double> amplitude;
#if 0
    ofstream fout_amp("D://summer_intern//Bearing_Diagnosis_Kurtogram_AccEnvelopeMethod_C++code_ver6//FFT_amplitude_index.txt");

    // ckeck whether the file reading is successful or not

    if (!fout_amp)
    {
        cout << "false to write data" << endl;
    }
#endif
    // #0 ~ #2*(2^17-1)+1 of data_test is the output of FFT
#if 0
    for (int i = 0; i < FFT.size; i++)
    {
        fout_amp.precision(15); //improve accuracy in fft 
        fout_amp << data_test[2 * i] << " ";
        fout_amp << data_test[2 * i + 1] << " ";
    }
#else
    for (int i = 0; i < (FFT.size + 1); i++)
    {
#if 0
        fout_amp.precision(15); //improve accuracy in fft 
        fout_amp << data_test[i] << " ";
#endif
        double m = data_test[i];
        amplitude.push_back(m);
    }
#endif
    //write data to file
#if 0
    fout_amp.close();   // close file
#endif


    //////////////////////////////*write data for frequency index*//////////////////////////////
    vector<double> hzindex;
#if 0
    ofstream fout_freq("D://summer_intern//Bearing_Diagnosis_Kurtogram_AccEnvelopeMethod_C++code_ver6//FFT_frequency_index.txt");

    // ckeck whether the file reading is successful or not

    if (!fout_freq)
    {
        cout << "false to write data" << endl;
    }
#endif

    for (int i = (FFT.size + 1); i < (2 * (FFT.size + 1)); i++)
    {
#if 0
        fout_freq.precision(15); //improve accuracy in fft 
        fout_freq << data_test[i] << " ";
#endif   
        double n = data_test[i];
        hzindex.push_back(n);
    }
    //write data to file

#if 0
    fout_freq.close();   // close file
#endif

    //////////////////////////////////////////////////////////////////////////////////
    ////////////////////////*identify fault freq. in spectrum*////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    double Characteristic_f[4]; // charactermatic f {BPFI , BPFO , BSF, FTF }
    int Characteristic_f_points[4];
    int window[4];


    Characteristic_f[0] = (pd / bd) * (omega / 60.0) * (1.0 - pow(bd / pd * cos(theta * PI / 180.0), 2));  //BSF
    Characteristic_f[1] = 0.5 * n * omega / 60.0 * (1.0 - bd * cos(theta * PI / 180.0) / pd);  //BPFO
    Characteristic_f[2] = 0.5 * n * omega / 60.0 * (1.0 + bd * cos(theta * PI / 180.0) / pd);  //BPFI
    Characteristic_f[3] = 0.5 * omega / 60.0 * (1.0 - bd * cos(theta * PI / 180.0) / pd);  //FTF

    Characteristic_f_points[0] = round(Characteristic_f[0] / (hzindex[1] - hzindex[0]));
    Characteristic_f_points[1] = round(Characteristic_f[1] / (hzindex[1] - hzindex[0]));
    Characteristic_f_points[2] = round(Characteristic_f[2] / (hzindex[1] - hzindex[0]));
    Characteristic_f_points[3] = round(Characteristic_f[3] / (hzindex[1] - hzindex[0]));

    window[0] = round(Characteristic_f[0] * TOLERANCE / (hzindex[1] - hzindex[0]));
    window[1] = round(Characteristic_f[1] * TOLERANCE / (hzindex[1] - hzindex[0]));
    window[2] = round(Characteristic_f[2] * TOLERANCE / (hzindex[1] - hzindex[0]));
    window[3] = round(Characteristic_f[3] * TOLERANCE / (hzindex[1] - hzindex[0]));


    cout << "BSF = " << Characteristic_f[0] << "  start from " << Characteristic_f_points[0] << "\n"
        << "BPFO = " << Characteristic_f[1] << "  start from " << Characteristic_f_points[1] << "\n"
        << "BPFI = " << Characteristic_f[2] << "  start from " << Characteristic_f_points[2] << "\n"
        << "FTF = " << Characteristic_f[3] << "  start from " << Characteristic_f_points[3] << "\n";

    failure_confirm(amplitude, Characteristic_f_points, window, FFT.size, sampling_rate, DOWN_SAMPLE_FACTOR);

    // extra engineering information
    string s0 = "ver6,";
    string s1 = to_string(max_value);
    string s2 = to_string(bandwidth);
    string s3 = to_string(center_frequency);
    s0.append(s1);
    s0.append(",");
    s0.append(s2);
    s0.append(",");
    s0.append(s3);
    log_text = s0.c_str();
    /*cout << log_text << "\n"
        << diagnosis_result_flag <<"\n";*/
}


