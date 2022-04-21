/* Example for using predefined function to diagnosis the bearing */

//extern predefined function
extern void set_bearing_spec(double bn, double ball_diameter, double pitch_diameter, double contact_angle, double motor_speed);
extern void set_exp_condition(int time_len, int samplerate);
extern void set_input_filepath(const char* argv_in);

extern void Kurtogram_and_Acceleraton_Envelope_Method();

                            /////////////////////////////////////////
                            ////////////* Main Function *////////////
                            /////////////////////////////////////////
int main()
{
    //////*parameters setting*//////

        //600RPM front:(10, 4.762, 26, 11.055, 600) back:(7, 4.762, 18, 11.055, 600); 25-1(6, 3.5, 12.9, 30, 5000)
    set_bearing_spec(6, 3.5, 12.9, 30, 5000); //1.bn 2.bd 3.pd 4.angle 5.speed
    set_exp_condition(20, 20000); //1.time_length(sec.) 2.sampling_rate(per sec.)
    set_input_filepath("D://summer_intern//Bearing_Diagnosis_Kurtogram_AccEnvelopeMethod_C++code_ver6//25-1_pos.txt");

    //////*calculation and outputfile generation*//////
    Kurtogram_and_Acceleraton_Envelope_Method();

    return 0;

}