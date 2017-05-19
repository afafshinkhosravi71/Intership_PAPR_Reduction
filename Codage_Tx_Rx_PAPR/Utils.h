//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
// IEEE 802.11 OFDM PHYSICAL LAYER
// Author : Boubacar DIALLO
// Date   : May 2017
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <itpp/itcomm.h>
#include <itpp/itstat.h>
#include <iostream>
#include <vector>
#include <string>
#include <iterator>
#include <fstream>

using namespace std;
using namespace itpp;
using std::cin;
using std::cout;
using std::endl;

void data_rate_choice(double* R, int* Nbpsc, int* data_rate,int *Ncbps, int *Ndbps);
bvec Punct_Conv_Code(const bvec &uncoded_bits, double &R);
bvec Punct_Conv_Decode(const vec &received_symbol, double &R);
cvec modulation(const bvec &input_bits, int Nbpsc,  ivec &deci, cvec &constl);
vec soft_demodulation(const cvec &received_symb, const int Nbpsc, const double &N0);
bvec Hard_demodulation(const cvec &received_symb, const int Nbpsc);
bvec scrambler(const bvec& input,const bvec& State);
bvec descrambler(const bvec& input,const bvec& State);
bvec interleaver(const bvec &coded_bits, int &Ncbps,int &Nbpsc);
vec deinterleaver(const vec &received_symbol, int &Ncbps,int &Nbpsc);
cvec preamble(void)	;
cvec Signal(int &len,int &data_rate);
bvec Data(const bvec& input_data, int &Ncbps, int &Ndbps);
int Pilot_Polarity_Sequence(const int index);
cvec ofdm_modulation(const  cvec &input, const int &Nfft ,const int &Ncp, const int & norm_factor_flag);
cvec ofdm_demodulation(const  cvec &input, const int &Nfft ,const int &Ncp, const int & norm_factor_flag);
cvec Pilot_insertion( const  cvec &input, const int &N_SD ,const int &Nfft,int &index);
cvec Pilot_Data_separation( const  cvec &input, const int &N_SD ,const int &Nfft,int &index);
cvec window(const cvec &input);			//windowing function
cvec OFDM_Data_Symbols_Assembler(const cvec &data_ofdm_transmited_symb);
cvec PPDU_Frame_Assembler(const cvec &PREAMBLE, const cvec &SIGNAL, const cvec &DATA);
void  PPDU_Frame_Deassembler(const cvec &ofdm_received_symb,  cvec &received_short_squence,  cvec &received_long_squence,  cvec &received_signal_field,  cvec &received_data_field);
ivec Receiver_data_rate(const cvec &received_signal_field,const double &N0,const cvec Hhat);
cvec LS_estimator(const cvec &rec_long_squence_CP);
cvec MMSE_estimator(const cvec &rec_long_squence_CP);
cmat Fmatrix(const int N1,const int N2);
double noise_variance(const cvec &rec_long_squence_CP);
cvec Equalization(const cvec &input,const cvec &Hhat);
bvec file2binary_converter(void);
void binary2file_converter(const bvec &binary_data);
cvec frames_nbr_field(int &len);
int frame_nbr_receiver(const cvec &received_frame_nbr_field,const double &N0,const cvec Hhat);
bvec text_example(void);
void BER_SNR(char* file_name, vec EbN0dB, vec bit_error_rate, int mode);
void Sent_Frame_Saver(cvec ofdm_transmited_frame);
int chain_length(char c[32]);
void header_display(void);

//********************************************************** Punctured Convolutional Coder and Decoder ****************************************************//
bvec Punct_Conv_Code(const bvec &uncoded_bits, double &R)
{
	  	Punctured_Convolutional_Code code;
		ivec generator(2);
		generator(0)=0133;				//polynomial used in 802.11
		generator(1)=0171;
		code.set_generator_polynomials(generator, 7);	// Constraint length K =7
		code.set_method(Tail);				//that the encoder starts in the zero state
		bvec coded_bits;				//coded bits
		//code.reset();

		// R  code rates 1/2 , 3/4 and 2/3
	  
		int sel;	//selection 
		if (R==1/2.0) sel = 1;	if (R==3/4.0) sel = 2;	if (R==2/3.0) sel = 3;

		int tail;		// varialble contains the tail size that has been added by the Punctured Convolutional Code
		bmat puncture_matrix;	//puncture matrix 

		switch (sel)
		{	
			case 1 :{	puncture_matrix = "1; 1";
					    tail=12;
					}
			break;
			case 2 : {	puncture_matrix = "1 1 0; 1 0 1";
					    tail=8;
					 }
			break;
			case 3 : {	puncture_matrix = "1 1 1 1 1 1; 1 0 1 0 1 0";
						tail=9;
					 }
			break;
			default:cout<<"Not valid Coding Rate"<<endl;
            exit(-1);
		}
		
		code.set_puncture_matrix(puncture_matrix);		//cout<<"puncture matrix"<<code.get_puncture_matrix();
		coded_bits=code.encode(uncoded_bits);
		return  coded_bits.mid(0,coded_bits.length()-tail); 	//removes the tail 
}
//****************************************************   Decoder  ********************************************************/
bvec Punct_Conv_Decode(const vec &received_symbol, double &R)
{
	  	Punctured_Convolutional_Code code;
		ivec generator(2);
		generator(0)=0133;//polynomial used in 802.11
		generator(1)=0171;
		code.set_generator_polynomials(generator, 7);	// Constraint length K =7
		code.set_method(Tail);							// encoder starts in the zero state
		bvec coded_bits;								//coded bits
		// cout<<code.get_encoder_state() ;
		// R  code rates 1/2 , 3/4 and 2/3

		int sel;	//selection 
		if (R==1/2.0) sel = 1;	if (R==3/4.0) sel = 2;	if (R==2/3.0) sel = 3;

		int tail;	// varialble contains the tail size that has been added by the Punctured Convolutional Code
		bmat puncture_matrix;	//puncture matrix 

		switch (sel)
		{	
			case 1 :{	puncture_matrix = "1; 1";
						tail=12;
					}
			break;
			case 2 : {	puncture_matrix = "1 1 0; 1 0 1";
						tail=8;
					 }
			break;
			case 3 : {	puncture_matrix = "1 1 1 1 1 1; 1 0 1 0 1 0";
						tail=9;
					 }
			break;
			default:cout<<"Not valid Coding Rate"<<endl;
            exit(-1);
		}
		vec tail_zeros(tail);
		tail_zeros.zeros();
		code.set_puncture_matrix(puncture_matrix);						//cout<<"punctured matrix"<<code.get_puncture_matrix();
		bvec received_decoded_bits=code.decode(concat(zeros(12),received_symbol,tail_zeros));	//adding the tail
		return  received_decoded_bits(int(12*R),-1); 						//added zeros just to make sure that the decoder begins in the zero state. 
}
//***************************************************  Modulation and Demodulation ********************************************************//

cvec modulation(const bvec &input_bits, int Nbpsc,  ivec &deci, cvec &constl)
{	
	cvec temp( input_bits.length() / Nbpsc ); 

	//BPSK
	if(Nbpsc == 1)	
	{
		for (int i = 0; i < input_bits.length(); i++)  
		{temp(i) = (input_bits(i) == 0 ? -1.0 : 1.0);
		}
		constl = "-1 1"; 
		deci = "0 1 "; 
	
	}
	//QPSK
	if(Nbpsc == 2)
	{
		vec qpsk_r = "-1 -1 1 1"; 
		vec qpsk_i = "-1 1 -1 1"; 
		ivec cont_d = "0:1:3";	//bits2symbols
		Modulator_2D qam; 
		deci = cont_d; 
		cvec comple;			//symbols
		comple.set_length(deci.length(), false); 
		for(int i = 0; i < deci.length(); i++)
			comple(i) = std::complex<double>( (double)qpsk_r(i), (double)qpsk_i(i))/sqrt((double)2 );
		qam.set(comple, cont_d);		//Set the constellation to use in the modulator. 
		constl = qam.get_symbols();		 //Get the symbol values used in the modulator.
		temp = qam.modulate_bits(input_bits); 
	}

	//16-QAM
	if(Nbpsc == 4)
	{
		cvec comple; 
		Modulator_2D qam; 
		//ivec cont_d = "0:1:15";  //bits2symbols
		ivec cont_d = "0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10";// to follow the mapping prposed by the standard
		deci = cont_d; 
		vec qam_r = "-3 -3 -3 -3 -1 -1 -1 -1  1  1  1  1  3  3  3  3";
		vec qam_i = "-3 -1  1  3 -3 -1  1  3 -3 -1  1  3 -3 -1  1  3";
		//vec qam_r = "1 -1  1 -1 3 -3  3 -3 1 -1  1 -1 3 -3  3 -3";  // mapping TRILION
		//vec qam_i = "1  1 -1 -1 1  1 -1 -1 3  3 -3 -3 3  3 -3 -3";  // mapping TRILION
		comple.set_length(qam_i.length(),false); 
		for(int i = 0; i < qam_i.length(); i++)
			comple(i) = std::complex<double>( (qam_r(i)/sqrt((double)10)) , (qam_i(i)/sqrt((double)10)) );
		qam.set(comple, cont_d);	//Set the constellation to use in the modulator. 
		constl = qam.get_symbols();  //Get the symbol values used in the modulator.
		temp = qam.modulate_bits( input_bits ); 
	}
	//64-QAM
	if(Nbpsc == 6)
	{
		cvec comple; 
		Modulator_2D qam; 
		ivec cont_d = "0:1:63";		//bits2symbols
		deci = cont_d; 
		vec qam_r = "-7 -7 -7 -7 -7 -7 -7 -7 -5 -5 -5 -5  -5 -5 -5 -5 -3 -3 -3 -3 -3 -3 -3 -3 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1 1  3  3  3  3  3  3  3  3  5  5  5  5  5  5  5 5  7  7  7  7  7  7  7  7";
		vec qam_i = "-7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1   1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7	 -5 -3 -1  1  3  5  7";
		comple.set_length(qam_i.length(),false);
		for(int i=0; i<qam_i.length(); i++)
			comple(i) = std::complex<double>(qam_r(i)/sqrt((double)42),qam_i(i)/sqrt((double)42));
		qam.set(comple, cont_d);	//Set the constellation to use in the modulator. 
		constl= qam.get_symbols();  //Get the symbol values used in the modulator.
		temp = qam.modulate_bits( input_bits ); 
	}
	return temp;
}

/*****************************************************     Demodulation    *****************************************************************/
vec soft_demodulation(const cvec &received_symb, const int Nbpsc, const double &N0)
{	
	//Soft demodulation 
	vec temp; 

	if(Nbpsc == 1)
	{
	BPSK_c cbpsk;//Complex cbpsk constructor
	temp= cbpsk.demodulate_soft_bits (received_symb,N0);
	}
	if(Nbpsc == 2)
	{
		vec qpsk_r = "-1 -1 1 1"; 
		vec qpsk_i = "-1 1 -1 1"; 
		Modulator_2D qam; 
		ivec cont_d = "0:1:3";  
		cvec comple; 
		comple.set_length(cont_d.length(), false); 
		for(int i = 0; i < cont_d.length(); i++)
			comple(i) = std::complex<double>( (double)qpsk_r(i), (double)qpsk_i(i))/sqrt((double)2 );		
		qam.set(comple, cont_d); 
		temp = qam.demodulate_soft_bits(received_symb, N0) ;
	}
	if(Nbpsc == 4)
	{
		cvec comple; 
		Modulator_2D qam; 
		//ivec cont_d = "0:1:15";  
		ivec cont_d = "0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10";// to follow the mapping prposed by the standard
		vec qam_r = "-3 -3 -3 -3 -1 -1 -1 -1  1  1  1  1  3  3  3  3";
		vec qam_i = "-3 -1  1  3 -3 -1  1  3 -3 -1  1  3 -3 -1  1  3";
		comple.set_length(qam_i.length(),false); 
		for(int i=0; i<qam_i.length(); i++)
			comple(i) = std::complex<double>( (qam_r(i)/sqrt((double)10)) , (qam_i(i)/sqrt((double)10)) );
		qam.set(comple, cont_d); 
		temp = qam.demodulate_soft_bits(received_symb, N0 ); 
	}
	if(Nbpsc == 6)
	{
		cvec comple; 
		Modulator_2D qam; 
		ivec cont_d = "0:1:63";
		vec qam_r = "-7 -7 -7 -7 -7 -7 -7 -7 -5 -5 -5 -5  -5 -5 -5 -5 -3 -3 -3 -3 -3 -3 -3 -3 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1 1  3  3  3  3  3  3  3  3  5  5  5  5  5  5  5 5  7  7  7  7  7  7  7  7";
		vec qam_i = "-7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1   1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7 -5 -3 -1  1  3  5  7";
		comple.set_length(qam_i.length(),false); 
		for(int i=0; i<qam_i.length(); i++)
			comple(i) = std::complex<double>(qam_r(i)/sqrt((double)42),qam_i(i)/sqrt((double)42));
		qam.set(comple, cont_d); 
		temp = qam.demodulate_soft_bits(received_symb, N0);  
	}
	return temp; 
}

//Hard Demodulation
bvec Hard_demodulation(const cvec &received_symb, const int Nbpsc)
{
	cvec comple; 
	Modulator_2D qam;
	bvec temp;

	if(Nbpsc == 1)	
	{	temp.set_length(received_symb.length());
		BPSK_c bpsk;
		temp=bpsk.demodulate_bits(received_symb);
	}

	if(Nbpsc == 2)
	  { 
		ivec cont_d = "0:1:3";
		vec qam_r = "-1 -1 1 1"; 
	 	vec qam_i = "-1 1 -1 1";  
		comple.set_length(qam_i.length(),false);
		for(int i = 0; i < qam_i.length(); i++)
	    		comple(i) = std::complex<double>( (double)qam_r(i), (double)qam_i(i))/sqrt((double)2 );
		qam.set(comple, cont_d); 
		temp = qam.demodulate_bits(received_symb); 
	  }

 	if(Nbpsc == 4)
	  {
		//ivec cont_d = "0:1:15"; 
		ivec cont_d = "0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10";// to follow the mapping prposed by the standard
		vec qam_r = "-3 -3 -3 -3 -1 -1 -1 -1  1  1  1  1  3  3  3  3";
		vec qam_i = "-3 -1  1  3 -3 -1  1  3 -3 -1  1  3 -3 -1  1  3";
		//vec qam_r = "1 -1  1 -1 3 -3  3 -3 1 -1  1 -1 3 -3  3 -3";  // mapping TRILION
		//vec qam_i = "1  1 -1 -1 1  1 -1 -1 3  3 -3 -3 3  3 -3 -3";  // mapping TRILION
		comple.set_length(qam_i.length(),false); 
		for(int i = 0; i < qam_i.length(); i++)
			comple(i) = std::complex<double>( (qam_r(i)/sqrt((double)10)) , (qam_i(i)/sqrt((double)10)) );
		qam.set(comple, cont_d); 
		temp = qam.demodulate_bits(received_symb); 
	  }
	
	if(Nbpsc == 6)
	  { 
		ivec cont_d = "0:1:63";
		vec qam_r = "-7 -7 -7 -7 -7 -7 -7 -7 -5 -5 -5 -5  -5 -5 -5 -5 -3 -3 -3 -3 -3 -3 -3 -3 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1  1  1  1 1  3  3  3  3  3  3  3  3  5  5  5  5  5  5  5 5  7  7  7  7  7  7  7  7";
		vec qam_i = "-7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1   1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7 -5 -3 -1  1  3  5  7 -7 -5 -3 -1  1  3  5 7 -7 -5 -3 -1  1  3  5  7";
		comple.set_length(qam_i.length(),false); 
		for(int i=0; i<qam_i.length(); i++)
			comple(i) = std::complex<double>(qam_r(i)/sqrt((double)42),qam_i(i)/sqrt((double)42));
		qam.set(comple, cont_d); 
		temp = qam.demodulate_bits(received_symb);  
	  }
	
	return temp; 
}

//************************************************* Interleaver and deinterleaver *****************************************************************//
bvec interleaver(const bvec &coded_bits, int &Ncbps,int &Nbpsc)
{		
	//Ncbps Coded bits per OFDM symbol
	//Nbpsc Coded bits per subcarrier 
	vec k;	//index before the first permutation shall be denoted by k;
	vec i;	//after the first and before the second permutation
	vec j;	//index after the second permutation, just prior to modulation mapping.
	
	k=linspace(0,Ncbps-1,Ncbps);//cout << k<<endl;
	//bvec First_output,Second_output;

	vec kmod(k.length()),jmod(k.length()),imod(k.length());//k modulus 16

	for (int p=0;p<k.length();p++)
		kmod(p)=mod ( k(p), 16);

	i = (Ncbps/16)*kmod + floor(k/16);//cout << i<<endl;
	
	vec s1(2); s1(0)= Nbpsc/2.0; s1(1)=1;
	double s = max(s1);

	for (int p=0;p<k.length();p++)
		jmod(p)= mod ( i(p)-floor(16 * i(p)/Ncbps)+Ncbps, s);

	j = s * floor(i/s) + jmod;//cout << j<<endl;

	bvec interleaved_bits(coded_bits.length());
	
	 for (int p = 0; p < coded_bits.length(); p=p+Ncbps)
		  {
			for (int m = 0; m < Ncbps; m++)
				interleaved_bits(j(m)+p) =coded_bits(k(m)+p);
	          	//cout <<interleaved_bits(p,p+Ncbps-1)<<endl;//ok //tested and compared L.1.4.3 Interleaving the SIGNAL field bits
		 }

	////first permutation 
	//Sequence_Interleaver<bin> FirstPermutation;//ok
	//FirstPermutation.set_interleaver_sequence(to_ivec(i));

	////second permutation 
	//Sequence_Interleaver<bin> SecondPermutation;//
	//SecondPermutation.set_interleaver_sequence(to_ivec(j));

	//First_output = FirstPermutation.interleave(coded_bits);
	//Second_output= SecondPermutation.interleave(First_output);

	return interleaved_bits;//interleaved data
}

//deinterlever
vec deinterleaver(const vec &received_symbol, int &Ncbps,int &Nbpsc)
{		
	//Ncbps Coded bits per OFDM symbol
	//Nbpsc Coded bits per subcarrier 
	vec k;	//index before the first permutation shall be denoted by k;
	vec i;	//after the first and before the second permutation
	vec j;	//index after the second permutation, just prior to modulation mapping.
	
	k=linspace(0,Ncbps-1,Ncbps);

	vec kmod(k.length()),jmod(k.length()),imod(k.length());

	for (int p=0;p<k.length();p++)
		kmod(p)=mod ( k(p), 16);

	i = (Ncbps/16)*kmod + floor(k/16);//cout << i;
	
	vec s1(2); s1(0)= Nbpsc/2; s1(1)=1;
	
	double s = max(s1);

	for (int p=0;p<k.length();p++)
		jmod(p)= mod ( i(p)-floor(16 * i(p)/Ncbps)+Ncbps, s);

	j = s * floor(i/s) + jmod;//cout << j;
	
	vec input_de, First_output_de,Second_output_de;

	for (int p=0;p<k.length();p++)
		imod(p)= mod ( j(p)+floor(16 * j(p)/Ncbps), s);

	i = s * floor(j/s) + imod;//cout << i;

	k = 16 * i - (Ncbps - 1)*floor(16 * i/Ncbps);//cout << k;

	vec deinterleaved_bits(received_symbol.length());

	for (int p = 0; p < received_symbol.length(); p=p+Ncbps)
		  {
			 for (int m = 0; m < Ncbps; m++)
			 deinterleaved_bits(k(m)+p) =received_symbol(j(m)+p);
	        // cout <<deinterleaved(p,p+Ncbps-1)<<endl;//ok //tested and compared with table Table L-8\97SIGNAL field bits after encoding
		 }

	////first permutation 
	//Sequence_Interleaver<double> FirstPermutation_de;//
	//FirstPermutation_de.set_interleaver_sequence(to_ivec(j));

	////second permutation 
	//Sequence_Interleaver<double> SecondPermutation_de;
	//SecondPermutation_de.set_interleaver_sequence(to_ivec(i));

	//First_output_de = FirstPermutation_de.deinterleave(received_symbol);
	//Second_output_de= SecondPermutation_de.deinterleave(First_output_de);

	return deinterleaved_bits;
}
//************************************************** Scrambler and Descrambler *********************************************************************//
bvec scrambler(const bvec& input,const bvec& State)
	{
        bvec Register(State); 		//cout<<Register;
        bvec output(input.length());	//cout<<size(input);
        bin FirstBit;
        for (int i = 0; i < size(input); i++) {
                FirstBit = Register(3) ^ Register(6);		//bitwise xOR (reg 4  xor reg 7)
                output(i) = FirstBit ^ input(i);		//bitwise xOR (reg 4  xor reg 7 xor input)
                Register.shift_right(FirstBit);			//shifting right by one bit and inserting reg 4 xor reg 7
         }
        return output;
	}

//Descrambler

bvec descrambler(const bvec& input,const bvec& State)
	{
      	return scrambler(input,State);				// the same scrambler is used for descrambling the data
	}
//************************************************* Preamble Generator *****************************************************************//
cvec preamble(void)	
{
	//cout << "/************* Preamble genered *******************/" << endl;
	// Generation of the short sequences
	cvec Short_training_seq(53);       //including the value 0 at dc
	Short_training_seq.zeros();

	//Short -26,26 = sqrt(13/6) * {0, 0, 1+j, 0, 0, 0, -1-j, 0, 0, 0, 1+j, 0, 0, 0, -1-j, 0, 0, 0, \961\96j, 0, 0, 0, 1+j, 0, 0, 0, 0, 0, 0, 0, \961\96j, 0, 0, 0, \961\96j, 0, 0, 0, 1+j, 0, 0, 0, 1+j, 0, 0, 0, 1+j, 0, 0, 0, 1+j, 0,0}
	Short_training_seq(2)=std::complex<double>(1, 1);
	Short_training_seq(6)=std::complex<double>(-1, -1);
	Short_training_seq(10)=std::complex<double>(1, 1);
	Short_training_seq(14)=std::complex<double>(-1, -1);
	Short_training_seq(18)=std::complex<double>(-1, -1);
	Short_training_seq(22)=std::complex<double>(1, 1);
	Short_training_seq(30)=std::complex<double>(-1, -1);
	Short_training_seq(34)=std::complex<double>(-1, -1);
	Short_training_seq(38)=std::complex<double>(1, 1);
	Short_training_seq(42)=std::complex<double>(1, 1);
	Short_training_seq(46)=std::complex<double>(1, 1);
	Short_training_seq(50)=std::complex<double>(1, 1);

	//Add the zeros to have 64-symbols 
	cvec a(6);a.zeros();
	cvec b(5);b.zeros();

	Short_training_seq=sqrt(13/6.0)*concat(a,Short_training_seq,b);
	cvec Short_training_seq_time;
	cvec Short_training_seq_time_extended;

	//the coefficients 1 to 26 are mapped to the same numbered IFFT inputs, while the coefficients \9626 to \961 are copied into IFFT inputs 38 to 63
	Short_training_seq_time=ifft(concat(Short_training_seq.right(32),Short_training_seq.left(32)),64);
	
	//the short training sequence is extended periodically for 161 samples
	Short_training_seq_time_extended=concat(Short_training_seq_time,Short_training_seq_time,Short_training_seq_time.left(33));
	
	Short_training_seq_time_extended=window(Short_training_seq_time_extended);	//window function

	cout << "Mean value of transmitted short training seq = " << mean(abs(Short_training_seq_time_extended)) << endl;
	//Random test : tested with Table L-4\97
	//cout << "Time domain representation of the short sequence " << endl;
	//for (int i=0; i<10;i++)
	//	cout << i*11+1 << " --> " << Short_training_seq_time_extended(i*11+1) << endl;
	

	//Generation of the long sequences
	cvec Long_training_seq(64);
	cvec Long_training_seq_time(64);
	cvec Long_training_seq_time_extended(161);

	Long_training_seq = " 0 0 0 0 0 0  1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0  "; 
	
	//the coefficients 1 to 26 are mapped to the same numbered IFFT inputs, while the coefficients \9626 to \961 are copied into IFFT inputs 38 to 63	
	Long_training_seq_time=ifft(concat(Long_training_seq.right(32),Long_training_seq.left(32)),64);
	Long_training_seq_time_extended=concat(Long_training_seq_time.right(32),Long_training_seq_time,Long_training_seq_time);	//CP + 2 long sequences

	//adding element 161 which is the first sample of the non CP part of the OFDM symbol
	Long_training_seq_time_extended=concat(Long_training_seq_time_extended,Long_training_seq_time(0));
	Long_training_seq_time_extended=window(Long_training_seq_time_extended);						//windowing
	
	//random test : tested with Table L-6\97
	//cout << "Time domain representation of the long sequence" << endl;
	//for (int i=0; i<10;i++)
	//	cout << i*12 << " --> " << Long_training_seq_time_extended(i*12) << endl;
	

	//Generation of the preamble
	//161+161-1 (overlapping and adding element 160 of the short sequence to element 0 ofthe long sequence)
	cvec preamble(321);
	cvec overlap(1);
	overlap=Short_training_seq_time_extended(160)+Long_training_seq_time_extended(0);
	preamble=concat(Short_training_seq_time_extended.left(160),overlap,Long_training_seq_time_extended.right(160));
	
	//tested ok
	//cout<< "Preamble generated " << preamble(159,163) << endl;
        //cout<< "Preamble length " << preamble.length() << endl;
	
	//cout << "/************* Preamble genered : Done *******************/" << endl;
	return preamble;
}
//******************************************************    windowing function     ****************************************************************//
cvec window(const cvec &input)//windowing function
	{	
		cvec output(size(input));
		output=input;
		output(0)=.5*input(0);
		output(input.length()-1)=.5*input(input.length()-1);
		return output;
	}
//******************************************************      Signal field       *******************************************************************//
cvec Signal(int &len,int &data_rate)
{
	//len is the length of bytes of the PSDU 
	
	//Signal field 
	bvec Rate(4);
	bvec Reserved(1);Reserved.zeros();
	bvec leng(12);leng.zeros();
	bvec parity(1);
	bvec tail(6);tail.zeros();
	bvec signal(24);signal.zeros();

	switch(data_rate)
	{
		case 6 :Rate="1 1 0 1";
		break;
		case 9 :Rate="1 1 1 1";
		break;
		case 12:Rate="0 1 0 1";
		break;
		case 18:Rate="0 1 1 1";
		break;
		case 24:Rate="1 0 0 1";
		break;
		case 36 :Rate="1 0 1 1";
		break;
		case 48 :Rate="0 0 0 1";
		break;
		case 54 :Rate="0 0 1 1";
		break;

	}

	//cout<<len<<"  "<<endl;
	leng=dec2bin (12, len);//convert the length of needed bytes of PSDU into the binary representation

	//transform the length (from LSB to MSB) to follow the order imposed by the standard
	bvec temp(leng.length());
	for( int i=0; i<leng.length();i++)
		temp(leng.length()-1-i)=leng(i);
	leng=temp;

	//Parity bit : even parity 
	if ((sum(to_ivec(Rate))+sum(to_ivec(Reserved))+sum(to_ivec(leng)))%2==0)
		parity="0";
	else
		parity="1";
	//Even parity functionality Test
	//cout<<sum(to_ivec(Rate))+sum(to_ivec(Reserved))+sum(to_ivec(leng));  

	signal = concat(Rate,Reserved,leng,parity,tail);
	//cout<<"signal" <<signal<<endl;//test the signal field // ok tested and compared with table L-7 Bit assignment for SIGNAL field
	
	//Coding the signal field
	bvec coded_signal_bits;
	double R =1/2.0;
	coded_signal_bits= Punct_Conv_Code(signal,R);		//The bits are encoded by the rate 1/2 convolutional encoder
	//cout <<coded_signal_bits;				// ok tested and compared with Table L-8\97 SIGNAL field bits after encoding

	//interleaving
	bvec interleaved_signal_bits;
	int Ncbps1=48;
	int Nbpsc1=1;
	/*bvec temp2(coded_signal_bits.length());
	for( int i=0; i<coded_signal_bits.length();i++)
		temp2(coded_signal_bits.length()-1-i)=coded_signal_bits(i);
	coded_signal_bits=temp2;*/

	interleaved_signal_bits=interleaver(coded_signal_bits, Ncbps1,Nbpsc1);
	//cout <<interleaved_signal_bits; // ok tested and compared with Table L-9\97SIGNAL field bits after interleaving
	
	//Modulation BPSK
	cvec modulated_signal_bits;
	ivec deci;cvec constl;
	modulated_signal_bits= modulation(interleaved_signal_bits, Nbpsc1,deci,constl);
	//for(int i = 0;  i < modulated_signal_bits.length(); i++)
	//cout <<interleaved_signal_bits(i)<<"\t"<<modulated_signal_bits(i)<<endl;
	//tested ok and compared with Table L-10\97Frequency domain representation of SIGNAL field
	//for (int i = 0; i < modulated_signal_bits.length(); i++)  cout<<i<<"  "<<modulated_signal_bits(i)<<endl;

	//Pilot insertion and OFDM modulation
	int N_SD=48;
	int Nfft=64;
	int Ncp=16;
	int index=0;		//to choose first element of the polarity sequence p0-126
	cvec signal_ofdm_freq;
    	signal_ofdm_freq= Pilot_insertion( modulated_signal_bits,N_SD ,Nfft,index);
	//for(int i = 0;  i < ofdm_signal_time.length(); i++)
	//cout <<-32+i<<"\t"<<ofdm_signal_time(i)<<endl;
	//tested ok and compared with Table L-11\97Frequency domain representation of SIGNAL field with pilots inserted
	//for (int i = 0; i < signal_ofdm_freq.length(); i++)  cout<<i<<"  "<<signal_ofdm_freq(i)<<endl;

	cvec signal_ofdm_modulated=ofdm_modulation(signal_ofdm_freq,Nfft ,Ncp,1);		//for real transmission flag = 1
	//cvec signal_ofdm_modulated=ofdm_modulation(signal_ofdm_freq,Nfft ,Ncp,0);
	//for testing purposes flag = 0 to have the results before the normalization in order to compare the results with tables in annex L
	//cout<<signal_ofdm_modulated.left(10)<<endl;

	signal_ofdm_modulated.ins(signal_ofdm_modulated.length(),signal_ofdm_modulated(16));   		// add the element 81 with the first symbol after the CP
	signal_ofdm_modulated= window(signal_ofdm_modulated);						//window function
	//for(int i = 0;  i < signal_ofdm_modulated.length(); i++)
	//cout <<i<<"\t"<<signal_ofdm_modulated(i)<<endl;
	//tested ok and compared with Table L-12\97 Time domain representation of SIGNAL field
	
	return signal_ofdm_modulated;
}
//******************************************************       Data field        *******************************************************************//
bvec Data(const bvec &input_data, int &Ncbps, int &Ndbps)
{
	//Ncbps, the number of coded bits in an OFDM symbol (48, 96, 192, or 288 bits).
	//Ncbps, the number of data bits in an OFDM symbol (24, 36, 48, 72, 96, 144, 192, or 216 bits).

	bvec service(16);//the first 7 bits (0\966) are used to synchronize the descrambler. The remaining 9 bits (7\9615) are reserved for future use.
	service.zeros();
	bvec tail(6);//to return the convolutional encoder to the zero state
	tail.zeros();//six scrambled zero bits should be rplaced by six nonscrambled zero bits.

	int Nsys;		//the number of OFDM symbols
	int Ndata;		//the number of bits in the DATA field
	int Npad;		//the number of pad bits
	double length;		//the length of the PSDU in bytes
	length=input_data.length()/8.0;				// to follow the convention used in the IEEE 802.11 document
	Nsys=ceil((16 + 8*length + 6)/(double)Ndbps);		//16 bit of the service field and 6 bits of the tail + length of the data
	Ndata=Nsys*Ndbps;
	Npad = Ndata - (16 + 8*length + 6);
	//cout<<Nsys<<"  "<<Ndata<<"  "<<Npad<<endl;//test ok 
	bvec Pad_Bits(Npad);
	Pad_Bits.zeros();
	bvec data(Ndata);
	data=concat(service,input_data,tail,Pad_Bits);
	return data;	
}
//************************************************* PPDU frame assembler ******************************************************//
cvec PPDU_Frame_Assembler(const cvec &PREAMBLE, const cvec &SIGNAL, const cvec &DATA)
{
	cvec PPDU(PREAMBLE.length());
	PPDU=PREAMBLE;
	PPDU(PPDU.length()-1)+=SIGNAL(0);			//overlapping the signal with the preamble
	PPDU.ins(PPDU.length(),SIGNAL(1,SIGNAL.length()-1));	//Inserting  the signal field
	PPDU(PPDU.length()-1)+=DATA(0);				//overlapping the Data
	PPDU.ins(PPDU.length(),DATA(1,DATA.length()-1));	//Inserting  the Data field
	return PPDU;//PPDU frame 
}
//************************************************* PPDU frame deassembler ******************************************************//
void  PPDU_Frame_Deassembler(const cvec &ofdm_received_symb,  cvec &received_short_squence,  cvec &received_long_squence,  cvec &received_signal_field,  cvec &received_data_field)
{
		 received_short_squence=ofdm_received_symb.left(160);
		 received_long_squence=ofdm_received_symb(160,319);
		 received_signal_field=ofdm_received_symb(320,399);
		 received_data_field=ofdm_received_symb(400,ofdm_received_symb.length()-2);//-2 to remove the extra overlapping symbol	
	return;
}
//******************************************** Data Rate Determination at the receiver side *****************************************************//
ivec Receiver_data_rate(const cvec &received_signal_field,const double &N0,const cvec Hhat)
{
	//this block determines and returns the data rate used at transmission side as well as the number of bytes in the current  PPDU frame
	int Ncbps=48;
	int Ncp=16;
	int Nbpsc=1;
	int index=0;
	int N_SD=48;
	int Nfft=64;
	double R=1/2.0;
	
	cvec ofdm_demodulated_signal=ofdm_demodulation(received_signal_field,Nfft ,Ncp,1);
	//for (int i = 0; i < ofdm_demodulated_signal.length(); i++)  cout<<i<<"  "<<ofdm_demodulated_signal(i)<<endl;
	cvec Equalized_output=Equalization(ofdm_demodulated_signal,Hhat);//equalization 

	cvec ofdm_demodulated_signal_data_carriers=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index);
	//for (int i = 0; i < ofdm_demodulated_signal_data_carriers.length(); i++)  cout<<i<<"  "<<ofdm_demodulated_signal_data_carriers(i)<<endl;
	vec soft_demodulated_symbols=soft_demodulation(ofdm_demodulated_signal_data_carriers, Nbpsc,N0);
	//for (int i = 0; i < soft_demodulated_symbols.length(); i++)  cout<<i<<"  "<<soft_demodulated_symbols(i)<<endl;
	vec received_symbols_deinterleaved=deinterleaver(soft_demodulated_symbols,Ncbps,Nbpsc);
	//for (int i = 0; i < received_symbols_deinterleaved.length(); i++)  cout<<i<<"   "<<received_symbols_deinterleaved(i)<<endl;
	bvec decoded_bits= Punct_Conv_Decode(received_symbols_deinterleaved, R);
	//signal field was not scrambled

	for (int i = 0; i < decoded_bits.length(); i++)  
		decoded_bits(i) = (decoded_bits(i) == 0 ? 1 : 0);	//flip all bits (that is, convert zeros to ones, and ones to zeros) // for BPSK
	//cout<<decoded_bits<<endl;
	
	//check parity bit
	bvec Parity_bit(1);
	Parity_bit=decoded_bits(17);	//received parity bit
	bvec RATE=decoded_bits.left(4);	//Rate
	bvec LENGTH=decoded_bits(5,16);	//Length

	//transform the length (from MSB to LSB) to inverse the order as at the sender side the standard 802.11 impose that length is represented (from LSB to MSB)
	bvec temp(LENGTH.length());
	for( int i=0; i<LENGTH.length();i++)
		temp(LENGTH.length()-1-i)=LENGTH(i);
	LENGTH=temp;

	int leng=bin2dec(LENGTH);//convert the length of bytes of PSDU into the decimal representation
	//cout<<leng<<endl;

	//Parity bit : even parity 
	bvec tested_parity(1);
	if ((sum(to_ivec(decoded_bits(0,16)))%2==0))
		tested_parity="0";
	else
		tested_parity="1";

	if (Parity_bit!=tested_parity)
		cout<< "ERROR IN THE PARITY BIT OF THE RECEIVED SIGNAL FIELD::Wrong or not valid data rate will be chosen!!!!"<<endl;

	//Receiver data rate determination
	int data_rate;

	if (RATE=="1 1 0 1") data_rate=6;
	else if (RATE=="1 1 1 1") data_rate=9;
		else if (RATE=="0 1 0 1") data_rate=12;
			else if (RATE=="0 1 1 1") data_rate=18;
				else if (RATE=="1 0 0 1") data_rate=24;
					else if (RATE=="1 0 1 1") data_rate=36;
						else if (RATE=="0 0 0 1") data_rate=48;
							else if (RATE=="0 0 1 1") data_rate=54;
								else cout<< "ERROR IN THE RECEIVED SIGNAL FIELD::not valid data rate !!!!"<<endl;
	ivec output(2);
	output(0)=data_rate;
	output(1)=leng;
	return output;
}
//************************************************* OFDM Data Symbols Assembler ******************************************************//
cvec OFDM_Data_Symbols_Assembler(const cvec &data_ofdm_transmited_symb)
{
	cvec data_ofdm_transmited_symb_final(81);		// size = 81 (Nfft+CP+1) 1 for the overlapping

	for(int i = 0;  i < data_ofdm_transmited_symb.length()/80; i++)
        {		
		// cout<<i<<" out of :"<<data_ofdm_transmited_symb.length()/80<<endl;
		if (data_ofdm_transmited_symb.length()%80!=0) cout<<"Error size exced 80 octets";	//normally there is no error at this level
		cvec block(80);
		block=data_ofdm_transmited_symb(i*80,(i+1)*80-1);
		block.ins(block.length(), block(16));		
		//Need to add the first element of the actual OFDM time sample (not the prefix)
		//at the end. This will be added as an overlap with the next time sample
                block= window(block);						// windowing function for block
				
                 //Consider the block overlap
                 if (i == 0) //First block 
				    data_ofdm_transmited_symb_final=block;
                 else {
                         //overlapping the ofdm symbol
			data_ofdm_transmited_symb_final(data_ofdm_transmited_symb_final.length() - 1) += block(0);
                        data_ofdm_transmited_symb_final.ins(data_ofdm_transmited_symb_final.length(), block(1, block.length() - 1));
						
                       }
		} 
	//cout<< data_ofdm_transmited_symb_final;cin>>x;
	return data_ofdm_transmited_symb_final;
}


/**************************************************   text example  ******************************************************************/
bvec text_example(void)//this function is used to generate the bits of the message that is used in the example in annex L
{	
	string text= "Joy, bright spark of divinity,\nDaughter of Elysium,\nFire-insired we trea";
	int text_i;
	bvec text_bin_t;	//temporary binary vector
	bvec text_bin; 		// binary representation of the text
	
	for (int i = 0; i < text.length(); i++) {
		text_i = text[i];
		text_bin_t = dec2bin (8, text_i);
		
		//transform the length (from LSB to MSB)
		bvec temp(text_bin_t.length());
		for( int i=0; i<text_bin_t.length();i++)
			temp(text_bin_t.length()-1-i)=text_bin_t(i);
		text_bin_t=temp;
		text_bin.ins(text_bin.length(),text_bin_t);//
		//cout<<text_bin<<endl;//ok tested
	}

	//MAC HEADER
	// MAC Header length() = 24 
	int MAC_header_int[] = { 4, 2, 0, 46, 0, 96, 8, 205, 55, 166, 0, 32, 214, 1,60, 241, 0, 96, 8, 173, 59, 175, 0, 0 };  
	bvec MAC_header_bin_t;		//tempo
	bvec MAC_header_bin;
	
	// transform MAC header to binary
	for (int i = 0; i < 24; i++) {
		MAC_header_bin_t = dec2bin (8, MAC_header_int[i]);
		//transform the length (from LSB to MSB)
		bvec temp(MAC_header_bin_t.length());
		for( int i=0; i<MAC_header_bin_t.length();i++)
			temp(MAC_header_bin_t.length()-1-i)=MAC_header_bin_t(i);
		//cout<<i<<" "<<MAC_header_bin_t<<"  "<<temp<<endl;//ok tested
		MAC_header_bin_t=temp;
		MAC_header_bin.ins(MAC_header_bin.length(),MAC_header_bin_t);//
	}

	//CRC 32 --> CRC length() = 4
	int CRC[] = { 103, 51, 33, 182 };//as per the example 
	bvec CRC_bin_t;//tempo
	bvec CRC_bin;

	for (int i = 0; i < 4; i++) {
		CRC_bin_t = dec2bin (8, CRC[i]);
		//transform the length (from LSB to MSB)
		bvec temp(CRC_bin_t.length());
		for( int i=0; i<CRC_bin_t.length();i++)
			temp(CRC_bin_t.length()-1-i)=CRC_bin_t(i);
		//cout<<i<<" "<<CRC_bin_t<<"  "<<temp<<endl;//ok tested
		CRC_bin_t=temp;
		CRC_bin.ins(CRC_bin.length(),CRC_bin_t);//
	}

	return concat(MAC_header_bin,text_bin,CRC_bin); // Frame format: | MAC Header | Data | CRC |

}
//****************************************************** Pilot_insertion and ofdm Symbol Generation  *************************************************//

cvec Pilot_insertion( const  cvec &input, const int &N_SD ,const int &Nfft,int &index)
{
	int nb_row;
	nb_row = input.length()/N_SD;
	cmat out_put(nb_row , Nfft);
	out_put.zeros(); 
	//index is the starting index value
	
	int polarity;
	//the four pilot subcarriers are added by taking the values the sequence p{1.0,1.0,1.0,-1.0}, multiplying them by the an element of sequence p0\85126,

	// insertion of pilot subcarriers and data 
	for(int i = 0;  i < nb_row; i++) 
		{	polarity=Pilot_Polarity_Sequence(index);	//cout<<index<<endl;
			int k = 0; 
			for (int j = 0 ; j < Nfft ; j++)
				{
					if( j == 11 || j == 25 || j == 39 || j == 53) //Pilot subcarriers  -21  -7  7  21
						if(j== 53)
							 out_put(i,j) = std::complex<double>(-1*polarity,(double)0); 		//conversion into complex
						else out_put(i,j) = std::complex<double>(1*polarity,(double)0);
					//Data subcarriers
					else if( j != 32 &&  j > 5 && j < 59 )
							{
								out_put(i,j) = input( i*N_SD + k); 			
								k++;
							}//end else if
				}//end for
			index=index+1;

		}//end for
	 
	// 
	cvec output_s( out_put.rows()*out_put.cols() ); 
	for(int i = 0; i < output_s.length(); i++)
		output_s(i) =  out_put(floor((long double)i/Nfft) , mod(i,Nfft)); 
	return output_s; 	
}
//************************************ Pilot Data separation ***************************************/
cvec Pilot_Data_separation( const  cvec &input, const int &N_SD ,const int &Nfft,int &index)
{
	cmat matrix(input.length()/Nfft , Nfft);
	matrix.zeros(); 

	for(int i = 0; i < input.length() ; i++)
		matrix(floor((float)i/Nfft), mod(i,Nfft)) = input(i); // transform into matrix form
	
	int nb_data_subcar = matrix.rows() * N_SD; // N_SD = Number of data subcarriers per ofdm symbol
	cvec data_subcar(nb_data_subcar);// data subcarriers

	int nb_pilot_subcar = matrix.rows() * 4; // total Number of pilot subcarriers = 4 * nbr of symbols
	cvec pilot_subcar(nb_pilot_subcar);// data subcarriers
	
	// Data and pilot recovery 
	for(int i = 0;  i < matrix.rows(); i++) 
	{
		int k = 0; //Data subcariers
		for (int j = 0 ; j < Nfft ; j++)
			if( j > 5 && j < 59 &&  j != 32 && j != 11 && j != 25 && j != 39 && j != 53 )
			{
				data_subcar( (i * N_SD) + k) = matrix( i, j );		
				k++;
			}

		int l = 0;  //Pilot subcariers
		for (int j = 0 ; j < Nfft ; j++)
			if( j == 11 || j == 25 || j == 39 || j == 53) 
			{
				pilot_subcar( (i * 4) + l) = matrix( i, j );
				l++;
			}
	}
	//cout << real(pilot_subcar);tested ok 1 1 1 -1

	return data_subcar; 	
}

//************************************** Pilot_Polarity_Sequence *****************************/
int Pilot_Polarity_Sequence(const int index)
{
	ivec Polarity_Sequence = "1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1,1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1,-1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1";
	return Polarity_Sequence(index % 126) ; 
}

//********************************* ofdm_modulation ******************************/
cvec ofdm_modulation(const  cvec &input, const int &Nfft ,const int &Ncp, const int & norm_factor_flag)
{
	double norm_factor;	//Normalizing factor

	if (norm_factor_flag==1)
		//norm_factor=sqrt(double(Nfft*Nfft)/double(Nfft+Ncp));
		norm_factor=sqrt(double(Nfft));
	else
		norm_factor=1;

	int nb_row;
	nb_row = input.length()/Nfft;
	cmat out_put(nb_row , Nfft+Ncp);
	out_put.zeros(); 
	cvec temp_t(Nfft);
	cvec temp_f(Nfft);
	
	for(int i = 0;  i < nb_row; i++) 
		{	
			temp_f=input(i*Nfft,(1+i)*Nfft-1);
			//temp_t=ifft (temp_f, Nfft);	//ifft
			temp_t=ifft (concat(temp_f.right(32),temp_f.left(32)), Nfft);	//the method used by this standard
			out_put.set_row(i,concat(temp_t.right(Ncp),temp_t));			//adding CP
		}//end for
	 
	// parallel to series conversion
	cvec ofdm_symbols( out_put.rows()*out_put.cols() ); 

	for(int i = 0; i < ofdm_symbols.length(); i++)
		ofdm_symbols(i) =  out_put(floor((long double)i/(Nfft+Ncp)) , mod(i,Nfft+Ncp)); //cout<<floor((long double)i/(Nfft+Ncp))<<"   "<<mod(i,Nfft+Ncp)<<endl;
	return norm_factor*ofdm_symbols; //tested ok
}
//************************************** ofdm_demodulation ***************************************/
cvec ofdm_demodulation(const  cvec &input, const int &Nfft ,const int &Ncp,const int & norm_factor_flag)

{	
	double norm_factor;//Normalizing factor

	if (norm_factor_flag==1)
		//norm_factor=sqrt(double(Nfft*Nfft)/double(Nfft+Ncp));
		norm_factor=sqrt(double(Nfft));
	else
		norm_factor=1;

	int nb_row;
	nb_row = input.length()/(Nfft+Ncp);
	cmat out_put(nb_row , Nfft);
	out_put.zeros(); 
	cvec temp_t1(Nfft+Ncp);
	cvec temp_t2(Nfft);
	cvec temp_f(Nfft);
	
	for(int i = 0;  i < nb_row; i++)
		{	temp_t1=input(i*(Nfft+Ncp),(1+i)*(Nfft+Ncp)-1);//cout<<i*(Nfft+Ncp)<<"  "<<(1+i)*(Nfft+Ncp)-1;
			temp_t2=temp_t1.right(Nfft);				  //removing CP
			temp_f=fft(temp_t2, Nfft);					 //fft
			out_put.set_row(i,concat(temp_f.right(32),temp_f.left(32)));//to compensate the inversion at the ifft stage		
		}//end for

	// parallel to series conversion
	cvec demodulated_ofdm_symbols( out_put.rows()*out_put.cols() ); 
	for(int i = 0; i < demodulated_ofdm_symbols.length(); i++)
		demodulated_ofdm_symbols(i) =  out_put(floor((long double)i/Nfft) , mod(i,Nfft)); 
	
	return demodulated_ofdm_symbols/norm_factor; 	//tested ok
}

//************************************************* Channel Estimation and Equalization **************************************************************//
	
//************************************************* LS estimator **************************************************************//
cvec LS_estimator(const cvec &rec_long_squence_CP)
{
	//rec_long_squence_CP is the received long training sequence with the CP
	int Nfft=64;
	double norm_factor=sqrt(double(Nfft));								// ofdm normalizing factor
	cvec rec_long_squence(128);
	rec_long_squence=rec_long_squence_CP(32,159);				// without the CP
	//cout << "rec_long_squence without CP = " << rec_long_squence.length()<<endl;// tested ok
	
	cvec rec_long_squence1=rec_long_squence.left(64);		    // 1st received symbol of the long training sequence without the CP
	cvec rec_long_squence2=rec_long_squence.right(64);			// 2nd received symbol of the long training sequence without the CP
	cvec dem_long_squence_temp;									// temporary demodulated long training sequence symbol		
	cvec dem_long_squence1,dem_long_squence2;					// received 1st and 2nd long sequences in the frequency domain

	// demodulated long training sequence symbol #1
	dem_long_squence_temp=fft(rec_long_squence1, Nfft);										
	dem_long_squence1=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence1=dem_long_squence1/norm_factor;
	//cout << "dem_long_squence1 = " << dem_long_squence1.left(10) <<" length "<<dem_long_squence1.length()<<endl;// tested ok
	
	// demodulated long training sequence symbol #2
	dem_long_squence_temp=fft(rec_long_squence2, Nfft);										
	dem_long_squence2=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence2=dem_long_squence2/norm_factor;
	//cout << "dem_long_squence2 = " << dem_long_squence2.left(10)<<" length "<<dem_long_squence2.length()<<endl;// tested ok
	
	// Long Training Sequence
	cvec Long_training_seq = " 0 0 0 0 0 0  1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0  "; 

	cvec Hh_LS1,Hh_LS2,Hh_LS;
	//cout <<dem_long_squence1.length()<<" "<<Long_training_seq.length();
	Hh_LS1=elem_div(dem_long_squence1,Long_training_seq);
	Hh_LS2=elem_div(dem_long_squence2,Long_training_seq);
	Hh_LS=(Hh_LS1+Hh_LS2)*0.5;
	Hh_LS.set_subvector(0,"0 0 0 0 0 0");	//cout<<Hh_LS(0,5)<<endl<<endl;
	Hh_LS(32)=0;							//cout<<Hh_LS(32)<<endl<<endl;
	Hh_LS.set_subvector(59,"0 0 0 0 0");

	return Hh_LS;
}
//************************************************* MMSE estimator **************************************************************//
cvec MMSE_estimator(const cvec &rec_long_squence_CP)
{
	//rec_long_squence_CP is the received long training sequence with the CP
	int Nfft=64;
	double norm_factor=sqrt(double(Nfft));						// ofdm normalizing factor
	cvec rec_long_squence(128);									
	rec_long_squence=rec_long_squence_CP(32,159);				// without the CP
	//cout << "rec_long_squence without CP = " << rec_long_squence.length()<<endl;// tested ok
	
	cvec rec_long_squence1=rec_long_squence.left(64);		    // 1st received symbol of the long training sequence without the CP
	cvec rec_long_squence2=rec_long_squence.right(64);			// 2nd received symbol of the long training sequence without the CP
	cvec dem_long_squence_temp;									// temporary demodulated long training sequence symbol		
	cvec dem_long_squence1,dem_long_squence2;					// received 1st and 2nd long sequences in the frequency domain

	// demodulated long training sequence symbol #1
	dem_long_squence_temp=fft(rec_long_squence1, Nfft);										
	dem_long_squence1=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence1=dem_long_squence1/norm_factor;
	//cout << "dem_long_squence1 = " << dem_long_squence1.left(10) <<" length "<<dem_long_squence1.length()<<endl;// tested ok
	
	// demodulated long training sequence symbol #2
	dem_long_squence_temp=fft(rec_long_squence2, Nfft);										
	dem_long_squence2=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence2=dem_long_squence2/norm_factor;
	//cout << "dem_long_squence2 = " << dem_long_squence2.left(10)<<" length "<<dem_long_squence2.length()<<endl;// tested ok
	
	//noise variance calculation
	double noise_var = 0;
	double temp;

	for (int i=0; i<Nfft; i++)
		{
			temp=(conj(dem_long_squence1(i)-dem_long_squence2(i))*(dem_long_squence1(i)-dem_long_squence2(i))).real();
			noise_var = noise_var +temp;
		}
	noise_var=noise_var/(Nfft*2);
	
	// Long Training Sequence (frequency domain)
	cvec Long_training_seq = " 0 0 0 0 0 0  1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0  "; 
	
	cmat X=diag(Long_training_seq);
	cvec Y=(dem_long_squence1+dem_long_squence2)/2;
	cmat F=Fmatrix(Nfft, Nfft);

	cvec Hh_LS=LS_estimator(rec_long_squence_CP);
	cvec cor=xcorr(Hh_LS,Hh_LS);//for(int i=0; i<cor.size(); i++) cout<<i+1<<" "<<cor(i)<<endl;
	cmat Rhh=toeplitz(cor(0,63));

	cmat Rhy = Rhh * hermitian_transpose(F) * hermitian_transpose(X);
	cmat Ryy = X * F * Rhh * hermitian_transpose(F) * hermitian_transpose(X) + noise_var * eye_c(64);
	cvec Hh_MMSE = F * Rhy * inv(Ryy) * Y;
	//for (int i=0; i<Hh_MMSE.length(); i++ ) cout <<i+1<<" "<<Hh_MMSE(i)<<endl;
	
	Hh_MMSE.set_subvector(0,"0 0 0 0 0 0");	//
	Hh_MMSE(32)=0;							//
	Hh_MMSE.set_subvector(59,"0 0 0 0 0");

	return Hh_MMSE;
}

cmat Fmatrix(const int N1,const int N2)
	{
	//N1 is the number of subcarriers
	//N2 is the number of number of columns
	mat k1k2(N1,N1);		// matrix contains the product of k1*k2
	cmat F;					// the DFT-matrix

	complex<double> j(0, 1.0);
	
	for (int k1=0; k1<N1; k1++)
		for (int k2=0; k2<N1; k2++)
			k1k2(k1,k2)=k2*k1;

	k1k2 = k1k2(0,N1-1,0,N2-1);
	F= exp(- j*2.0 * pi * k1k2/ N1) / sqrt(1.0 * N1);
	return F;
	}
//***************************************** noise_variance ************************************************************//
double noise_variance(const cvec &rec_long_squence_CP)
{
	//rec_long_squence_CP is the received long training sequence with the CP
	int Nfft=64;
	double norm_factor=sqrt(double(Nfft));				// ofdm normalizing factor
	cvec rec_long_squence(128);									
	rec_long_squence=rec_long_squence_CP(32,159);			// without the CP
	
	cvec rec_long_squence1=rec_long_squence.left(64);		// 1st received symbol of the long training sequence without the CP
	cvec rec_long_squence2=rec_long_squence.right(64);		// 2nd received symbol of the long training sequence without the CP
	cvec dem_long_squence_temp;					// temporary demodulated long training sequence symbol		
	cvec dem_long_squence1,dem_long_squence2;			// received 1st and 2nd long sequences in the frequency domain

	// demodulated long training sequence symbol #1
	dem_long_squence_temp=fft(rec_long_squence1, Nfft);										
	dem_long_squence1=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence1=dem_long_squence1/norm_factor;
	
	// demodulated long training sequence symbol #2
	dem_long_squence_temp=fft(rec_long_squence2, Nfft);										
	dem_long_squence2=concat(dem_long_squence_temp.right(32),dem_long_squence_temp.left(32));		//to compensate the inversion at the ifft stage	
	dem_long_squence2=dem_long_squence2/norm_factor;
	
	//noise variance calculation
	double noise_var = 0;
	double temp;

	for (int i=0; i<Nfft; i++)
		{
			temp=(conj(dem_long_squence1(i)-dem_long_squence2(i))*(dem_long_squence1(i)-dem_long_squence2(i))).real();
			noise_var = noise_var +temp;
		}
	noise_var=noise_var/(Nfft*2);
	return noise_var;
}
	
//*********************************************************   Equalization  ***********************************************************************//
cvec Equalization(const cvec &input,const cvec &Hhat)
{
	int nb_row;
	int Nfft=64;
	nb_row = input.length()/Nfft;
	cmat out_put(nb_row , Nfft);
	out_put.zeros(); 
	cvec temp(Nfft);
	
	for(int i = 0;  i < nb_row; i++) 
		{	
			temp=elem_div(input(i*Nfft,(1+i)*Nfft-1), Hhat);	//Applying the equalization
			out_put.set_row(i,temp);
		}//end for
	 
	// parallel to series conversion
	cvec equalized_output( out_put.rows()*out_put.cols() ); 
	for(int i = 0; i < equalized_output.length(); i++)
		equalized_output(i) =  out_put(floor((long double)i/Nfft) , mod(i,Nfft));
	return equalized_output;
	
}
//*********************************************     file to binary vector conversion    ************************************************//
bvec file2binary_converter(void)
{
	cout<<endl;
	cout<<"Please enter the file name you want to transmit with its extension eg image.jpg"<<endl;
	cout<<"Please note that this file must be in the same directory"<<endl;
	char file_name[20];
	cin>>file_name;

	unsigned char Buf;
	int i, Databin;
	char buf[20];
	sprintf(buf,"%s",file_name);
	FILE* fileIn=fopen(buf,"rb");
	FILE* fileOut=fopen("file.txt","wb");
	
	//ifstream fileIn(buf);
   	//ofstream fileOut("out.bin", std::ios::binary);

	double max_file_size = 1e10; 		//maximum number of bits that can be stored
	bvec binary_data(max_file_size);

	int k=0;
	cout<<"Now, file is being converted into a binary vector.... ";
	while (!feof(fileIn))
	{
	
		size_t result = fread(&Buf,sizeof(char),1,fileIn);
		// Test en cas de non lecture du fichier
		//if (result == 0) {  cout << " File reading error " << endl; }  //fputs ("Reading error",stderr); exit (3);}
	
		for(i=0; i<8; i++)

			{
				Databin=((int)Buf>>i)&1;
				//printf("%d",Databin);
				binary_data(k)=Databin;
				fprintf(fileOut,"%d\n",Databin); // to save the bits in a text file
				k++;
			}
		}
	binary_data.set_size(k,true);
	cout<<"Done"<<endl;
	fclose(fileIn);
	fclose(fileOut);
	return binary_data;
}
//*********************************************     binary vector to file conversion ************************************************//
void binary2file_converter(const bvec &binary_data)
{
	cout<<endl;
	cout<<"Please enter a file name of the received file with its extension eg. image.jpg"<<endl;
	cout<<"this file will be saved in the same directory"<<endl;

	char file_name[20];
	cin>>file_name;
	char buf[20];

	sprintf(buf,"%s",file_name);
	FILE* fileOut=fopen(buf,"wb");

	int bit;
	int cpt=0,tmp=0;  
	char c;

	cout<<"Now, file is being converted into its orignal format .... "<<endl;
	for(int i=0; i<binary_data.size(); i++)
		{
			bit=binary_data[i];
			tmp += bit*pow2(cpt);
			cpt++;  
		
			if(cpt == 8)
			{ 
				 fprintf(fileOut,"%c",((char)tmp));
				  tmp = 0;
				  cpt = 0;
			}
		}

	fclose(fileOut);
	cout<<"Done"<<endl;
}
//***************************************************************************************************************************************//
//****************************************************** Saving Results to files *******************************************************************//
void BER_SNR(char* file_name, vec EbN0dB, vec bit_error_rate, int mode)
{
	FILE *f;
	int i;
	char c[32]; //becuase double is 8 bytes in GCC compiler but take 10 for safety 
	int leng;
	char buf[32];

	if (mode == 0) //Save the results to it file:
		 {
		   cout <<"saving data to  " <<file_name<<".it and " <<file_name<<".txt"<<endl;//
		   sprintf(buf, "%s.it",file_name);
		   cout << endl;

		  it_file ff;
		  ff.open(buf);
		  ff << Name("EbN0dB") << EbN0dB;
		  ff << Name("ber") << bit_error_rate;
		  ff.close();
	}
		else  //Save the results to text file:
			{
			sprintf(buf, "%s.txt",file_name);
			f = fopen (buf, "w+");

			if (f == NULL)
    				cout << "Unable to open the file !" << endl;
  			else
   			{	fwrite ("SNR\t\tBER\n",sizeof(char),9,f);
				for (i=0; i<EbN0dB.length(); i++)
				{
					sprintf(c , "%f" , EbN0dB(i));
					leng = chain_length(c);
					fwrite (c,sizeof(char),leng,f);
					fwrite ("\t",sizeof(char),1,f);			
					sprintf(c , "%f" ,bit_error_rate(i));
					leng = chain_length(c);
					fwrite (c,sizeof(char),leng,f);
					fwrite ("\n",sizeof(char),1,f);	
				}		
				fclose (f);	
			}
	}
}
int chain_length(char c[32])
{
	int i,l=0;	
	
	for(i=0;c[i]!='\0';i++)
		l++;
	
	return l;
}

//*****************************   Sent_Frame_Saver ******************************************//
void Sent_Frame_Saver(cvec PPDU_Frame_Final)
{
	  char file_name[11]="Sent_Frame";
	  //file_name
	  char buf[14];
	  cout << endl;
   
	   
	   sprintf(buf, "%s.it",file_name);
	   cout << endl;

		it_file ff;
		ff.open(buf);
		ff << Name("PPDU_Frame_Final") << PPDU_Frame_Final;
		ff.close();
		cout <<"saving data to  " <<file_name<<".it" << " ... Done"<< endl;//
}
//************************************************************************************************************************************************/

//**************************************Simulation Parameters*********************************************/
void header_display(void)
{
	cout << " ****************************************************************************\n";	
	cout << " ************************** 802.11 OFDM Physical Layer **********************\n";
	cout << " ****************************************************************************\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " Modulation\tCoding\tCoded bits/ \tCoded bits/ \tData bits/ \tData\n";
	cout << " \t\t rate\tsubcarrier\tOFDM symbol\tOFDM symbol\trate\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " BPSK \t\t 1/2 \t     1 \t\t 48 \t\t 24\t\t 6\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " BPSK \t\t 3/4 \t     1 \t\t 48 \t\t 36\t\t 9\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " QPSK \t\t 1/2 \t     2 \t\t 96 \t\t 48\t\t 12\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " QPSK \t\t 3/4 \t     2 \t\t 96 \t\t 72\t\t 18\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " 16-QAM\t\t 1/2 \t     4 \t\t 192 \t\t 96\t\t 24\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " 16-QAM\t\t 3/4 \t     4 \t\t 192 \t\t 144\t\t 36\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " 64-QAM\t\t 2/3 \t     6 \t\t 288 \t\t 192\t\t 48\n";
	cout << " ----------------------------------------------------------------------------\n";
	cout << " 64-QAM\t\t 3/4 \t     6 \t\t 288 \t\t 216\t\t 54\n";
	cout << " ----------------------------------------------------------------------------\n";

}
void data_rate_choice(double* R, int* Nbpsc, int* data_rate,int *Ncbps,int *Ndbps)
{	//function to determine the data rate and set all the other parameters accordingly
	
	if (*data_rate==0)	//if data_rate==0 the data rate value will be taken from the user
		{
		cout << "Setting Simulation Parameters :\n";	
		cout << "Data rate(Mb/s):\n";
		cout << "\t- Enter one of these values (6, 9, 12, 18, 24, 36 ,48,54 )(Mb/s): ";
		cin  >> *data_rate;
		cout << endl;
		cout << endl;
		}
		
	switch(*data_rate)
	{
		case 6 :
				*R=1/2.0;
				*Nbpsc=1;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				/*cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : BPSK\n";
				cout << "\t\t- Coding rate : 1/2\n";
				cout << "\t\t- Coded bits per subcarrier : "<< *Nbpsc << endl;
				cout << "\t\t- Coded bits per OFDM symbol :"<< *Ncbps << endl;
				cout << "\t\t- Coded data bits per OFDM symbol :"<< *Ndbps << endl;*/
		break;
		case 9 :
				*R=3/4.0;
				*Nbpsc=1;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : BPSK\n";
				cout << "\t\t- Coding rate : 3/4\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		case 12 :
				*R=1/2.0;
				*Nbpsc=2;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : QPSK\n";
				cout << "\t\t- Coding rate : 1/2\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		case 18 :
				*R=3/4.0;
				*Nbpsc=2;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : QPSK\n";
				cout << "\t\t- Coding rate : 3/4\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";

		break;
		case 24 :
				*R=1/2.0;
				*Nbpsc=4;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : 16-QAM\n";
				cout << "\t\t- Coding rate : 1/2\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		case 36 :
				*R=3/4.0;
				*Nbpsc=4;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : 16-QAM\n";
				cout << "\t\t- Coding rate : 3/4\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		case 48 :
				*R=2/3.0;
				*Nbpsc=6;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : 64-QAM\n";
				cout << "\t\t- Coding rate : 2/3\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		case 54 :
				*R=3/4.0;
				*Nbpsc=6;
				*Ncbps=48**Nbpsc;
				*Ndbps=(int)*Ncbps**R;
				cout << "Simulation parameters :\n\n";
				cout << "\t\t- Modulation : 64-QAM\n";
				cout << "\t\t- Coding rate : 3/4\n";
				cout << "\t\t- Coded bits per subcarrier :"<<*Nbpsc<<"\n";
				cout << "\t\t- Coded bits per OFDM symbol :"<<*Ncbps<<"\n";	
				cout << "\t\t- Coded data bits per OFDM symbol :"<<*Ndbps<<"\n";
		break;
		default:{
			     cout << "Not a valid choice"<<endl;
			     cout<< "!!! data rate must be either =0 to take the value from the user as an input \
						  value or one of these values (6, 9, 12, 18, 24, 36 ,48,54 )(Mb/s):";
				}
		break;
	}
}
//*********************************************** number of frames **********************************************************/
cvec frames_nbr_field(int &len)
{
	// len is number of the frames 
	
	bvec leng(24);
	leng.zeros();
	if (len>pow2(24))
		cout<<"Number of frames is extremely large !!! "<<endl;	//normally number of frames is much less than this value 

	leng=dec2bin (24, len);//convert the length of needed bytes of PSDU into the binary representation

	//transform the length (from LSB to MSB) to follow the order imposed by the standard
	bvec temp(leng.length());
	for( int i=0; i<leng.length();i++)
		temp(leng.length()-1-i)=leng(i);
	leng=temp;

	//frame number field is not scrambled

	//Coding the frame number field
	bvec coded_frame_nbr_bits;
	double R =1/2.0;
	coded_frame_nbr_bits= Punct_Conv_Code(leng,R);		//The bits are encoded by the rate 1/2 convolutional encoder		

	//interleaving
	bvec interleaved_frame_nbr_bits;
	int Ncbps=48;
	int Nbpsc=1;
	interleaved_frame_nbr_bits=interleaver(coded_frame_nbr_bits, Ncbps,Nbpsc);
	
	//Modulation BPSK
	cvec modulated_frame_nbr_bits;
	ivec deci;cvec constl;
	modulated_frame_nbr_bits= modulation(interleaved_frame_nbr_bits, Nbpsc,deci,constl);
	
	//Pilot insertion and OFDM modulation
	int N_SD=48;
	int Nfft=64;
	int Ncp=16;
	int index=0;//to choose first element of the polarity sequence (like the signal field)
	cvec frame_nbr_ofdm_freq;
    	frame_nbr_ofdm_freq= Pilot_insertion( modulated_frame_nbr_bits,N_SD ,Nfft,index);

	//OFDM Modulation	
	cvec frame_nbr_ofdm_modulated=ofdm_modulation(frame_nbr_ofdm_freq,Nfft ,Ncp,1);			//for real transmission flag = 1
	frame_nbr_ofdm_modulated.ins(frame_nbr_ofdm_modulated.length(),frame_nbr_ofdm_modulated(16));	// add the element 81 with the first symbol after the CP
	frame_nbr_ofdm_modulated= window(frame_nbr_ofdm_modulated);					//window function
	
	return frame_nbr_ofdm_modulated(0,79);
}

//***********************************   determination of number of frames ****************************************//
int frame_nbr_receiver(const cvec &received_frame_nbr_field,const double &N0,const cvec Hhat)
{
	// this block determines and returns the number of frames that was sent at transmission phase 

	int Ncbps=48;
	int Ncp=16;
	int Nbpsc=1;
	int index=0;
	int N_SD=48;
	int Nfft=64;
	double R=1/2.0;
	// OFDM demodulation
	cvec ofdm_demodulated_frame_nbr=ofdm_demodulation(received_frame_nbr_field,Nfft ,Ncp,1);
	//equalization
	cvec Equalized_output=Equalization(ofdm_demodulated_frame_nbr,Hhat);	 
	// Pilote Data separation
	cvec ofdm_demodulated_frame_nbr_data_carriers=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index);
	// soft demodulation
	vec soft_demodulated_symbols=soft_demodulation(ofdm_demodulated_frame_nbr_data_carriers, Nbpsc,N0);
	//deinterleaver
	vec received_symbols_deinterleaved=deinterleaver(soft_demodulated_symbols,Ncbps,Nbpsc);
	//decoder
	bvec decoded_bits= Punct_Conv_Decode(received_symbols_deinterleaved, R);

	//frame number field was not scrambled

	//flip all bits (that is, convert zeros to ones, and ones to zeros) // for BPSK
	for (int i = 0; i < decoded_bits.length(); i++)  
		decoded_bits(i) = (decoded_bits(i) == 0 ? 1 : 0);	
	
	bvec LENGTH=decoded_bits;	//Length

	//transform the length (from MSB to LSB) to inverse the order
	bvec temp(LENGTH.length());
	for( int i=0; i<LENGTH.length();i++)
		temp(LENGTH.length()-1-i)=LENGTH(i);
	LENGTH=temp;

	//convert the number of frames  into the decimal representation
	int leng=bin2dec(LENGTH);	
	
	return leng;
}
