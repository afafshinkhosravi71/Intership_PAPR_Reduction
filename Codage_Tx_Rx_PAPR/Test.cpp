/************************************************************************/  
// IEEE 802.11 OFDM PHYSICAL LAYER
// Author : Boubacar DIALLO 
// Date   : May 2017
/************************************************************************/

#include "itpp/itcomm.h"
#include "itpp/itstat.h"
#include "Utils.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <iterator>

// g++ `pkg-config --cflags itpp` -o Test Utils.h Test.cpp `pkg-config --libs itpp`

using namespace std;
using namespace itpp;
using std::cin;
using std::cout;
using std::endl;

int main()
{
	vec EbN0_dB = "0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0";
    	double R = 1.0/2.0;    					//coding rate (non punctured PCCC)
    	double Ec = 1.0;      					//coded bit energy
	cout << "Now simulating EbN0db = " << EbN0_dB << endl;

	vec sigma2 = (0.5*Ec/R)*pow(inv_dB(EbN0_dB), -1.0); 	// N0/2
	cout << "Noise variance Sigma = " << sigma2 << endl;

	return 0;
}

/*
int main()
{
	
	string text= "Joy, bright spark of divinity,\nDaughter of Elysium,\nFire-insired";
	int text_i;
	bvec text_bin_t;	//temporary binary vector
	bvec text_bin; 		// binary representation of the text
	
	for (int i = 0; i < text.length(); i++) {
		text_i = text[i];
		text_bin_t = dec2bin (8, text_i);
		//cout<<"text_bin_t 0" <<text_bin_t<<endl;//ok tested
		cout<<"text_i = " << text[i] <<" == " << static_cast<int>(text[i]) <<endl;//ok tested
		//transform the length (from LSB to MSB)
		bvec temp(text_bin_t.length());
		for( int i=0; i<text_bin_t.length();i++){
			temp(text_bin_t.length()-1-i)=text_bin_t(i);
			//cout << "temp" << temp <<endl;
			}
		text_bin_t=temp;
		text_bin.ins(text_bin.length(),text_bin_t);//

		
		//cout<<"text_bin_t" <<text_bin_t<<endl;//ok tested
		//cout<<"text_bin" <<text_bin<<endl;//ok tested
		//cout<<"text_bin" <<text_bin.length()<<endl;//ok tested
		//cout<<"text size" <<text.length()<<endl;//ok tested
	}
	return 0;
} 
*/
/*
cvec window(const cvec &input)//windowing function
	{	
		cvec output(size(input));
		output=input;
		output(0)=.5*input(0);
		output(input.length()-1)=.5*input(input.length()-1);
		return output;
	}

int main()
{
	// Generation of the short sequences
	cvec Short_training_seq(53);       //including the value 0 at dc
	Short_training_seq.zeros();

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

	cout << "Short_training 1" << Short_training_seq << endl;

	//Add the zeros to have 64-symbols 
	cvec a(6);a.zeros();
	cvec b(5);b.zeros();

	Short_training_seq=sqrt(13/6.0)*concat(a,Short_training_seq,b);
	cout << "Short_training 2" << Short_training_seq << endl;

	cvec Short_training_seq_time;
	cvec Short_training_seq_time_extended;

	//the coefficients 1 to 26 are mapped to the same numbered IFFT inputs, while the coefficients \9626 to \961 are copied into IFFT inputs 38 to 63
	Short_training_seq_time=ifft(concat(Short_training_seq.right(32),Short_training_seq.left(32)),64);
	cout << "Short_training after IFFT" << Short_training_seq_time << endl;

	//the short training sequence is extended periodically for 161 samples
	Short_training_seq_time_extended=concat(Short_training_seq_time,Short_training_seq_time,Short_training_seq_time.left(33));
	cout << "Short_training 4 extended 161 samples" << Short_training_seq_time_extended << endl;

	Short_training_seq_time_extended=window(Short_training_seq_time_extended);	//window function

	//Random test : tested with Table L-4\97
	cout << "Time domain representation of the short sequence " << endl;
	for (int i=0; i<10;i++)
		cout << i*11+1 << " --> " << Short_training_seq_time_extended(i*11+1) << endl;

	
	int index = 4;
	ivec Polarity_Sequence = "1,1,1,1, -1,-1,-1,1, -1,-1,-1,-1, 1,1,-1,1, -1,-1,1,1, -1,1,1,-1, 1,1,1,1, 1,1,-1,1,1,1,-1,1, 1,-1,-1,1, 1,1,-1,1, -1,-1,-1,1, -1,1,-1,-1, 1,-1,-1,1, 1,1,1,1, -1,-1,1,1,-1,-1,1,-1, 1,-1,1,1, -1,-1,-1,1, 1,-1,-1,-1, -1,1,-1,-1, 1,-1,1,1, 1,1,-1,1, -1,1,-1,1, -1,-1,-1,-1, -1,1,-1,1, 1,-1,1,-1, 1,1,1,-1, -1,1,-1,-1, -1,1,1,1, -1,-1,-1,-1, -1,-1,-1";
	cout << "Polarity " << Polarity_Sequence(index % 126) << endl;
	cout << "Polarity " << index % 126 << endl;
	
	return 0;
} 
*/

