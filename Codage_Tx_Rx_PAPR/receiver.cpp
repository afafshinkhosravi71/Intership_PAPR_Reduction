//#include "stdafx.h"
#include "itpp/itcomm.h"
#include "802.11_Phy.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main()
{


	double R;		//Coding rate
   
	//ofdm parameters
	int Ncbps;			//Coded bits per OFDM symbol 
	int Nbpsc;			//Coded bits per subcarrier 
	int Ndbps;			//Coded data bits OFDM symbol
	int N_SD=48;			//Number of data subcarriers
	int N_SP=4;			//Number of pilot subcarriers
	int N_ST=N_SD+N_SP;		//Number of subcarriers,total
	int Nfft=64;
	int Ncp=16;
	bvec State = " 1 0 1 1 1 0 1 ";		//scrambler intial state 

	
	cvec ofdm_received_symb;
	// Declare the it_file class
	it_file ff;
        it_file ff1;

	// Open the file
	ff.open("IQDATA_matlab.it");

	// Read the variable freq_cor_IQData from the file. Put result in vector ofdm_received_symb.
	ff >> Name("freq_cor_IQData") >> ofdm_received_symb;
	ff.close();
	
	//intially is set to zero then the real number of frames will be determined from the received frames_nbr_field that was added at the end of the transmission
	int nbr_frames=1;
	cvec frame_nbr_field=ofdm_received_symb(400,400+80-1);
	ofdm_received_symb.del(400,400+80-1);
	
	bvec received_bits_final;
        cvec received_pilot_symbols_final;

	for (int p=0;p<nbr_frames;p++)
	{
//***** Receiver Side ****//
		cout<<"Receiving frame number : "<<p+1<<endl;
		cvec received_short_squence=ofdm_received_symb.left(160);

		//AGC for every frame
		//Average value criterion
		double mean1= 0.8663;		//This is the mean value of the transmitted Short Training Seq.
		double mean2=mean(abs(received_short_squence));
		double AGC=mean1/mean2;
		ofdm_received_symb=AGC*ofdm_received_symb;
		//cout<<AGC<<endl;

		cvec received_long_squence=ofdm_received_symb(160,319);
		cvec received_signal_field=ofdm_received_symb(320,399);


		//Channel estimation
		cvec Hh_LS=LS_estimator(received_long_squence);
		cvec Hh_MMSE=MMSE_estimator(received_long_squence);

		//choose one of these methods
		cvec Hhat=Hh_LS;
		//cvec Hhat=Hh_MMSE;

		double N0=noise_variance(received_long_squence);//cout<<N0<<endl;

		if (N0==0)		//a zero value poses problems for the channel estimation
			N0=1e-6;
		
		//update the number of frames and get real number of frames 
		if (p==0)
		nbr_frames=frame_nbr_receiver(frame_nbr_field,N0,Hhat);
		cout<<"Number of frames: "<<nbr_frames<<endl;

		ivec output=Receiver_data_rate(received_signal_field,N0,Hhat);
		int rec_data_rate=output(0);			// received  data rate must be one of these values (6, 9, 12, 18, 24, 36 ,48,54 )(Mb/s)
		int rec_length=output(1);				// received length 
		//cout <<"data rate : "<<rec_data_rate<<endl<<"length : "<<rec_length<<endl;//tested ok

		data_rate_choice( &R , &Nbpsc, &rec_data_rate, &Ncbps, &Ndbps);	//this function will set all the parameters according to the data rate value determined by the receiver

		int Nsys_data;					//the number of OFDM symbols in data field
		int Ndata;						//the number of bits in the DATA field
		
		Nsys_data=ceil((16 + 8*rec_length + 6)/(double)Ndbps);		// 16 bit of the service field and 6 bits of the tail + length of the data
		
		Ndata=Nsys_data*Ndbps;
		
		cout<<"Nsys_data : "<<Nsys_data+5<<endl;
		//cout<<"Ndata : "<<Ndata<<endl;

		//AGC for every frame
		//Average value criterion
		//double mean1= 0.8663;		//This is the mean value of the transmitted Short Training Seq.
		//double mean2=mean(abs(received_short_squence));
		//double AGC=mean1/mean2;
		//cout<<AGC<<endl;
		//int y;cin>>y;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////				
		int nbr_zeros=0;														// added to separate two consecutive frames you have to use the same value used at the receiver side
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		cvec received_data_field=ofdm_received_symb(5*80,5*80+Nsys_data*80-1);	// data field for this frame (without the five first symbols)
		//cout<<ofdm_received_symb(1000+400+Nsys_data*80-1,1000+400+Nsys_data*80+1);// in case of there were 1000 zeros had been added
		ofdm_received_symb.del(0,nbr_zeros+80*5+Nsys_data*80-1);				// 80*5 is the first five ofdm symbols ; nbr_zeros sample were added to separate two consecutive frames

		cvec ofdm_demodulated_symb=ofdm_demodulation(received_data_field,Nfft ,Ncp,1);
		//Equalization
		cvec Equalized_output=Equalization(ofdm_demodulated_symb,Hhat);

		int index=1;														//index starts from 1 for the data field
		cvec received_data_symbols=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index,1);
                cvec received_pilot_symbols=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index,0);

		vec soft_demodulated_symbols=soft_demodulation(received_data_symbols, Nbpsc,N0);

		vec	received_symbols_deinterleaved=deinterleaver(soft_demodulated_symbols, Ncbps,Nbpsc);
	
		bvec decoded_bits= Punct_Conv_Decode(received_symbols_deinterleaved, R);

		bvec descrambled_bits= descrambler(decoded_bits,State);
	
		// flip all bits (that is, convert zeros to ones, and ones to zeros) for the BPSK
		if (Nbpsc==1)
			for (int i = 0; i < descrambled_bits.length(); i++)  
				descrambled_bits(i) = (descrambled_bits(i) == 0 ? 1 : 0);	
	
		// Removing the service bits, the tail, and pad bits		
		bvec received_bits=descrambled_bits(16,16+8*rec_length-1);
	
		// Regrouping the received bits 
		if (p == 0) {
			received_bits_final=received_bits;
                        received_pilot_symbols_final=received_pilot_symbols;
                }
		else {
			received_bits_final.ins(received_bits_final.length(),received_bits);
		        received_pilot_symbols_final.ins(received_pilot_symbols_final.length(),received_pilot_symbols);
		}
        cout<<"end of for loop for iteration number"<<p<<endl;
	
	 }//end for 
	ff1.open("received_pilots_file.it");
        ff1 << Name("rec_pilots") << received_pilot_symbols_final;
	ff1 << Name("number_of_frames") << nbr_frames;
        //ff1 << Name("Nsys_data") << Nsys_data;
        ff1.close();
	// Converting the file into its original form
	binary2file_converter(received_bits_final);		
	

	
//  int x; cin>> x;

	return 0;
}
