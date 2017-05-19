/************************************************************************/  
// IEEE 802.11 OFDM PHYSICAL LAYER
// Author : Boubacar DIALLO 
// Date   : May 2017
/************************************************************************/

#include "itpp/itcomm.h"
#include "Utils.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main()
{
	// Data Rate (6, 9, 12, 18, 24, 36 ,48,54 )(Mb/s)
	int data_rate;	
        
	//Coding rate
	double R;		
   
	// OFDM parameters
	int Ncbps;			// Coded bits per OFDM symbol 
	int Nbpsc;			// Coded bits per subcarrier 
	int Ndbps;			// Coded data bits OFDM symbol
	int N_SD = 48;		        // Number of data subcarriers
	int N_SP = 4;			// Number of pilot subcarriers
	int N_ST = N_SD + N_SP;	        // Number of subcarriers total = 52
	int Nfft = 64; 			// FFT size
	int Ncp = 16;			// Cyclic prefix size

	// Initialisation
	cvec data_ofdm_freq_symb;
	cvec data_ofdm_transmited_symb;
	cvec ofdm_received_symb;
	cvec ofdm_demodulated_symb;
	cvec received_data_symbols;
        cvec sent_pilots_final;

	/*************************** START of programm ***********************************/
   	header_display();
   	
	// Intially set data rate=0 to get the desired data rate value from the user
	data_rate = 0;
	// Selon le debit, il donne les valleurs des parameters				
	data_rate_choice( &R ,&Nbpsc, &data_rate, &Ncbps, &Ndbps);

  
	//Declarations of scalars and vectors: 
  	int i,x;
  	//double noise_variance;
	
	//bvec is a vector containing bits
  	bvec transmitted_bits;                  

  	Real_Timer tt;                   // The timer used to measure the execution time
  	tt.tic();			 // Reset and start the timer:

	/************* Declarations of parameters: ****************************************/
  	//Modulation Parameters
  	ivec deci;
  	cvec constl;				//constellation 
  	cvec modulated_symbols;			//modulated symbols
  	vec soft_demodulated_symbols;		//demodulated symbols

  	//scrambler parameters
  	bvec scrambled_bits,descrambled_bits;
  	bvec State = " 1 0 1 1 1 0 1 ";		//intial state 
  	//bvec State = randb(7);		//Random intial state
	//cout << "State Randb" << State << endl; 
  
  	//Code parameters
  	bvec coded_bits;
  	bvec decoded_bits;
  
  	//interleaver
  	bvec interleaved_bits;
  	vec received_symbols_deinterleaved;

	/****************************************************************************************/
	//Generating the preamble
  	cvec PREAMBLE = preamble();	
  	PREAMBLE = PREAMBLE*sqrt((double)(Nfft));	//multiplied by square root of Nfft to compensate for the normalizing factor of the OFDM 	

  	//Generating the Data field
   	int N_bits;					//number of bits being transmitted
	
	// from a file 
	bvec input_binary_vector;
	input_binary_vector=file2binary_converter();	//read data from a file 
	N_bits=input_binary_vector.length();		//number of bits in the file being transmitted 
	transmitted_bits=input_binary_vector;	

	/*****************************************************************************/	
	int max_bit_nbr;
	// max number of bits per frame. An ofdm frame in 802.11a is 4095 octets in the data field + the preamble and signal fields
	// this value needs a very sophisticated phase tracking algo so adapt another soution which is to estimate the channel more often 
	// max_bit_nbr=4095*8;							
	max_bit_nbr=256*8;	// here you can specify the max number of bits per frame (must be divisible by 8) to do the channel estimation for a shorter frame

	//number of ofdm frames needed to transmit this file
	int nbr_frames=ceil(N_bits/(double)(max_bit_nbr)); 	

	cout<<"Number of bits being transmitted : "<< N_bits << endl;
	cout<<"Number of frames being transmitted : "<< nbr_frames << endl;

	cvec PPDU_Frame_Final;

	for (int p=0;p<nbr_frames;p++)		//this for loop is for segementing data greater than 4095 octets into frames
	{	
		//cout << "Now simulating frame # :" << p+1 << "/" << nbr_frames << "..." << endl;
		
		bvec transmitted_bits_seg;		// vector contains the segmented bits  

		if (p!=nbr_frames-1)			// all except for the last one
			transmitted_bits_seg=transmitted_bits(p*(max_bit_nbr),(p+1)*(max_bit_nbr)-1);
			
		if (p==nbr_frames-1)
			transmitted_bits_seg=transmitted_bits(p*(max_bit_nbr),N_bits-1);

		//cout << "transmitted_bits_seg.length()" << transmitted_bits_seg.length() << endl;

        	int Nsys_data;							//the number of OFDM symbols in data field
		int Ndata;							//the number of bits in the DATA field
		int length;							// the length of the PSDU in bytes in the DATA field
		length=transmitted_bits_seg.length()/8.0;			// to follow the convention used in the IEEE 802.11 document (length in octets)
		Nsys_data=ceil((16 + 8*length + 6)/(double)Ndbps);		// 16 bit of the service field and 6 bits of the tail + length of the data
		Ndata=Nsys_data*Ndbps;
		
		//cout << "Number of bits in frame : "<< p+1 <<"  is : "<< length*8 << endl;
		// 5 = (2 short training sequence symobols + 2 long training sequence symobols + 1 signal symbol)
		//cout << "Number of OFDM symbols in frame : "<< p+1 <<"  is : "<< Nsys_data+5 << endl;
		//cout << "Size (nbr of samples) of frame : "<< p+1 <<"  is : "<< (Nsys_data+5)*80 << endl;
		//cout << "Number of bytes in the DATA field in frame : "<< p+1 <<"  is : "<< length << endl;


		//Generating the Signal field
		//cout << "/********** Generating the Signal field ************************/" << endl;
		if (length>4095)
			cout << "ERORR !!! Number of bytes of the PSDU is greater than 4095 octets "<<endl;	// normally there should not be any error here
    		
		cvec SIGNAL = Signal(length,data_rate);							// will be coded with R=1/2, BPSK-OFDM modulated
		//cout<< " Done " << endl;		

		//cout << "/***** Generating the data filed for the current frame ************/" << endl;
		//Generating the data filed for the current frame
		bvec DATA;
		DATA = Data(transmitted_bits_seg,Ncbps,Ndbps);			 
		//cout<< " Done " << endl;		

		//cout << "/*************** Scrambling ******************/" << endl;
		//***** Scrambling *****//
		//cout<<"Scrambling ... ";
		scrambled_bits= scrambler(DATA,State);
		//cout<< " Done " << endl;
 		
		//cout << "/*************** Coding ******************/" << endl;
		//***** Coding *****//
		//cout<<"Coding ... ";
		coded_bits= Punct_Conv_Code(scrambled_bits, R);
		//cout << "Done" << endl;

		//***** Interleaving *****//
		//cout << "/*************** Interleaving ******************/" << endl;
		interleaved_bits=interleaver(coded_bits, Ncbps,Nbpsc);
		//cout << "Done" << endl;
	
		//***** Modulation *****//
		//cout << "/*************** Modulation ******************/" << endl;
		modulated_symbols= modulation(interleaved_bits,Nbpsc, deci, constl);
		//cout << "Done" << endl;
	 
		//***** Pilot insertion and OFDM modulation *****//
		//cout << "/********* Pilot insertion **************/" << endl;
		int index=1;			
		// first index to choose the second element of the polarity sequence p0-126 and on since zero was chosen for the signal field
        	data_ofdm_freq_symb= Pilot_insertion( modulated_symbols,N_SD ,Nfft,index);
        
        	int polarity;
        	cvec sent_pilots;
        	sent_pilots.set_size(4*Nsys_data,false);
        	index=1;
        	
		for (int i=0;i<Nsys_data*4;i+=4) 
		{ 
             		polarity=Pilot_Polarity_Sequence(index); 
             		sent_pilots(i)=std::complex<double>(1*polarity,(double)0);
             		sent_pilots(i+1)=std::complex<double>(1*polarity,(double)0);
             		sent_pilots(i+2)=std::complex<double>(1*polarity,(double)0);
             		sent_pilots(i+3)=std::complex<double>(-1*polarity,(double)0);
             		index +=1;
        	}
		//cout << "Done" << endl;
	
		//cout << "/********* OFDM modulation **************/" << endl;
		data_ofdm_transmited_symb=ofdm_modulation(data_ofdm_freq_symb,Nfft ,Ncp,1);	
		//cout << "Done" << endl;

		//cout << "/********* Data Assembler **************/" << endl;
		cvec DATA_final=OFDM_Data_Symbols_Assembler(data_ofdm_transmited_symb);
		//cout<<"Done"<<endl;

		//cout<<"Assembling Frame ... ";
		//cout << "/********* Assembling Frame **************/" << endl;
		cvec PPDU_Frame= PPDU_Frame_Assembler(PREAMBLE,SIGNAL,DATA_final);
		//cout << "Done" << endl;

		//cout<<"PREAMBLE:  " << PREAMBLE.length() << endl;
		//cout<<"SIGNAL:  " << SIGNAL.length() << endl;
		//cout<<"DATA_final:  " << DATA_final.length() << endl;
		//cout<<"PPDU_Frame:  " << PPDU_Frame.length() << endl;
		//overlap the segmented frames : this method was not proposed by the standard 802.11
	
		/***************************************************************************/
		int nbr_zeros=0;		
		// number of zeros between the frames used if you want to separate between the ofdm frames
		cvec zerovec=zeros_c(nbr_zeros);
     
		if (p == 0) 		//First PPDU_Frame 
		{
			PPDU_Frame_Final=PPDU_Frame(0,PPDU_Frame.length() - 2);//-2 because no need for the last sample since it's used for the overlapping
			PPDU_Frame_Final.ins(PPDU_Frame_Final.length(), zerovec);
                        sent_pilots_final = sent_pilots;
		}
		else{
                    	//overlapping the frames
			PPDU_Frame_Final.ins(PPDU_Frame_Final.length(),PPDU_Frame(0,PPDU_Frame.length() - 2));	//-2 because no need for the last sample since it's used for the overlapping
			PPDU_Frame_Final.ins(PPDU_Frame_Final.length(), zerovec);
                        sent_pilots_final.ins(sent_pilots_final.length(),sent_pilots);
                 }
		//cout << "PPDU_Frame_Final:  " << PPDU_Frame_Final.length() << endl;
  	}
	

	//cout << sent_pilots_final << endl;
	// save into a file
        it_file ff1;
        ff1.open("sent_pilots_file.it");
        ff1 << Name("sent_pilots") << sent_pilots_final;
	ff1 << Name("number_of_frames") << nbr_frames;
        //ff1 << Name("Nsys_data") << Nsys_data;
        ff1.close();	

        // zeros to be to be added at the end of the transmission in order to be used by the MXA to start triggering in the RFBurst mode and in order to be used by the synch alog in the Emidate 
	// 30000 was chosen carefully to be compatabile with the synchronyzation algo used at the receiption side and the interframe distance in case of 6 Mbps 
	PPDU_Frame_Final.ins(PPDU_Frame_Final.length(), zeros_c(30000));
	
	// now we will add an additional ofdm symbol that contains the number of frames 
	// this is not part of the 802.11a standard but it is necessary for our system because it's a half duplex system 
	// this will be used by the receiver to know how many frames to take into account and ignore all the other frames.
	// this was necessary because the MXA aquires multiple copies of the same data
	// this field looks like the signal field but with no data rate or parity bit subfields
	// this field will be added only once after the signal field
	 
	cvec frame_nbr_field = frames_nbr_field(nbr_frames);
	PPDU_Frame_Final.ins(400,frame_nbr_field);
	//cout<<frame_nbr_field.length()<<endl;
	
	//Saving PPDU_Frame into Sent_Frame.it
	//Sent_Frame_Saver(PPDU_Frame_Final);

	cout << "/********* Assembling Frame **************/" << endl;
	cout << "PPDU_Frame_Final length() :  " << PPDU_Frame_Final.length() << endl;

	/********************************************************************************************************************/

	int op = 13;
	vec EbN0_dB = "0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0";
    	//double R = 1.0/3.0;    				//coding rate (non punctured PCCC)
    	double Ec = 1.0;      					//coded bit energy
	cout << "Now simulating EbN0db = " << EbN0_dB << endl;

	vec sigma2 = (0.5*Ec/R)*pow(inv_dB(EbN0_dB), -1.0); 	// N0/2
	cout << "Noise variance Sigma = " << sigma2 << endl;
	cout << "Noise variance Sigma used = " << sigma2(op) << endl;
	
	//AWGN channel
    	AWGN_Channel channel;
        channel.set_noise(sigma2(op));

	//Signal + AWGN channel
        PPDU_Frame_Final = channel(PPDU_Frame_Final);

	cout << " Transmission time : " << tt.toc() << " seconds" << endl;

	//Saving PPDU_Frame into Sent_Frame.it
	Sent_Frame_Saver(PPDU_Frame_Final);


	/********************************************************************************************************************/

	tt.tic();

	// Declare the it_file class
	it_file ff;
        //it_file ff1;
	bvec received_bits_final;
        cvec received_pilot_symbols_final;
	//cvec ofdm_received_symb;

	// File permition
	int res = system("sudo chmod 777 Sent_Frame.it");
	//cout << "result system(sudo chmod 777 Sent_Frame.it) = " << res << endl;

	// Open the file
	ff.open("Sent_Frame.it");

	// Read the variable freq_cor_IQData from the file. Put result in vector ofdm_received_symb.
	ff >> Name("PPDU_Frame_Final") >> ofdm_received_symb;
	ff.close();
	
	//intially is set to zero then the real number of frames will be determined from the received frames_nbr_field that was added at the end of the transmission
	nbr_frames = 1;
	frame_nbr_field = ofdm_received_symb(400,400+80-1);
	// suppress this part
	ofdm_received_symb.del(400,400+80-1);

	/****************** Receiver Side ********************************/
	for (int p=0;p<nbr_frames;p++)
	{
		//cout<<"Receiving frame number : "<<p+1<<endl;
		// short sequence recuperation length = 160
		cvec received_short_squence = ofdm_received_symb.left(160);

		//AGC for every frame: (Average value criterion)
		double mean1 = 0.8663;					//This is the mean value of the transmitted Short Training Seq.
		double mean2 = mean(abs(received_short_squence));
		double AGC = mean1/mean2;
		ofdm_received_symb=AGC*ofdm_received_symb;
		//cout<< "AGC for every frame = " << AGC<<endl;
		
		// long training sequence recuperation length = 160
		cvec received_long_squence=ofdm_received_symb(160,319);

		// Signal field recuperation length = 80
		cvec received_signal_field = ofdm_received_symb(320,399);


		//Channel estimation
		cvec Hh_LS = LS_estimator(received_long_squence);
		//cvec Hh_MMSE=MMSE_estimator(received_long_squence);

		//choose one of these methods
		cvec Hhat = Hh_LS;
		//cvec Hhat=Hh_MMSE;

		// Noise variance of the received long training seq
		double N0 = noise_variance(received_long_squence); 
		//cout << N0 << endl;

		//a zero value poses problems for the channel estimation
		if (N0==0)		
			N0=1e-6;
		
		//update the number of frames and get real number of frames 
		if (p==0)
			nbr_frames = frame_nbr_receiver(frame_nbr_field, N0, Hhat);
		//cout<<"Number of frames: "<< nbr_frames << endl;

		// received data rate recuperationnoise variance vs SNR
		ivec output = Receiver_data_rate(received_signal_field, N0, Hhat);

		int rec_data_rate=output(0);			// received  data rate must be one of these values (6, 9, 12, 18, 24, 36 ,48,54 )(Mb/s)
		int rec_length=output(1);			// received length 
		//cout <<"data rate : "<< rec_data_rate << endl << "length : "<< rec_length << endl; //tested ok

		//this function will set all the parameters according to the data rate value determined by the receiver
		data_rate_choice( &R , &Nbpsc, &rec_data_rate, &Ncbps, &Ndbps);	

		int Nsys_data;					//the number of OFDM symbols in data field
		int Ndata;					//the number of bits in the DATA field
		
		// 16 bit of the service field and 6 bits of the tail + length of the data
		Nsys_data=ceil((16 + 8*rec_length + 6)/(double)Ndbps);		
		
		Ndata=Nsys_data*Ndbps;
		
		//cout<<"Nsys_data : "<<Nsys_data+5<<endl;
		//cout<<"Ndata : "<<Ndata<<endl;

		// added to separate two consecutive frames you have to use the same value used at the receiver side				
		int nbr_zeros=0;														

		// data field for this frame (without the five first symbols)
		cvec received_data_field = ofdm_received_symb(5*80,5*80+Nsys_data*80-1);	
		//cout<<ofdm_received_symb(1000+400+Nsys_data*80-1,1000+400+Nsys_data*80+1);// in case of there were 1000 zeros had been added

		// 80*5 is the first five ofdm symbols ; nbr_zeros sample were added to separate two consecutive frames
		ofdm_received_symb.del(0,nbr_zeros+80*5+Nsys_data*80-1);				

		// OFDM demodulation
		cvec ofdm_demodulated_symb=ofdm_demodulation(received_data_field,Nfft ,Ncp,1);
		
		//Equalization
		cvec Equalized_output=Equalization(ofdm_demodulated_symb,Hhat);

		//index starts from 1 for the data field
		int index=1;														
		cvec received_data_symbols=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index);
		
		//index starts from 1 for the received pilote field
		int index2 = 0;
                cvec received_pilot_symbols=Pilot_Data_separation(Equalized_output, N_SD ,Nfft, index2);

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
        	//cout<<"end of for loop for iteration number"<<p+1<<endl;
	
	 }//end for 
	ff1.open("received_pilots_file.it");
        ff1 << Name("rec_pilots") << received_pilot_symbols_final;
	ff1 << Name("number_of_frames") << nbr_frames;
        //ff1 << Name("Nsys_data") << Nsys_data;
        ff1.close();
	
	cout << " Reception time : " << tt.toc() << " seconds" << endl;

	// Converting the file into its original form
	binary2file_converter(received_bits_final);

	/****************************************************************************************/
	/****************************************************************************************/

	//Calculate the bit error rate:
	BERC berc;    
    	berc.clear();						//Clear the bit error rate counter
    	berc.count(input_binary_vector,received_bits_final);	//Count the bit errors
    	double bit_error_rate = berc.get_errorrate();
	double bit_error = berc.get_errors();
	cout << "Bit Error Rate : " << bit_error_rate << endl;
	cout << "Number of bit Errors : " << bit_error << endl;
	cout << "Number of bits : " << input_binary_vector.length() << endl;

	// save the bit vectors into an it file
	it_file a;
	a.open("bits.it");
	a << Name("sent_bits") << input_binary_vector;
	a << Name("rec_bits") << received_bits_final;
	a.close();

	//**************************************** END ******************************************************//	
	return 0;
}

