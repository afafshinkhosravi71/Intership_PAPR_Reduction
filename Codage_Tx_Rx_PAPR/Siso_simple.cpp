#include <itpp/itcomm.h>
#include <itpp/itstat.h>
#include <iostream>
#include "itpp/itcomm.h"
#include "Utils.h"

using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

// g++ `pkg-config --cflags itpp` -o Siso Siso_simple.cpp `pkg-config --libs itpp`

#include <string>
#include <fstream>

int main(int argc, char *argv[])
{		

		cout << " Let's start with finale program " << endl;
        
		// generation d'une source binaire
        	double msg_size = 4096;
        	Bernoulli_RNG genb(0.7);
        
   		bvec data;
   		genb.sample_vector(msg_size, data);
	
		cout << "/*********************** Message Source *******************************/" << endl;
   		cout << "Message source Taille ==> "<< data.length() << endl;
		
		
		//Modulation Parameters
		int Nbpsc = 4;
  		ivec deci;
  		cvec constl;				//constellation 
  		cvec modulated_symbols;			//modulated symbols
		cvec data_ofdm_freq_symb;
		int N_SD = 48; 
		int Nfft = 64;
		int Ncp = 16;				//Cyclic prefix size

  		cout << "/*************** Modulation ******************/" << endl;

		modulated_symbols = modulation(data, Nbpsc, deci, constl);
   		cout << "Message modulé ==> "<< modulated_symbols.length() << endl;
		//cout << "Done" << endl;

	 	it_file ff1;
       		ff1.open("modulated.it");
        	ff1 << Name("modulated_symbols") << modulated_symbols;
        	ff1.close();

		
		cout << "/********* Insertion de zeros **************/" << endl;
		
		int NSymb = ceil(modulated_symbols.length()/N_SD) + 1; 
		cout << "NSymb " << NSymb << endl;                     
		
		cvec seros = zeros_c(N_SD * NSymb - modulated_symbols.length());
		cout << "Nbre zeros " << seros.length() << endl;

		modulated_symbols = concat( modulated_symbols, seros); 
		cout << "Insertion zeros " << modulated_symbols.length() << endl;

		cout << "/********* Pilot insertion **************/" << endl;

		int index=0;			
        	data_ofdm_freq_symb = Pilot_insertion( modulated_symbols, N_SD , Nfft, index);
		cout << "Insertion pilote " << data_ofdm_freq_symb.length() << endl;
		//cout << "Done" << endl;
		
       		ff1.open("Pilot.it");
        	ff1 << Name("modulated_symbols_pilots") << data_ofdm_freq_symb;
        	ff1.close();
		
		//cout << "/********* Tone Reservation **************/" << endl;
		
		/*data_ofdm_freq_symb = tone_reserv(data_ofdm_freq_symb);
		cout << "Tone Reservation " << data_ofdm_freq_symb.length() << endl;
		//cout << "Done" << endl;

		ff1.open("Tone_reser.it");
        	ff1 << Name("Tone_reserv") << data_ofdm_freq_symb;
        	ff1.close();
		*/

		cout << "/********* OFDM modulation and add GI **************/" << endl;

		cvec data_ofdm_transmited_symb = ofdm_modulation(data_ofdm_freq_symb ,Nfft ,Ncp, 1);	
		cout << "Modulation OFDM and add GI " << data_ofdm_transmited_symb.length() << endl;
		//cout << "Done "<< endl;

		ff1.open("data_ofdm.it");
        	ff1 << Name("data_ofdm") << data_ofdm_transmited_symb;
        	ff1.close();
		

		cout << "/********* Clipping method **************/" << endl;
		
 		double A_clip = 1.65;
		cvec data_ofdm_transmited_symb_clip = clipping(data_ofdm_transmited_symb, A_clip);
		cout << "Clipping " << data_ofdm_transmited_symb_clip.length() << endl;
		//cout << "Done "<< endl;

		ff1.open("Clipping.it");
        	ff1 << Name("Clipping") << data_ofdm_transmited_symb_clip;
		
        	ff1.close();
		
		//cout<<"Done"<<endl;
	
}

