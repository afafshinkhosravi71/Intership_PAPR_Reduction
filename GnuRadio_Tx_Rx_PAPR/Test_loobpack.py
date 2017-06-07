#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Test Loobpack
# Generated: Tue Jun  6 16:54:10 2017
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

import os
import sys
sys.path.append(os.environ.get('GRC_HIER_PATH', os.path.expanduser('~/.grc_gnuradio')))

from Error_Rate2 import Error_Rate2  # grc-generated hier_block
from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import channels
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
#from gnuradio.qtgui import Range, RangeWidget
from optparse import OptionParser
import numpy ; import math;
from gnuradio import qtgui

#################################################################################
try:
    from scipy.special import erfc
except ImportError:
    print "Error: could not import scipy (http://www.scipy.org/)"
    sys.exit(1)

try:
    import pylab
except ImportError:
    print "Error: could not import pylab (http://matplotlib.sourceforge.net/)"
    sys.exit(1)

# to choose
RAND_SEED = 42
N_BITS = 1e7

def berawgn(EbN0):
    """ Calculates theoretical bit error rate in AWGN (for BPSK and given Eb/N0) """
    return 0.5 * erfc(math.sqrt(10**(float(EbN0)/10)))

######################################################################################

class Test_loobpack(gr.top_block, Qt.QWidget):

    def __init__(self, EbN0):


        ##################################################
        # Variables
        ##################################################
        self.bits_per_symb = bits_per_symb = 2
        #self.EbN0 = EbN0 = 10.0
        self.samp_rate = samp_rate = 100000
        self.packet_len = packet_len = 50
        self.len_tag_key = len_tag_key = "packet_len"
        self.freq_offset = freq_offset = 0
        self.fft_len = fft_len = 64
        self.Noise_voltage = Noise_voltage = 1.0 / math.sqrt(bits_per_symb * (10**(EbN0/10.0)) )
	self.mod = mod = digital.constellation_qpsk()

	#-----------------------------------------------------------------------#
	# Source is N_BITS bits, non-repeated
        data = map(int, numpy.random.randint(0, self.mod.arity(), N_BITS/self.mod.bits_per_symbol()))
	#-----------------------------------------------------------------------#

        ##################################################
        # Blocks
        ##################################################
      
        self.digital_ofdm_tx_0 = digital.ofdm_tx(
        	  fft_len=fft_len, cp_len=fft_len/4,
        	  packet_length_tag_key=len_tag_key,
        	  bps_header=1,
        	  bps_payload=2,
        	  rolloff=0,
        	  debug_log=False,
        	  scramble_bits=False
        	 )
        self.digital_ofdm_rx_0 = digital.ofdm_rx(
        	  fft_len=fft_len, cp_len=fft_len/4,
        	  frame_length_tag_key='frame_'+"rx_len",
        	  packet_length_tag_key="rx_len",
        	  bps_header=1,
        	  bps_payload=2,
        	  debug_log=False,
        	  scramble_bits=False
        	 )
        self.channels_channel_model_0 = channels.channel_model(
        	noise_voltage=Noise_voltage,
        	frequency_offset=freq_offset * 1.0/fft_len,
        	epsilon=1.0,
        	taps=(1.0 + 1.0j, ),
        	noise_seed=0,
        	block_tags=True
        )
        self.blocks_vector_source_x_0 = blocks.vector_source_b(data, False, 1, ())
        self.blocks_vector_sink_x_0 = blocks.vector_sink_f(1)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_tag_debug_0 = blocks.tag_debug(gr.sizeof_char*1, 'Rx Packets', ""); self.blocks_tag_debug_0.set_display(True)
        self.blocks_stream_to_tagged_stream_0 = blocks.stream_to_tagged_stream(gr.sizeof_char, 1, packet_len, len_tag_key)
        self.Error_Rate2_0 = Error_Rate2()

        ##################################################
        # Connections
        ##################################################
        self.connect((self.Error_Rate2_0, 0), (self.blocks_vector_sink_x_0, 0))
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.digital_ofdm_tx_0, 0))
        self.connect((self.blocks_throttle_0, 0), (self.digital_ofdm_rx_0, 0))
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.Error_Rate2_0, 1))
        self.connect((self.blocks_vector_source_x_0, 0), (self.blocks_stream_to_tagged_stream_0, 0))
        self.connect((self.channels_channel_model_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.digital_ofdm_rx_0, 0), (self.Error_Rate2_0, 0))
        self.connect((self.digital_ofdm_rx_0, 0), (self.blocks_tag_debug_0, 0))
        self.connect((self.digital_ofdm_tx_0, 0), (self.channels_channel_model_0, 0))


    def get_bits_per_symb(self):
        return self.bits_per_symb

    def set_bits_per_symb(self, bits_per_symb):
        self.bits_per_symb = bits_per_symb
        self.set_Noise_voltage(1.0 / math.sqrt(self.bits_per_symb * (10**(self.EbN0/10.0)) ))

    def get_EbN0(self):
        return self.EbN0

    def set_EbN0(self, EbN0):
        self.EbN0 = EbN0
        self.set_Noise_voltage(1.0 / math.sqrt(self.bits_per_symb * (10**(self.EbN0/10.0)) ))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len
        self.blocks_vector_source_x_0.set_data(range(self.packet_len), ())
        self.blocks_stream_to_tagged_stream_0.set_packet_len(self.packet_len)
        self.blocks_stream_to_tagged_stream_0.set_packet_len_pmt(self.packet_len)

    def get_len_tag_key(self):
        return self.len_tag_key

    def set_len_tag_key(self, len_tag_key):
        self.len_tag_key = len_tag_key

    def get_freq_offset(self):
        return self.freq_offset

    def set_freq_offset(self, freq_offset):
        self.freq_offset = freq_offset
        self.channels_channel_model_0.set_frequency_offset(self.freq_offset * 1.0/self.fft_len)

    def get_fft_len(self):
        return self.fft_len

    def set_fft_len(self, fft_len):
        self.fft_len = fft_len
        self.channels_channel_model_0.set_frequency_offset(self.freq_offset * 1.0/self.fft_len)

    def get_Noise_voltage(self):
        return self.Noise_voltage

    def set_Noise_voltage(self, Noise_voltage):
        self.Noise_voltage = Noise_voltage
        self.channels_channel_model_0.set_noise_voltage(self.Noise_voltage)


############################################################################################

def simulate_ber(EbN0):
    """ All the work's done here: create flow graph, run, read out BER """
    print "Eb/N0 = %d dB" % EbN0
    fg = Test_loobpack(EbN0)
    fg.run()
    print " BER = %.4E " % numpy.sum(fg.blocks_vector_sink_x_0.data())
    return numpy.sum(fg.blocks_vector_sink_x_0.data())

if __name__ == "__main__":
    
    EbN0_min = 0
    EbN0_max = 15
    EbN0_range = range(EbN0_min, EbN0_max+1)
    ber_theory = [berawgn(x) for x in EbN0_range]
    print ber_theory
    print "Simulating..."
    ber_simu   = [simulate_ber(x) for x in EbN0_range]
    print ber_simu

    f = pylab.figure()
    s = f.add_subplot(1,1,1)
    s.semilogy(EbN0_range, ber_theory, 'g-.', label="Theoretical")
    s.semilogy(EbN0_range, ber_simu, 'b-o', label="Simulated")
    s.set_title('BER Simulation')
    s.set_xlabel('Eb/N0 (dB)')
    s.set_ylabel('BER')
    s.legend()
    s.grid()
    pylab.show()

##########################################################################################
