#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Dvbt Rx Demo 8K
# Generated: Tue Apr 18 22:11:36 2017
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

from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import dtv
from gnuradio import eng_notation
from gnuradio import fft
from gnuradio import gr
from gnuradio import qtgui
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from optparse import OptionParser
import sip
import sys
import time
from gnuradio import qtgui


class dvbt_rx_demo_8k(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Dvbt Rx Demo 8K")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Dvbt Rx Demo 8K")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "dvbt_rx_demo_8k")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = (8000000.0 * 8) / 7
        self.center_freq = center_freq = 429e6

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("num_recv_frames=128,master_clock_rate=" + str(samp_rate*4), "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_0.set_gain(30, 0)
        self.uhd_usrp_source_0.set_antenna('RX2', 0)
        self.qtgui_freq_sink_x_0 = qtgui.freq_sink_c(
        	1024, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	center_freq, #fc
        	samp_rate, #bw
        	"", #name
        	1 #number of inputs
        )
        self.qtgui_freq_sink_x_0.set_update_time(0.10)
        self.qtgui_freq_sink_x_0.set_y_axis(-170, -10)
        self.qtgui_freq_sink_x_0.set_y_label('Relative Gain', 'dB')
        self.qtgui_freq_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_0.enable_grid(True)
        self.qtgui_freq_sink_x_0.set_fft_average(0.2)
        self.qtgui_freq_sink_x_0.enable_axis_labels(True)
        self.qtgui_freq_sink_x_0.enable_control_panel(False)

        if not True:
          self.qtgui_freq_sink_x_0.disable_legend()

        if "complex" == "float" or "complex" == "msg_float":
          self.qtgui_freq_sink_x_0.set_plot_pos_half(not True)

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_freq_sink_x_0_win)
        self.qtgui_const_sink_x_0_0 = qtgui.const_sink_c(
        	1024, #size
        	"", #name
        	1 #number of inputs
        )
        self.qtgui_const_sink_x_0_0.set_update_time(0.10)
        self.qtgui_const_sink_x_0_0.set_y_axis(-1.5, 1.5)
        self.qtgui_const_sink_x_0_0.set_x_axis(-1.5, 1.5)
        self.qtgui_const_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, "")
        self.qtgui_const_sink_x_0_0.enable_autoscale(True)
        self.qtgui_const_sink_x_0_0.enable_grid(False)
        self.qtgui_const_sink_x_0_0.enable_axis_labels(True)

        if not True:
          self.qtgui_const_sink_x_0_0.disable_legend()

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "red", "red", "red",
                  "red", "red", "red", "red", "red"]
        styles = [0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0]
        markers = [0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_const_sink_x_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_const_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_const_sink_x_0_0.set_line_width(i, widths[i])
            self.qtgui_const_sink_x_0_0.set_line_color(i, colors[i])
            self.qtgui_const_sink_x_0_0.set_line_style(i, styles[i])
            self.qtgui_const_sink_x_0_0.set_line_marker(i, markers[i])
            self.qtgui_const_sink_x_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_const_sink_x_0_0_win = sip.wrapinstance(self.qtgui_const_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_const_sink_x_0_0_win)
        self.fft_vxx_0 = fft.fft_vcc(8192, True, (window.rectangular(8192)), True, 1)
        self.dtv_dvbt_viterbi_decoder_0 = dtv.dvbt_viterbi_decoder(dtv.MOD_16QAM, dtv.NH, dtv.C2_3, 768)
        self.dtv_dvbt_symbol_inner_interleaver_0 = dtv.dvbt_symbol_inner_interleaver(6048, dtv.T8k, 0)
        self.dtv_dvbt_reed_solomon_dec_0 = dtv.dvbt_reed_solomon_dec(2, 8, 0x11d, 255, 239, 8, 51, 8)
        self.dtv_dvbt_ofdm_sym_acquisition_0 = dtv.dvbt_ofdm_sym_acquisition(1, 8192, 6817, 256, 30)
        self.dtv_dvbt_energy_descramble_0 = dtv.dvbt_energy_descramble(8)
        self.dtv_dvbt_demod_reference_signals_0 = dtv.dvbt_demod_reference_signals(gr.sizeof_gr_complex, 8192, 6048, dtv.MOD_16QAM, dtv.NH, dtv.C2_3, dtv.C2_3, dtv.GI_1_32, dtv.T8k, 1, 0)
        self.dtv_dvbt_demap_0 = dtv.dvbt_demap(6048, dtv.MOD_16QAM, dtv.NH, dtv.T8k, 1)
        self.dtv_dvbt_convolutional_deinterleaver_0 = dtv.dvbt_convolutional_deinterleaver(136, 12, 17)
        self.dtv_dvbt_bit_inner_deinterleaver_0 = dtv.dvbt_bit_inner_deinterleaver(6048, dtv.MOD_16QAM, dtv.NH, dtv.T8k)
        self.blocks_vector_to_stream_1 = blocks.vector_to_stream(gr.sizeof_char*1, 6048)
        self.blocks_vector_to_stream_0_0 = blocks.vector_to_stream(gr.sizeof_gr_complex*1, 6048)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, 'advdvbt.ts', False)
        self.blocks_file_sink_0.set_unbuffered(False)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_vector_to_stream_0_0, 0), (self.qtgui_const_sink_x_0_0, 0))
        self.connect((self.blocks_vector_to_stream_1, 0), (self.dtv_dvbt_viterbi_decoder_0, 0))
        self.connect((self.dtv_dvbt_bit_inner_deinterleaver_0, 0), (self.blocks_vector_to_stream_1, 0))
        self.connect((self.dtv_dvbt_convolutional_deinterleaver_0, 0), (self.dtv_dvbt_reed_solomon_dec_0, 0))
        self.connect((self.dtv_dvbt_demap_0, 0), (self.dtv_dvbt_symbol_inner_interleaver_0, 0))
        self.connect((self.dtv_dvbt_demod_reference_signals_0, 0), (self.blocks_vector_to_stream_0_0, 0))
        self.connect((self.dtv_dvbt_demod_reference_signals_0, 0), (self.dtv_dvbt_demap_0, 0))
        self.connect((self.dtv_dvbt_energy_descramble_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.dtv_dvbt_ofdm_sym_acquisition_0, 0), (self.fft_vxx_0, 0))
        self.connect((self.dtv_dvbt_reed_solomon_dec_0, 0), (self.dtv_dvbt_energy_descramble_0, 0))
        self.connect((self.dtv_dvbt_symbol_inner_interleaver_0, 0), (self.dtv_dvbt_bit_inner_deinterleaver_0, 0))
        self.connect((self.dtv_dvbt_viterbi_decoder_0, 0), (self.dtv_dvbt_convolutional_deinterleaver_0, 0))
        self.connect((self.fft_vxx_0, 0), (self.dtv_dvbt_demod_reference_signals_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.dtv_dvbt_ofdm_sym_acquisition_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.qtgui_freq_sink_x_0, 0))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "dvbt_rx_demo_8k")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.qtgui_freq_sink_x_0.set_frequency_range(self.center_freq, self.samp_rate)

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_source_0.set_center_freq(self.center_freq, 0)
        self.qtgui_freq_sink_x_0.set_frequency_range(self.center_freq, self.samp_rate)


def main(top_block_cls=dvbt_rx_demo_8k, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
