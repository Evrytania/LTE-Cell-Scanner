.bin files are captured by rtlsdr, hackrf, usrp, or bladerf.
The content of file is in raw captured binary format. No any header!

rtlsdr capture example:
rtl_sdr -f 1860e6 -s 1.92e6 -n 1.92e6 f1860_s1.92_1s_rtlsdr.bin

hackrf capture example:
hackrf_transfer -s 19200000 -f 1860000000 -a 0 -l 32 -g 34 -b 20000000 -n 19200000 -r f1860_s19.2_bw20_1s_hackrf.bin

naming rule:

begin with fxxxx : xxxx represents carrier frequency in MHz
...
_syyyy : yyyy represents sampling rate in MHz
...
board string will also be contained in file name: rtlsdr/hackrf/usrp/bladerf : represent board which is used for signal capturing.

You may insert other segment for more information, but program only parses fxxxx, _sxxxx, and rtlsdr/hackrf/usrp/bladerf



