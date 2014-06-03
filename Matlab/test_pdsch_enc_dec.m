% function test_pdsch_enc_dec
clear all;
close all;

enb.DuplexMode = 'TDD';
enb.TDDConfig = 2;

enb.PDSCH.Modulation = 'QPSK';
enb.PDSCH.NLayers = 1;
enb.PDSCH.RV = 0;
enb.PDSCH.NLayers = 1;
enb.PDSCH.TxScheme = 'Port0';
enb.PDSCH.NTurboDecIts = 5;

tbsize = 32;
trblkData = round(rand(tbsize,1));
[cw, chinfo] = lteDLSCH(enb, enb.PDSCH, 3*tbsize+12+24, trblkData);

cw(cw==0) = -1;

[bits, blkcrc, decState] = lteDLSCHDecode(enb, enb.PDSCH, tbsize, cw,[]);

tdd_max_num_HARQ_table = [4 7 10 9 12 15 6];
max_num_HARQ = tdd_max_num_HARQ_table(enb.TDDConfig+1);
[bits1, blkcrc1] = pdsch_bit_level_proc( cw.', tbsize, max_num_HARQ, enb.PDSCH.RV);
plot(xor(bits.', bits1(1:tbsize)));
