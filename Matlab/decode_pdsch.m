function [sib_info, syms] = decode_pdsch(peak, reg_info, dci_info, subframe_idx, tfg, ce, np_ce)

n_ports = peak.n_ports;

[syms, ce] = pdsch_extract(peak, reg_info, dci_info, subframe_idx, tfg, ce);

% % mmse_flag = 0;
% [syms, ~] = LTE_SFBC_demod(syms, ce, np_ce, n_ports, mmse_flag);

if dci_info.MCS >=0 && dci_info.MCS <= 9
    enb.PDSCH.Modulation = {'QPSK'};
    mod_type = 'QAM';
elseif dci_info.MCS >=10 && dci_info.MCS <= 16
    enb.PDSCH.Modulation = {'16QAM'};
    mod_type = 'QAM16';
elseif dci_info.MCS >=17 && dci_info.MCS <= 28
    enb.PDSCH.Modulation = {'64QAM'};
    mod_type = 'QAM24'; 
end

enb.PDSCH.RNTI = hex2dec('FFFF'); % SI-RNTI

if peak.n_ports == 1
    enb.PDSCH.TxScheme = 'Port0';
else
    enb.PDSCH.TxScheme = 'TxDiversity';
end
enb.PDSCH.NLayers = 1;
enb.PDSCH.RV = dci_info.RV;
enb.PDSCH.NTurboDecIts = 6;

enb.NCellID = peak.n_id_cell;
enb.NSubframe = subframe_idx;
enb.CellRefP = peak.n_ports;
if peak.duplex_mode == 1
    enb.DuplexMode = 'TDD';
else
    enb.DuplexMode = 'FDD';
end
enb.TDDConfig = peak.uldl_cfg;
enb.NDLRB = peak.n_rb_dl;
enb.CFI = reg_info.CFI;
if peak.cp_type_val == 0
    enb.CyclicPrefix = 'Normal';
else
    enb.CyclicPrefix = 'Extended';
end

% mmse_flag = 0;
% [syms, np] = LTE_SFBC_demod(syms, ce, np_ce, n_ports, mmse_flag);
% 
% % ------Extract the bits---------
% e_est = deqam(syms, np, mod_type, 'LTE');
% 
% % ----------Unscramble---------------
% c_init = subframe_idx*(2^9) + peak.n_id_cell + enb.PDSCH.RNTI*(2^14);
% scr = lte_pn(c_init, length(e_est));
% e_est(scr==1)=1-e_est(scr==1);
% llr = log(e_est./(1-e_est)).';

% ce_tmp = zeros(length(syms),1,n_ports);
% for i=1:n_ports
%     ce_tmp(:,1,i) = ce(i,:).';
% end

[llr, syms] = ltePDSCHDecode(enb, enb.PDSCH, syms(:));

if mod(dci_info.TPC, 2)
    SIB_nRB = 3;
else
    SIB_nRB = 2;
end
tbsize = get_tbsize(dci_info.MCS, SIB_nRB);

[bits, blkcrc, ~] = lteDLSCHDecode(enb, enb.PDSCH, tbsize, llr{1},[]);
sib_info.bits = bits;
sib_info.blkcrc = blkcrc;
