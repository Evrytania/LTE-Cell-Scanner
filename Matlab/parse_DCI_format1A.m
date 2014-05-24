function dci_str = parse_DCI_format1A(peak, space_ind, bits)
n_rb_dl = peak.n_rb_dl;
% n_rb_ul = n_rb_dl;
duplex_mode = peak.duplex_mode;
% uldl_cfg = peak.uldl_cfg;

if space_ind ~= 0
    disp('Now only common space is supported (space_ind = 0)');
    return;
end

FLAG_len = 1;
VRB_len = 1;
RBassign_len = ceil( log2( n_rb_dl*(n_rb_dl+1)/2 ) );

MCS_len = 5;
if duplex_mode == 1
    HARQ_len = 4;
else
    HARQ_len = 3;
end
NEW_len = 1;
RV_len = 2;
TPC_len = 2;
if duplex_mode == 1
    DAI_len = 2;
else
    DAI_len = 0;
end

sp = 1;

FLAG = bi2de( bits(sp : sp+FLAG_len-1), 'left-msb');
sp = sp + FLAG_len;

VRB = bi2de( bits(sp : sp+VRB_len-1), 'left-msb');
sp = sp + VRB_len;

RBassign = bi2de( bits(sp : sp+RBassign_len-1), 'left-msb');
sp = sp + RBassign_len;

MCS = bi2de( bits(sp : sp+MCS_len-1), 'left-msb');
sp = sp + MCS_len;

HARQ = bi2de( bits(sp : sp+HARQ_len-1), 'left-msb');
sp = sp + HARQ_len;

NEW = bi2de( bits(sp : sp+NEW_len-1), 'left-msb');
sp = sp + NEW_len;

RV = bi2de( bits(sp : sp+RV_len-1), 'left-msb');
sp = sp + RV_len;

TPC = bi2de( bits(sp : sp+TPC_len-1), 'left-msb');
sp = sp + TPC_len;

DAI = bi2de( bits(sp : sp+DAI_len-1), 'left-msb');

dci_str = ['FLAG-' num2str(FLAG) ' VRB-' num2str(VRB) ...
     ' RBassign-' num2str(RBassign)  ' MCS-' num2str(MCS) ...
      ' HARQ-' num2str(HARQ)  ' NEWind-' num2str(NEW) ...
       ' RV-' num2str(RV)  ' TPC-' num2str(TPC) ...
        ' DAI-' num2str(DAI)];
