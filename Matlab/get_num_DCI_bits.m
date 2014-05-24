function num_bits = get_num_DCI_bits(peak, format_str, space_ind)
% space_ind: 0 common space; 1 UE specific space
% 
% number of DCI bits in 36.212
% now only 0/1A/1C/3/3A are supported.

if space_ind ~= 0
    disp('Now only common space is supported (space_ind = 0)');
    return;
end
% n_rb_dl = peak.n_rb_dl;
% n_rb_ul = n_rb_dl;
% duplex_mode = peak.duplex_mode;
% uldl_cfg = peak.uldl_cfg;
% 
% abgs_bits = [12, 14, 16 ,20, 24, 26, 32, 40, 44, 56];
% 
% num_bits = 0;
if strcmpi(format_str, '0-with-1a')
    num_bits = get_num_DCI_format0_bits(peak, '1a', space_ind);
elseif strcmpi(format_str, '0-with-1a-ra-c-rnti')
    num_bits = get_num_DCI_format0_bits(peak, '1a-ra-c-rnti', space_ind);
elseif strcmpi(format_str, '1a')
    num_bits = get_num_DCI_format1A_bits(peak, space_ind);
elseif strcmpi(format_str, '1a-ra-c-rnti')
    num_bits = get_num_DCI_format1A_RA_C_RNTI_bits(peak, space_ind);
elseif strcmpi(format_str, '1c')
    disp('format1C not supported currently!');
    num_bits = -inf;
    return;
elseif strcmpi(format_str, '3-with-1a')
    num_bits_format0 = get_num_DCI_format0_bits(peak, '1a', space_ind);
    N = floor(num_bits_format0/2);
    num_bits = 2*N;
    if N < (num_bits_format0/2)
        num_bits = num_bits + 1;
    end
elseif strcmpi(format_str, '3-with-1a-ra-c-rnti')
    num_bits_format0 = get_num_DCI_format0_bits(peak, '1a-ra-c-rnti', space_ind);
    N = floor(num_bits_format0/2);
    num_bits = 2*N;
    if N < (num_bits_format0/2)
        num_bits = num_bits + 1;
    end
elseif strcmpi(format_str, '3a-with-1a')
    num_bits_format0 = get_num_DCI_format0_bits(peak, '1a', space_ind);
    M = num_bits_format0;
    num_bits = M;
elseif strcmpi(format_str, '3a-with-1a-ra-c-rnti')
    num_bits_format0 = get_num_DCI_format0_bits(peak, '1a-ra-c-rnti', space_ind);
    M = num_bits_format0;
    num_bits = M;
else
    disp('Now only 0-with-1A/0-with-1A-RA-C-RNTI/1A/1A-RA-C-RNTI/3-with-1A/3-with-1A-RA-C-RNTI/3A-with-1A/3A-with-1A-RA-C-RNTI are supported!');
    num_bits = -inf;
    return;
end

function num_bits = get_num_DCI_format0_bits(peak, format_str, space_ind)
n_rb_dl = peak.n_rb_dl;
n_rb_ul = n_rb_dl;
duplex_mode = peak.duplex_mode;
uldl_cfg = peak.uldl_cfg;

% abgs_bits = [12, 14, 16 ,20, 24, 26, 32, 40, 44, 56];

num_bits = 0;

if space_ind == 1
    num_bits = num_bits + 3;
end
num_bits = num_bits + 1;
num_bits = num_bits + ceil( log2( n_rb_ul*(n_rb_ul+1)/2 ) );
num_bits = num_bits + 5;
num_bits = num_bits + 1;
num_bits = num_bits + 2;
num_bits = num_bits + 3;
if duplex_mode == 1 && uldl_cfg == 0
    num_bits = num_bits + 2;
end
if duplex_mode == 1 && ( uldl_cfg == 1 || uldl_cfg == 2 || uldl_cfg == 3 || uldl_cfg == 4 || uldl_cfg == 5 || uldl_cfg == 6 )
    num_bits = num_bits + 2;
end
num_bits = num_bits + 1; % careful
num_bits = num_bits + 0; % careful
num_bits = num_bits + 0; % careful

if strcmpi(format_str, '1a')
    num_bits_1A = get_num_DCI_format1A_bits(peak, space_ind);
elseif strcmpi(format_str, '1a-ra-c-rnti')
    num_bits_1A = get_num_DCI_format1A_RA_C_RNTI_bits(peak, space_ind);
end

if num_bits < num_bits_1A
    num_bits = num_bits_1A;
end

function num_bits = get_num_DCI_format1A_bits(peak, space_ind)
n_rb_dl = peak.n_rb_dl;
% n_rb_ul = n_rb_dl;
duplex_mode = peak.duplex_mode;
% uldl_cfg = peak.uldl_cfg;

abgs_bits = [12, 14, 16 ,20, 24, 26, 32, 40, 44, 56];

num_bits = 0;
if space_ind == 1
    num_bits = num_bits + 3;
end
num_bits = num_bits + 1;
num_bits = num_bits + 1;
num_bits = num_bits + ceil( log2( n_rb_dl*(n_rb_dl+1)/2 ) );

num_bits = num_bits + 5;
if duplex_mode == 1
    num_bits = num_bits + 4;
else
    num_bits = num_bits + 3;
end
num_bits = num_bits + 1;
num_bits = num_bits + 2;
num_bits = num_bits + 2;
if duplex_mode == 1
    num_bits = num_bits + 2;
end
num_bits = num_bits + 0;

if prod(num_bits - abgs_bits) == 0
    num_bits = num_bits + 1;
end

function num_bits = get_num_DCI_format1A_RA_C_RNTI_bits(peak, space_ind)
n_rb_dl = peak.n_rb_dl;
% n_rb_ul = n_rb_dl;
duplex_mode = peak.duplex_mode;
% uldl_cfg = peak.uldl_cfg;

abgs_bits = [12, 14, 16 ,20, 24, 26, 32, 40, 44, 56];

num_bits = 0;

if space_ind == 1
    num_bits = num_bits + 3;
end
num_bits = num_bits + 1;
num_bits = num_bits + 1;
num_bits = num_bits + ceil( log2( n_rb_dl*(n_rb_dl+1)/2 ) );
num_bits = num_bits + 6;
num_bits = num_bits + 4;

num_bits = num_bits + 5;
if duplex_mode == 1
    num_bits = num_bits + 4;
else
    num_bits = num_bits + 3;
end
num_bits = num_bits + 1;
num_bits = num_bits + 2;
num_bits = num_bits + 2;
if duplex_mode == 1
    num_bits = num_bits + 2;
end
num_bits = num_bits + 0;

if prod(num_bits - abgs_bits) == 0
    num_bits = num_bits + 1;
end
