function gather_format_signal_file_for_future
clear all; close all;

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-1850-1880MHz/f1860_s1.92_g0_1s_strong.bin', 'f1860_s1.92_g0_1s_strong_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-1850-1880MHz/f1860_s1.92_g0_1s.bin', 'f1860_s1.92_g0_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-1880-1900MHz/f1890_s1.92_g0_1s.bin',  'f1890_s1.92_g0_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2555-2575MHz/f2564.9_s1.92_g0_0.8s.bin',  'f2564.9_s1.92_g0_0.8s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2555-2575MHz/f2564.9_s1.92_g0_1s.bin',  'f2564.9_s1.92_g0_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2555-2575MHz/f2565_s1.92_g0_1s.bin',  'f2565_s1.92_g0_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2575-2595MHz/f2585_s1.92_g0_0.8s.bin',  'f2585_s1.92_g0_0.8s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2575-2595MHz/f2585_s1.92_g0_1s.bin',  'f2585_s1.92_g0_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2575-2595MHz/f2584.9_s1.92_g0_1s.bin', 'f2584.9_s1.92_g0_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2595-2615MHz/f2604.9_s1.92_g0_0.8s.bin', 'f2604.9_s1.92_g0_0.8s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2595-2615MHz/f2604.9_s1.92_g0_1s.bin', 'f2604.9_s1.92_g0_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2595-2615MHz/f2605_s1.92_g0_1s.bin', 'f2605_s1.92_g0_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_1s.bin', 'f2645_s1.92_g0_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz/f2645_s1.92_g0_0.8s.bin', 'f2645_s1.92_g0_0.8s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz/f2644.9_s1.92_g0_1s.bin', 'f2644.9_s1.92_g0_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g0_SamplingPPM26.2_1s.bin', 'f2645_s1.92_g0_SamplingPPM26.2_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz-know-PPM/f2645_s1.92_g20_SamplingPPM26.2_1s.bin', 'f2645_s1.92_g20_SamplingPPM26.2_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz-know-PPM/f2644.9_s1.92_g0_SamplingPPM26.2_1s.bin', 'f2644.9_s1.92_g0_SamplingPPM26.2_1s_rtlsdr.bin');
rtlsdr_remove_header_and_save_as('~/git/rtl-sdr-LTE/scan-capture/frequency-2635-2655MHz-know-PPM/f2644.9_s1.92_g20_SamplingPPM26.2_1s.bin', 'f2644.9_s1.92_g20_SamplingPPM26.2_1s_rtlsdr.bin');

rtlsdr_remove_header_and_save_as('../test/f1860_s1.92_g0_1s.bin', 'f1860_s1.92_g0_1s_rtlsdr1.bin');
rtlsdr_remove_header_and_save_as('../test/f1815_s1.92_from_dimitri.bin', 'f1815_s1.92_from_dimitri_rtlsdr.bin');

function [fc_requested, fc_programmed, fs_requested, fs_programmed] = read_header_from_bin(bin_filename)

fc_requested = inf;
fc_programmed = inf;
fs_requested = inf;
fs_programmed = inf;

fid = fopen(bin_filename);

if fid==-1
    disp('read_header_from_bin: Can not open file for read!');
    return;
end

magic1 = fread(fid, 1, 'double');
tmp1 = fread(fid, 1, 'uint64');

magic2 = fread(fid, 1, 'double');
tmp2 = fread(fid, 1, 'uint64');

magic3 = fread(fid, 1, 'double');
tmp3 = fread(fid, 1, 'uint64');

magic4 = fread(fid, 1, 'double');
tmp4 = fread(fid, 1, 'uint64');

magic5 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic6 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic7 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic8 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

fclose(fid);

fc_requested_magic = 73492.215;
fc_programmed_magic = -0.7923597;
fs_requested_magic = -189978508;
fs_programmed_magic = 93.126712;

reserve1_magic = -53243.129;
reserve2_magic = 0.0008123898;
reserve3_magic = -6.0098321;
reserve4_magic = 237.09983;

if magic1 == fc_requested_magic && ...
   magic2 == fc_programmed_magic && ...
   magic3 == fs_requested_magic && ...
   magic4 == fs_programmed_magic && ...
   magic5 == reserve1_magic && ...
   magic6 == reserve2_magic && ...
   magic7 == reserve3_magic && ...
   magic8 == reserve4_magic

    fc_requested = tmp1;
    fc_programmed = tmp2;
    fs_requested = tmp3;
    fs_programmed = tmp4;
end

function [fc_requested, fc_programmed, fs_requested, fs_programmed] = read_header_from_bin_fid(fid)

fc_requested = inf;
fc_programmed = inf;
fs_requested = inf;
fs_programmed = inf;

% fid = fopen(bin_filename);
% 
% if fid==-1
%     disp('read_header_from_bin: Can not open file for read!');
%     return;
% end

magic1 = fread(fid, 1, 'double');
tmp1 = fread(fid, 1, 'uint64');

magic2 = fread(fid, 1, 'double');
tmp2 = fread(fid, 1, 'uint64');

magic3 = fread(fid, 1, 'double');
tmp3 = fread(fid, 1, 'uint64');

magic4 = fread(fid, 1, 'double');
tmp4 = fread(fid, 1, 'uint64');

magic5 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic6 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic7 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

magic8 = fread(fid, 1, 'double');
fread(fid, 1, 'uint64');

% fclose(fid);

fc_requested_magic = 73492.215;
fc_programmed_magic = -0.7923597;
fs_requested_magic = -189978508;
fs_programmed_magic = 93.126712;

reserve1_magic = -53243.129;
reserve2_magic = 0.0008123898;
reserve3_magic = -6.0098321;
reserve4_magic = 237.09983;

if magic1 == fc_requested_magic && ...
   magic2 == fc_programmed_magic && ...
   magic3 == fs_requested_magic && ...
   magic4 == fs_programmed_magic && ...
   magic5 == reserve1_magic && ...
   magic6 == reserve2_magic && ...
   magic7 == reserve3_magic && ...
   magic8 == reserve4_magic

    fc_requested = tmp1;
    fc_programmed = tmp2;
    fs_requested = tmp3;
    fs_programmed = tmp4;
end

function s = rtlsdr_remove_header_and_save_as(rtl_sdr_bin_filename, output_filename)
disp([rtl_sdr_bin_filename ' --> ' output_filename]);

header_exist = false;

[fc_requested_exist, fc_programmed_exist, fs_requested_exist, fs_programmed_exist] = read_header_from_bin(rtl_sdr_bin_filename);

if fc_requested_exist~=inf
    header_exist = true;
    disp('There is already a header!');
    disp(['fc_requested ' num2str(fc_requested_exist) ' fc_programmed ' num2str(fc_programmed_exist)  ' fs_requested ' num2str(fs_requested_exist)  ' fs_programmed ' num2str(fs_programmed_exist) ]);
end

fid = fopen(rtl_sdr_bin_filename);

if fid==-1
    disp('rtlsdr_remove_header_and_save_as: Can not open input file!');
    return;
end

if header_exist
    read_header_from_bin_fid(fid);
end

[s, count] = fread(fid, inf, 'uint8');
fclose(fid);

plot(s);

fid = fopen(output_filename, 'w');

if fid==-1
    disp('rtlsdr_remove_header_and_save_as: Can not open output file!');
    return;
end

fwrite(fid, s, 'uint8');
fclose(fid);
