% Jiao Xianjun (putaoshu@gmail.com)
% lte_receiver.m
% HACKRF LTE Receiver

clear all;
close all;

% filename = '../test/f1860_s15.36_bw10_1s.bin';
filename = '../test/f1890_s15.36_bw10_1s.bin';

fid = fopen(filename);
a = fread(fid, inf, 'int8');
plot(a);
