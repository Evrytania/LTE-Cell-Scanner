function [ant_gain_new, lna_gain_new, vga_gain_new] = hackrf_gain_regulation(ant_gain, lna_gain, vga_gain)
% ant_gain = 0; % 0 turn off, 1 turn on
% lna_gain = 20; %0-40dB, 8dB steps
% vga_gain = 20; %0-62dB, 2dB steps

if ant_gain > 0
    ant_gain_new = 1;
else
    ant_gain_new = 0;
end

lna_gain_new = round(lna_gain/8);
lna_gain_new = lna_gain_new*8;
if lna_gain_new > 40
    lna_gain_new = 40;
elseif lna_gain_new < 0
    lna_gain_new = 0;
end

vga_gain_new = round(vga_gain/2);
vga_gain_new = vga_gain_new*2;
if vga_gain_new > 62
    vga_gain_new = 62;
elseif vga_gain_new < 0
    vga_gain_new = 0;
end
