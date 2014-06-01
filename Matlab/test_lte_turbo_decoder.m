% function test_lte_turbo_decoder
clear all;
close all;

bits = round(rand(1,6144));
bits_encoded = lte_turbo_encoder(bits);

bits_encoded(bits_encoded==0) = -1;

bits_rx = lte_turbo_decoder(bits_encoded, 6);

plot(bits - bits_rx);
