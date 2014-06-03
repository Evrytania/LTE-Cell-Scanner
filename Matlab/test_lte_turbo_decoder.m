% function test_lte_turbo_decoder
clear all;
close all;

code_length = 6144;

EbN0 = [0.6];
snr = EbN0 - 10.*log10((3*code_length+12)/code_length) + 10.*log10(2);

ber = zeros(1, length(EbN0));
bler = zeros(1, length(EbN0));
for i = 1 : length(EbN0)
    snr_tmp = snr(i);
    sigma2 = 10^(-snr_tmp/10);
    
    blk_err = 0;
    bit_err = 0;
    for j = 1 : 10000
        bits = round(rand(1,code_length));
        bits_encoded = lte_turbo_encoder(bits);

        bits_encoded(bits_encoded==0) = -1;

        noise_vec = sqrt(sigma2).*randn(1, length(bits_encoded));
        bits_rx = lte_turbo_decoder(bits_encoded + noise_vec, 6);
        
        bit_err_tmp = sum(xor(bits_rx, bits));
        blk_err = blk_err + sign(bit_err_tmp);
        bit_err = bit_err + bit_err_tmp;
        
        if blk_err > 20
            break;
        end
        
        if mod(j, 20) == 0
            disp(num2str([j, EbN0(i), bit_err/(j*code_length), blk_err/j]));
        end
    end
    ber(i) = bit_err/(j*code_length);
    bler(i) = blk_err/j;
    disp(['Eb/N0 ' num2str(EbN0(i)) 'dB ber ' num2str(ber(i)) ' bler ' num2str(bler(i))]);
end

semilogy(EbN0, ber); hold on;
semilogy(EbN0, bler, 'r'); grid on;

