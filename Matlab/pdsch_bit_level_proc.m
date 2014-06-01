function [bits, blkcrc] = pdsch_bit_level_proc( llr, tbsize)
min_val = 1;
llr(llr>min_val) = min_val;
llr(llr<-min_val) = -min_val;

% blkcrc = 0;

currentCBS = length(llr);
niter_turbocode = 6;
L_total = tbsize+24;

punct_pattern = LTE_intlv_punct(1:(3*L_total + 12), currentCBS);
% alphaLTE = intlvLTE(L_total);
% cl = 4; cg = [13 15]; fg = 13;
% trellis_dec = ctrellis_gen(cl, cg, fg);

LLR_eq_e = LTE_MCS_depunct_LLRcombine(llr, punct_pattern, L_total);

% hTDec  = comm.TurboDecoder('TrellisStructure', poly2trellis(4, ...
%              [13 15 ], 13), 'InterleaverIndices', alphaLTE.', ...
%              'NumIterations', 1, 'Algorithm', 'Max*');
% bits  = step(hTDec, LLR_eq_e.');

% LLR_eq_e = LLR_eq_e( 1 : (end-12) );
% rec_s = TurboDePunct_x(LLR_eq_e, alphaLTE); % demultiplex to get input for decoder 1 and 2
% bits = Turbo13Dec_x(rec_s, trellis_dec, alphaLTE, niter_turbocode);

bits = lte_turbo_decoder(LLR_eq_e, niter_turbocode);

crc_bits = lte_calc_crc(bits(1:tbsize), '24A');
blkcrc = sum(xor(crc_bits, bits(tbsize+1:end)));
