function peak_out = get_uldl_cfg(peak, pcfich_info)

peak_out = peak;

if peak.duplex_mode == 0
    peak_out.uldl_cfg = -2; % default UL DL configuration: all DL (FDD) or invalid TDD
    peak_out.max_num_HARQ = 8;
elseif peak.duplex_mode == 1 % TDD
    uldl_table = [ ...
        1 1 0 0 0 1 1 0 0 0; ...
        1 1 0 0 1 1 1 0 0 1; ...
        1 1 0 1 1 1 1 0 1 1; ...
        1 1 0 0 0 1 1 1 1 1; ...
        1 1 0 0 1 1 1 1 1 1;
        1 1 0 1 1 1 1 1 1 1;
        1 1 0 0 0 1 1 0 0 1
        ];
    pcfich_info = sign(pcfich_info);
    
    diff_flag = sum(abs(kron(ones(7,1), pcfich_info) - uldl_table), 2);
    
    [~, min_idx] = min(diff_flag);
    peak_out.uldl_cfg = min_idx - 1;
    
    tdd_max_num_HARQ_table = [4 7 10 9 12 15 6];
    peak_out.max_num_HARQ = tdd_max_num_HARQ_table(min_idx);
%     if isempty(peak_out.uldl_cfg)
%         peak_out.uldl_cfg = -1;
%     end
else
    disp('duplex_mode has to be 0(FDD) or 1(TDD)');
end
