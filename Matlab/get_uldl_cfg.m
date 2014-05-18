function peak_out = get_uldl_cfg(peak, pcfich_info)

if peak.duplex_mode == 0
    peak_out.uldl_cfg = -2; % default UL DL configuration: all DL (FDD) or invalid TDD
elseif peak.duplex_mode == 1 % TDD
    peak_out = peak;
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
    
%     if isempty(peak_out.uldl_cfg)
%         peak_out.uldl_cfg = -1;
%     end
else
    disp('duplex_mode has to be 0(FDD) or 1(TDD)');
end
