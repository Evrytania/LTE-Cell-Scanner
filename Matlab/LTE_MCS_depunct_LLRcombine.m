function outbit = LTE_MCS_depunct_LLRcombine(inbit, punct_pattern, L_total)
if punct_pattern ~= -1
    len_inbit = length(inbit);
    numCbit_unpunct = L_total*3 + 12;
    len_outbit = numCbit_unpunct;
    outbit = zeros(1, numCbit_unpunct);
    if len_inbit <= len_outbit
        outbit(punct_pattern) = inbit;
    else
%         tmp_val = zeros(1, length(punct_pattern) );
        for idx = 1 : len_inbit
%             disp(num2str([outbit(punct_pattern(idx)) inbit(idx)]));
%             tmp_val(idx) = sign(outbit(punct_pattern(idx))*inbit(idx));
            outbit(punct_pattern(idx)) = outbit(punct_pattern(idx)) + inbit(idx);
        end
%         plot(tmp_val);
    end
else
    outbit = inbit;
end
