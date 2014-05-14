function [pcfich_info, pcfich_corr] = decode_pcfich(peak, subframe_idx, tfg, ce_tfg)

% nRB = peak.n_rb_dl;
% nSC = nRB*12;

% Local shortcuts
% cp_type = peak.cp_type;
n_ports = peak.n_ports;

n_id_cell = peak.n_id_cell;

[pcfich_sym, pcfich_ce] = pcfich_extract(peak, tfg, ce_tfg);

num_pcfich = size(pcfich_sym, 1);

pcfich_cw = [0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1
    1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0
    1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1].';
pcfich_cw = 1-2.*pcfich_cw;

pcfich_info = zeros(1, num_pcfich);
pcfich_corr = zeros(1, num_pcfich);

subframe_idx = mod(subframe_idx, 10);
start_slot_idx = subframe_idx*2;

slot_idx = start_slot_idx;
% slot_idx = 1;
for i=1:num_pcfich
    scr=lte_pn((floor(slot_idx/2) + 1)*(2*n_id_cell+1)*(2^9) + n_id_cell, 32);
    
%     % MMSE equalizer --- seems worse than ZF! let's use ZF
%     if (n_ports==1)
%         h = pcfich_ce(i,:,1);
%         hconj = conj(h);
%         habs2 = real(pcfich_ce(i,:,1).*conj(pcfich_ce(i,:,1)));
%         syms = hconj.*pcfich_sym(i,:)./(habs2 + np_ce);
%     elseif (n_ports==2)
%         np_mean=mean([np_ce(1) np_ce(2)]);
%         syms=NaN(1,length(pcfich_sym(i,:)));
% 
%         for t=1:2:length(syms)
%             % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
%             h1=mean(pcfich_ce(i,t:t+1, 1));
%             h2=mean(pcfich_ce(i,t:t+1, 2));
%             r1=pcfich_sym(i, t);
%             r2=pcfich_sym(i, t+1);
%             r = [r1; conj(r2)];
% 
%             h = [h1, -h2; conj(h2), conj(h1)]./sqrt(2);
%             hconj = h';
%             habs2 = h*h';
%             x = hconj*(habs2 + diag(ones(1,2).*np_mean))\r;
%             
%             syms(t) = x(1);
%             syms(t+1) = conj( x(2) );
%         end
%     elseif (n_ports==4)
%         np_mean=mean([np_ce(1) np_ce(2) np_ce(3) np_ce(4)]);
%         syms=NaN(1,length(pcfich_sym(i,:)));
%         for t=1:2:length(syms)
%             % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
%             if (mod(t,4)==1)
%               h1=mean(pcfich_ce(i,t:t+1, 1));
%               h2=mean(pcfich_ce(i,t:t+1, 3));
%             else
%               h1=mean(pcfich_ce(i,t:t+1, 2));
%               h2=mean(pcfich_ce(i,t:t+1, 4));
%             end
%             r1=pcfich_sym(i, t);
%             r2=pcfich_sym(i, t+1);
%             r = [r1; conj(r2)];
% 
%             h = [h1, -h2; conj(h2), conj(h1)]./sqrt(2);
%             hconj = h';
%             habs2 = h*h';
%             x = hconj*(habs2 + diag(ones(1,2).*np_mean))\r;
%             
%             syms(t) = x(1);
%             syms(t+1) = conj( x(2) );
%         end
%     else
%         error('Check code...');
%     end
    
    if (n_ports==1)
        syms=pcfich_sym(i,:)./pcfich_ce(i,:,1);
    elseif (n_ports==2)
        syms=NaN(1,length(pcfich_sym(i,:)));
        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            h1=mean(pcfich_ce(i,t:t+1, 1));
            h2=mean(pcfich_ce(i,t:t+1, 2));
            x1=pcfich_sym(i, t);
            x2=pcfich_sym(i, t+1);
            scale=sum(absx2([h1 h2]));
            syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
            syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
        end
    elseif (n_ports==4)
        syms=NaN(1,length(pcfich_sym(i,:)));
        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            if (mod(t,4)==1)
              h1=mean(pcfich_ce(i,t:t+1, 1));
              h2=mean(pcfich_ce(i,t:t+1, 3));
            else
              h1=mean(pcfich_ce(i,t:t+1, 2));
              h2=mean(pcfich_ce(i,t:t+1, 4));
            end
            x1=pcfich_sym(i, t);
            x2=pcfich_sym(i, t+1);
            scale=sum(absx2([h1 h2]));
            syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
            syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
        end
    else
        error('Check code...');
    end
    
    % Extract signal of sequence
    e_est = zeros(1,32);
    e_est(1:2:end) = sign(real(syms));
    e_est(2:2:end) = sign(imag(syms));

    % Unscramble
    e_est(scr==1) = -1.*e_est(scr==1);
    
    corr_val = abs( e_est*pcfich_cw );
    [max_val, max_idx] = max(corr_val);
    
    pcfich_corr(i) = max_val;
    if max_val>=20
        pcfich_info(i) = max_idx;
    end
    
    slot_idx = mod(slot_idx+2, 20);
end
