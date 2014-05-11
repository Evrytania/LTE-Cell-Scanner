function pcfich_info = decode_pcfich(peak,tfg)

nRB = peak.n_rb_dl;
nSC = nRB*12;

peak_out=peak;

% Local shortcuts
cp_type = peak.cp_type;
n_ports = peak.n_ports;

% Derive some values
n_ofdm=size(tfg,1);
if (strcmpi(cp_type,'normal'))
  n_symb_dl=7;
elseif (strcmpi(cp_type,'extended'))
  n_symb_dl=6;
else
 error('Check code...');
end

n_id_cell = peak.n_id_cell;

% Channel estimation
ce_tfg = NaN(n_ofdm,nSC,n_ports);
np_ce = zeros(1, n_ports);
for i=1:n_ports
    [ce_tfg(:,:,i), np_ce(i)]=chan_est(peak, tfg, i-1, nRB);
end

[pcfich_sym, pcfich_ce] = pcfich_extract(peak, tfg, ce_tfg);

num_pcfich = size(pcfich_sym, 1);

pcfich_cw = [0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1
    1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0
    1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1].';

pcfich_info = zeros(1, num_pcfich);

scr=lte_pn((floor(0/2) + 1)*(2*n_id_cell+1)*(2^9) + n_id_cell, 32);

for i=1:num_pcfich
    if (n_ports==1)

%         np=np_ce(1);
%         gain=conj(pcfich_ce(i,:,1))./absx2(pcfich_ce(i,:,1));
%         syms=pcfich_sym(i,:).*gain;
%         np=np*absx2(gain);

        syms=pcfich_sym(i,:)./pcfich_ce(i,:,1);
%         np = 0.001;

    elseif (n_ports==2)

        np_mean=mean([np_ce(1) np_ce(2)]);
        syms=NaN(1,length(pcfich_sym(i,:)));
        np=NaN(1,length(pcfich_sym(i,:)));

        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            h1=mean(pcfich_ce(i,t:t+1, 1));
            h2=mean(pcfich_ce(i,t:t+1, 2));
            x1=pcfich_sym(i, t);
            x2=pcfich_sym(i, t+1);
            scale=sum(absx2([h1 h2]));
            syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
            syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
            np(t)=(abs(h1)/scale)^2*np_mean+(abs(h2)/scale)^2*np_mean;
            np(t+1)=np(t);
        end
        % 3dB factor comes from precoding for transmit diversity
        syms=syms*sqrt(2);

    elseif (n_ports==4)

        np_mean=mean([np_ce(1) np_ce(2) np_ce(3) np_ce(4)]);
        syms=NaN(1,length(pcfich_sym(i,:)));
        np=NaN(1,length(pcfich_sym(i,:)));

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
            np(t)=(abs(h1)/scale)^2*np_mean+(abs(h2)/scale)^2*np_mean;
            np(t+1)=np(t);
        end
        % 3dB factor comes from precoding for transmit diversity
        syms=syms*sqrt(2);

    else
        error('Check code...');
    end
    
    % Extract the bits
%     e_est=deqam(syms,np,'QAM','LTE');
    
    a = zeros(1,32);
    a(1:2:end) = (1-sign(real(syms)))./2;
    a(2:2:end) = (1-sign(imag(syms)))./2;
    e_est = a;
    
    % Unscramble
    e_est(scr==1)=1-e_est(scr==1);
    
    corr_val = abs( ((2.*e_est)-1)*((2.*pcfich_cw)-1) );
    [max_val, max_idx] = max(corr_val);
    if max_val>24
        pcfich_info(i) = max_idx;
    else
        pcfich_info(i) = -max_val;
    end
end
