function peak_out=decode_pcfich(peak,tfg)

nRB = 100;
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
np = zeros(1, n_ports);
for i=1:n_ports
    [ce_tfg(:,:,i), np(i)]=chan_est(peak, tfg, i, nRB);
end

[pcfich_sym, pcfich_ce] = pcfich_extract(peak, tfg, ce_tfg);

if (n_ports==1)

  gain=conj(pcfich_ce(1,:))./absx2(pcfich_ce(1,:));
  syms=pcfich_sym.*gain;

elseif (n_ports==2)

  syms=NaN(1,length(pcfich_sym));

  for t=1:2:length(syms)
    % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
    h1=mean(pcfich_ce(1,t:t+1));
    h2=mean(pcfich_ce(2,t:t+1));
    x1=pcfich_sym(t);
    x2=pcfich_sym(t+1);
    scale=sum(absx2([h1 h2]));
    syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
    syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);

  end
  % 3dB factor comes from precoding for transmit diversity
  syms=syms*sqrt(2);

elseif (n_ports==4)

  syms=NaN(1,length(pcfich_sym));

  for t=1:2:length(syms)
    % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
    if (mod(t,4)==1)
      h1=mean(pcfich_ce(1,t:t+1));
      h2=mean(pcfich_ce(3,t:t+1));
    else
      h1=mean(pcfich_ce(2,t:t+1));
      h2=mean(pcfich_ce(4,t:t+1));
    end
    x1=pcfich_sym(t);
    x2=pcfich_sym(t+1);
    scale=sum(absx2([h1 h2]));
    syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
    syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);

  end
  % 3dB factor comes from precoding for transmit diversity
  syms=syms*sqrt(2);

else
  error('Check code...');
end
