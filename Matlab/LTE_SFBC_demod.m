function [syms, np] = LTE_SFBC_demod(sym_in, ce_in, np_ce, n_ports, mmse_flag)

np_mean = mean(np_ce);
syms=NaN(1,length(sym_in));
np=NaN(1,length(sym_in));

if mmse_flag == 0 % ZF
    if (n_ports==1)

        gain=conj(ce_in(1,:))./absx2(ce_in(1,:));
        syms=sym_in.*gain;
        np=np_mean*absx2(gain);

    elseif (n_ports==2)

        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            h1=mean(ce_in(1,t:t+1));
            h2=mean(ce_in(2,t:t+1));
            x1=sym_in(t);
            x2=sym_in(t+1);
            scale=sum(absx2([h1 h2]));
            syms(t)=(conj(h1)*x1+h2*conj(x2))/scale;
            syms(t+1)=conj((-conj(h2)*x1+h1*conj(x2))/scale);
            np(t)=(abs(h1)/scale)^2*np_mean+(abs(h2)/scale)^2*np_mean;
            np(t+1)=np(t);
        end
        % 3dB factor comes from precoding for transmit diversity
        syms=syms*sqrt(2);

    elseif (n_ports==4)

        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            if (mod(t,4)==1)
              h1=mean(ce_in(1,t:t+1));
              h2=mean(ce_in(3,t:t+1));
            else
              h1=mean(ce_in(2,t:t+1));
              h2=mean(ce_in(4,t:t+1));
            end
            x1=sym_in(t);
            x2=sym_in(t+1);
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

else % MMSE
    if (n_ports==1)

        h = ce_in(1,:);
        hconj = conj(h);
        habs2 = real(ce_in(1,:).*conj(ce_in(1,:)));
        syms = hconj.*sym_in./(habs2 + np_mean);

    elseif (n_ports==2)

        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            h1=mean(ce_in(1,t:t+1));
            h2=mean(ce_in(2,t:t+1));
            r1=sym_in(t);
            r2=sym_in(t+1);
            r = [r1; conj(r2)];

            h = [h1, -h2; conj(h2), conj(h1)]./sqrt(2);
            hconj = h';
            habs2 = h*h';
            x = hconj*inv(habs2 + diag(ones(1,2).*np_mean))*r;

            syms(t) = x(1);
            syms(t+1) = conj( x(2) );
        end

    elseif (n_ports==4)

        for t=1:2:length(syms)
            % http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
            if (mod(t,4)==1)
              h1=mean(ce_in(1,t:t+1));
              h2=mean(ce_in(3,t:t+1));
            else
              h1=mean(ce_in(2,t:t+1));
              h2=mean(ce_in(4,t:t+1));
            end
            r1=sym_in(t);
            r2=sym_in(t+1);
            r = [r1; conj(r2)];

            h = [h1, -h2; conj(h2), conj(h1)]./sqrt(2);
            hconj = h';
            habs2 = h*h';
            x = hconj*inv(habs2 + diag(ones(1,2).*np_mean))*r;
            
            syms(t) = x(1);
            syms(t+1) = conj( x(2) );
        end

    else
        error('Check code...');
    end
end
