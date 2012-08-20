% Full PSS search algorithm including searching for various frequency offsets.
%
% Final, cleaned up version.

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Design considerations for this algorithm:
%
% Algorithm must work in a capture-then-process mode. All data will be
%   captured, and then analyzed to determine cell information.
%
% Want this algorithm to operate in low SNR regions.
%   Therefore, want to capture the full BCH which spans 4 frames (40ms).
%   To guarantee that a full BCH is captured, 80ms of data must be captured.
%
% Use all available data to detect PSS/SSS.
%   From frame to frame, different antennas can be used to transmit PSS.
%     Coherent combining not possible, only incohrent combining.
%   SSS can be coherently combined based on channel estimates received
%     from the nearby PSS.

n_dft=128;
n_rb=6;
%cp_type='normal';
%cp_type='extended';
n_subframes=10*8;
%slot_start=0;
%n_id_1=18;
%n_id_2=0;
%scenario='UMa';
scenario='flat';
%scenario='EVA';
fc=739e6;
% Noise added in the time domain
% SNR over the synchrnization signals is db10(1/rx_noise_power)
rx_noise_power=udb10(-6.2+10-300);
% In Hz
n_trials=1;
thresh1_n_nines=12;
%thresh1_n_nines=6;
thresh2_n_sigma=3.0;
ds_comb_arm=2;
%f_search_set=[-10000 -5000 0 5000 10000];
%f_search_set=20e3:5e3:60e3;
f_search_set=-35e3:5e3:35e3;
%f_search_set=-5e3:5e3:5e3;
%sig_source='cpp';
%sig_source='recorded';
sig_source='generated';

n_f=length(f_search_set);

% Create the time domain pss signals
pss_td=NaN(3,128+9);
for t=0:2
  temp=pss(t);
  %temp=temp(6:67);
  temp=[0 temp(32:end) zeros(1,65) temp(1:31)];
  temp_td=idft(temp)*sqrt(128/62);
  pss_td(t+1,:)=[temp_td(end-8:end) temp_td];
end

log_success=NaN(1,n_trials);
log_success_mib=NaN(1,n_trials);
log_thresh1_fail=NaN(1,n_trials);
log_thresh2_fail=NaN(1,n_trials);
log_false_detection=NaN(1,n_trials);
log_foe=NaN(1,n_trials);
state=etc;
%if (strcmp(sig_source,'cpp'))
%  itload('log_capbuf.it');
%end
for nt=1:n_trials
  if (rand<.5)
    cp_type='normal';
  else
    cp_type='extended';
  end
  slot_start=floor(rand*20);
  n_id_1=floor(rand*168);
  n_id_2=floor(rand*3);
  freq_offset=0000;
  load_factor=max([abs(randn/9) .5]);
  if (rand>.5)
    load_factor=1-load_factor;
  end
  %load_factor=0

  %cp_type='normal'
  %slot_start=0;
  %n_id_1=18
  %n_id_2=0
  %freq_offset=2000;
  %load_factor=0

  % Derive some values
  n_sc=n_rb*12;
  if strcmpi(cp_type,'normal')
    n_ofdm=7;
  elseif strcmpi(cp_type,'extended')
    n_ofdm=6;
  else
    error('Check code...');
  end

  if (nt~=1)
    disp(sprintf('Success rate (excluding MIB): %.2f%%',sum(log_success(1:nt-1))/(nt-1)*100));
    disp(sprintf('Success rate (including MIB): %.2f%%',sum(log_success_mib(1:nt-1))/(nt-1)*100));
    disp(sprintf('Thresh1 rejection rate: %.2f%%',sum(log_thresh1_fail(1:nt-1))/(nt-1)*100));
    disp(sprintf('Thresh2 rejection rate: %.2f%%',sum(log_thresh2_fail(1:nt-1))/(nt-1-sum(log_thresh1_fail(1:nt-1)))*100));
    disp(sprintf('False detection rate: %.2f%%',sum(log_false_detection(1:nt-1))/(nt-1)*100));
  end
  log_success(nt)=0;
  log_thresh1_fail(nt)=0;
  log_thresh2_fail(nt)=0;
  log_false_detection(nt)=0;

  % Create the wideband transmitted signal
  if (strcmpi(sig_source,'generated'))
    sig_tx=create_dl_sig(cp_type,n_subframes,slot_start,n_id_1,n_id_2,load_factor);
    assert(length(sig_tx)==n_subframes*.001*(fs_lte/16));
  elseif (strcmpi(sig_source,'recorded'))
    %data_read=transpose(read_complex_binary('LTE_data_cap.dat'));
    %data_read=rtl_cap_read('LTE_cap_commandline_sync.dat');
    %data_read=rtl_cap_read('LTE_739.000.dat');
    data_read=rtl_cap_read('LTE_739.000_twocells.dat');
    %data_read=conj(data_read);
    %sdata_read=complex(imag(data_read),-real(data_read));
    sig_tx=data_read(3e6:3e6+n_subframes*.001*(fs_lte/16)-1);
    clear data_read;
  elseif (strcmpi(sig_source,'cpp'))
    itload('capbuf_0000.it');
    sig_tx=capbuf;
  else
    error('Check code...');
  end

  % Add the channel
  [h h_timestamp delays]=channel_gen(0,1,1,0,fc,100,0,scenario);
  % Artificial channel with very wide delay spread
  %h=[rayleigh(1,.5,2) rayleigh(1,.5,2)];
  %h=repmat(rayleigh(1,.5,2),1,2);
  %h=[1 exp(j*rand*2*pi)];
  %h=shiftdim(h,-2);
  %delays=[0 3*16*(1/fs_lte)];
  % Normalize to 0dB power on the occupied subcarriers.
  H=channel2tf(h,delays,[-31:-1 1:31]*15e3);
  H=shiftdim(H,1);
  h=h/sqrt(sigpower(H(:)));
  % Remove any time offset that might be present in the frequency
  % domain.
  to=angle(sum([conj(H(1:30)).*H(2:31) conj(H(32:61)).*H(33:62)]));
  delays=delays+to/(2*pi)*2048*1/fs_lte;
  H=channel2tf(h,delays,[-31:-1 1:31]*15e3);
  H=shiftdim(H,1);
  %H_temp=[0 H(32:end) zeros(1,500), H(1:31)];
  %figure(3);
  %plot(tshift(abs(idft(H_temp)),10));
  %xlim([0 20]);
  %angle(sum([conj(H(1:30)).*H(2:31) conj(H(32:61)).*H(33:62)]))
  chan_out=zeros(1,length(sig_tx));
  for t=1:length(delays)
    chan_out=chan_out+h(1,1,t)*tshift(sig_tx,delays(t)/(1/(fs_lte/16)));
  end

  % Add the frequency offset
  sig_rx=fshift(chan_out,freq_offset,fs_lte/16);

  % Add the receiver noise
  sig_rx=sig_rx+sqrt(rx_noise_power)*blnoise(length(sig_rx));

  % Filter out the portion that the UE is looking at.
  % No need. The interpolated time domain PSS sequence automatically filters
  % out energy outside of the bad of interest.
  ue_rx=sig_rx;
  rx_cutoff=(n_rb*12*15e3/2+4*15e3)/(fs_lte/16/2);
  %disp('cyc_filt bypassed');
  %ue_rx=cyc_filt(ue_rx,@(f)(abs(f)<=rx_cutoff));

  n_samp=length(ue_rx);

  % Extract samples for the capture buffer
  %n_cap=(2*(2048+512)+307200*7+4*(2048+512))/16;
  n_cap=307200*8/16;
  capbuf=ue_rx(1:n_cap);
  %disp('look here!');
  %capbuf=ue_rx_sav(1:n_cap);

  % Correlate against the PSS
  [xc_incoherent_collapsed_pow xc_incoherent_collapsed_frq n_comb_xc n_comb_sp xc_incoherent_single xc_incoherent sp_incoherent xc sp]=xcorr_pss(capbuf,f_search_set,ds_comb_arm,fc);

  % The following analysis is accurate when ds_comb_arm is 0. It loses
  % accuracy as ds_comb_arm increases because adjacent despread values are
  % correlated in the time domain. Hence, if ds_comb_arm is > 0, the thresholds
  % must be set higher than what is predicted below.
  %
  % The received signal had a power level of sp_incoherent in a bandwidth of
  % rx_cutoff. If the cutoff frequency was actually 1, the noise power would
  % be:
  % sp_incoherent/rx_cutoff.
  % Assuming that the received signal is actually just gaussian noise with a
  % power level of sp_incoherent/rx_cutoff, after despreading, we would
  % have a gaussian noise signal Y with a power level of:
  % sp_incoherent/(rx_cutoff*137)
  % In other words:
  % E[Y]=0
  % E[Y*conj(Y)]=E[abs(Y)^2]=sp_incoherent/(rx_cutoff*137)
  % Y=xc(1,:)
  %
  % abs(Y)^2 is composed of two uncorrelated components: Y_real^2+Y_imag^2
  % Each component is gaussian:
  % E[Y_real]=E[Y_imag]=0
  % Var[Y_real]=Var[Y_imag]=sp_incoherent/(rx_cutoff*137*2)
  % Var[abs(Y)^2]=(sp_incoherent/(rx_cutoff*137))^2
  %
  % Let Q=Y/sqrt(sp_incoherent/(rx_cutoff*137*2));
  % E[Q]=0
  % E[Q*conj(Q)]=2;
  % Var[Q_real]=Var[Q_imag]=1;
  %
  % Suppose R=abs(Q1)^2+abs(Q2)^2+...abs(QN)^2:
  % Then, R is chi-squared distributed with 2N degrees of freedom.
  R99=chi2inv(.99,2*n_comb_xc*(2*ds_comb_arm+1));
  R999=chi2inv(.999,2*n_comb_xc*(2*ds_comb_arm+1));
  R9999=chi2inv(.9999,2*n_comb_xc*(2*ds_comb_arm+1));
  R99999=chi2inv(.99999,2*n_comb_xc*(2*ds_comb_arm+1));
  R999999=chi2inv(.999999,2*n_comb_xc*(2*ds_comb_arm+1));
  R_mean=2*n_comb_xc*(2*ds_comb_arm+1);
  % Suppose Z=(abs(Y1)^2+abs(Y2)^2+...abs(YN)^2)/N:
  % Z=R*(sp_incoherent/rx_cutoff/137/2)/N
  Z99=R99*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  Z999=R999*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  Z9999=R9999*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  Z99999=R99999*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  Z999999=R999999*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  Z_mean=R_mean*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);

  % Find the candidate peaks
  R_th1=chi2inv(1-10^-thresh1_n_nines,2*n_comb_xc*(2*ds_comb_arm+1));
  Z_th1=R_th1*sp_incoherent/rx_cutoff/137/2/n_comb_xc/(2*ds_comb_arm+1);
  %Z_th1=Z_mean;
  peaks=peak_search(xc_incoherent_collapsed_pow,xc_incoherent_collapsed_frq,Z_th1,f_search_set);
  %peaks(1).ind=829;

    %1:9600,db10(despread_power_expected), ...
    %1:9600,db10(despread_power_expected+despread_power_sigma*3), ...
    %1:9600,db10(thresh), ...
  figure(1);
  plot( ...
    1:9600,db10([transpose(xc_incoherent_collapsed_pow) transpose(sp_incoherent)]), ...
    1:9600,db10(Z_mean), ...
    1:9600,db10(Z99), ...
    1:9600,db10(Z999), ...
    1:9600,db10(Z9999), ...
    1:9600,db10(Z99999), ...
    1:9600,db10(Z999999), ...
    1:9600,db10(Z_th1),'--', ...
    NaN,NaN ...
  );
  zgo;
  %axis([818 830 -11 0]);
  %%pause
  %drawnow;

  %%%%%%
  % Search through the detected PSS peaks
  %%%%%%
  %for pk=1:length(peaks)
  disp('full peak display disabled...');
  for pk=1:1
    peak_out=sss_detect(peaks(pk),capbuf,thresh2_n_sigma,fc);
    peaks(pk)=peak_out;
    if (isnan(peaks(pk).n_id_1))
      disp('Thresh2 fail');
      log_thresh2_fail(nt)=1;
      continue;
    end
    %peaks(1).frame_start=1

    peak_out=pss_sss_foe(peaks(pk),capbuf,fc);
    peaks(pk)=peak_out;

    disp(sprintf('%.2f dB',db10(peaks(pk).pow)));
    disp(sprintf('n_id_1/2 = %i/%i',peaks(pk).n_id_1,peaks(pk).n_id_2));
    disp(sprintf('n_id_cell = %i',3*peaks(pk).n_id_1+peaks(pk).n_id_2));
    disp(peaks(pk).cp_type)
    disp(sprintf('%i',peaks(pk).frame_start));
    disp(sprintf('%.3f Hz',peaks(pk).freq_fine));

    % Create the time/frequency grid
    [tfg tfg_timestamp]=extract_tfg(peaks(pk),capbuf,fc);

    % Perform FOE-C and TOE-C
    [tfg_comp tfg_comp_timestamp peaks_out]=tfoec(peaks(pk),tfg,tfg_timestamp,fc);
    peaks(pk)=peaks_out;

    % Decode the MIB!!! Finally!!!
    peaks_out=decode_mib(peaks(pk),tfg_comp);
    peaks(pk)=peaks_out;
  end

  %log_success(nt)=0;
  %for pk=1:length(peaks)
  %  if ( ...
  %    (length(peaks)>=1)&& ...
  %    (peaks(pk).n_id_1==n_id_1)&& ...
  %    (peaks(pk).n_id_2==n_id_2)&& ...
  %    ((abs(peaks(pk).frame_start+9600*2-(wrap(-slot_start*960-2,0,9600*2)+1))<=5)|| ...
  %    (abs(peaks(pk).frame_start-(wrap(-slot_start*960-2,0,9600*2)+1))<=5))&& ...
  %    strcmpi(peaks(pk).cp_type,cp_type) ...
  %  )
  %    log_success(nt)=1;
  %    log_foe(nt)=peaks(pk).freq_fine;
  %    break;
  %  end
  %end
  %mean(log_foe(~isnan(log_foe)))
  %sqrt(var(log_foe(~isnan(log_foe))))
  %sum(log_success(1:nt))/nt
  %log_false_detection(nt)=1-log_success(nt);
  %if (~log_success(nt))
  %  return
  %end
  % Checking to see whethere peak was properly identified in captured
  % data.
  log_success(nt)=0;
  log_success_mib(nt)=0;
  for pk=1:length(peaks)
    if ( ...
      (peaks(pk).n_id_1==92)&& ...
      (peaks(pk).n_id_2==1)&& ...
      (abs(peaks(pk).frame_start-17449.5250338295)<4.5)&& ...
      strcmpi(peaks(pk).cp_type,'normal') ...
    )
      peaks(pk)
      log_success(nt)=1;
      log_foe(nt)=peaks(pk).freq_fine;
      if ( ...
        (peaks(pk).n_rb_dl==50)&& ...
        (peaks(pk).sfn==649) ...
      )
        log_success_mib(nt)=1;
      end
      break;
    end
  end
  %sum(log_success(1:nt)/nt)
  sum(log_success_mib(1:nt)/nt)

  state=etc(state,nt/n_trials);
end

