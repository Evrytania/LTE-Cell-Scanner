% Measure PSS correlations in the presence of a frequency offset.
%n_id_2=0;

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

freq_off_set=linspace(-30e3,30e3,1001);
time_off_set=-3000:3000;

n_freqs=length(freq_off_set);
n_times=length(time_off_set);

pss_td=NaN(3,144+2048+0*2048);
for n_id_2=0:2
  pss_freq=pss(n_id_2);
  pss_freq=[0 pss_freq(32:end) zeros(1,2048-73+5+5) pss_freq(1:31)];
  pss_td(n_id_2+1,145:end)=[idft(pss_freq)*sqrt(2048/62) zeros(1,0*2048)];
  pss_td(n_id_2+1,1:144)=pss_td(n_id_2+1,end-143:end);
  pss_td(n_id_2+1,:)=pss_td(n_id_2+1,:)/sqrt(sigpower(pss_td(n_id_2+1,:)));
end

% Frequency domain correlations
figure(1);
log_xc_freq=NaN(3,n_freqs);
for n_id_2=0:2
  for t=1:n_freqs
    freq=freq_off_set(t);
    log_xc_freq(n_id_2+1,t)=sum(conj(pss_td(n_id_2+1,:)).*fshift(pss_td(n_id_2+1,:),freq/(fs_lte/2)))/(2048+144);
  end
end
plot(freq_off_set,db20(abs(transpose(log_xc_freq))));
%ylim([-20 5]);
drawnow;
zgo;

% Time domain auto correlations
pss_td_ext=[pss_td zeros(3,(2048+144)*2)];
log_xc_td=NaN(3,n_times);
figure(2);
for n_id_2=0:2
  for t=1:n_times
    to=time_off_set(t);
    log_xc_td(n_id_2+1,t)=sum(conj(pss_td_ext(n_id_2+1,:)).*tshift(pss_td_ext(n_id_2+1,:),to))/(2048+144);
  end
end
plot(time_off_set,db20(abs(transpose(log_xc_td))));
%ylim([-50 5]);
drawnow;
zgo;

% Time domain cross correlations
log_xc_td_cross=NaN(3,n_times);
log_xc_td_cross_max16=NaN(3,n_times);
figure(3);
xc_set=[1 2;1 3;2 3]-1;
for k=1:3
  n_id_2_1=xc_set(k,1);
  n_id_2_2=xc_set(k,2);
  for t=1:n_times
    to=time_off_set(t);
    log_xc_td_cross(k,t)=sum(conj(pss_td_ext(n_id_2_1+1,:)).*tshift(pss_td_ext(n_id_2_2+1,:),to))/(2048+144);
  end
  for t=1:n_times
    log_xc_td_cross_max16(k,t)=max(log_xc_td_cross(k,max([t-32 1]):min([t+80 n_times])));
  end
end
plot(time_off_set,db20(abs(transpose(log_xc_td_cross))), ...
time_off_set,db20(abs(transpose(log_xc_td_cross_max16))));
legend('0 1','0 2','2 3','location','se');
%ylim([-50 5]);
drawnow;
zgo;


