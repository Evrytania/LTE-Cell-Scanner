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

%n_trials=100*50;
n_trials=200;

mib=bitstream(40);
mib_ce=lte_conv_encode(mib);
mib_rm=lte_conv_ratematch(mib_ce,1920);
mib_sym=qam(mib_rm,'QAM','LTE');

%np_set_db=linspace(8,20,13);
np_set_db=[14];
n_np=length(np_set_db);

log_errors=NaN(1,n_np);
log_mib_drm=NaN(n_trials,120);
for k=1:n_np
  log_errors(k)=0;
  np=udb10(np_set_db(k));

  state=etc;
  for t=1:n_trials
    added_noise=blnoise(960)*sqrt(np);
    mib_rx=mib_sym+added_noise;
    mib_bp=deqam(mib_rx,np,'QAM','LTE');
    mib_drm=lte_conv_deratematch(mib_bp,40);
    log_mib_drm(t,:)=log((1-mib_drm(:))./mib_drm(:));
    mib_bits_rx=lte_conv_decode(mib_drm);
    log_errors(k)=log_errors(k)+any(mib_bits_rx~=mib);
    if (t==n_trials)||(mod(t,20)==0)
      t
      log_errors/t
      temp=log_errors;
      temp(1:k-1)=temp(1:k-1)/n_trials;
      temp(k)=temp(k)/t;
      plot(np_set_db,temp,'x-');
      drawnow;
    end
    state=etc(state,t/n_trials);
  end
end

