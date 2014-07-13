% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% convert original uint8 data captured by rtl_sdr to double complex data.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

% all variable are assumed to be column vector
function b = raw2iq(a)
c = a(1:2:end,:) + 1i.*a(2:2:end,:);
% c = a(2:2:end,:) + 1i.*a(1:2:end,:); % swap

% b = c- kron(ones(size(c,1),1), ( sum(c, 1)./size(c,1) ));

% b = c - 128 - 1i.*128;

b = (c - 128 - 1i.*128)./128; % to comform to C code of LTE-Cell-Scanner

