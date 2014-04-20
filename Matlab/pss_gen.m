% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% generate frequency domain PSS and time domain PSS signal
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function [fd_pss, td_pss] = pss_gen
fd_pss = zeros(128,3);
td_pss = zeros(128+9,3);
for i=1:3
    temp=pss(i-1);
    fd_pss(:,i)=[0 temp(32:end) zeros(1,65) temp(1:31)].';
    temp_td=idft(fd_pss(:,i))*sqrt(128/62);
    td_pss(:,i)=[temp_td(end-8:end); temp_td];
end
