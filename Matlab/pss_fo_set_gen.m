% Jiao Xianjun (putaoshu@msn.com; putaoshu@gmail.com)
% Generate pss time domain sequences under different frequency offsets.
% A script of project: https://github.com/JiaoXianjun/rtl-sdr-LTE

function pss_fo_set = pss_fo_set_gen(pss, fo_search_set)
num_pss = size(pss, 2);
len_pss = size(pss, 1);

sampling_rate = 1.92e6; % LTE spec
num_fo = length(fo_search_set);

pss_fo_set = zeros(len_pss, num_fo*num_pss);
for i=1:num_pss
    sp = (i-1)*num_fo + 1;
    ep = sp + num_fo - 1;
    pss_fo_set(:,sp:ep) = kron(ones(1, num_fo), pss(:,i)).*exp(1i.*2.*pi.*(1./sampling_rate).*(0:(len_pss-1)).'*fo_search_set);
    pss_fo_set(:,sp:ep) = conj(pss_fo_set(:,sp:ep))./len_pss;
end

% % normalize along with column
% pss_fo_set = sqrt(len_pss).*pss_fo_set./kron( ones(len_pss,1), sqrt( sum(abs(pss_fo_set).^2,1) ) ); 
