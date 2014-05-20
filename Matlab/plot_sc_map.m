function plot_sc_map(sc_map, tfg)
% sc_map is 1200*n_pdcch_symb matrix
% values definition:
% 1 -- rs
% 2 -- pcfich
% 3 -- phich
% 4 -- pdcch

num_sub = size(sc_map, 2);
nSC = size(sc_map, 1);
sc_all = 0 : (nSC-1);

for i = 1 : num_sub
    subplot(num_sub,1,i);
    sub_sc_map = sc_map(:,i);
    legend_str_count = 1;
    legend_str = cell(1, 4);
    if sum(sub_sc_map==1)>0
        bar(sc_all(sub_sc_map==1), abs( tfg(i, sub_sc_map==1) ), 'r', 'EdgeColor', 'none', 'BarWidth', 0.001); axis tight; hold on;
%         bar(sc_all(sub_sc_map==1), abs( tfg(i, sub_sc_map==1) ), 'r', 'linewidth', 0.05); axis tight; hold on;
        legend_str{legend_str_count} = 'RS';
        legend_str_count = legend_str_count + 1;
    end
    if sum(sub_sc_map==2)>0
        bar(sc_all(sub_sc_map==2), abs( tfg(i, sub_sc_map==2) ), 'g', 'EdgeColor', 'none', 'BarWidth', 0.001); axis tight; hold on;
%         bar(sc_all(sub_sc_map==2), abs( tfg(i, sub_sc_map==2) ), 'g', 'linewidth', 0.05); axis tight; hold on;
        legend_str{legend_str_count} = 'PCFICH';
        legend_str_count = legend_str_count + 1;
    end
    if sum(sub_sc_map==3)>0
        bar(sc_all(sub_sc_map==3), abs( tfg(i, sub_sc_map==3) ), 'b', 'EdgeColor', 'none', 'BarWidth', 0.001); axis tight; hold on;
%         bar(sc_all(sub_sc_map==3), abs( tfg(i, sub_sc_map==3) ), 'b', 'linewidth', 0.05); axis tight; hold on;      
        legend_str{legend_str_count} = 'PHICH';
        legend_str_count = legend_str_count + 1;
    end
    if sum(sub_sc_map==4)>0
        bar(sc_all(sub_sc_map==4), abs( tfg(i, sub_sc_map==4) ), 'k', 'EdgeColor', 'none', 'BarWidth', 0.001); axis tight; hold on;
%         bar(sc_all(sub_sc_map==4), abs( tfg(i, sub_sc_map==4) ), 'k', 'linewidth', 0.05); axis tight; hold on;
        legend_str{legend_str_count} = 'PDCCH';
        legend_str_count = legend_str_count + 1;
    end
    legend_str = legend_str(1 : (legend_str_count-1));
    legend(legend_str);
end
