% function merge_result
clear all; close all;

load CellSearch_test1to8_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;
cell_info_save = cell_info;

load CellSearch_test9to16_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;
cell_info_save(test_sp:test_ep) = cell_info(test_sp:test_ep);

load CellSearch_test17to24_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;
cell_info_save(test_sp:test_ep) = cell_info(test_sp:test_ep);

load CellSearch_test25to33_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat;
cell_info_save(test_sp:test_ep) = cell_info(test_sp:test_ep);

cell_info = cell_info_save;

test_sp = 1;
test_ep = 33;

clear cell_info_save;

save CellSearch_test1to33_twist0_fo-140to140_resv2_numPtry1_Prange-1_parTh8_numPth0.5.mat

