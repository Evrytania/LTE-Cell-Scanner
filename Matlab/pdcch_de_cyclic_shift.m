function b = pdcch_de_cyclic_shift(a, n_id_cell)
b = a;

for i = 1 : size(a,3);
    b(:,:,i) = circshift(a(:,:,i), [n_id_cell, 0]);
end
