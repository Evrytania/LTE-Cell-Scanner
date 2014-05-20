function v = lte_subblock_interleave(d)

% Create the matrix before permutation.
n_c=32;
n_r=ceil(size(d,2)/n_c);
num_stream = size(d, 1);

% Subblock interleaving
perm_pattern=1+[1 17 9 25 5 21 13 29 3 19 11 27 7 23 15 31 0 16 8 24 4 20 12 28 2 18 10 26 6 22 14 30];
v=NaN(num_stream, n_r*n_c);
for t=1:num_stream
  y=[NaN(1,n_r*n_c-size(d,2)) d(t,:)];
  y=transpose(reshape(y,n_c,n_r));

  % Permute the columns
  y_perm=y(:,perm_pattern);
  v(t,:)=y_perm(:);
end
