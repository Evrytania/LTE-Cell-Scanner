function a = ex_odd_even(a)
temp = a(:, 2:2:end);
a(:, 2:2:end) = a(:, 1:2:end);
a(:, 1:2:end) = temp;
