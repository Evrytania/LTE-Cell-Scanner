function show_time_frequency_grid_raw(s, sampling_rate)

grid_size = 8;
space_SC = 15e3*grid_size;
nSC = 1200/grid_size;

fft_size = sampling_rate/space_SC;

a = abs(fft( vec2mat(s, fft_size).', fft_size, 1)).^2;

b = [a(end-(nSC/2)+1 : end, :); a(1 : (nSC/2), :)];

pcolor(b); shading flat; drawnow;
