
numSignals = 10;
iterations = 1000;

phases = 2*pi*rand(numSignals, iterations);
amplitude = sum(exp(1i*phases), 1);

power = 20*log10(abs(amplitude))';