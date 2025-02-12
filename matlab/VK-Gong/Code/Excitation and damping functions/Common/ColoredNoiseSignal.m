function [ fex ] = ColoredNoiseSignal( Color, T0, Amplitude, fmin, fmax, deltaf, TimeLength, fs, Tn )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computes a vector of Tn samples of colored noise and maximum amplitude fm. 



% Even number of samples is forced
N = round(TimeLength*fs);
if mod(N,2)
    N = N + 1;
end

% Generate limited band white noise
f = (fmin : deltaf : fmax);
phase = rand(1, length(f));

n = linspace(0, TimeLength, N);

y = 0;

for i = 1: length(f)
    y = y + cos(2*pi*f(i)*n + phase(i)*2*pi);
end

Y = fft(y); 

Nm = N/2 + 1;
n = 1:Nm;



switch Color
    case 'White'
        
    case 'Pink'
        n = sqrt(n);
        Y(1:Nm) = Y(1:Nm)./n;
    case 'Blue'
        n = sqrt(n);
        Y(1:Nm) = Y(1:Nm).*n;

    case 'Red'
        Y(1:Nm) = Y(1:Nm)./n;
        
    case 'Purple'
        Y(1:Nm) = Y(1:Nm).*n;
    otherwise
        disp('Color not recognised');
end

Y(Nm+1:N) = real(Y(N/2:-1:2)) -1i*imag(Y(N/2:-1:2));

y = ifft(Y);
y = real(y(1, 1:round(TimeLength*fs)));

y = y - mean(y);

y = y/max(y)*Amplitude;

fex = zeros(Tn,1);

n = (1:length(y));
N0 = floor(fs*T0);

n = n + N0;

if n(end)>Tn
    disp('Warning: Excitation lasts longer than the simulation time. Samples after Ts will be omitted.');
    disp(['Deleted samples: ', num2str(sum(n>Tn))]);
    n(n>Tn) = [];
    y = y(1:length(n));
end


fex(n) = y;

end

