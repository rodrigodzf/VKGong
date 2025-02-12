function [ fex ] = HarmonicSignal(T0, f, Amplitude, phase, Times, fs, Tn )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a harmonic wave signal.

%% Input arguments
% Amplitude: Amplitude of the signal
% Ts: Time length
% fs: Sampling rate

if ((length(Amplitude) ~= length(Times))||(length(Amplitude) ~= length(f)))
    disp('Harmonic excitation: Amplitude, Times and f must have the same length');
    return
end

if ((length(Times)>1)&&(Times(1) > 0))
    disp('Warning! Harmonic excitation: Times(1) should be 0');
end

fex = zeros(Tn,1);

t = (0:1/fs:Times(end));

ff = ones(size(t))*f(1);
aa = ones(size(t))*Amplitude(1);

for n = 1:length(Times)-1
    t0 = floor(Times(n)*fs);
    if t0 == 0
        t0 = t0+1;
    end
    tend = floor(Times(n+1)*fs);
    T = tend - t0 + 1 ;
    ff(t0:tend) = linspace(f(n), f(n+1), T);
    aa(t0:tend) = linspace(Amplitude(n), Amplitude(n+1), T);
end

%y = Amplitude*sin(2*pi*f*t + phase);

y = aa.*sin(2*pi*ff.*t + phase);

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

