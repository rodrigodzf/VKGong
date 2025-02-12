function [ h ] = DisplayResults(signal, fsd, Fc, TimeSignal, FFT, Spectrogram, E, hd, nu, Rd, rho   )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that displays the results of the time integration.

k = 1;

if TimeSignal == 1
    h(k) = figure;    
    plot((1:length(signal))/fsd, signal); 
    xlabel('Time [s]');
    k = k+1;
end

if FFT >0 
    %% Define values por FFT representation
    n = floor(log(fsd)/log(2));
    NFFT = 2^(n+2);
    win = 2^(n+1);
    sample_0 = length(signal)-NFFT;
    
    while sample_0 < 0
        n = n-1;
        NFFT = 2^(n+2);
        win = 2^(n+1);
        sample_0 = length(signal)-NFFT;
    end
    
    sig = signal(sample_0:sample_0+NFFT-1);
    
    zero_padding_factor = NFFT/length(sig);
    zp=zero_padding_factor;

    F=(0:fsd/NFFT:(((NFFT/2)-1)*fsd)/NFFT);
    

    sig=sig'.*(win);
    spectre = fft(sig,NFFT);
    sp = abs(spectre(1:NFFT/2)).^2;
    splog=20*log(sp)/log(10);
    MAXI =  max(splog);
    splog = splog - MAXI;
    h(k) = figure;
    
    
    
    if FFT == 1 % plot in terms of frequency
        plot(F,splog);
        xlim([0 Fc]);
        xlabel('Frequency [Hz]')

    else % plot in terms of non-dimensional angular frequency
        D = E*hd^3/(12*(1-nu^2));
        tnd = Rd^2*sqrt(rho*hd/D);
        
        omega=F*tnd*2*pi;
        plot(omega,splog);
        xlim([0 Fc*tnd*2*pi]);
        xlabel('Dimensionless \omega')
    end

    k = k+1;
end

if Spectrogram == 1
    n = floor(log(fsd)/log(2));
    N = 2^(n-2);
    HS = 2^(n-3);

    Nwin = floor((max(size(signal))-N)/HS);
    hamming = 0.5*(1-cos(2*pi*[0:N-1]'/N));
    X = zeros(N,Nwin);

    for qq = 1:Nwin
        vec = hamming.*signal((qq-1)*HS+1:(qq-1)*HS+N);
        X(:,qq) = 10*log10(abs(fft(vec)));
    end

    Nc = floor(Fc*N/fsd);
    
    h(k) = figure;
    surf((0:Nwin-1)*HS/fsd, (0:Nc-1)*Fc/Nc, X(1:Nc,:));

    colormap jet
    shading interp
    view(2)
    axis tight
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
end



savefig(h,'OutputPlots.fig');
save('OutputPlots.mat','signal','fsd','Fc');
end

