function [ f_time, Tn ] = score_circ( ScoreFileName, Rd, hd, E, BC, e, Nphi, scheme, C, k_t, c_t, xkn, JJ, II, Kkn, tnd, fs, Tsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parser of the excitation file. 


%% Score_cell structure
% score_cell = {
%     'Strike' [T0 fm Twid fp_th fp_r ];
%     'Harmonic' {T0 [Frequency] [Amplitude] Phase [Times] fp_th fp_r};
%     'ColoredNoise' {Color [T0 Amplitude fmin fmax deltaf TimeLength fp_th fp_r]};
% };

if exist(ScoreFileName, 'file') == 2
    load(ScoreFileName);
else
    disp('Score file not found. Using preset instead.');

    score_cell = {
        'Strike' [.01 500 1e-3 pi 0.5];
         'Harmonic' {70e-3 [1 10 1] [1 50 1] 0 [0 0.5 1] pi 0.7};
    %     'ColoredNoise' {'Blue' [1.1 50 20 5000 1 1 pi 0.7]};
    %     'ColoredNoise' {'White' [2.1 50 20 5000 1 1 pi 0.7]};
    %     'ColoredNoise' {'Purple' [3.1 50 20 5000 1 1 pi 0.7]};
    %     'ColoredNoise' {'Pink' [4.1 50 20 5000 1 1 pi 0.7]};
    %    'ColoredNoise' {'Red' [5.1 50 20 5000 1 1 pi 0.7]};
        };
end

Ts = Tsd/tnd;

Tn = round(Ts*fs);


Ns = size(score_cell,1);

fex = zeros(Tn,Ns);
fp_th = zeros(Ns,1);
fp_r = zeros(Ns,1);

T0 = 0;

for n = 1:Ns
    switch score_cell{n,1}
        case 'Strike'
            score_mat = cell2mat(score_cell(n,2));
            
            T0 = T0 + score_mat(1)/tnd;
            fm = score_mat(2)*Rd^4/(E*hd^4);
            Twid = score_mat(3)/tnd;             
            fp_th(n) = score_mat(4);
            fp_r(n) = score_mat(5);
            
            fex(:,n) = StrikeExcitation( T0, fm, Twid,  fs, Tn );
            
        case 'Harmonic'
            score_mat = score_cell{n,2};
            T0 = T0 + cell2mat(score_mat(1))/tnd;
            f = cell2mat(score_mat(2))*tnd;
            Amplitude = cell2mat(score_mat(3))*Rd^4/(E*hd^4);
            Phase = cell2mat(score_mat(4));
            Times = cell2mat(score_mat(5))/tnd; 
            fp_th(n) = cell2mat(score_mat(6));
            fp_r(n) = cell2mat(score_mat(7));
            
            fex(:,n) = HarmonicSignal(T0, f, Amplitude, Phase, Times, fs, Tn );
            
        case 'ColoredNoise'
            Color = score_cell{n,2}{1};
            score_mat = cell2mat(score_cell{n,2}(2));
            T0 = T0 + score_mat(1)/tnd;
            Amplitude = score_mat(2)*Rd^4/(E*hd^4);
            fmin = score_mat(3)*tnd;
            fmax = score_mat(4)*tnd;
            deltaf = score_mat(5)*tnd;
            TimeLength = score_mat(6)/tnd; 
            fp_th(n) = score_mat(7);
            fp_r(n) = score_mat(8);
            
            fex(:,n) = ColoredNoiseSignal( Color, T0, Amplitude, fmin, fmax, deltaf, TimeLength, fs, Tn );
    end
end

%% Compute rp and P, factors applied to input and output
P =  zeros(Nphi,Ns);
f_time = zeros(Nphi,Tn);


for n = 1 : Ns
    switch BC
        case {'free', 'elastic'}
            Jtild = JJ;
            Itild = II;
            for ii = 1 : Nphi

                JJ0p = besselj(k_t(ii),xkn(ii)*fp_r(n));
                II0p = besseli(k_t(ii),xkn(ii)*fp_r(n));


                P(ii,n) = (JJ0p - (Jtild(ii)*II0p/(Itild(ii))))*cos(k_t(ii)*fp_th(n)-(c_t(ii)-1)/2*pi);   % Force


                P(ii,n) = P(ii,n) * Kkn(ii);

            end
        case 'clamped'
            Jkn = JJ;
            Ikn = II;
            
            for ii = 1 : Nphi
                Jp=besselj(k_t(ii),xkn(ii)*fp_r(n));
                Ip=besseli(k_t(ii),xkn(ii)*fp_r(n));
                
                P(ii,n) = (Ikn(ii)*Jp-Jkn(ii)*Ip)*cos(k_t(ii)*fp_th(n)-(c_t(ii)-1)/2*pi); % Jtild = Jkn, Itild = Ikn
                
                P(ii,n) = P(ii,n) * Kkn(ii); 
                
            end
            
    end
    P(:,n) = e*P(:,n);
    P(:,n) = P(:,n)/Rd^2;    
    switch scheme
        case 'ECS'
            P(:,n) = P(:,n);
        case 'verlet'
            P(:,n) = P(:,n)./C;
    end
    f_time = f_time + P(:,n)*fex(:,n)';
end
    
disp ('Score file parsed');
        
end

