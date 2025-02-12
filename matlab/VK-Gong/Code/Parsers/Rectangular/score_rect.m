function [ f_time, Tn ] = score_rect( ScoreFileName, Lx, Ly, hd, rho, kx, ky, BC, Nphi, scheme, C, fsd, Tsd)
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
%     'Strike' [T0 fm Twid fp_x fp_y ];
%     'Harmonic' {T0 [Frequency] [Amplitude] Phase [Times] fp_x fp_y};
%     'ColoredNoise' {Color [T0 Amplitude fmin fmax deltaf TimeLength fp_x fp_y]};
% };

if exist(ScoreFileName, 'file') == 2
    load(ScoreFileName);
else
    disp('Score file not found. Using preset instead.');

    score_cell = {
        'Strike' [.01 500 1e-3 0.3 0.5];
         'Harmonic' {70e-3 [1 10 1] [1 50 1] 0 [0 0.5 1]  0.3 0.7};
    %     'ColoredNoise' {'Blue' [1.1 50 20 5000 1 1 0.3 0.7]};
    %     'ColoredNoise' {'White' [2.1 50 20 5000 1 1 0.3 0.7]};
    %     'ColoredNoise' {'Purple' [3.1 50 20 5000 1 1 0.3 0.7]};
    %     'ColoredNoise' {'Pink' [4.1 50 20 5000 1 1 0.3 0.7]};
    %    'ColoredNoise' {'Red' [5.1 50 20 5000 1 1 0.3 0.7]};
        };
end

Tn = round(Tsd*fsd);


Ns = size(score_cell,1);

fex = zeros(Tn,Ns);
fp_x = zeros(Ns,1);
fp_y = zeros(Ns,1);

T0 = 0;

for n = 1:Ns
    switch score_cell{n,1}
        case 'Strike'
            score_mat = cell2mat(score_cell(n,2));
            
            T0 = T0 + score_mat(1);
            fm = score_mat(2);
            Twid = score_mat(3);             
            fp_x(n) = score_mat(4);
            fp_y(n) = score_mat(5);
            
            fex(:,n) = StrikeExcitation( T0, fm, Twid,  fsd, Tn );
            
        case 'Harmonic'
            score_mat = score_cell{n,2};
            T0 = T0 + cell2mat(score_mat(1));
            f = cell2mat(score_mat(2));
            Amplitude = cell2mat(score_mat(3));
            Phase = cell2mat(score_mat(4));
            Times = cell2mat(score_mat(5)); 
            fp_x(n) = cell2mat(score_mat(6));
            fp_y(n) = cell2mat(score_mat(7));
            
            fex(:,n) = HarmonicSignal(T0, f, Amplitude, Phase, Times, fsd, Tn );
            
        case 'ColoredNoise'
            Color = score_cell{n,2}{1};
            score_mat = cell2mat(score_cell{n,2}(2));
            T0 = T0 + score_mat(1);
            Amplitude = score_mat(2);
            fmin = score_mat(3)*tnd;
            fmax = score_mat(4)*tnd;
            deltaf = score_mat(5)*tnd;
            TimeLength = score_mat(6); 
            fp_x(n) = score_mat(7);
            fp_y(n) = score_mat(8);
            
            fex(:,n) = ColoredNoiseSignal( Color, T0, Amplitude, fmin, fmax, deltaf, TimeLength, fs, Tn );
    end
end

%% Compute rp and P, factors applied to input and output
P =  zeros(Nphi,Ns);
f_time = zeros(Nphi,Tn);


for n = 1 : Ns
    switch BC
        case {'SimplySupported'}                        
            P(:,n) = sin(fp_x(n)*pi*kx/Lx).*sin(fp_y(n)*pi*ky/Ly);       
    end

    P(:,n) = P(:,n)/rho/(Lx*Ly/4)/hd;    
    
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

