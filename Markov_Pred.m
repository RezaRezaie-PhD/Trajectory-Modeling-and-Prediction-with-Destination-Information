

% Markov trajectory prediction without destination

clear all
% close all
clc

%--------------------

N = 100;

T = 15;  %7

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];

q = 0.01;

Q = [q*T^3/3 q*T^2/2 0 0;...
    q*T^2/2 q*T 0 0;...
    0 0 q*T^3/3 q*T^2/2;...
    0 0 q*T^2/2 q*T];

D = chol(Q,'lower');

sigma = [10 10];
R = diag(sigma.^2);

H = [1 0 0 0;...
    0 0 1 0];

%-------------------Truth

X0_m_t = [2000;70;5000;0];
X0_Cov_t = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];
        
XN_m_t = [130000;70;2000;0];
XN_Cov_t = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];
        
%         C_N0_t = zeros(4,4);
        C_N0_t = [800 20 0 0;
                  20 7 0 0;
                  0 0 800 20;
                  0 0 20 7];



X0_m_m = [2000;70;5000;0];
X0_Cov_m = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];
%------------------------------

kf = 9;

iteration = 1000;

%..................

EEp = zeros(iteration,N+1);
EEv = zeros(iteration,N+1);

Xh = zeros(4,N+1);

for iter=1:iteration

%     iter

    Xr = zeros(4,N+1);

    %       %--------------------
    %
    %       Dxn = chol(X1_Cov_m_r);
    %       Xr(:,1) = X1_M_m + Dxn*randn(4,1);
    %
    %       for k=2:N
    %
    %           v = randn(4,1);
    %           Xr(:,k) = F*Xr(:,k-1) + D*v;
    %
    %       end
    %
    %       %--------------------

    %++++++++++++++++++++++++++++
    %---------------

     Dxn = chol(X0_Cov_t,'lower');
    Xr(:,1) = X0_m_t + Dxn*randn(4,1);
    
    %----------
    
    % k = N
    
     Cn = XN_Cov_t - C_N0_t/X0_Cov_t*C_N0_t';
     Dxn = chol(Cn,'lower');
     % Xr(:,N+1) = XN_m_t + Dxn*randn(4,1);
      Xr(:,N+1) = XN_m_t + C_N0_t/X0_Cov_t*(Xr(:,1) - X0_m_t) + Dxn*randn(4,1);
      
    %----------
    %+++++++++++++++++++++++++++++++++++++


    for k=1:N-1

        %           k

        %-----------CN|k

        CNk = zeros(4,4);

        for ii=0:N-k-1

            CNk = CNk + F^ii*Q*(F^ii)';

        end

        %.................
        Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;  
        Gk_N = Gk*F^(N-k)'/(CNk);
        Gk_km1 = F - Gk_N*F^(N-k+1);
        %...............
        
        DG = chol(Gk,'lower');
        
        Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N+1) + DG*randn(4,1);
        
        %..........

    end % for k=2:N-1 measurements

    %--------------------

    %-------------------------------------------- Measurement

   Dxn = chol(R,'lower');
    Z = [Xr(1,2:end);Xr(3,2:end)] + Dxn*randn(2,N);


    %--------------------------------------------

    %     figure(1)
    %     hold on
    %     plot(Xr(1,:),Xr(3,:),'.')
    %     %   axis([2000 4150 1800 2200])


   
   %----------------------------- time k = 1

    Xs = X0_m_m;
    Ps = X0_Cov_m;
    
    Xh(:,1) = Xs;
    Ph = Ps;

    %-----------------------------

    for k=1:N

        if k<kf+1

            %-----------Prediction

            Xp = F*Xs;
            Pp = F*Ps*F' + Q;

            %---------- Update

            S = H*Pp*H' + R;
            K = Pp*H'/S;

            Xh(:,k+1) = Xp + K*(Z(:,k) - H*Xp);
            Ph = Pp - K*S*K';

            %----------

            Xs = Xh(:,k+1);
            Ps = Ph;

        else

            Xp = F*Xs;
            Pp = F*Ps*F' + Q;

            Xh(:,k+1) = Xp;
            Ph = Pp;

            Xs = Xp;
            Ps = Pp;

        end

    end % for k=1:N

    EEp(iter,:) = sqrt((Xr(1,:) - Xh(1,:)).^2 + (Xr(3,:) - Xh(3,:)).^2);
    EEv(iter,:) = sqrt((Xr(2,:) - Xh(2,:)).^2 + (Xr(4,:) - Xh(4,:)).^2);
   


end  % for iter=1:iteration


%-----------------------------

AEEp = mean(EEp,1);
AEEv = mean(EEv,1);


% ----------------------------

kk=0:N;

figure(1)
hold on
plot(kk,log10(AEEp),'or')

figure(2)
hold on
plot(kk,log10(AEEv),'or')

% figure(3)
% hold on
% plot(Xr(2,:))







