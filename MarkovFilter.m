

% Tracker based on a Markov Model without destination

clear all
% close all
clc

%--------------------

N = 100;

T = 5;  %7

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];

q = 1;

Q = [q*T^3/3 q*T^2/2 0 0;...
    q*T^2/2 q*T 0 0;...
    0 0 q*T^3/3 q*T^2/2;...
    0 0 q*T^2/2 q*T];

D = chol(Q,'lower');

sigma = [10 10];
R = diag(sigma.^2);

H = [1 0 0 0;...
    0 0 1 0];

%================================= 

%-------------------Truth

X0_m_t = [2000;70;5000;0];
X0_Cov_t = diag([1000;10;1000;10]);

% XN_m_t = [30000;20;5000;-100];
% XN_Cov_t = diag([1000;10;1000;10]);


X0_m_m = [2500;95;5500;25];
X0_Cov_m = diag([1000;10;1000;10]);


%-----------------------------

iteration = 1000;

%..................

SEp = zeros(iteration,N+1);
SEv = zeros(iteration,N+1);

Xh = zeros(4,N+1);

for iter=1:iteration

%     iter

    Xr = zeros(4,N+1);

    %--------------------

          Dxn = chol(X0_Cov_t,'lower');
          Xr(:,1) = X0_m_t + Dxn*randn(4,1);
    
          for k=1:N
    
              v = randn(4,1);
              Xr(:,k+1) = F*Xr(:,k) + D*v;
    
          end
    
    %-------------------- 

    %-------------------------------------------- Measurement

   Dxn = chol(R,'lower');
    Z = [Xr(1,2:end);Xr(3,2:end)] + Dxn*randn(2,N);


    %--------------------------------------------

    %     figure(1)
    %     hold on
    %     plot(Xr(1,:),Xr(3,:),'.')
    %     %   axis([2000 4150 1800 2200])

    % %        %----------------------------- time k = 1

    Xs = X0_m_m;
    Ps = X0_Cov_m;
    
    Xh(:,1) = Xs;
    Ph = Ps;

    %-----------------------------

    for k=1:N


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

    end % for k=1:N

    SEp(iter,:) = ((Xr(1,:) - Xh(1,:)).^2 + (Xr(3,:) - Xh(3,:)).^2);
    SEv(iter,:) = ((Xr(2,:) - Xh(2,:)).^2 + (Xr(4,:) - Xh(4,:)).^2);
    
    

end  % for iter=1:iteration


%-----------------------------

RMSEp = sqrt(mean(SEp,1));
RMSEv = sqrt(mean(SEv,1));


% ----------------------------


figure(12)
hold on
plot((RMSEp),'k')

figure(13)
hold on
plot((RMSEv),'k')









