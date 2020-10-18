

% Tracker based on a CML model with destination information


clear all
% close all
clc

%--------------------

N = 100;

T = 15;

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];


q = 0.01;
%  Q = q*eye(4);
Q = [q*T^3/3 q*T^2/2 0 0;...
    q*T^2/2 q*T 0 0;...
    0 0 q*T^3/3 q*T^2/2;...
    0 0 q*T^2/2 q*T];


D = chol(Q,'lower');

%------------------------------

%=================================

  
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
  
%============================Known

%-------------------------
 X0_m_m = [2000;70;5000;0];
X0_Cov_m = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];
        
XN_m_m = [130000;70;2000;0];
XN_Cov_m = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];
        
%         C_N0_t = zeros(4,4);
        C_N0_m = [800 20 0 0;
                  20 7 0 0;
                  0 0 800 20;
                  0 0 20 7];

%-------------------------


%============================================
  

%------------------------------
sigma_x = 10;
sigma_y = 10;
R = diag([sigma_x sigma_y].^2);

H = [1 0 0 0;0 0 1 0];
%------------------------------

iteration = 1000;

%..................

EEp = zeros(iteration,N+1);
EEv = zeros(iteration,N+1);
% EEz = zeros(iteration,N);

EENp = zeros(iteration,N+1);
EENv = zeros(iteration,N+1);

% EEN = zeros(iteration,N+1);

tic

for iter=1:iteration
    
    %     iter
    
    Xr = zeros(4,N+1);
    
    
    %++++++++++++++++++++++++++++++++++
    %---------------
    
    % k = 0
    
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
        
        %-----------CN|k
        
        CNk = zeros(4,4);
        
        for ii=0:N-k-1
            
            CNk = CNk + F^ii*Q*(F^ii)';
            
        end
        
        %============
        %----------
         %..............
%         Gk = inv((inv(Q) + F^(N-k)'/(CNk)*F^(N-k)));
%         Gk_N = ((inv(Q) + F^(N-k)'/(CNk)*F^(N-k)))\F^(N-k)'/(CNk);
%         Gk_km1 = F - Gk_N*F^(N-k+1);
        
        %..............
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
    
    %==========================================================TRACKING
    
    %--------------------------------Initializaation
    %....................
    
    Ys = [X0_m_m;XN_m_m];
    Ps0 = X0_Cov_m;
    PsN = XN_Cov_m;
    PsN0 = C_N0_m;
    
    Pys = [Ps0 PsN0';
           PsN0 PsN];
    
    
    % k: 0 to 1
    
    %--------------------
    
    Yh = zeros(8,N);
    Xh = zeros(4,N+1);
    
    Yh(:,1) = Ys;
    Xh(:,1) = X0_m_m;
    
    XhN = zeros(4,N+1);
    XhN(:,1) = XN_m_m;
    
    %     SDhN11 = zeros(1,N+1);
    %     SDhN22 = zeros(1,N+1);
    %     SDhN33 = zeros(1,N+1);
    %     SDhN44 = zeros(1,N+1);
    %
    %     SDhN11(1,1) = sqrt(XN_Cov_m(1,1));
    %     SDhN22(1,1) = sqrt(XN_Cov_m(2,2));
    %     SDhN33(1,1) = sqrt(XN_Cov_m(3,3));
    %     SDhN44(1,1) = sqrt(XN_Cov_m(4,4));
    
    
    Hy = H*[eye(4) zeros(4,4)];
    %======================================= Main Loop
    
    for k = 1:N
        
        if k<N
            
            
            %---------------------------Prediction from k-1 to k
            
            %............................Dynamic model coefficients
            %............. CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
%             Gk = inv(inv(Q)+F^(N-k)'/CNk*F^(N-k));
%             
%             Gk_N = (inv(Q)+F^(N-k)'/CNk*F^(N-k))\F^(N-k)'/CNk;
%             Gk_km1 = F - Gk_N*F^(N-k+1);
            
            %..............
            
            Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;  
            
            Gk_N = Gk*F^(N-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N-k+1);
            
            %---------------------------
           
            
            Gy1 = horzcat(Gk_km1,Gk_N);
            Gy2 = horzcat(zeros(4,4),eye(4));
            Gy = vertcat(Gy1,Gy2);
            
            GCy = blkdiag(Gk,zeros(4,4));
            
            %.............
            
            Yp = Gy*Ys;
            Pyp = Gy*Pys*Gy' + GCy;
            
            
            %--------------------------------------------------
            
            Pyh = Pyp - (Pyp*Hy')/(Hy*Pyp*Hy'+R)*(Pyp*Hy')';
            
            Yh(:,k+1) = Yp + (Pyp*Hy')/(Hy*Pyp*Hy'+R)*(Z(:,k) - Hy*Yp);
            
            Xh(:,k+1) = [eye(4) zeros(4,4)]*Yh(:,k+1);
            %             Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
            
            %..................
            %.................. For the next time
            
            Ys = Yh(:,k+1);
            Pys = Pyh;
            
            %..................
            
            %--------------------------------------------------
            
            %+++++++++++++++++++++
            
            XhN(:,k+1) = [zeros(4,4) eye(4)]*Yh(:,k+1);
            %             PhN = [zeros(4,4) eye(4)]*Pyh*[zeros(4,4) eye(4)]';
            
            %             SDhN11(1,k+1) = sqrt(PhN(1,1));
            %             SDhN22(1,k+1) = sqrt(PhN(2,2));
            %             SDhN33(1,k+1) = sqrt(PhN(3,3));
            %             SDhN44(1,k+1) = sqrt(PhN(4,4));
            
            %             PN11(1,k) = PN(1,1);
            %             PN22(1,k) = PN(2,2);
            %             PN33(1,k) = PN(3,3);
            %             PN44(1,k) = PN(4,4);
            
            %             %+++++++++++++++++++++ Correlation Xk XN
            %
            %             PkN = Pyh(1:4,5:8);
            %
            %             CXkN = Pyh(1:4,5:8);
            %             CXN = Pyh(5:8,5:8);
            %
            %             B = CXkN/CXN;
            %
            %             B11(1,k) = B(1,1);
            %             B12(1,k) = B(1,2);
            %             B21(1,k) = B(2,1);
            %             B22(1,k) = B(2,2);
            
            %+++++++++++++++++++++
            
        elseif k == N
            
            
            Xp = [zeros(4,4) eye(4)]*Ys;
            Pp = [zeros(4,4) eye(4)]*Pys*[zeros(4,4) eye(4)]';
            
            Ph = Pp - (Pp*H')/(H*Pp*H'+R)*(Pp*H')';
            
            Xh(:,k+1) = Xp + (Pp*H')/(H*Pp*H'+R)*(Z(:,k) - H*Xp);
            
            
            %+++++++++++++++++++++
            
            XhN(:,k+1) = Xh(:,k+1);
            %             PhN = Ph;
            
            %             SDhN11(1,k+1) = sqrt(PhN(1,1));
            %             SDhN22(1,k+1) = sqrt(PhN(2,2));
            %             SDhN33(1,k+1) = sqrt(PhN(3,3));
            %             SDhN44(1,k+1) = sqrt(PhN(4,4));
            
            %             PN11(1,k) = PN(1,1);
            %             PN22(1,k) = PN(2,2);
            %             PN33(1,k) = PN(3,3);
            %             PN44(1,k) = PN(4,4);
            %
            %             %+++++++++++++++++++++ Correlation Xk XN
            %
            %             PkN = Ph;
            %
            %             B = Ph/Ph;
            %
            %             B11(1,k) = B(1,1);
            %             B12(1,k) = B(1,2);
            %             B21(1,k) = B(2,1);
            %             B22(1,k) = B(2,2);
            
            
            %+++++++++++++++++++++
            
        end
        
    end
    
    %-----------------------
    
    
    %-----------------------
    
    EEp(iter,:) = sqrt((Xr(1,:) - Xh(1,:)).^2 + (Xr(3,:) - Xh(3,:)).^2);
    EEv(iter,:) = sqrt((Xr(2,:) - Xh(2,:)).^2 + (Xr(4,:) - Xh(4,:)).^2);
    %     EEz(iter,:) = sqrt((Z(1,:) - Xr(1,2:end)).^2 + (Z(2,:) - Xr(3,2:end)).^2);
    
    EENp(iter,:) = sqrt((Xr(1,N+1) - XhN(1,:)).^2 + (Xr(3,N+1) - XhN(3,:)).^2);
    EENv(iter,:) = sqrt((Xr(2,N+1) - XhN(2,:)).^2 + (Xr(4,N+1) - XhN(4,:)).^2);
    
%     EEN(iter,:) = sqrt((Xr(1,N+1) - XhN(1,:)).^2 + (Xr(3,N+1) - XhN(3,:)).^2 + ... 
%     (Xr(2,N+1) - XhN(2,:)).^2 + (Xr(4,N+1) - XhN(4,:)).^2);
    
    
    %========================================================== END OF TRACKING
    
end  % for iter=1:iteration

TT = toc
%-----------------------------

AEEp = mean(EEp,1);
AEEv = mean(EEv,1);
% AEEz = mean(EEz,1);

AEENp = mean(EENp,1);
AEENv = mean(EENv,1);

% AEEN = mean(EEN,1);
%----------------------------


% figure(1)
% hold on
% plot(Xr(1,:),Xr(3,:),'b')

kk = 0:N;

figure(1)
hold on
plot(kk,AEEp,'--k')

figure(2)
hold on
plot(kk,AEEv,'--k')

figure(3)
hold on
plot(kk,log10(AEENp),'--k')

figure(4)
hold on
plot(kk,log10(AEENv),'--k')


% figure(5)
% hold on
% plot(kk,AEEN,'r')
% 
% figure(6)
% hold on
% plot(kk,log10(AEEN),'r')


% figure(6)
% hold on
% plot(XhN(1,:),'b')
%
% figure(7)
% hold on
% plot(SDhN11(1,:),'b')
%
% figure(8)
% hold on
% plot(XhN(2,:),'b')
%
% figure(9)
% hold on
% plot(SDhN22(1,:),'b')















