

% (CML) trajectory prediction with destination information

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
  
%============================

%-------------------------
 X0_m_m = [2000;70;5000;0];
X0_Cov_m = [10000 0 0 0;
            0 100 0 0;
            0 0 10000 0;
            0 0 0 100];
        
XN_m_m = [130000;70;2000;0];
XN_Cov_m = [10000 0 0 0;
            0 100 0 0;
            0 0 10000 0;
            0 0 0 100];
        
%         C_N0_t = zeros(4,4);
        C_N0_m = [7000 0 0 0;
                  0 60 0 0;
                  0 0 7000 0;
                  0 0 0 60];

%-------------------------
%===========================================


kf = 9;

iteration = 1000;

%..................

EEp = zeros(iteration,N+1);
EEv = zeros(iteration,N+1);

% Px11 = zeros(1,N+1);
% Px22 = zeros(1,N+1);
% Py11 = zeros(1,N+1);
% Py22 = zeros(1,N+1);
% 
% PPx11 = zeros(1,N+1);
% PPx22 = zeros(1,N+1);
% PPy11 = zeros(1,N+1);
% PPy22 = zeros(1,N+1);
% 
%   W11=zeros(1,N+1);
%   W22=zeros(1,N+1);  
%   
%   W33=zeros(1,N+1);
%   W44=zeros(1,N+1);


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
    
   
    
    %+++++++++++++++++++++++++++++++
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
        
        %-----------CN|k
        
        CNk = zeros(4,4);
        
        for ii=0:N-k-1
            
            CNk = CNk + F^ii*Q*(F^ii)';
            
        end
        
        %============
        %----------
         %.................
        Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;  
        Gk_N = Gk*F^(N-k)'/(CNk);
        Gk_km1 = F - Gk_N*F^(N-k+1);
        %...............
        
        DG = chol(Gk,'lower');
        
        Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N+1) + DG*randn(4,1);
        
        %..........
        
    end % for k=2:N-1 measurements
    
    
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
    
    %--------------------
    
    Yh = zeros(8,N);
    Xh = zeros(4,N+1);
    
    Yh(:,1) = Ys;
    Xh(:,1) = X0_m_m;
    
    %------------------------
    
    
    Hy = H*[eye(4) zeros(4,4)];
    %======================================= Main Loop
    
    for k = 1:N
        
        
        if k<kf+1
            
            
            
            %---------------------------Prediction from k-1 to k
            
            %............................Dynamic model coefficients
            %............. CN|k
            
            CNk = zeros(4,4);
            
            for ii=0:N-k-1
                
                CNk = CNk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;
            
            Gk_N = Gk*F^(N-k)'/(CNk);
            Gk_km1 = F - Gk_N*F^(N-k+1);
            
            
            
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
            Ph = [eye(4) zeros(4,4)]*Pyh*[eye(4) zeros(4,4)]';
            
            %..................
            %.................. For the next time
            
            Ys = Yh(:,k+1);
            Pys = Pyh;
            
            
            
            %+++++++++++++++++++++
            
            
            
        elseif k>kf && k<N
            
            %===================================== Prediction
            %=======================================
            
            %--------------------------------------W preparation
            %............. CN|k
            
            n=k-kf;
            
            Cknk = zeros(4,4);
            
            for ii=0:n-1
                
                Cknk = Cknk + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            CNkn = zeros(4,4);
            
            for ii=0:N-(kf+n)-1
                
                CNkn = CNkn + F^ii*Q*(F^ii)';
                
            end
            %.............
            
            W = Cknk - Cknk*F^(N-(kf+n))'/(CNkn + F^(N-(kf+n))*Cknk*F^(N-(kf+n))')*F^(N-(kf+n))*Cknk;
            
            w2 = W*F^(N-(kf+n))'/(CNkn);
            w1 = F^n - w2*F^(N-kf);
            
            %............ To know the Velocity behavior
%             W11(1,k)=W(1,1);
%             W22(1,k)=W(2,2);
%             
%             w1_11(1,k)=w1(1,1);
%             w1_22(1,k)=w1(2,2);
% %             w1_12(1,k)=w1(1,1);
%             w1_21(1,k)=w1(2,1);
%             
%             w2_11(1,k)=w2(1,1);
%             w2_22(1,k)=w2(2,2);
%             w2_21(1,k)=w2(2,1);

%=============================Gk
%            W11(1,k+1)=W(1,1);
%            W22(1,k+1)=W(2,2);
%            
%            W33(1,k+1)=W(3,3);
%            W44(1,k+1)=W(4,4);
           
           %==============================
            
            %........................................
            
            E = [w1 w2];
            
            Xh(:,k+1) = E*Yh(:,kf+1);
            Pkn = W + E*Pyh*E';
            
%             %================================================
%             PP=E*Pyh*E';
%             PPx11(1,k+1)=PP(1,1);
%             PPx22(1,k+1)=PP(2,2);
%             PPy11(1,k+1)=PP(3,3);
%             PPy22(1,k+1)=PP(4,4);
%             
%             Px11(1,k+1)=Pkn(1,1);
%             Py11(1,k+1)=Pkn(3,3);
%             
%             Px22(1,k+1)=Pkn(2,2);
%             Py22(1,k+1)=Pkn(4,4);
%             
%             %================================================
            
            
        elseif k==N
            
            Xh(:,N+1) = [zeros(4,4) eye(4)]*Yh(:,kf+1);
            Pkn = [zeros(4,4) eye(4)]*Pyh*[zeros(4,4) eye(4)]';
            
            
            
        end
        
        
    end
    
    %-----------------------
    
    
    %-----------------------
    
    EEp(iter,:) = sqrt((Xr(1,:) - Xh(1,:)).^2 + (Xr(3,:) - Xh(3,:)).^2);
    EEv(iter,:) = sqrt((Xr(2,:) - Xh(2,:)).^2 + (Xr(4,:) - Xh(4,:)).^2);
    %     EEz(iter,:) = sqrt((Z(1,:) - Xr(1,:)).^2 + (Z(2,:) - Xr(3,:)).^2);
    
    
    %========================================================== 
    
end  % for iter=1:iteration

%-----------------------------

AEEp = mean(EEp,1);
AEEv = mean(EEv,1);


%----------------------------




kk=0:N;


figure(1)
hold on
plot(kk,log10(AEEp),'.r')

figure(2)
hold on
plot(kk,log10(AEEv),'.r')

% figure(3)
% hold on
% subplot(2,1,1)
% plot(kk,Px11)
% subplot(2,1,2)
% plot(kk,Py11)
% 
% figure(4)
% hold on
% subplot(2,1,1)
% plot(kk,Px22)
% subplot(2,1,2)
% plot(kk,Py22)
% 
% 
% figure(5)
% subplot(2,1,1)
% plot(kk,W11)
% hold on
% subplot(2,1,2)
% plot(kk,W33)
% 
% 
% figure(6)
% subplot(2,1,1)
% plot(kk,W22)
% hold on
% subplot(2,1,2)
% plot(kk,W44)
% 
% figure(7)
% hold on
% subplot(2,1,1)
% plot(kk,PPx11)
% subplot(2,1,2)
% plot(kk,PPy11)
% 
% figure(8)
% hold on
% subplot(2,1,1)
% plot(kk,PPx22)
% subplot(2,1,2)
% plot(kk,PPy22)
% 
% % figure(4)
% % subplot(3,1,1)
% % plot(w1_11)
% % hold on
% % subplot(3,1,2)
% % plot(w1_22)
% % hold on
% % subplot(3,1,3)
% % plot(w1_21)
% % 
% % figure(5)
% % subplot(3,1,1)
% % plot(w2_11)
% % hold on
% % subplot(3,1,2)
% % plot(w2_22)
% % hold on
% % subplot(3,1,3)
% % plot(w2_21)
% 
% 
% 
% 
% 
% 
% 
% 
% 
