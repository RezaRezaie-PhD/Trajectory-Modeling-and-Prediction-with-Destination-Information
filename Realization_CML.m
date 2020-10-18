

% This code generates trajectories (CML trajectories) with a given destination 

clear all
close all
clc

%--------------------

N = 100;

T = 15;

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
  

  %======================================
   
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
        
C_N0 = [800 20 0 0;
        20 7 0 0;
        0 0 800 20;
        0 0 20 7];
  
  
  
  iteration = 50;
  
  Xr = zeros(4,N+1);

      
  for iter=1:iteration
      
%       iter
     
      
      Dxn = chol(X0_Cov_m,'lower');

      Xr(:,1) = X0_m_m + Dxn*randn(4,1);
      
      
      Cn = XN_Cov_m - C_N0/X0_Cov_m*C_N0';
     Dxn = chol(Cn,'lower');

      Xr(:,N+1) = XN_m_m + C_N0*X0_Cov_m\(Xr(:,1) - X0_m_m) + Dxn*randn(4,1);
          
      %----------
          
      for k=1:N-1
          
              
          %----------CN|k
          
          CNk = zeros(4,4);
          
          for ii=0:N-k-1
              
              CNk = CNk + F^ii*Q*(F^ii)';
              
          end
          
        %.....
        Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;
        DG = chol(Gk,'lower');
        
        Gk_N = Gk*F^(N-k)'/(CNk);
        Gk_km1 = F - Gk_N*F^(N-k+1);
 
           %==============================
          
          Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N+1) + DG*randn(4,1);
          
       %--------------------

      end
      

  
figure(1)
hold on
plot(Xr(1,:),Xr(3,:),'b')
grid


kk=0:1:100;

figure(2)
hold on
plot(kk,Xr(2,:),'b')
grid

figure(3)
hold on
plot(kk,Xr(4,:),'b')
grid


  end  
      

  
  
  
  