

% This code generates trajectories (Markov trajectories) without
% destination

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
 

Q = [q*T^3/3 q*T^2/2 0 0;...
    q*T^2/2 q*T 0 0;...
    0 0 q*T^3/3 q*T^2/2;...
    0 0 q*T^2/2 q*T];

 
D = chol(Q,'lower');
  
  
  %=================================
  

X0_m_m = [2000;70;5000;10];

X0_Cov_m = [1000 40 0 0;
            40 10 0 0;
            0 0 1000 40;
            0 0 40 10];

  %------------------------

XN_m_m = F^(N)*X0_m_m;

CN = zeros(4,4);
for ii=0:N-1

    CN = CN + F^ii*Q*(F^ii)';

end

XN_Cov_m = CN + F^(N)*X0_Cov_m*F^(N)';

   
  %------------------------------
  
  iteration = 50;
  
      
  for iter=1:iteration
      
      iter
      
      Xr = zeros(4,N+1);
      
       %--------------------
      
      Dx = chol(X0_Cov_m,'lower');
      Xr(:,1) = X0_m_m + Dx*randn(4,1);
      
      for k=2:N+1
          
          w = randn(4,1);
          Xr(:,k) = F*Xr(:,k-1) + D*w;
          
      end

  
figure(1)
hold on
plot(Xr(1,:),Xr(3,:),'--r')
grid
  

kk=0:1:100;

figure(2)
hold on
plot(kk,Xr(2,:),'--r')
grid

figure(3)
hold on
plot(kk,Xr(4,:),'--r')
grid


  end  


  
  
