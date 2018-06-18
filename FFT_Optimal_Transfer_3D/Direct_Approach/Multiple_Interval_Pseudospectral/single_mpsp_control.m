%% (C) Yajur Kumar, 2017. All rights reserved.

function U_k = single_mpsp_control(dFdX,dFkdUk,dYndXn,U0,k,n)
  
% Warning: n is time-step and may not be equal to 1 second.

  B0(n-1) = dYndXn;
  B(n-1)  = B0(n-1)*dFdU(n-1);
  
  for k=(n-2):-1:1
      B0(k) = B0(k+1)*(dFdX(i+1));
      B(k)  = B0(k)*dFdU(k);
  end
  
  A_lambda = 0;
  b_lambda = 0;
  
  R = 1;
  
  for k=1:1:(n-1)
      A_lambda = A_lambda - B(k)*inv(R)*B(k)';
      b_lambda = b_lambda + B(k)*U0(k);
  end
  
  if det(A_lambda)==0
      error('A_lambda is singular. ERROR!');
  end
  
  lambda = inv(A_lambda)*(b_lambda-dY(n));
   
  U_k   = -inv(R)*B(k)'*lambda;
