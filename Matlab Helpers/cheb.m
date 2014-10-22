function [D,x] = cheb(N)
%% CHEB computes the Chebyshev differentiation
%    matrix on N+1 points (i.e. N intervals)
%  D = differentiation matrix
%  x = Chebyshev grid

if N == 0
    D = 0; 
    x = 1;
else
    x  = cos(pi*(0:N)/N)';
    c  = [2; ones(N-1,1); 2].*(-1).^(0:N)';
    X  = repmat(x,1,N+1);
    dX = X-X';
    D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
    D  = D - diag(sum(D'));                 % diagonal entries
end

end
