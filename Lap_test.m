%% A script to test the Impose_Conditions function
%  on a 1D and 2D Poisson problem.

%% 1D
% clear all;
% 
% Nx = 100;
% 
% x = linspace(0,1,Nx);
% 
% Dx = FiniteDiff(x, 4, false, true);
% 
% A = Dx^2;
% B = speye(size(A));
% 
% % bs = [Dx(x==0,:); Dx(x==1,:)]; % Neumann conditions
% % bs = [B(x==0,:); B(x==1,:)]; % Dirichlet conditions
% bs = [Dx(x==0,:); B(x==1,:)];
% 
% [As,P] = Impose_Conditions({A,B}, bs);
% A = As{1};
% B = As{2};
% 
% [vec, val] = eig(full(A),full(B));
% 
% vec = P*vec; % Re-introduce the removed points
% 
% [~,ind] = sort(diag(real(val)), 'descend');
% vec = vec(:,ind);
% 
% for ii = 1:10
%   plot(x,real(vec(:,ii)), '-b', x, imag(vec(:,ii)), '--r')
%   drawnow()
%   pause(1)
% end

%% 2D

clear all;

Nx = 300;
Ny = 300;

x = linspace(0,1,Nx);
y = linspace(0,sqrt(2),Ny);

[xx,yy] = meshgrid(x,y);

[Dx, Dy] = FiniteDiff({x,y}, 4, true, true);

A = Dx^2 + Dy^2;
B = speye(size(A));

% bs = [B(xx(:)==min(x),:); ...
%       B(xx(:)==max(x),:); ...
%       Dy(yy(:)==min(y),:); ...
%       Dy(yy(:)==max(y),:)];
% bs = [Dx(xx(:)==min(x),:); ...
%       Dx(xx(:)==max(x),:); ...
%       B(yy(:)==min(y),:); ...
%       B(yy(:)==max(y),:)];
% bs = [Dx(xx(:)==min(x),:); ...
%       Dx(xx(:)==max(x),:); ...
%       Dy(yy(:)==min(y),:); ...
%       Dy(yy(:)==max(y),:)]; % Neumann conditions
bs = [B(xx(:)==min(x),:); ...
      B(xx(:)==max(x),:); ...
      B(yy(:)==min(y),:); ...
      B(yy(:)==max(y),:)]; % Dirichlet conditions

% [~,P] = Impose_Conditions({}, bs);
P = Impose_Conditions2(bs);
Pinv = P\speye(size(A));
% A = As{1};
% B = As{2};
A = Pinv*A*P;
B = Pinv*B*P;

% [vec, val] = eig(full(A),full(B));
% [vec, val] = eigs(A,B,5,'lr');
% val = diag(val);
% 
% vec = P*vec; % Re-introduce the removed points
% 
% vec = vec(:,real(val)<0);
% val = val(real(val)<0);
% 
% [val,ind] = sort(real(val), 'descend');
% vec = vec(:,ind);
% 
% for ii = 1:5
%     pcolor(xx,yy,reshape(vec(:,ii),[Nx,Ny]))
%     shading flat
%     drawnow()
%     pause(1)
% end


