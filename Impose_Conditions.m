function [As, P] = Impose_Conditions(As, bs)
% IMPOSE_CONDITIONS imposes arbitrary linear and homogenous conditions
% on the matrices given by As.
% > As is a cell array whose entries are NxN matrices.
% > bs is an MxN matrix, where each row corresponds to a condition of the form
%   bs(ii,:) /dot /vec{x} = 0.
% > M must be strictly less than N.
% > A warning message will be printed for any redundant conditions.
% > Output As contains, in the same order as input As, the input matrices
%   with the conditions of bs applied.
% > P is the recovery matrix. That is, if [A2,P] = Impose_Conditions({A},bs),
%   and [eigvec, eigval] = eig(A2), then P*eigvec re-introduces any points that
%   may have been removed in order to impose the conditions.

  %% Step 1: Make sure we don't have too many conditions
  if size(bs,1) > size(bs,2)
      error('Impose_Conditions: System is over-defined.')
  elseif size(bs,1) == size(bs,2)
      error(['Impose_Conditions: System is critically defined and can' ...
            'be found uniquely.'])
  end
  
  %% Step 2: Row reduce the problem to make sure that we don't have
  %          linearly dependent conditions.
  
  num_conds = size(bs,1);
  bs = rref(bs);
  tmp = sum(abs(bs),2);
  bs = bs(tmp~=0,:);
  
  if num_conds > size(bs,1)
      null_conds = num_conds - size(bs,1);
      num_conds = size(bs,1);
      warning(sprintf(['Impose_Conditions: %d Conditions were ' ...
                           'linear combinations of the other conditions.',...
                      ' There are %d remaining conditions.'],...
                          null_conds, num_conds))
  end
  
  %% Step 3: Impose the conditions
  
  % Initialize the transformation matrix P to recovered the deleted
  % portions.
  P = speye(size(As{1}));
  
  counter = 0;
  while size(bs,1) > 0
      counter = counter + 1;
      %[m,index] = max(abs(bs(1,:)));
      index = find(bs(1,:), 1, 'first');
      m = bs(1,index);
      
      if m == 0 || isempty(index)
        % If the maximum value was zero, or no non-zero entries were found, 
        % then one of our conditions is a linear combination of the others.
        % Step 2 above should have removed the necessity for this case,
        % but I've kept for the sake of being cautious.
          
          warning(sprintf(['Impose_Conditions: Conditions %d is a ' ...
                           'linear combination of the other conditions.'],...
                          counter))
          
          if size(bs,1) > 1
              bs = bs(2:end,:);
          else
              bs = [];
          end
        
      else
      
          b = bs(1,:); b(index) = [];
          b = -b/bs(1,index);
      
          % Modify the matrices of As
          for ii = 1:length(As)
              A = As{ii};
              
              aA = A(:,index); aA(index) = [];
              
              A(index,:) = [];
              A(:,index) = [];
              
              A = A + aA*b;
              
              As{ii} = A;
          end
          
          % Modify the Transformation Matrix P
          aP =  P(:,index);
          P(:,index) = [];
          P  = P  + aP*b;
          
          
          % Modify the remaining conditions to account
          % for the curent condition.
          % Note: Since the condition matrix bs was put
          % into RREF, this isn't necessary and so was
          % removed.
          % ab = bs(:,index);
          % bs = bs + ab*b;
          
          
          if size(bs,1) > 1
              bs = bs(2:end,:);
              bs(:,index) = [];
          else
              bs = [];
          end
      end
  end
end
