function [P] = Impose_Conditions2(bs)
% IMPOSE_CONDITIONS imposes arbitrary linear and homogenous conditions
% on the matrices given by As.
% > bs is an MxN matrix, where each row corresponds to a condition of the form
%   bs(ii,:) /dot /vec{x} = 0.
% > M must be strictly less than N.
% > A warning message will be printed for any redundant conditions.
% > Output As contains, in the same order as input As, the input matrices
%   with the conditions of bs applied.
% > P is the transfromation matrix. That is, if P = Impose_Conditions(bs),
%   then solving [eigvec, eigval] = eig(P\A*P) solves the constrained
%   problem. P*eigvec re-introduces any points that may have been removed 
%   in order to impose the constraints.

  %% Step 1: Check inputs    
  if size(bs,1) > size(bs,2)
      error('Impose_Conditions: System is over-defined.')
  elseif size(bs,1) == size(bs,2)
      error(['Impose_Conditions: System is critically defined and could' ...
            'be found uniquely.'])
  end
  
  %% Step 2: Row reduce the problem to make sure that we don't have
  %          linearly dependent conditions.
  
  num_conds = size(bs,1);
  num_elems = size(bs,2);
  
  bs = frref(bs);
  tmp = sum(abs(bs),2);
  bs = bs(tmp~=0,:);
  
  if num_conds > size(bs,1)
      null_conds = num_conds - size(bs,1);
      num_conds = size(bs,1);
      warning(['Impose_Conditions: %d Conditions were ' ...
                           'linear combinations of the other conditions.',...
                      ' There are %d conditions remaining.'],...
                          null_conds, num_conds)
  end
  
  %% Step 3: Impose the conditions
  
  % First, find the indices of the removed points
  removed_indices = zeros(num_conds,1);
  for ii = 1:num_conds
      removed_indices(ii) = find(bs(ii,:),1,'first');
  end
  [removed_indices, ind] = sort(removed_indices,'ascend');
  bs = bs(ind,:);
  
  % Now reduce the bs matrix
  bs(:,removed_indices) = [];
  cnt = sum(bs(:)~=0);
  meta_ind = 1:size(bs,2);
  
  % Initialize the transformation matrix P
  ind1 = zeros(cnt+num_elems-num_conds,1);
  ind2 = zeros(cnt+num_elems-num_conds,1);
  vals = zeros(cnt+num_elems-num_conds,1); 
  
  % Now find the index-value arrangement for the transformation matrix P
  ind = 1; % To keep track of how many removed points we've reconstructed.
  tracker = 1; % To keep track of where we are in ind1, ind2, and vals.
  for ii = 1:num_elems
      % We loop through, computing the reconstruction vector for each
      % point. That is, x_i in x_original = v .* x_reduced. Here we
      % determine v.
      if ii == removed_indices(ind)
          % If we're reconstructing one of the removed points...
%           tmp_ind = find(bs(ind,:));
          tmp_ind = meta_ind(bs(ind,:)~=0);
          len = length(tmp_ind); % How many other points are needed for the 
                                % reconstruction
          ind1(tracker:tracker+len-1) = ii;
          ind2(tracker:tracker+len-1) = tmp_ind;
          vals(tracker:tracker+len-1) = bs(ind,tmp_ind);
          tracker = tracker + len;
          ind = ind + 1;
          
      else
          % If we're reconstructing a point that wasn't removed, then it's
          % really easy.
          ind1(tracker) = ii; % We're returning the ii^th element
          ind2(tracker) = ii - ind + 1; % But in the reduced vector, it's 
                                        % (ii-ind+1)^st element
          vals(tracker) = 1;
          tracker = tracker + 1;
      end
  end
  
  P = sparse(ind1,ind2,vals,num_elems,num_elems-num_conds);
  
end