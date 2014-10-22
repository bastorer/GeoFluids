function proj_matr = ProjectionMatrix(y1, y2)
%% ProjectionMatrix creates a matrix A such that
%    if Y1 = fn(y1) and Y2 = fn(y2) are two vectors
%    representing functions evaluated at the grids 
%    y1 and y2, then A*Y1 = Y2
%  The inputs y1 and y2 must be vectors.

  proj_matr = zeros(length(y2),length(y1));
    
  for ii = 1:length(y1)
      vec = zeros(size(y1));
      vec(ii) = 1;
      proj_matr(:,ii) = interp1(y1,vec,y2,'spline');
  end

end
