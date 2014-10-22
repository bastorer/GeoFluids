function proj_matr = ProjectionMatrix(y1, y2)

    proj_matr = zeros(length(y2),length(y1));
    
    for ii = 1:length(y1)
        vec = zeros(size(y1));
        vec(ii) = 1;
        proj_matr(:,ii) = interp1(y1,vec,y2,'spline');
    end

end
