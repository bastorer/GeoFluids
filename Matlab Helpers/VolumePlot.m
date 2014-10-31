function VolumePlot(X,Y,Z,D,varargin)
%% VolumePlot creates a collection of 3D isosurfaces.
% X, Y, and Z correspond to the grids produced by meshgrid
% D is the 3D matrix for the scalar field to be displayed
% The optional output is a list of percentiles. For each percentile
%     an isosurface shell is created. e.g. [0.2 0.8] would plot the
%     20th and 80th percentiles. 
%   The default is [0.25, 0.5, 0.75].
%   If the data has both positive and negative value, then the
%     percentile value is based on the maximum absolute value.
%     A shell is made for both signs of the percentile value.


if isempty(varargin)
    p = [0.25, 0.5, 0.75];
elseif length(varargin) == 1
    p = varargin{1};
end

% fh = figure;

cm = jet(201);

M = max(abs(D(:)));

%% Plot the faces
hold on
for ii = 1:length(p)
    surf = patch(isosurface(X,Y,Z,D,p(ii)*M));
    set(surf, 'FaceAlpha', 1/(length(p)+1),...
        'FaceColor', cm(round(101+100*p(ii)),:),...
        'EdgeAlpha',0)
    
    verts = surf.vertices;
    faces = surf.faces;
    
    
    
    surf = patch(isosurface(X,Y,Z,D,-p(ii)*M));
    set(surf, 'FaceAlpha', 1/(length(p)+1),...
        'FaceColor', cm(round(101-100*p(ii)),:),...
        'EdgeAlpha',0)
end

%% Plot the edges


hold off

end
