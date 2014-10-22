function VolumePlot(X,Y,Z,D,varargin)

if isempty(varargin)
    p = [0.25, 0.5, 0.75];
elseif length(varargin) == 1
    p = varargin{1};
end

% fh = figure;

cm = darkjet(201);

M = max(abs(D(:)));

hold on
for ii = 1:length(p)
    surf = patch(isosurface(X,Y,Z,D,p(ii)*M));
    set(surf, 'FaceAlpha', 1/(length(p)+1),...
        'FaceColor', cm(round(101+100*p(ii)),:),...
        'EdgeAlpha', 1/(length(p)+1))
    
    surf = patch(isosurface(X,Y,Z,D,-p(ii)*M));
    set(surf, 'FaceAlpha', 1/(length(p)+1),...
        'FaceColor', cm(round(101-100*p(ii)),:),...
        'EdgeAlpha', 1/(length(p)+1))
end
hold off

end