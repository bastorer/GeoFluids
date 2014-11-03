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
    
    surf = isosurface(X,Y,Z,D,p(ii)*M);
    verts = surf.vertices;
    faces = surf.faces;

    patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', [0,0,1],...
        'FaceColor', cm(round(101+100*p(ii)),:), 'EdgeColor', [0,0,0],...
        'EdgeAlpha', 0, 'FaceAlpha', 1/(length(p)+1))
    
    PlotEdges(X,Y,Z,D,0.5*cm(round(101+100*p(ii)),:),1/(length(p)+1),p(ii)*M)

    surf = isosurface(X,Y,Z,D,-p(ii)*M);
    verts = surf.vertices;
    faces = surf.faces;
    
    patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', [0,0,1],...
        'FaceColor', cm(round(101-100*p(ii)),:), 'EdgeColor', [0,0,0],...
        'EdgeAlpha', 0, 'FaceAlpha', 1/(length(p)+1))
    
    PlotEdges(X,Y,Z,D,0.5*cm(round(101-100*p(ii)),:),1/(length(p)+1),-p(ii)*M)

end

hold off
box on

end

function PlotEdges(x,y,z,D,col,alp,val)
    
    for ii = linspace(1, size(x,2), 20)

        yS = squeeze(y(:,round(ii),1));
        zS = squeeze(z(1,round(ii),:));
        dS = squeeze(D(:,round(ii),:));

        pts = contourc(yS,zS,dS,[val,val]);

        yS = pts(1,:);
        zS = pts(2,:);
        
        while ~isempty(zS)
            num_pts = zS(1);
            Y = yS(2:1+num_pts);
            Z = zS(2:1+num_pts);
            X = x(1,round(ii),1)*ones(size(Y));
            
            patchline(X,Y,Z,'EdgeColor',col,'EdgeAlpha',alp)
            
            zS = zS(2+num_pts:end);
            yS = yS(2+num_pts:end);
            
        end

    end

    for ii = linspace(1, size(x,1), 20)

        xS = squeeze(x(round(ii),:,1));
        zS = squeeze(z(round(ii),1,:));
        dS = squeeze(D(round(ii),:,:));

        pts = contourc(xS,zS,dS,[val,val]);

        xS = pts(1,:);
        zS = pts(2,:);
        
        while ~isempty(zS)
            num_pts = zS(1);
            X = xS(2:1+num_pts);
            Y = y(round(ii),1,1)*ones(size(X));
            Z = zS(2:1+num_pts);
            
            patchline(X,Y,Z,'EdgeColor',col,'EdgeAlpha',alp)
            
            zS = zS(2+num_pts:end);
            xS = xS(2+num_pts:end);
            
        end

    end
    
    for ii = linspace(1, size(x,3), 20)

        xS = squeeze(x(1,:,round(ii)));
        yS = squeeze(y(:,1,round(ii)));
        dS = squeeze(D(:,:,round(ii)));

        pts = contourc(xS,yS,dS,[val,val]);

        xS = pts(1,:);
        yS = pts(2,:);
        
        while ~isempty(yS)
            num_pts = yS(1);
            X = xS(2:1+num_pts);
            Y = yS(2:1+num_pts);
            Z = z(1,1,round(ii))*ones(size(Y));
            
            patchline(X,Y,Z,'EdgeColor',col,'EdgeAlpha',alp)
            
            yS = yS(2+num_pts:end);
            xS = xS(2+num_pts:end);
            
        end

    end
    
end
