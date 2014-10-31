[x,y,z] = meshgrid(0:0.01:1, 0:0.01:1, 0:0.01:1);
C = sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

fc = isosurface(x,y,z,C,0.5);
verts = fc.vertices;
faces = fc.faces;

hold on

patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', [0,0,1],...
    'FaceColor', [0,0,0], 'EdgeColor', [0,0,0], 'EdgeAlpha', 0)


for x0 = 0.2
    test1 = verts(faces(:,1),1) == x0;
    test2 = verts(faces(:,2),1) == x0;
    test3 = verts(faces(:,3),1) == x0;
    
    T1 = test2 & test3;
    T2 = test1 & test3;
    T3 = test1 & test2;
    
    test = (test1 & test2) | (test1 & test3) | (test2 & test3);
    
    tmp = faces;
    tmp()
    
    face2 = tmp(test,:);
    
    
    

    patch('Vertices', verts, 'Faces', face2, 'FaceVertexCData', [0,0,0],...
        'FaceColor', [1,1,1], 'EdgeColor', [1,0,0], 'FaceAlpha', 0)
end

% for y0 = 0:0.05:1
%     face2 = faces(verts(faces(:,2),2)==y0,:);
% 
%     patch('Vertices', verts, 'Faces', face2, 'FaceVertexCData', [0,0,0],...
%         'FaceColor', [1,1,1], 'EdgeColor', [0,1,0], 'FaceAlpha', 0)
% end
% 
% for z0 = 0:0.05:1
%     face2 = faces(verts(faces(:,3),3)==z0,:);
% 
%     patch('Vertices', verts, 'Faces', face2, 'FaceVertexCData', [0,0,0],...
%         'FaceColor', [1,1,1], 'EdgeColor', [0,0,1], 'FaceAlpha', 0)
% end

hold off