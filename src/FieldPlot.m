function fig = FieldPlot(mesh,solution,titel,fieldTitle)

fig = figure('units','normalized','outerposition',[0 0 1 1])
title(titel);
patch('Faces',mesh.connectivity,'Vertices',mesh.nodes,'FaceVertexCData',...
    solution,'FaceColor','interp','EdgeColor','k')
h = colorbar('FontSize',22,'TickLabelInterpreter','latex')%,'Interpreter',...
    %'latex');
h.Label.Interpreter = 'latex';
h.Label.String = fieldTitle;
alpha(0.95);
colormap(jet(256));
end

