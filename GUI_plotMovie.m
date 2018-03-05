function GUI_plotMovie(pp, plot_h, EzPlot_h)

xx = pp.dimensionX;
yy = pp.dimensionY;

x0 = pp.targetPointY;
y0 = pp.targetPointX;

filename = 'output.bin';
fileID = fopen(filename, 'rb');

dim_x = fread(fileID, 1, 'int32');
dim_y = fread(fileID, 1, 'int32');
time_step = fread(fileID, 1, 'int32');

mat0 = fread(fileID, [dim_x*dim_y time_step], 'double');

maxEz = max(mat0(:));

line_h = animatedline(EzPlot_h, 'Marker','.');
Ez_time = zeros(1,time_step);
EzPlot_h.YLim = [-maxEz, +maxEz];
caxis(plot_h, [-maxEz, +maxEz]);

for i=1:time_step

    mat = reshape(mat0(:, i), [dim_x, dim_y]);  
    
    imagesc(plot_h, mat);
    
    colorbar(plot_h);
    axis(plot_h, 'image');
    plot_h.Title.String = sprintf('FDTD grid layout @ time = %d dt', i);
    plot_h.XLabel.String = 'x / grid index';
    plot_h.YLabel.String = 'y / grid index';
    drawnow limitrate;

    Ez_time(i) = mat(x0, y0);
    addpoints(line_h,i,Ez_time(i));

end

fclose(fileID)


%{
mat = zeros(xx, yy);

if pp.PMLFlag
    dd = pp.PMLthickness;

    xx = pp.dimensionX;
    yy = pp.dimensionY;
    mat(1:dd, :) = -1;
    mat(xx-dd+1:xx, :) = -1;
    mat(:, 1:dd) = -1;
    mat(:, yy-dd+1:yy) = -1;
end

if pp.PECScattFlag
    xx = pp.PECScattLocationX; lx = pp.PECScattDimensionX;
    yy = pp.PECScattLocationY; ly = pp.PECScattDimensionY;

    mat(xx:xx+lx-1, yy:yy+ly-1) = +1;
end

cla(plot_h,'reset');
plot_h.NextPlot = 'add';
xx = pp.dimensionX;
yy = pp.dimensionY;
imagesc(plot_h, 1:xx, 1:yy, mat, [-1, +1]);
axis(plot_h, 'image');
title('')

xx = pp.sourceLocationX;
yy = pp.sourceLocationY;
plot(plot_h, xx, yy, 'ro');
text(plot_h, xx+2, yy, 'source');

xx = pp.targetPointY;
yy = pp.targetPointX;
plot(plot_h, xx, yy, 'bx');
text(plot_h, xx+2, yy, 'observation point');

% https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#property_d119e54679
plot_h.Title.String = 'FDTD grid layout';
plot_h.XLabel.String = 'x / grid index';
plot_h.YLabel.String = 'y / grid index';
%}

end
