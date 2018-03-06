%
%  GUI_plotMovie.m
%  FDTD 2D TMz simulation
%
%  Created by Yifan Wang on 3/1/18.
%
function Ez_time = proj1_FDTDgradientPML_GUI(pp, plot_h, EzPlot_h)

c0=3e+8;
mu0=4*pi*1e-7;
epsilon0=(1/(36*pi))*1e-9;
eta0=sqrt(mu0/epsilon0);

dim_x = pp.dimensionX;
dim_y = pp.dimensionY;
source_x = pp.sourceLocationX;
source_y = pp.sourceLocationY;
total_time_steps = pp.time_step;
%Boundary width of PML in all directions
bound_width=pp.PMLthickness;


PEC_x0 = 0; PEC_x1 = 0;
PEC_y0 = 0; PEC_y1 = 0;

if pp.PECScattFlag
    PEC_x0 = pp.PECScattLocationX; PEC_x1 = pp.PECScattLocationX + pp.PECScattDimensionX -1;
    PEC_y0 = pp.PECScattLocationY; PEC_y1 = pp.PECScattLocationY + pp.PECScattDimensionY -1;
end

target_x = pp.targetPointX;
target_y = pp.targetPointY;

dx = 1e-6; dy = dx;
dt = dx/(sqrt(2)*c0); % (c0*dt)/dx < 1/sqrt(2), Courant, Courant-Friedrichs-Lewy_condition
Ez_pos_x = ([1, dim_x] - source_x)*dx;
Ez_pos_y = ([1, dim_y] - source_y)*dy;

lambda = 10*dx;
freq = c0/lambda;
period_dt = lambda/(c0*dt);

fprintf('FDTD Simulation on [%d,%d] x [%d,%d]\n', 1-source_x, dim_x-source_x, 1-source_y, dim_y-source_y);
fprintf('Wave length = %d dx, Period = %0.1f dt\n', lambda/dx, period_dt);

Ez = zeros(dim_x, dim_y);
Ezx = zeros(dim_x, dim_y);
Ezy = zeros(dim_x, dim_y);
Hx = zeros(dim_x, dim_y); % Hx(ii, jj) = Hx(ii+1/2, jj+1/2)
Hy = zeros(dim_x, dim_y);

% init permittivity and permeability matrix
epsilon=    epsilon0*ones(dim_x, dim_y);
mu     =    mu0     *ones(dim_x, dim_y);

% init electric conductivity matrices of x and y directions
sigma_e_x=zeros(dim_x, dim_y);
sigma_e_y=zeros(dim_x, dim_y);
sigma_m_x=zeros(dim_x, dim_y);
sigma_m_y=zeros(dim_x, dim_y);

if pp.PMLFlag
%polynomial grading sigma(x) = (x/d)^m * sigma_max ; x=0: 0, x=d: sigma_max
mm = 3;

%Required reflection co-efficient
R0_coeff=1e-6;
%{
R0(\phi) = exp( - 2*\eta*cos(\phi) \int_0^d \sigma_e(x) dx)
\sigma_e(x) = \sigma_emax * (x/m)^d
=> \int_0^d \sigma_e(x) dx = \sigma_emax * d / (m+1)

=> \sigma_emax = (-log( R0(\phi = 0) * (m+1) / (2*\eta*d)))
%}

sigma_max    =   ( -log(R0_coeff) * (mm+1) ) / (2*eta0*(bound_width*dx));

for i=0:bound_width-1
    xx0 = (i-0)*dx;
    xx1 = (i+1)*dx;
    sigma_e_value(i+1) = sigma_max/(bound_width*dx)^mm/(mm+1) * (xx1^(mm+1) - xx0^(mm+1)) /dx;

    xx0 = (i-0.5)*dx;
    xx1 = (i+0.5)*dx;
    if i ~= 0
        sigma_m_value(i+1) = (sigma_max*mu0/epsilon0)/(bound_width*dx)^mm/(mm+1) * (xx1^(mm+1) - xx0^(mm+1)) /dx; % HERE we assume everywhere we have mu0/epsilon0
    else
        sigma_m_value(i+1) = (sigma_max*mu0/epsilon0)/(bound_width*dx)^mm/(mm+1) * ( xx1^(mm+1) ) /dx;
    end
end

sigma_e_y(1:dim_x, bound_width        :-1:1    ) = repmat(sigma_e_value, [dim_x,1]);
sigma_e_y(1:dim_x, dim_y-bound_width+1: 1:dim_y) = repmat(sigma_e_value, [dim_x,1]);

sigma_e_x(bound_width        :-1:1   ,  1:dim_y) = repmat(sigma_e_value.', [1, dim_y]);
sigma_e_x(dim_x-bound_width+1: 1:dim_x, 1:dim_y) = repmat(sigma_e_value.', [1, dim_y]);

sigma_m_y(1:dim_x, bound_width+1      :-1:2    ) = repmat(sigma_m_value, [dim_x,1]);
sigma_m_y(1:dim_x, dim_y-bound_width+1: 1:dim_y) = repmat(sigma_m_value, [dim_x,1]);

sigma_m_x(bound_width+1      :-1:2   ,  1:dim_y) = repmat(sigma_m_value.', [1, dim_y]);
sigma_m_x(dim_x-bound_width+1: 1:dim_x, 1:dim_y) = repmat(sigma_m_value.', [1, dim_y]);

end

plotmovie = true;

% @ Ez
Ax = (epsilon - (0.5*dt)*sigma_e_x)./(epsilon + (0.5*dt)*sigma_e_x);
Ay = (epsilon - (0.5*dt)*sigma_e_y)./(epsilon + (0.5*dt)*sigma_e_y);

Bx = (dt/dy)./(epsilon + (0.5*dt)*sigma_e_x);
By = (dt/dx)./(epsilon + (0.5*dt)*sigma_e_y);

Cx = (mu - (0.5*dt)*sigma_m_y)./(mu + (0.5*dt)*sigma_m_y); % @ Hx , y-0.5
Cy = (mu - (0.5*dt)*sigma_m_x)./(mu + (0.5*dt)*sigma_m_x); % @ Hy , x-0.5

Dx = (dt/dy)./(mu + (0.5*dt)*sigma_m_y); % @ Hx
Dy = (dt/dx)./(mu + (0.5*dt)*sigma_m_x); % @ Hy

Ez_time = zeros(1, 1+total_time_steps);
line_h = animatedline(EzPlot_h, 'Marker','.');

cla(plot_h,'reset');
EzPlot_h.YLim = [-1, 1];
caxis(plot_h, [-1, 1]);

for time_i = 0:total_time_steps

    % : with lower operator precedence, i.e. 2:4-1 => [2, 3]
    % @ t = l, from t = l-0.5 -> l+0.5
    % @ (x,y) = (i, j) * dx
    xx = 2:dim_x-1;
    yy = 2:dim_y-1;
    xx_1 = 3:dim_x;
    yy_1 = 3:dim_y;


    Ezx(xx, yy) = Ax(xx, yy).*Ezx(xx, yy) - Bx(xx, yy).*(Hx(xx, yy_1) - Hx(xx, yy));
    Ezy(xx, yy) = Ay(xx, yy).*Ezy(xx, yy) + By(xx, yy).*(Hy(xx_1, yy) - Hy(xx, yy));

    % l = 0.5, t = 0;
    % hard source
    % Ez(source_x, source_y) = exp(1j*( 2*pi*freq*time_i*dt - pi/2));
    Ezx(source_x, source_y) =   0.5*sin( 2*pi*freq*time_i*dt       );
    Ezy(source_x, source_y) =   0.5*sin( 2*pi*freq*time_i*dt       );

    if pp.PECScattFlag
        Ezx(PEC_x0:PEC_x1, PEC_y0:PEC_y1) = 0;
        Ezy(PEC_x0:PEC_x1, PEC_y0:PEC_y1) = 0;
    end

    Ez = Ezx + Ezy;

    % @ t = l+0.5, from t = l -> l+1
    xx0 = 2:dim_x;
    xx0_0 = 1:dim_x-1;
    yy0 = 2:dim_y;
    yy0_0 = 1:dim_y-1;

    % @ (x,y) = (i, j-0.5) * dx
    Hx(xx, yy0) = Cx(xx, yy0).*Hx(xx, yy0) - Dx(xx, yy0).*(Ez(xx, yy0) - Ez(xx, yy0_0));

    % @ (x,y) = (i-0.5, j) * dx
    Hy(xx0, yy) = Cy(xx0, yy).*Hy(xx0, yy) + Dy(xx0, yy).*(Ez(xx0, yy) - Ez(xx0_0, yy));

    if plotmovie
        imagesc(plot_h, Ez_pos_x, Ez_pos_y, real(Ez.'));

        colorbar(plot_h);
        axis(plot_h, 'image');
        plot_h.Title.String = sprintf('FDTD grid layout @ time = %d dt', time_i);
        plot_h.XLabel.String = 'x / m';
        plot_h.YLabel.String = 'y / m';
        plot_h.YDir = 'normal';
        drawnow limitrate;

        Ez_time(time_i+1) = real(Ez(target_x, target_y));
        addpoints(line_h,time_i,Ez_time(time_i+1));
        drawnow limitrate;
        % Actions Equivalent to drawnow : https://www.mathworks.com/help/matlab/ref/drawnow.html#burd3gs-3
    end
end

end


