function [hList, xList, yList, shapeList, calAList] = runMesoCellExtension(NV,NPINS,calA0Init,Kl,Kb,cL,cB,plotIt,varargin)
%% FUNCTION to run single-cell extension code to test cell shape
% Code takes in cell parameters and extends equally-separated vertices
% to a set maximum extension, records shape data along the way
% 
% see README on github page for function information

% always Ka = 1
Ka = 1.0;

% -- Check Inputs

% check NV
if NV < 3
    error('runMesoCellExtension:inputNV','\nINPUT NV is < 3, cannot construct DPM with this input. Ending.');
else
    fprintf('** STARTING SINGLE CELL EXTENSION SIM WITH NV = %d',NV);
end

% check NPINS
if NPINS > NV
    error('runMesoCellExtension:inputNPINS','\nINPUT NPINS is > NV, cannot have more pins that vertices. Ending');
else
    fprintf(', NPINS = %d',NPINS);
end

% check mechanical parameters
if Kl < 0
    error('runMesoCellExtension:inputKl','\nINPUT Kl is < 0, energy would be incorrect. Ending.');
else
    fprintf(', Kl = %0.4g',Kl);
end

if Kb < 0
    error('runMesoCellExtension:inputKb','\nINPUT Kb is < 0, energy would be incorrect. Ending.');
else
    fprintf(', Kb = %0.4g',Kb);
end

% -- Plot information
if plotIt == 1
    % vertex color
    clr = [0 0.3 1];
end


% --  Check for movie string
makeAMovie = 0;
if ~isempty(varargin)
    if nargin == 9
        movieFileStr = varargin{1};
        if ischar(movieFileStr)
            fprintf('.\n** Also saving animation to location %s, setting plotIt -> 1.\n** Beginning simulation!',movieFileStr);
            plotIt      = 1;
            clr         = [0 0.3 1];
            makeAMovie  = 1;
            vobj        = VideoWriter(movieFileStr,'MPEG-4');
            open(vobj);
        else
            error('runMesoCellExtension:movieFileStringWrong','\n8th argument must be string for movie file name. Ending');
        end
    else
        error('runMesoCellExtension:tooManyArgs','\nToo many additional arguments. Must either have 7 or 8 input, last must be string. Ending');
    end
else
    fprintf('.\n** Beginning simulation!\n');
end

%% Initialization

% simulation parameters
phi0    = 0.1;                              % cell packing fraction
calA0   = calA0Init*NV*tan(pi/NV)/pi;       % initial shape parameter (1.01 * minimum)
dt0     = 0.001;                            % time step magnitude
a0      = 1.0;                              % initial preferred area
Ftol    = 1e-8;                             % force tolerance for minimization

% seed random number generator
rng('shuffle');

% initial coordinates of deformable polygon
x = zeros(NV,1);
y = zeros(NV,1);
r = zeros(NV,1);

% polygon vertex-based radii
r0 = sqrt((2.0*a0)/(NV*sin((2.0*pi)/NV)));
l0 = 2.0*sqrt(pi*calA0*a0)/NV;
th0 = ((2.0*pi)/NV).*ones(NV,1);

% vertex positions
L = sqrt(a0/phi0);
for vv = 1:NV
    x(vv) = r0*cos(2.0*pi*vv/NV) + 0.5*L + 0.01*l0*randn;
    y(vv) = r0*sin(2.0*pi*vv/NV) + 0.5*L + 0.01*l0*randn;
    r(vv) = 0.5*l0;
end

% vertex indexing
ip1 = [2:NV 1];

% PIN INFO
pstep   = round(NV/NPINS);
pinds   = 1:pstep:NV;
if length(pinds) <= NPINS
    NPINS   = length(pinds);
elseif length(pinds) > NPINS
    pinds   = pinds(1:NPINS);
end

% use FIRE to relax coordinates to energy minimum without pins
[x, y, ~] = vertexFIREPinnedVerts(x,y,[],a0,l0,th0,Ka,Kl,Kb,dt0,Ftol);

% draw initial cell
calA = computeShapeParam(x,y);
if plotIt == 1
    figure(1), clf, hold on, box on;
    for vv = 1:NV
        xp = x(vv) - r(vv);
        yp = y(vv) - r(vv);
        wp = 2.0*r(vv);
        rectangle('Position',[xp yp wp wp],'Curvature',[1 1],'EdgeColor','k','FaceColor','b');
    end
    axis equal;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    ax.XLim = L.*[-0.25 1.25];
    ax.YLim = ax.XLim;
    title(['Initial condition, $\mathcal{A} = ' sprintf('%0.5g',calA) '$'],'Interpreter','latex','fontsize',18);
end

% save initial locations
x0 = x;
y0 = y;

%% Extension sim setup

% overlap list
hStep       = -0.02;
hI          = 1.0;
hF          = -1.5;
hList       = hI:hStep:hF;
NSTEPS      = length(hList);

% -- get initial pin indices, loop will move radially outward

% choose vertex positions based on radial location from centroid
cx = mean(x);
cy = mean(y);

rx = x - cx;
ry = y - cy;
rs = sqrt(rx.*rx + ry.*ry);

ux = rx./rs;
uy = ry./rs;

% save data
xList = cell(NSTEPS,1);
yList = cell(NSTEPS,1);
shapeList = cell(NSTEPS,3);
calAList = zeros(NSTEPS,1);

% Loop over pinned positions, plot / animate
for hh = 1:NSTEPS
    % overlap value
    h = hList(hh);
    
    % -- Assign next pin locations
    x(pinds) = x0(pinds) + (1 - h)*l0*ux(pinds);
    y(pinds) = y0(pinds) + (1 - h)*l0*uy(pinds);
    
    % relax particle shape based on external force
    [x, y, thi] = vertexFIREPinnedVerts(x,y,pinds,a0,l0,th0,Ka,Kl,Kb,dt0,Ftol);
    thi = -thi;
    
    % age perimeter and bending
    lx = x(ip1) - x;
    ly = y(ip1) - y;
    l = sqrt(lx.^2 + ly.^2);
    lmean = mean(l);
    
    % lag toward preferred
    l0 = l0 + cL*(lmean - l0);
    th0 = th0 + cB*(thi - th0);
    
    % save state data
    xList{hh}           = x;
    yList{hh}           = y;
    shapeList{hh,1}     = calA0;
    shapeList{hh,2}     = Kb;
    shapeList{hh,3}     = th0;
    
    calA                = computeShapeParam(x,y);
    calAList(hh)        = calA;
    
    % print to console
    fprintf('-- On cumulative h = %0.4g, calA0 = %0.4g, calA = %0.4g\n',h,calA0,calA);
   
    % draw overlap of initial and final
    if plotIt == 1
        figure(2), clf, hold on, box on;

        for vv = 1:NV
            xv = x(vv) - r(vv);
            yv = y(vv) - r(vv);
            wv = 2.0*r(vv);
            if sum(vv == pinds) == 0
                rectangle('Position',[xv, yv, wv, wv],'Curvature',[1 1],'EdgeColor','k','FaceColor',clr);
            else
                rectangle('Position',[xv, yv, wv, wv],'Curvature',[1 1],'EdgeColor','k','FaceColor',1 - clr);
            end

            xv = x0(vv) - r(vv);
            yv = y0(vv) - r(vv);
            rectangle('Position',[xv, yv, wv, wv],'Curvature',[1 1],'EdgeColor',clr);
        end

        axis equal;
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.XLim = [-0.25*L 1.25*L];
        ax.YLim = [-0.25*L 1.25*L];

        tstr = ['$h = ' sprintf('%0.4g',h) '$, $\mathcal{A}_0 = ' sprintf('%0.4g',calA0) '$, $\mathcal{A} = ' sprintf('%0.4g',calA) '$'];
        title(tstr,'Interpreter','latex','FontSize',18);
        
        % make movie if writing
        if makeAMovie == 1
            set(gcf,'color','white');
            currframe = getframe(gcf);
            writeVideo(vobj,currframe);
        end
    end
end

% close video if writing
if makeAMovie == 1
    close(vobj);
end


end




function calA = computeShapeParam(x,y)
%% Small function to compute shape parameter

% perimeter
lx = x([2:end 1]) - x;
ly = y([2:end 1]) - y;
l = sqrt(lx.^2 + ly.^2);
p = sum(l);

% area
a = polyarea(x,y);

% shape parameter
calA = (p*p)/(4.0*pi*a);

end