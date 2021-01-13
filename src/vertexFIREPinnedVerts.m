function [x, y] = vertexFIREPinnedVerts(x,y,pinds,a0,l0,Kl,Kb,dt0,Ftol)
%% Function to relax shape energy with pinned vertices

% get number of vertices (check inputs)
NVx = size(x,1);
NVy = size(y,1);
if (NVx ~= NVy)
    fprintf('lengths of vx and vy do not match in dpVDOS function...\n');
    fprintf('size(vx) = %d %d, size(vy) = %d %d, ending\n',size(x,1),size(x,2),size(y,1),size(y,2));
    error('size of inputs is incorrect');
else
    NV = NVx;
end

% indexing
im1 = [NV 1:NV-1];
ip1 = [2:NV 1];

% force parameters
rho0            = sqrt(a0);                 % units: length
fa              = 1/rho0;                   % units: inv length, because of grad a
fl              = Kl*(rho0/l0);             % units: dim. less
fb              = Kb*(rho0/(l0*l0));        % units: inv length, because of length in force expression

% initialize velocites to zero
vx = zeros(NV,1);
vy = zeros(NV,1);

% initialize forces
fx = zeros(NV,1);
fy = zeros(NV,1);

fxold = zeros(NV,1);
fyold = zeros(NV,1);

% FIRE VARIABLES (hard code in)
alpha0      = 0.1;
finc        = 1.01;
fdec        = 0.5;
falpha      = 0.99;
dtmax       = 10*dt0;
dtmin       = 1e-8*dt0;
dt          = dt0;
NMIN        = 20;
NNEGMAX     = 5000;
NDELAY      = 100;
npPos       = 0;
npNeg       = 0;
npPMin      = 0;
alpha       = alpha0;
plotskip    = 2e4;

% force check
fcheck = 10*Ftol;
kcheck = 0;

% get dynamic indices
di = false(NV,1);
for ii = 1:NV
    if sum(ii == pinds) == 0
        di(ii) = true;
    end
end
NVMOB = NV - length(pinds);

% USE FIRE to relax forces
it = 0;
itmax = 5e6;
% loop over FIRE protocol while fcheck and kcheck are above minimum
while((fcheck > Ftol || npPMin < NMIN) && it < itmax)
    % update iterate
    it = it + 1;
    
    % Step 1. calculate P, fnorm, vnorm
    P = sum(fx(di).*vx(di)) + sum(fy(di).*vy(di));
    vnorm = sqrt(sum(vx(di).*vx(di)) + sum(vy(di).*vy(di)));
    fnorm = sqrt(sum(fx(di).*fx(di)) + sum(fy(di).*fy(di)));
    
    % plot FIRE information 
    if mod(it,plotskip) == 0 || it == 1
        fprintf('\nOn FIRE (pinned vertices) step %d\n',it);
        fprintf('\t ** F = %0.5g\n',fcheck);
        fprintf('\t ** K = %0.5g\n',kcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** P = %0.5g\n',P);
        fprintf('\t ** Pdir = %0.5g\n',P/(vnorm*fnorm));
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** vnorm = %0.5g\n',vnorm);
        fprintf('\t ** fnorm = %0.5g\n\n',fnorm);
    end
    
    % Step 2. adjust simulation based on net motion of system
    if (P > 0)
        % increase positive counter
        npPos = npPos + 1;
        
        % reset negative counter
        npNeg = 0;
        
        % alter simulation if enough positive steps have been taken
        if (npPos > NMIN)
            % change time step
            if (dt*finc < dtmax)
                dt = dt*finc;
            end
        end
        
        % decrease alpha
        alpha = alpha*falpha;
    else
        % reset positive counter
        npPos = 0;
        
        % increase negative counter
        npNeg = npNeg + 1;
        
        % check for stuck simulation
        if (npNeg > NNEGMAX)
            fprintf('Simulation negative for too long, ending program here.\n')
            error('FIRE did not converge');
        end
        
        % decrease time step if past initial delay
        if (it > NMIN)
            % decrease time step
            if (dt*fdec > dtmin)
                dt = dt*fdec;
            end
            
            % reset alpha
            alpha = alpha0;
        end
        
        % take a half step backwards
        x(di) = x(di) - 0.5*dt*vx(di) - 0.25*dt*dt*fxold(di);
        y(di) = y(di) - 0.5*dt*vy(di) - 0.25*dt*dt*fyold(di);
        
        % reset velocities to 0
        vx = zeros(NV,1);
        vy = zeros(NV,1);
    end
    
    % update velocities if forces are acting
    if fnorm > 0
        vx(di) = (1 - alpha).*vx(di) + alpha.*(fx(di)./fnorm)*vnorm;
        vy(di) = (1 - alpha).*vy(di) + alpha.*(fy(di)./fnorm)*vnorm;
    end
    
    % do first verlet update for vertices (assume unit mass)
    x(di) = x(di) + dt*vx(di) + 0.5*dt*dt*fxold(di);
    y(di) = y(di) + dt*vy(di) + 0.5*dt*dt*fyold(di);
    
    % update area
    a = polyarea(x, y);
    
    % * * * * * * * * * * * * * * * * * *
    % calculate forces based on positions
    % * * * * * * * * * * * * * * * * * *
    
    % reset forces
    fx = zeros(NV,1);
    fy = zeros(NV,1);
    
    % -- perimeter force
    
    % segment vectors
    lvx = x(ip1) - x;
    lvy = y(ip1) - y;
    
    % update perimeter segment lengths
    l = sqrt(lvx.^2 + lvy.^2);
    
    % segment unit vectos
    ulvx = lvx./l;
    ulvy = lvy./l;
    
    % segment strain
    dli = (l./l0) - 1.0;
    dlim1 = (l(im1)./l0) - 1.0;
    
    % perimeter force at this iteration
    flx = fl*(dli.*ulvx - dlim1.*ulvx(im1));
    fly = fl*(dli.*ulvy - dlim1.*ulvy(im1));
    
    % add to total force
    fx(di) = fx(di) + flx(di);
    fy(di) = fy(di) + fly(di);
    
    
    % -- area force
    areaStrain = (a/a0) - 1.0;
    
    % area force at this iteration
    fax = fa*0.5*areaStrain.*(y(im1) - y(ip1));
    fay = fa*0.5*areaStrain.*(x(ip1) - x(im1));
    
    % add to total force
    fx(di) = fx(di) + fax(di);
    fy(di) = fy(di) + fay(di);
    
    % -- bending force
    
    % s vectors
    six = lvx - lvx(im1);
    siy = lvy - lvy(im1);
    
    % bending force at this iteration
    fbx = fb*(2.0*six - six(im1) - six(ip1));
    fby = fb*(2.0*siy - siy(im1) - siy(ip1));
    
    % add to force
    fx(di) = fx(di) + fbx(di);
    fy(di) = fy(di) + fby(di);
    
    % do second verlet update for vertices
    vx(di) = vx(di) + dt*0.5*(fx(di) + fxold(di));
    vy(di) = vy(di) + dt*0.5*(fy(di) + fyold(di));
    
    fxold = fx;
    fyold = fy;
    
    % update Kcheck and Fcheck
    fcheck = sqrt(sum(fx(di).^2 + fy(di).^2))/NVMOB;
    kcheck = 0.5*sum(vx(di).^2 + vy(di).^2)/NVMOB;
    if fcheck < Ftol && it > NDELAY
        npPMin = npPMin + 1;
    else
        npPMin = 0;
    end
end
if (it == itmax)
    error('     ERROR: FIRE did not converge in itmax iterations, ending.\n');
else
    fprintf('\t FIRE converged, F = %0.5g, K = %0.5g\n',fcheck,kcheck);
end

end