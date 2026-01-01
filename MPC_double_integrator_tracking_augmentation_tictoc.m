clc
clear
close all

% =========================
% Continuous time model
% =========================
Ac = [0,1;0,0];
Bc = [0;1];
Cc = [1,0;0,1];
Dc = [0;0];

% =========================
% Discrete-time model    
% =========================
Ts = 0.1;  % sampling period, sec
sysd = c2d( ss(Ac, Bc, Cc, Dc), Ts );

Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

% =========================
% Augmented model
% =========================
Ax = [Ad,Bd,0*Ad(:,1);...   
      0*Ad(1,:),1,0;...
      0*Ad(1,:),0,1];
Bx = [Bd;1;0];                 
Cx = [1,0,0,0;...
      0,1,0,0;...
      0,0,1,0;...
      0,0,0,1];
Dx = [0;0;0;0];

Ex = [Cd(1,:),0,-1]; % tracking error

model = LTISystem('A',Ax,'B',Bx,'C',Cx,'D',Dx,'Ts',Ts);

% =========================
% Constraints
% =========================
model.u.min = -1;
model.u.max =  1;

model.x.min = [-3, -3, -1, -3];
model.x.max = [ 3,  3,  1,  3];

% =========================
% Cost function
% =========================
model.u.penalty = QuadFunction(diag(0.1));
model.x.penalty = QuadFunction(Ex'*Ex);

% =========================
% Horizon
% =========================
N = 12;

% =========================
% Measure MPC construction time (tic/toc)
% =========================
tic
ctrl = MPCController(model, N);
t_build = toc;

fprintf('MPC controller build time: %.6f s\n', t_build);

% =========================
% Closed-loop simulation
% =========================
x0   = [0;0];
u    = 0;
r0   = -2;
Nsim = 200;

x = x0;

data.X     = zeros(length(x0), Nsim);
data.U     = zeros(size(Bd,2), Nsim);
data.R     = zeros(1, Nsim);
data.W     = zeros(1, Nsim);
data.Tcomp = zeros(1, Nsim);   % MPC computation time

for i = 1:Nsim

    % Reference
    if (i-1)*Ts < 1
        r = r0; 
        w = 0;
    elseif (i-1)*Ts < 10
        r = r0; 
        w = 0.0; 
    elseif (i-1)*Ts < 25
        r = r0;  
        w = 0;
    else
        r = r0;
        w = 0;
    end

    % =========================
    % MPC evaluation (tic/toc)
    % =========================
    z = [x(:); u; r];

    try
        tic
        [du, feasible, openloop] = ctrl.evaluate(z);
        data.Tcomp(i) = toc;
    catch
        du = 0;
        data.Tcomp(i) = NaN;
    end

    % =========================
    % System update
    % =========================
    u = u + du;

    data.X(:,i) = x;
    data.U(:,i) = u;
    data.R(:,i) = r;
    data.W(:,i) = w;

    x = Ad*x + Bd*(u + w);
end

% =========================
% Timing statistics
% =========================
validTimes = data.Tcomp(~isnan(data.Tcomp));

fprintf('Average MPC evaluation time: %.6f s\n', mean(validTimes));
% fprintf('Maximum MPC evaluation time: %.6f s\n', max(validTimes));
% fprintf('Minimum MPC evaluation time: %.6f s\n', min(validTimes));
fprintf('MPC evaluation time: %.6f s\n', max(validTimes));
fprintf('Total MPC controller sim time: %.6f s\n', t_build + mean(validTimes));

% =========================
% Plot responses
% =========================
figure;
set(gcf,'Color','white')

subplot(3,1,1)
plot(((1:Nsim)-1)*Ts, data.X(1,:), ((1:Nsim)-1)*Ts, data.R, 'k--','LineWidth',2);
ylabel('$x_1$','Interpreter','latex','FontSize',14)
grid on

subplot(3,1,2)
plot(((1:Nsim)-1)*Ts, data.X(2,:), 'LineWidth',2)
ylabel('$x_2$','Interpreter','latex','FontSize',14)
grid on

subplot(3,1,3)
plot(((1:Nsim)-1)*Ts, data.U(1,:), 'LineWidth',2)
ylabel('$u$','Interpreter','latex','FontSize',14)
grid on

