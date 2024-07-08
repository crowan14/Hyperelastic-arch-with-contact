close all; 

global N1;
global N2;

global L1;
global L2;
global L3;

global C1; global D1;

global pts;
global tx; global dxt; 
global cx; global dxc;
global r2x; global r2y; global dx2; global dy2;
global r3x; global r3y; global dx3; global dy3;
global r4x; global r4y; global dx4; global dy4;

global traction;

global pen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of neurons in hidden layer
N1 = 8; 

%geometric parameters
L1 = 1; 
L2 = .2;
L3 = 0.1;

%material parameters
C1 = 1; D1 = 1;

%contact penalty parameter
pen = 1E9;

%side length of integration element
sz = 5E-3;

%for plotting
pts = 250;

%define traction vector (applied on upper horizontal bar)
traction = @(x1,x2) ( 2E0*[ .5 ; -2 ] );

%defines location of traction boundary
tx_pts = fix(L1/sz);
tx = linspace(L1/3,2*L1/3,tx_pts);
dxt = tx(2) - tx(1);

%defines contact boundary
cx_pts = fix((L1-2*L3)/sz);
cx = linspace(L3,L1-L3,tx_pts);
dxc = cx(2) - cx(1);

%region 2, left vertical bar
r2x_pts = fix(L3/sz); r2y_pts = fix((L2-L3)/sz);
r2x = linspace(0,L3,r2x_pts); r2y = linspace(0,L2-L3,r2y_pts);
dx2 = r2x(2) - r2x(1); dy2 = r2y(2) - r2y(1);

%region 3, top horizontal bar
r3x_pts = fix(L1/sz); r3y_pts = fix((L3)/sz);
r3x = linspace(0,L1,r3x_pts); r3y = linspace(L2-L3,L2,r3y_pts);
dx3 = r3x(2) - r3x(1); dy3 = r3y(2) - r3y(1);

%region 4, right vertical bar
r4x_pts = fix(L3/sz); r4y_pts = fix((L2-L3)/sz);
r4x = linspace(L1-L3,L1,r4x_pts); r4y = linspace(0,L2-L3,r4y_pts);
dx4 = r4x(2) - r4x(1); dy4 = r4y(2) - r4y(1);

int_pts = r2x_pts*r2y_pts + r3x_pts*r3y_pts + r4x_pts*r4y_pts;

fprintf('%d integration points \n\n',int_pts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1 = yes , 0 = no

%run symbolic computations?
symbolic = 0;

%perform optimization?
opt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%symbolic computations
[num_params] = one_hidden_hyperelasticity(symbolic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of steps in fmincon
evals = 300;

%optimization
p0 = .1*randn(num_params,1);

options = optimoptions('fmincon', ...
    'OptimalityTolerance', 0, ...
    'GradObj' , 'on' , ...
    'StepTolerance', 0, ...
    'MaxFunctionEvaluations', evals,...
    'MaxIterations', evals);

if opt == 1
    param = fmincon( @energy , p0 , [] , [] , [] , [] , [] , [] , [] , options );
end

%plot deformed structure
u = @(x1,x2) ( N( [x1;x2] , param ) );

pgridx = linspace(0,L1,pts); pgridy = linspace(0,L2,pts);
pos1 = zeros(pts,pts); pos2 = zeros(pts,pts); 

[ X , Y ] = meshgrid( pgridx,pgridy );

scale = 1;

for i=1:pts
    for j=1:pts
        if or( or( X(i,j) < L3 , X(i,j) > (L1-L3) ) , Y(i,j) > (L2-L3) )
            pos =  [ X(i,j) ; Y(i,j) ] + scale*u( X(i,j) , Y(i,j) );
            pos1(i,j) = pos(1);
            pos2(i,j) = pos(2);
     
        end
    end
end

pos1(pos1==0) = nan; 
pos2(pos2==0) = nan;

x_pos = reshape( pos1 , [pts^2,1] );
y_pos = reshape( pos2 , [pts^2,1] );

figure(2)
scatter(x_pos,y_pos)
grid on
axis equal
title('Deformed Structure')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






