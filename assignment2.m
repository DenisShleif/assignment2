function [] = assignment2()
clear
clc
close all

board.h = 0.05; %mm
board.L = 3; %mm
board.W = 2; %mm
board.Vo = 5; %V
board.Lb = 1; %mm
board.Wb = 0.75; %mm
board.sigBox = 1e-2; %S
board.sigOut = 1; %s
board = initializeBoard(board,true);

question1A(board);
question1b(board);
question1bTheory(board);
question2a(board,true);
board.h = 0.1; %mm
board = initializeBoard(board);
question2b(board);
question2c(board);
question2d(board);

end

function [] = question2b(board)
N = 100;
x = linspace(0.05,1,N);
y = zeros(1,N);
for n = 1:N
    board.h = x(n);
    board = initializeBoard(board);    
    y(n) = question2a(board);
end
h = figure;
plot(x,y,'*-');
xlabel('Mesh Size (mm)');
ylabel('Current (I)');
title('Question 2b: Current vs Mesh Size');
grid on;
saveas(h,'2b.png');

end

function [] = question2c(board)
N = 100;
x = linspace(0,board.W,N);
y = zeros(1,N);
for n = 1:N
    board.Wb = (board.W - x(n))/2;
    board = initializeBoard(board);    
    y(n) = question2a(board);
end
h = figure;
plot(x,y,'*-');
xlabel('Bottle Neck Width (mm)');
ylabel('Current (I)');
title('Question 2c: Current vs Bottle Neck Width');
grid on;
saveas(h,'2c.png');
end

function [] = question2d(board)
N = 100000;
x = zeros(1,N);
y = zeros(1,N);

currentPercent = 0;
tic;
for n = 1:N
    board.sigBox = generateRand(0,1);
    board.Wb = generateRand(0,board.W/2);
    board.Lb = generateRand(0,board.L);
    board = initializeBoard(board); 
    x(n) = board.sigAv;
    y(n) = question2a(board);
    
    tmp = floor(n*100/N);
    if tmp > currentPercent
        currDuration = toc;
        ETR = currDuration*(N-n)/n;
        fprintf('%d %% Complete, ETR: %f min\n',tmp,ETR/60);
        currentPercent = tmp;
    end
end
h = figure;
plot(x,y,'.');
xlabel('Conductivity (S)');
ylabel('Current (I)');
title('Question 2d: Current vs Average Conductivity');
grid on;
saveas(h,'2d.png');
end

function [out] = generateRand(a,b)
    out = rand()*(b-a) + a;
end

function [board] = initializeBoard(board,createPlot)
if nargin == 1
   createPlot = false; 
end

board.x = 0:board.h:board.L;
board.Nx = length(board.x);
board.y = 0:board.h:board.W;
board.Ny = length(board.y);
[board.X,board.Y] = meshgrid(board.x,board.y);
board.N = board.Nx*board.Ny;
board.sig = zeros(board.Ny,board.Nx);
board.xlim = [(board.L - board.Lb)/2,(board.L + board.Lb)/2];
board.ylim = [board.Wb,board.W - board.Wb];
I = board.X >= board.xlim(1) & board.X <= board.xlim(2) & ~(board.Y > board.ylim(1) & board.Y < board.ylim(2));
board.sig(I) = board.sigBox;
board.sig(~I) = board.sigOut;
board.sigAv = mean(board.sig,'all');

if ~createPlot
   return; 
end

h = figure;
surf(board.X,board.Y,board.sig);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Conductance (S)');
title('Question 2a Plot');
view(0,90)
saveas(h,'2aS.png');
end

function [] = question1A(board)
G = spalloc(board.Nx,board.Nx,board.Nx*5);
B = zeros(board.Nx,1);
for n = 1:board.Nx
    if n == 1
        B(n) = board.Vo;
        G(n,n) = 1;
    elseif n == board.Nx
        B(n) = 0;
        G(n,n) = 1;
    else
        B(n) = 0;
        G(n,n) = -2;
        G(n,n + 1) = 1;
        G(n,n - 1) = 1;
    end
end

V = mldivide(G,B);
h = figure;
plot(board.x,V,'*-')
xlabel('x (mm)');
ylabel('Voltage (V)');
title('Question 1a: Voltage Plot');
saveas(h,'1a.png');
end

function [V] = question1b(board)
G = spalloc(board.N,board.N,board.N*5);
B = zeros(board.N,1);
for n = 1:board.N
    curr.n = n;
    [curr.y,curr.x] = getPosSpace(board.Nx,curr.n);
    if curr.x == 1 || curr.x == board.Nx
        B(curr.n) = board.Vo;
        G = assignG(board,curr,G,[0,0],1);
    elseif curr.y == 1 || curr.y == board.Ny
        G = assignG(board,curr,G,[0,0],1);
        B(curr.n) = 0;
    else
        G = assignG(board,curr,G,[0,0],-4);
        G = assignG(board,curr,G,[-1,0],1);
        G = assignG(board,curr,G,[+1,0],1);
        G = assignG(board,curr,G,[0,-1],1);
        G = assignG(board,curr,G,[0,+1],1);
        B(n) = 0;
    end
end

V = mldivide(G,B);
V = reshape(V,[board.Nx,board.Ny])';
h = figure;
surf(board.X,board.Y,V);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Voltage (V)');
title('Question 1b Plot');
saveas(h,'1b.png');
end

function [] = question1bTheory(board)
V = zeros(board.Ny,board.Nx);
h = figure;
ax = subplot(1,1,1,'parent',h);
plt = surf(ax,board.X,board.Y,V);
xlabel(ax,'x (mm)');
ylabel(ax,'y (mm)');
zlabel(ax,'Voltage (V)');
title(ax,'Question 1b Theory Plot');

a = board.W;
b = board.L/2;
x = board.X - board.L/2;
y = board.Y;
for n = 1:2:79
    V = V + (4.*board.Vo./pi).*(cosh(n.*pi.*x./a)).*(sin(n.*pi.*y./a))./(n.*(cosh(n.*pi.*b./a)));
end
set(plt,'ZData',V);
saveas(h,'1bT.png');
end

function [I] = question2a(board,createPlot)
if nargin == 1
    createPlot = false;
end

G = spalloc(board.N,board.N,board.N*5);
B = zeros(board.N,1);
val = zeros(1,4);
for n = 1:board.N
    curr.n = n;
    [curr.y,curr.x] = getPosSpace(board.Nx,curr.n);
    if curr.x == 1
        B(curr.n) = board.Vo;
        G = assignG(board,curr,G,[0,0],1);
    elseif curr.x == board.Nx
        B(curr.n) = 0;
        G = assignG(board,curr,G,[0,0],1);
    elseif curr.y == 1
        [G,val(1)] = assignG(board,curr,G,[1,0]);
        [G,val(2)] = assignG(board,curr,G,[0,-1]);
        [G,val(3)] = assignG(board,curr,G,[0,+1]);
         G = assignG(board,curr,G,[0,0],-sum(val(1:3)));
        B(curr.n) = 0;
    elseif curr.y == board.Ny
        [G,val(1)] = assignG(board,curr,G,[-1,0]);
        [G,val(2)] = assignG(board,curr,G,[0,-1]);
        [G,val(3)] = assignG(board,curr,G,[0,+1]);
        G = assignG(board,curr,G,[0,0],-sum(val(1:3)));
        B(curr.n) = 0;
    else
        [G,val(1)] = assignG(board,curr,G,[-1,0]);
        [G,val(2)] = assignG(board,curr,G,[1,0]);
        [G,val(3)] = assignG(board,curr,G,[0,-1]);
        [G,val(4)] = assignG(board,curr,G,[0,+1]);
        G = assignG(board,curr,G,[0,0],-sum(val));
        B(curr.n) = 0;
    end
end

V = mldivide(G,B);
V = reshape(V,[board.Nx,board.Ny])';
[Ex,Ey] = gradient(V,board.h);
Ex = -Ex;
Ey = -Ey;
Jx = Ex.*board.sig;
Jy = Ey.*board.sig;
I = mean(Jx,'all')*board.W;

if ~createPlot
   return; 
end

tmp = zeros(board.Ny,board.Nx);

h = figure;
surf(board.X,board.Y,V);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Voltage (V)');
title('Question 2a: Voltage Plot');
view(25,20)
saveas(h,'2aV.png');

h = figure;
quiver(board.X,board.Y,Ex,tmp);
hold on
drawLines(board);
axis([0,board.L,0,board.W]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Question 2a: Ex Plot');
grid on
saveas(h,'2aEx.png');

h = figure;
quiver(board.X,board.Y,tmp,Ey);
hold on
drawLines(board);
axis([0,board.L,0,board.W]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Question 2a: Ey Plot');
grid on
saveas(h,'2aEy.png');

h = figure;
quiver(board.X,board.Y,Jx,Jy);
hold on
drawLines(board);
axis([0,board.L,0,board.W]);
xlabel('x (mm)');
ylabel('y (mm)');
title('Question 2a: Current Density Plot');
grid on
saveas(h,'2aJ.png');
end

function [] = drawLines(board)
x = [0,board.xlim(1),board.xlim(1),board.xlim(2),board.xlim(2),board.L,board.L,board.xlim(2),board.xlim(2),board.xlim(1),board.xlim(1),0,0];
y = [0,0,board.ylim(1),board.ylim(1),0,0,board.W,board.W,board.ylim(2),board.ylim(2),board.W,board.W,0];
line(x,y,'Color','k');
end

function [G,val] = assignG(board,curr,G,delta,val)
eq = getEqSpace(board.Nx,curr.y + delta(1),curr.x + delta(2));
if eq <= 0 || eq > board.N
    return;
end
if nargin == 4
    val = mean(board.sig(curr.y + [0,delta(1)],curr.x + [0,delta(2)]),'all');
end
G(curr.n,eq) = val;
end

function [eq] = getEqSpace(nx,yPos,xPos)
eq = xPos + (yPos-1)*nx;
end

function [yPos,xPos] = getPosSpace(nx,eq)
xPos = mod(eq-1,nx)+1;
yPos = floor((eq-1)/nx) + 1;
end