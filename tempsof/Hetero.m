%% Pre-requisite
clear all;
close all;
% set(0,'DefaultFigureWindowStyle','docked');

%% Parameters
n = 0.3;
A = 1.4;
MaxL = 2*(A-n) / ((1+n)^2);
L = MaxL*0.3;

%% Constants
m = 0;
D = 1+2*A-m-n;
a = (2+2*A-n) / D   + (2+2*A)     / D * L;
b = (1+m)     / D   + (1+m+n)     / D * L;
c = (2+2*m)   / D   + (2+2*m+2*n) / D * L;

%% Equilibrium points M0, M1
% M0
p0 = 0;
q0 = 0;
r0_= c;
r0 = r0_^(1/(1+n));

% M1
p1 = 0;
q1 = 1;
r1_= c - (1+n)/(A-n)*L;
r1 = r1_^(1/(1+n));

% Eigenvectors at M0 
V1 = [0 ; (-n + (A-n)/L*r0_ )/r0; -1 ];
V2 = [1 ; b*r0 ; (L+b)*r0^2 / ( 2*n - (A-n)/L*r0_ ) ];
V3 = [0 ; 0 ; -1];
V1 = V1/(norm(V1,2));
V2 = V2/(norm(V2,2));
V3 = V3/(norm(V3,2));

% Eigenvectors at M1 for L\ne b 
W1 = [ ( -n*(1+n)/(A-n) - (A-n)/L*r1_ )*( 1-(1+n)/(A-n) )/( (b-L)*r1^2 ) ; -...
       ( -n*(1+n)/(A-n) - (A-n)/L*r1_ )/r1 ; ...
       ( ( 1-(1+n)/(A-n) )*L + (b-L) )/( b-L ) ];
W2 = [0 ; ( -n - (A-n)/L*r1_ )/r1 ; 1];
W3 = [0 ; 0 ; 1];
W1 = W1/(norm(W1,2));
W2 = W2/(norm(W2,2));
W3 = W3/(norm(W3,2));

%% Solve the ode for n NOT equal to 0 

% Gram Schmidt orthogonalization in the order V2 -> V3 -> V1
B2 = V2;
B3 = V3 - (V3'*B2)/(B2'*B2) * B2; B3 = B3/norm(B3,2);
B1 = V1 - (V1'*B2)/(B2'*B2) * B2 - (V1'*B3)/(B3'*B3) * B3; B1 = B1/norm(B1,2);

figure(10); view(90,20);  hold on;
figure(2);  view(99,58); hold on;
flag = 0;

for w = [32.2157:0.00001:32.2159]
    % initial position of ther orbit near M1
    Z0    = [p1;q1;r1] + 0.00001*(50*W2 + w*W1);

    % solve in negative time
    myodefuncM = @(t,Z,L,A,n) myodefunc(t,Z,L,A,n);
    etaspan = [50 0];      % fast time variable t = eta/n;
    tspan   = etaspan/n;
    options = odeset('AbsTol',1e-12,'RelTol',1e-10);
    [t, Z] = ode45(myodefuncM,tspan,Z0,options,L,A,n);
    eta = t*n;
    
  
    % stop if the last vector contains B1 component
    Zdot = myodefunc(0,Z(end,:),L,A,n);
    Zdot = Zdot / norm(Zdot,2);
    % calculate angle
    [w (B1')*Zdot (B2')*Zdot (B3')*Zdot]
    [Z(end,:)' Zdot]
    if ((B1')*Zdot < 0 && flag==0) 
        % draw orbit
        figure(10); 
        plot3(Z(:,1),Z(:,2),Z(:,3),'r--','LineWidth',2); 
        flag = 1;
        break;
    else
        figure(10); 
        plot3(Z(:,1),Z(:,2),Z(:,3),'b');
%         plot3(Z(1:250,1),Z(1:250,2),Z(1:250,3),'b'); 
    end
    
end
figure(10);hold off;

%% Prepare plots

% draw the last selected orbit
figure(1); 
xlabel('$p$','interpreter','latex'); 
ylabel('$q$','interpreter', 'latex'); 
zlabel('$r$','interpreter', 'latex'); 
set(get(gca,'ZLabel'),'Rotation',0);

hold on;
plot3(Z(:,1),Z(:,2),Z(:,3));
plot3(0,      0,      r0,'rx','MarkerSize',15,'LineWidth',2);
plot3(0,      q1,     r1,'rx','MarkerSize',15,'LineWidth',2);
%plot3(0,      0,      0, 'rx','MarkerSize',15,'LineWidth',2);

set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view(69,12); box on; grid on;
% axis([-0.05 0.2 0 1 0 3.5]); 
hold off;
% ylabh = get(gca,'yLabel');
% set(ylabh, 'Units', 'Normalized')
% set(ylabh, 'Position',get(ylabh,'Position').*[1,1,0.8]);
% axis([0 max(Z(:,1)) -0.1 qs 0 max(Z(:,3))]);

% draw eigenvectors at M0
figure(2); hold on;
xlabel('p');ylabel('q'); zlabel('r');
xlabel('$p$','interpreter','latex'); 
ylabel('$q$','interpreter', 'latex'); 
zlabel('$r$','interpreter', 'latex'); 
set(get(gca,'ZLabel'),'Rotation',0);
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
% set(gca,'XTick',[],'YTick',[],'ZTick',[]);
grid on; box on;
view(99,58);
for w=0:0.00001:0.00010
    myV = [p0;q0;r0] + w*V1;
    plot3(myV(1),myV(2),myV(3),'k.');
    myV = [p0;q0;r0] - w*V1;
    plot3(myV(1),myV(2),myV(3),'k.');
    myV = [p0;q0;r0] + w*V2;
    plot3(myV(1),myV(2),myV(3),'r.');
    myV = [p0;q0;r0] - w*V2;
    plot3(myV(1),myV(2),myV(3),'r.');
    myV = [p0;q0;r0] + w*V3;
    plot3(myV(1),myV(2),myV(3),'b.');
    myV = [p0;q0;r0] - w*V3;
    plot3(myV(1),myV(2),myV(3),'b.');
end
CUT=500;
plot3(Z(end-CUT:end,1),Z(end-CUT:end,2),Z(end-CUT:end,3)); 
plot3(0,      0,      r0,'ro','MarkerSize',15,'LineWidth',2);
% axis([-1e-7 1e-7 -0.25e-6 1e-6 11.2563 11.2565])
figure(2); hold off;

figure(10); 
xlabel('p');ylabel('q'); zlabel('r');
xlabel('$p$','interpreter','latex'); 
ylabel('$q$','interpreter', 'latex'); 
zlabel('$r$','interpreter', 'latex'); 
set(get(gca,'ZLabel'),'Rotation',0);
hold on;
%plot3(0,      qs,     0 ,'rx','MarkerSize',15,'LineWidth',2);
plot3(0,      0,      r0,'ro','MarkerSize',4,'LineWidth',2);
plot3(0,      q1,     r1,'ro','MarkerSize',4,'LineWidth',2);
%plot3(0,      0,      0, 'rx','MarkerSize',15,'LineWidth',2);
hold off;
view(71,76); grid on; box on; 
% % axis([-0.05 0.35 -0.1 1 0 6]);
% % axis([-0.05 0.15 -0.1 0.1 0 6]);

%% Save Data
save pqr.mat A n B1 B2 B3 V1 V2 V3 W1 W2 W3 p0 q0 r0 p1 q1 r1 eta Z








