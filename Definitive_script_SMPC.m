clear;clc;close; clear all;
animated_plot = false;

MPCtype='prob'; % 'robust', 'prob', 'nominal'
disturb_flag= 0;  %0 no, 1 yes
constant_dist= 0; %0 no, 1 yes

%% Reference position and velocity for Com and ZMP
omega = sqrt(9.8/0.231);   % it first was /0.33   
 % center of feet
foot_distance_x = 0.06;  % we set it as 0.07 
foot_distance_y = 0.098/2; % we set it as 0.09 
% time discretization
delta = 0.008;    % we set it as 0.01 
% foot half width
w = 0.05/2; 

number_of_steps=12;
t=0:0.01:20;

%% Define walk
mult=1;
g=0;             %contatore debugging
superlista=[];   %lista debugging
counter_list=[]; %lista debugging counter
break_flag=0;
counter=0;
%new part for disturbance simulations
rng(1);
disturb_mult=0.4;   %It fisrt was 0.4   Prob: 4.5 is the maximum achievable. Nominal: 4.5 oscillates much more

if disturb_flag == 0
    disturb_mult=1;
    constant_dist=0;
end
S = 50;       % we set it as 30*mult %%%%30, improved with 50, much improved with 70, best value until explosion: 90
D = 30;       % we set it as 20*mult %%%%20, improved with 40, much improved with 60, best value until explosion: 60
N = 80*mult;  % it was 50 %%%%It was 100. 30 not ok. 50 ok.     100Explodes if S=90, D=60, N=50 or 80 or 90. With N=120 it is worse 


balance = false;

if balance == true
    fs_matrix = zeros (26,2);
else
    fs_matrix = [0,-foot_distance_y; %middle of rhe foot point
        foot_distance_x,foot_distance_y;
        2*foot_distance_x,-foot_distance_y;
        3*foot_distance_x,foot_distance_y;
        4*foot_distance_x,-foot_distance_y;
        5*foot_distance_x,foot_distance_y;
        6*foot_distance_x,-foot_distance_y;
        7*foot_distance_x,foot_distance_y;
        8*foot_distance_x,-foot_distance_y;
        9*foot_distance_x,+foot_distance_y;
        10*foot_distance_x,-foot_distance_y;
        11*foot_distance_x,+foot_distance_y;
        12*foot_distance_x,-foot_distance_y
        13*foot_distance_x,+foot_distance_y;
        14*foot_distance_x,-foot_distance_y;
        15*foot_distance_x,+foot_distance_y
        16*foot_distance_x,-foot_distance_y;
        17*foot_distance_x,+foot_distance_y;
        18*foot_distance_x,-foot_distance_y;
        19*foot_distance_x,+foot_distance_y;
        20*foot_distance_x,-foot_distance_y;
        21*foot_distance_x,+foot_distance_y;
        22*foot_distance_x,-foot_distance_y;
        23*foot_distance_x,+foot_distance_y;
        24*foot_distance_x,-foot_distance_y;
        25*foot_distance_x,+foot_distance_y;
        26*foot_distance_x,-foot_distance_y];
end

fs_matrix = [     0      0;     
    0  0.049;   
    0 -0.049;    
    0  0.049;    
    0 -0.049;  
 0.06  0.049;     
 0.12 -0.049;   
 0.18  0.049;    
 0.24 -0.049;    
  0.3  0.049;    
 0.36 -0.049;    
 0.36      0;   
 0.36      0;     
  0.36      0;    
  0.36      0;     
  0.36      0;     
  0.36      0;    
  0.36      0;    
  0.36      0];

%% General parameters
w = 0.05/2;

x(1) = 0.0;
xd(1) = 0.0;
zx(1) = 0.0;

y(1) = 0;
yd(1) = 0;
zy(1) = 0;


%% Compute constraints
f1_y = -foot_distance_y;
f2_y = foot_distance_y;
bound_y_com=0.06;   % It was 0.1

additionalFirstStepDuration = 20*mult*0; 

fs_sequence_x = zeros(S+D+additionalFirstStepDuration,1);
fs_sequence_y = zeros(S+D+additionalFirstStepDuration,1);

for i = 1:size(fs_matrix,1)-2
    f1_x = fs_matrix(i,1);
    f2_x = fs_matrix(i+1,1);
    
    f1_y = fs_matrix(i,2);
    f2_y = fs_matrix(i+1,2);
    
    fs_sequence_x = [fs_sequence_x; ones(S,1) * f1_x; f1_x + (1:D)'*(f2_x-f1_x)/D];
    fs_sequence_y = [fs_sequence_y; ones(S,1) * f1_y; f1_y + (1:D)'*(f2_y-f1_y)/D];
end

fs_sequence_x(1) = [];
fs_sequence_y(1) = [];

zx_min = fs_sequence_x - w;
zx_max = fs_sequence_x + w;
if balance == true
    zy_min = fs_sequence_y - 5*w;
    zy_max = fs_sequence_y + 5*w;
else
    zy_min = fs_sequence_y - w;
    zy_max = fs_sequence_y + w;
end
zy_min(1:S+D+additionalFirstStepDuration-1) = zy_min(1:S+D+additionalFirstStepDuration-1) - foot_distance_y;
zy_max(1:S+D+additionalFirstStepDuration-1) = zy_max(1:S+D+additionalFirstStepDuration-1) + foot_distance_y;

%% compute references 
%ZMP reference----------------------------------------------------
x_speed=0.02;    %it was 0.02  %%%%%constant speed x direction motion : no influence on the spikes at all
zx_ref=fs_sequence_x;
zy_ref=fs_sequence_y;


%COM reference----------------------------------------------------
cx_delta=foot_distance_x/(S+D);
cx_ref=zeros(size(fs_sequence_x));    
cy_ref=zeros(size(fs_sequence_x));
for i=2:size(fs_sequence_x)
   %if i<=(S+D+additionalFirstStepDuration)
   if fs_sequence_x(i)==0
       cx_ref(i)=0;
       cy_ref(i)=0;
   else
       cx_ref(i)=cx_ref(i-1)+cx_delta;
       cy_ref(i)=0;
   end
end


%% Cost index matrix H
H=zeros(3*N,3*N);
alpha=0.5; %%%%1, 0 is better slightly
beta=1;    %%%%1, 0 explodes, 0.5 worsens, 0.8 improves
gamma=7;    % was 7%%%%1, and 0 worsens, 1.5 improves, 2.5 worsens


for i=1:N
    H(2*i-1: 2*i,2*i-1:2*i)=[alpha,0;
                            0 ,  beta];
end
H(2*N+1:end,2*N+1:end)=eye(N)*gamma;

%% State dynamics constraints

A=[0         1 ;
   omega^2   0];

B=[0;
   -omega^2];

%% Compute discretized state space 

A_discrete=expm(A*delta);

B_discrete=[0;0];
syms x

% The first time we set it as:
%fun1= @(x) [2.72*exp(-5.45*x) - 2.72*exp(5.45*x)];
%fun2= @(x) [- 14.8*exp(-5.45*x) - 14.8*exp(5.45*x)];

% The second time we set it as:
fun1= @(x) [3.26*exp(-6.51*x) - 3.26*exp(6.51*x)];
fun2= @(x) [- 21.2*exp(-6.51*x) - 21.2*exp(6.51*x)];


B_discrete(1,1)=integral(fun1,0,delta);
B_discrete(2,1)=integral(fun2,0,delta);



s=zeros((N)*2,N);
for i=1:N
    for j=1:i
        s(i*2-1:i*2,j)=A_discrete^(i-j)*B_discrete;
    end
end

T=zeros(2*N,2);
for i=1:N
    T(2*i-1:2*i,1:2)=A_discrete^i;
end

A_eq=[eye(2*N,2*N),-s];


%% Robust control constraints
k_rob=exp(omega*delta)/(exp(omega*delta)-1);
K_robust=k_rob*[1 1/omega];
% Ahmed W_bounded=[0.0016;0.016];
W_bounded=disturb_mult*[0.0005;0.005]; % expressed in meters
eta_rob=K_robust*W_bounded             % a scalar
eta_rob=eta_rob*ones(1,N);




%% Probabilistic contraints
K=-place(A_discrete,B_discrete,[0.9 0.8]); 
if isequal(MPCtype,'robust')
    K=K_robust;
end
Ak=A_discrete+B_discrete*K;

% The sd in stochastic case is nearly 1/3 wrt the robust case for coherence

Sigma_w=disturb_mult^2*[ 0.00015^2 , 0   ;
           0     , 0.0015^2];
variances=[];
variances_y=[]; % for the bound on the y_com 
eta_list=[];
eta_list_y=[];  % for the bound on the y_com 

b=0.05;
b_y=0.05;

for i=0:N-1
        % compute variance
        temporary=zeros(2,2);
        for j=0:i
            temporary=temporary+Ak^j*Sigma_w*(Ak^j)';
        end

        variances(i+1)=K*temporary*K'; 
        variances_y(i+1)=temporary(1,1);
        eta_list(i+1)=norminv(1-b,0,sqrt(K*temporary*K'));
        eta_list_y(i+1)=norminv(1-b_y,0,sqrt(temporary(1,1)));
end

%this eta list is valid for both components x and y
figure
plot(eta_list_y)
figure
plot(eta_list)


%% Robust y compute omega 

om_A=W_bounded;
for i=1:50
    om_A=Ak^i*W_bounded+om_A;
end
OMEGA=abs(om_A(2));

%% Span time 

% initial conditions
cx=[0];
cxd=[0];
zmp_x=[0];
cy=[0];
cyd=[0];
zmp_y=[0];

bx=[];
reference_stack_x=[];
ZMP_reference_stack_x=[];

by=[];
reference_stack_y=[];
ZMP_reference_stack_y=[];

%disturb
w_dist=[];


for i=1:N
    alt_id(2*i-1:2*i,2*i-1:2*i)=[1,0;
                            0 ,  0];
end
    A_zmp_inequality_y=[zeros(N,2*N),eye(N);
                       zeros(N,2*N),-eye(N);
                       alt_id, zeros(2*N,N);
                       -alt_id, zeros(2*N,N)];
 
    A_zmp_inequality=[zeros(N,2*N),eye(N);
                       zeros(N,2*N),-eye(N);];


time=10*(S+D);   %It was 6.



% Plot disturbance
w_vector_1 = [];
w_vector_2 = [];




for k=0:time
 %constraint on zmp-----------------------------
 
 %disturb realization
 w_dist=[randn*sqrt(Sigma_w(1,1)); randn*sqrt(Sigma_w(2,2))]*disturb_flag; %hp good knowledge of disturbe variance
 if w_dist(1)>W_bounded(1)
     w_dist(1)=W_bounded(1);
 end
 if w_dist(2)>W_bounded(2)
     w_dist(2)=W_bounded(2);
 end
 
 if constant_dist == 1
     w_dist(1)=0.9*W_bounded(1);
     w_dist(2)=0.9*W_bounded(2);
     
 end
 
 w_vector_1(k+1) = w_dist(1);  %position
 w_vector_2(k+1) = w_dist(2);  %speed
 
  % X direction
    bx=[];
    reference_stack_x=[];
    ZMP_reference_stack_x=[];
    
    
               
    
                   
                   
   for i=1:N
        reference_stack_x=[reference_stack_x,cx_ref(k+i),(cx_ref(k+i+1)-cx_ref(k+i))/delta]; 
        ZMP_reference_stack_x=[ZMP_reference_stack_x,zx_ref(k+i)];
        % ZMP constraints
        if isequal(MPCtype,'nominal')
            % NOMINAL CASE
            bx(i)=zx_max(k+i);
            bx(i+N)=-zx_min(k+i);
        end
        
        if isequal(MPCtype,'prob')
            % PROBABILISTIC CASE
            bx(i)=zx_ref(k+i)+w-eta_list(i);
            bx(i+N)=-zx_ref(k+i)+w-eta_list(i);
        end
        if isequal(MPCtype,'robust')
            % ROBUST CASE CASE
            bx(i)=zx_ref(k+i)+w-eta_rob(i);
            bx(i+N)=-zx_ref(k+i)+w-eta_rob(i);
        end
    end
    
    bx=bx';        
    
    %constraint on model---------------------------
     b_eq_x=T*[cx(k+1);cxd(k+1)];
     
     %cost function f----------------------------------
     f=-2*[reference_stack_x*H(1:2*N,1:2*N),gamma*ZMP_reference_stack_x]';
     % X optimization-------------------------
     
    %  X=quadprog(2*H,f,A_zmp_inequality,bx,A_eq,b_eq_x);
    X=quadprog(2*H,f,A_zmp_inequality,bx,A_eq,b_eq_x);
    % X(2*N+1) is the zmp value that we want
    
        % update discrtet time model
    update=A_discrete*[cx(k+1);cxd(k+1)]+B_discrete*X(2*N+1)+0*w_dist; %no disturbance on x 
    cx(k+2)=update(1);
    cxd(k+2)=update(2);
    zmp_x(k+1)=X(2*N+1);
    
      %-----------------------------
    %-----------------------------
    %-----------------------------
    %-----------------------------
    
    
    % Y direction
        by=[];
    reference_stack_y=[];
    ZMP_reference_stack_y=[];
    
%      w_dist=[randn*sqrt(Sigma_w(1,1)); randn*sqrt(Sigma_w(2,2))]*disturb_flag; %hp good knowledge of disturbe variance
%     if w_dist(1)>W_bounded(1)
%      w_dist(1)=W_bounded(1);
%  end
%  if w_dist(2)>W_bounded(2)
%      w_dist(2)=W_bounded(2);
%  end
        for i=1:N
        reference_stack_y=[reference_stack_y,cy_ref(k+i),(cy_ref(k+i+1)-cy_ref(k+i))/delta];
        ZMP_reference_stack_y=[ZMP_reference_stack_y,zy_ref(k+i)];
        % ZMP constraints

        if isequal(MPCtype,'nominal')
            % NOMINAL CASE
            by(i)=zy_max(k+i);
            by(i+N)=-zy_min(k+i);
            by(2*i+2*N-1:2*i+2*N)=[bound_y_com, 0];
            by(2*i+4*N-1:2*i+4*N )=[bound_y_com,0];
        end
        if isequal(MPCtype,'prob')
            % PROBABILISTIC CASE
            by(i)=zy_ref(k+i)+w-eta_list(i);
            by(i+N)=-zy_ref(k+i)+w-eta_list(i);
            
            by(2*i+2*N-1:2*i+2*N)=[bound_y_com-eta_list_y(i), 0];
            by(2*i+4*N-1:2*i+4*N )=[bound_y_com-eta_list_y(i),0];
        end
            
        if isequal(MPCtype,'robust')
            % ROBUST CASE CASE
            by(i)=zy_ref(k+i)+w-eta_rob(i);
            by(i+N)=-zy_ref(k+i)+w-eta_rob(i);
            by(2*i+2*N-1:2*i+2*N)=[bound_y_com-OMEGA, 0];
            by(2*i+4*N-1:2*i+4*N )=[bound_y_com-OMEGA,0];
        
        end
    end
    
    by=by';
    
     f=-2*[reference_stack_y*H(1:2*N,1:2*N),gamma*ZMP_reference_stack_y]';
    
    
    % update s_t 
%     b_eq=T*s_t
    b_eq_y=T*[cy(k+1);cyd(k+1)];
    
    
    
    Y=quadprog(2*H,f,A_zmp_inequality_y,by,A_eq,b_eq_y);
    % X(2*N+1) is the zmp value that we want
    if length(Y)~=0 %recursive feasibility
        Last_Y=Y;
        
        % update discrtet time model
        zmp_y(k+1)=Y(2*N+1);
        %zmp_y(k+2)=Y(2*N+2)+K*([cy(k+1);cyd(k+1)]-[; %COMPLETA
    
        update=A_discrete*[cy(k+1);cyd(k+1)]+B_discrete*zmp_y(k+1)+w_dist;
        cy(k+2)=update(1);
        cyd(k+2)=update(2);
        
        if counter >0
            counter_list=[counter_list, counter];
        end
        
        counter=0;
    
    else %fix no feasibility
        g=g+1;
        superlista=[superlista, k];
        zmp_y(k+1)=Last_Y(2*N+2+counter);
        update=A_discrete*[cy(k+1);cyd(k+1)]+...
            B_discrete*(zmp_y(k+1)+K*([cy(k+1);cyd(k+1)]-[Last_Y(3+2*counter);Last_Y(4+2*counter)]))+w_dist;
        cy(k+2)=update(1);
        cyd(k+2)=update(2);
        counter=counter+1;
        if counter == N-1     %% It was 99
            break_flag=1;
        end
    end
    
    
    
    
    
    
    
    
     if break_flag==1
         break
     end
end


% figure()
% hold on
% plot(cx_ref(1:time+1),cy_ref(1:time+1),'m',zx_ref(1:time+1),zy_ref(1:time+1),'b','lineWidth',1)
% plot(cx(1:time+1),cy(1:time+1),'g','lineWidth',2)
% plot(zmp_x(1:time+1),zmp_y(1:time+1),'r','lineWidth',2)

figure()
hold on
plot(cx_ref(1:time+1),cy_ref(1:time+1),'m',zx_ref(1:time+1),zy_ref(1:time+1),'b','lineWidth',1)
plot(cx,cy,'g','lineWidth',2)
plot(zmp_x,zmp_y,'r','lineWidth',2)



  rect_x = [w,w,-w,-w,w];
    rect_y = [foot_distance_y+w,-foot_distance_y-w,-foot_distance_y-w,foot_distance_y+w,foot_distance_y+w];
    plot(rect_x,rect_y,'m','lineWidth',2,'HandleVisibility','off');
    
    rect_x = [w,w,-w,-w,w];
    rect_y = [w,-w,-w,w,w];
    
    nPlottedFootsteps = 8;
    
    for j = 1:nPlottedFootsteps
        rect_x = [w,w,-w,-w,w];
        rect_y = [w,-w,-w,w,w];
        h1 = plot(fs_matrix(j,1)+rect_x,fs_matrix(j,2)+rect_y,'m','lineWidth',2,'HandleVisibility','off');
    end

    
    
% Compute velocities on x and y
xdot = (gradient(cx)/delta).';
ydot = (gradient(cy)/delta).';

% Transpose cx and cy
cx_transpose=cx.';
cy_transpose=cy.';

% Compute un' altra footstep plan matrix 
%f_sequence_matrix = [fs_sequence_x.' , fs_sequence_y.'];




%Plot disturbance
figure()
hold on
%plot(w_vector_1,'g','lineWidth',2) %posizione
plot(w_vector_2,'r','lineWidth',2)
axis([0 180 -0.0025 0.0025]);
xlabel('t');
ylabel('m/s^2');
legend("velocity");

figure
plot(cy)


