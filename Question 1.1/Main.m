%% Solve heat equation 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
Alpha=1/16;
T=0.5;
N=5;% number of mesh points of first mesh
N_k=10; %number of subintervals of time interval [0,T]
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);

%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    number_mesh_point(inumber_mesh)=N;
    h=(bx-ax)/N;
    k=T/N_k; %divide the interval [0,T] into N_k sub-intervals of constant length
    r=(Alpha*k)/(h^2);
%% Create the mesh point x(i+1/2)
    x=zeros(N+1,1);
    for i_iter=1:N+1
        x(i_iter)=ax+(i_iter-1)*h; %luoi deu
    end
    
%% create control point x(i)
    x_cp=zeros(N+2,1);
    for i_iter=1:N+2
        if(i_iter==1)
            x_cp(i_iter)=x(i_iter); % voi  x_1/2=x_0
        else
            if(i_iter==N+2)
                x_cp(i_iter)=x(i_iter-1); % voi x_N+1/2=x_N+1
            else
                x_cp(i_iter)=(x(i_iter-1)+x(i_iter))/2.0; %do x_cp(i_iter) la trung diem cua x(i_iter-1)va x(i_iter)
            end
        end
    end  
    
%% create control point t(i)   
    t=zeros(N_k+2,1);
    for i_iter=1:N_k+2
        t(i_iter)=(i_iter-1)*k;
    end
    
 %% Create matrix I
    I=eye(N);
 %% Create matrix A
    A=sparse(N);
    for i_iter=1:N
        if (i_iter==1)
            A(i_iter,i_iter)=-2*r;
            A(i_iter,i_iter+1)=r;
        elseif(i_iter==N)
            A(i_iter,i_iter)=-2*r;           
            A(i_iter,i_iter-1)=r;
        else
            A(i_iter,i_iter)=-2*r;
            A(i_iter,i_iter+1)=r;
            A(i_iter,i_iter-1)=r;
        end
    end

 %% Create discrete solution with boundary 
    u_dis=zeros(N+2,N_k+2);
    for i_iter=2:N+1
        u_dis(i_iter,1)=sin(2*pi*x(i_iter));
    end
    u_dis(1,:)=zeros(1,N_k+2);
    u_dis(N+2,:)=zeros(1,N_k+2);
    G=I+k*A;
    for i_iter=2:N_k+2
       u_dis(2:N+1,i_iter)=(I+k*A)*u_dis(2:N+1,i_iter-1);  %Forward Euler method
    end

 %% Get exact solution
    u_ex=zeros(N+2,N_k+2);
    for i_iter=1:N+2
        for j_iter=1:N_k+2
            if i_iter==1||i_iter==N+2
                u_ex(i_iter,j_iter)=0;
            else
            u_ex(i_iter,j_iter)=exact_solution(x_cp(i_iter),t(j_iter));
            end
        end
    end

%%  Calculate the error on L^2 
    norm_l2(inumber_mesh)=0;
    for i_iter=1:N+2
        norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i_iter) - u_ex(i_iter))^2*h;
    end
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    L=norm_l2(inumber_mesh)
%% Calculate the error on H1
    norm_h1(inumber_mesh)=0;
    for i_iter=1:N+1
        if i_iter==1
           norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i_iter+1)-u_ex(i_iter+1))-(u_dis(i_iter)-u_ex(i_iter)))/h)^2*h*0.5;
        elseif i_iter==N+1
            norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i_iter+1)-u_ex(i_iter+1))-(u_dis(i_iter)-u_ex(i_iter)))/h)^2*h*0.5;
        else
            norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i_iter+1)-u_ex(i_iter+1))-(u_dis(i_iter)-u_ex(i_iter)))/h)^2*h;
        end
    end
  norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
  H=norm_h1(inumber_mesh)
    %% Figure exact and dicrete solutions    
    figure
    for i_iter=1:5
    pause(0.25)
    plot(x_cp,u_dis(:,i_iter),'ro',x_cp,u_ex(:,i_iter),'b');
    title('Discrete solution and Exact solution');
    legend('Discrete solution','Exact solution');
    end
    
 %% Refine mesh (increse mesh point)    
     N=2*N;
end

%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point),-log(norm_l2),'r', log(number_mesh_point), -log(norm_h1),'blue', log(number_mesh_point),log(number_mesh_point), 'black', log(number_mesh_point), 2*log(number_mesh_point),'green');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Error');
legend('L^2 Norm', 'H^1 norm', 'x', '2x','Location','NorthEastOutside');