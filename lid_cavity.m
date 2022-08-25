clc;
nx=51; %no of grid points in x direction
ny=51; %no of grid points in y direction
Re=100; %Reynolds number
a=1;    %length of the square cavity
dx=a/(nx-1);
dy=a/(ny-1);
no_of_iterations=500;
tol=0.001;  %parameter for the point gauss siedel iteration
opengl software
%--------------------------------------------------------
%initializing all the variables
omega=zeros(nx,ny);
omega_new=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);
phi_new=zeros(nx,ny);
p=zeros(nx,ny);
p_new=zeros(nx,ny);
X=zeros(nx,ny);
Y=zeros(nx,ny);
%----------------------------------------------------------
%calculation of stream function
for b=1:no_of_iterations
    k=1;
    error=100;
    phi=zeros(nx,ny);
    while error>tol
        for i=2:nx-1
            for j=2:ny-1
                phi_new(i,j)=0.25*(phi(i+1,j)+phi_new(i-1,j)+phi(i,j+1)+phi_new(i,j-1)+(dx*dx)*omega(i,j));
            end
        end
        error=sum(abs(phi-phi_new),'all');
        phi=phi_new;
        k=k+1;
    end
    %calculation of velocities
    for i=2:nx-1
        for j=2:ny-1
            u(i,j)=(phi_new(i,j+1)-phi_new(i,j-1))/(2*dy);
            v(i,j)=(phi_new(i-1,j)-phi_new(i+1,j))/(2*dx);
        end
    end
    u(:,1)=0;
    u(:,ny)=1;
    v(1,:)=0;
    v(nx,:)=0;
    %boundary condition for vorticity
    for i=2:nx-1
        for j=2:ny-1
            omega(i,1)=-2.0*phi_new(i,2)/(dy*dy);                 % bottom wall
            omega(i,ny)=-2.0*phi_new(i,ny-1)/(dy*dy)-2.0/dy;      % top wall
            omega(1,j)=-2.0*phi_new(2,j)/(dx*dx);                 % left wall
            omega(nx,j)=-2.0*phi_new(nx-1,j)/(dx*dx);             %right wall
        end
    end
    %calculation of vorticity
    for i=2:nx-1
        for j=2:ny-1
            omega_new(i,j)=0.25*((1-(dx*Re/2)*u(i,j))*omega(i+1,j)+(1+(dx*Re/2)*u(i,j))*omega(i-1,j)+(1-(dx*Re/2)*v(i,j)) ...
                *omega(i,j+1)+(1+(dx*Re/2)*v(i,j))*omega(i,j-1));
            
        end
    end
    for i=2:nx-1
        for j=2:ny-1
            omega(i,j)=omega_new(i,j);
        end
    end

    
end
%calculating pressure
error1=100;
tol1=0.001;
m=1;
gradientp=zeros(nx,ny);
while error1>tol1
    for i=2:nx-1
        for j=2:ny-1
            gradientp(i,j)=(((phi_new(i-1,j)-2*phi_new(i,j)+phi_new(i+1,j))/(dx*dx))...
            *((phi_new(i,j-1)-2*phi_new(i,j)+phi_new(i,j+1))/(dy*dy)))...
            - ((phi_new(i+1,j+1)-phi_new(i+1,j-1)-phi_new(i-1,j+1)+phi_new(i-1,j-1))/(4*(dx*dy)))^2;
          
              p_new(i,j)= (0.25*(p(i+1,j)+p_new(i-1,j) + p(i,j+1)+p_new(i,j-1))- 0.5*((gradientp(i,j)*dx^2)));
        end
    end
    error1=sum(abs(p-p_new),'all');
    p=p_new;
    m=m+1;
end


%visualizations
x=linspace(0,1,nx);
y=linspace(0,1,ny);
for i=1:nx
    for j=1:ny
        X(i,j)=x(i);
        Y(i,j)=y(j);
    end
end
figure(1)
contourf(X,Y,omega_new),xlabel('X direction'),ylabel('Y direction'),title('Vorticity');
figure(2)
contourf(X,Y,phi_new),xlabel('X direction'),ylabel('Y direction'),title('Stream function');
figure(3)
contourf(X,Y,p_new),xlabel('X direction'),ylabel('Y direction'),title('Pressure Distribution');
figure(4)
plot(u(:,26),y),xlabel('X direction'),ylabel('Y direction'),title('u velocity plot');
figure(5)
plot(x,v(:,26)),xlabel('X direction'),ylabel('Y direction'),title('v velocity plot');


