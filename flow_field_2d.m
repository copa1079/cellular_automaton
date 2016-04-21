% 2D Lattice Boltzmann (BGK) model of a fluid. By Ian Haslam.
%  c4  c3   c2  D2Q9 model. At each timestep, particle densities propagate
%    \  |  /    outwards in the directions indicated in the figure. An
%  c5 -c9 - c1  equivalent 'equilibrium' density is found, and the densities
%    /  |  \    relax towards that state, in a proportion governed by omega.
%  c6  c7   c8  Iain Haslam, March 2006.
% Comments added by GT, April 2007
% Edited by Cole Pazar, April 2016, for playing around with 2D models
%
% go to: https://github.com/copa1079/cellular_automaton  for more outputs etc.
%
%% initialize code
% constants
font = 18;
omega=1.0;           % = 1/tau where tau=relaxation time scale
factor = 10;         % factor to skew the density
density=factor*1.0;  % initial density
% this sneakily incorporates the factor of of 3 that shows up in the
% equilibrium distributions, so c_squ actually equals c^2/3.
c_squ=1/3;  
% these are the weights for equilibrium distribution
t1=4/9;   % Rest particle
t2=1/9;   % "Straight" directions
t3=1/36;  % Diagonal directions
% lattice size:
size = 50; nx=size; ny=size;
% set up matrix:
matrix_length = 9;
F=repmat(density/matrix_length,[nx ny matrix_length]); FEQ=F; 
msize=nx*ny;
CI=(0:msize:msize*matrix_length);
% set up flow domain
% initial fluid flow is in the positive x direction, flowing from left to
% right, which can be seen by setting the space to zero:
space = 0.02;   % very minimal obstructions = 0.01, max debris = 0.5
gaps = 1-space; % = a number between 0.5?0.99 for varying porosities
BOUND=rand(nx,ny)>gaps; % porous (depending on gaps) random domain
ON=find(BOUND); % matrix offset of each Occupied Node
TO_REFLECT=[ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4) ...
            ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8)];
REFLECTED= [ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8) ...
            ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4)];
% Set up initial flow velocities
avu=2.0;         % Average x-directed velocity
prevavu=2.0;     % Previous "
ts=0;            % Time-step counter
deltaU=1e-8;     % Velocity added to first row (left side)
numactivenodes=sum(sum(1-BOUND));

%% run the main loop, stopping at equilibrium or 4000 time steps:
while (ts<4000 && 1e-10<abs((prevavu-avu)/avu)) || ts<100
    
    % Propagate (be careful about directions here ...)
    F(:,:,4)=F([2:nx 1],[ny 1:ny-1],4);     % up and right
    F(:,:,3)=F(:,[ny 1:ny-1],3);            % right
    F(:,:,2)=F([nx 1:nx-1],[ny 1:ny-1],2);  % down and right
    F(:,:,5)=F([2:nx 1],:,5);               % up
    F(:,:,1)=F([nx 1:nx-1],:,1);            % propagate down
    F(:,:,6)=F([2:nx 1],[2:ny 1],6);        % up and left
    F(:,:,7)=F(:,[2:ny 1],7);               % left
    F(:,:,8)=F([nx 1:nx-1],[2:ny 1],8);     % down and left
    
    BOUNCEDBACK=F(TO_REFLECT); % densities bouncing back at next timestep
    
    % calculate local densities as the sum of all velocity components
    DENSITY=sum(F,3);
    
    % calculate velocities:
    % sum of "rightward" minus sum of "leftward"
    UX=(sum(F(:,:,[1 2 8]),3)-sum(F(:,:,[4 5 6]),3))./DENSITY; 
    % sum of "upward" minus sum of "downward"
    UY=(sum(F(:,:,[2 3 4]),3)-sum(F(:,:,[6 7 8]),3))./DENSITY;
    
    % apply boundary conditions:
    UX(1,1:ny)=UX(1,1:ny)+deltaU; % increase inlet pressure
    UX(ON)=0; UY(ON)=0; DENSITY(ON)=0;
    U_SQU=UX.^2+UY.^2; U_C2=UX+UY; U_C4=-UX+UY; U_C6=-U_C2; U_C8=-U_C4;
    
    % calculate equilibrium distribution: stationary
    FEQ(:,:,9)=t1*DENSITY.*(1-U_SQU/(2*c_squ));
    % nearest-neighbours
    FEQ(:,:,1)=t2*DENSITY.*(1+UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,3)=t2*DENSITY.*(1+UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,5)=t2*DENSITY.*(1-UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,7)=t2*DENSITY.*(1-UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    % next-nearest neighbours
    FEQ(:,:,2)=t3*DENSITY.*(1+U_C2/c_squ+0.5*(U_C2/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,4)=t3*DENSITY.*(1+U_C4/c_squ+0.5*(U_C4/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,6)=t3*DENSITY.*(1+U_C6/c_squ+0.5*(U_C6/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,8)=t3*DENSITY.*(1+U_C8/c_squ+0.5*(U_C8/c_squ).^2-U_SQU/(2*c_squ));
    
    % compute new particle velocities
    F=omega*FEQ+(1-omega)*F;
    
    % bounce-back from boundaries
    F(REFLECTED)=BOUNCEDBACK;
    
    % update time step and check for equilibrium
    prevavu=avu;avu=sum(sum(UX))/numactivenodes; ts=ts+1;
end
%% finalize plotting
figure(1);clf;colormap(gray(2));image(2-BOUND');hold on;
quiver(2:nx,1:ny,UX(2:nx,:)',UY(2:nx,:)');
title(['Flow field after \deltat = ',num2str(ts)],'fontsize',font);
xlabel('x-distance','fontsize',font);
ylabel('y-distance','fontsize',font);
set(gca,'fontsize',font)
% end of code
