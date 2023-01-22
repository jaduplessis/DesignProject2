
clc
clear
tic

% set problem constants
k=0.033; %thermal conductivity of polystyrene
h=1; %convectivity
temp_amb=30; % temperature at external boundary
temp_sink=0; % temperature at internal boundary (sink)
scalar = 0.1; % lengths are scaled down to reduce size of program
length_x = 4; % corresponds to 40 cm
length_y = 3; % 30cm
length_z = 2; % 20cm
d=0.05; % Element size. Resolution. Smaller is better
tol=1e-10; % to stop iteration

% box properties
t = 0.4; % thickness

% set location of sink - two values are [min, max]
% will be rounded to the nearest whole element
% values parametrised by thickness and outer dimensions
source_x=[t,length_x-t];
source_y=[t,length_y-t];
source_z=[t,length_z-t];

% calculate some values for the sink
nx=round(length_x/d); % number of nodes must be an integer
ny=round(length_y/d);
nz=round(length_z/d);
sx=round(source_x/d); % nodes for sink
sy=round(source_y/d);
sz=round(source_z/d);
sx=sx(1)+1:sx(2); % turn into ranges
sy=sy(1)+1:sy(2);
sz=sz(1)+1:sz(2);
if ismember(1,sx) || ismember(nx,sx) || ismember(1,sy) || ismember(ny,sy) || ismember(1,sz) || ismember(nz, sz)
    error('Sink elements cannot be on the boundary.')
end

% fill matrices with initial values (sink temperature)
T_new=zeros(nx,ny,nz)+temp_sink;
T_old=T_new;

% set up constants for convection and conduction
c2=k/(6*k + 4*h*d); %c2,c3 are for edges
c3=(2*h*d)/(3*k+2*h*d);
c4=k /(3*k+3*h*d); %c4,c5 for corners
c5=h*d/(k+h*d);
c6=k/(6*k+2*h*d); %c6,c7 will be for surfaces
c7=h*d/(3*k+h*d);

% set up neighbour kernel for convolution. Outputs a 3x3x3 matrix
neighbours(:,:,1)=[[0,0,0];[0,1/6,0];[0,0,0]];
neighbours(:,:,2)=[[0,1/6, 0];[1/6,0,1/6];[0,1/6,0]];
neighbours(:,:,3)=[[0,0,0];[0,1/6,0];[0,0,0]];

% set up main iteration loop
errc=1; erro=errc; it=0;
while errc/erro>tol

    % calculate updated temperatures for internal nodes of problem matrix using convolution
    T_new(2:nx-1,2:ny-1,2:nz-1)= convn(T_old,neighbours,'valid');
    
    % equations for the edges
    T_new(2:nx-1, ny, nz) = c2*(T_old(1:nx-2, ny, nz) + T_old(3:nx, ny, nz) + 2*T_old(2:nx-1, ny-1,nz) + 2*T_old(2:nx-1, ny, nz-1)) + c3*temp_amb;
    T_new(2:nx-1, 1, nz) = c2*(T_old(1:nx-2, 1, nz) + T_old(3:nx, 1, nz) + 2*T_old(2:nx-1, 2,nz) + 2*T_old(2:nx-1, 1, nz-1)) + c3*temp_amb;
    T_new(2:nx-1, ny, 1) = c2*(T_old(1:nx-2, ny, 1) + T_old(3:nx, ny, 1) + 2*T_old(2:nx-1, ny-1,1) + 2*T_old(2:nx-1, ny, 2)) + c3*temp_amb;
    T_new(2:nx-1, 1, 1) = c2*(T_old(1:nx-2, 1, 1) + T_old(3:nx, 1, 1) + 2*T_old(2:nx-1, 2 , 1) + 2*T_old(2:nx-1, 1, 2)) + c3*temp_amb;
    
    T_new(nx, 2:ny-1, nz) = c2*(T_old(nx, 1:ny-2, nz) + T_old(nx, 3:ny, nz) + 2*T_old(nx-1, 2:ny-1,nz) + 2*T_old(nx, 2:ny-1, nz-1)) + c3*temp_amb;
    T_new(1, 2:ny-1, nz) = c2*(T_old(1, 1:ny-2, nz) + T_old(1, 3:ny, nz) + 2*T_old(2, 2:ny-1,nz) + 2*T_old(1, 2:ny-1, nz-1)) + c3*temp_amb;
    T_new(nx, 2:ny-1, 1) = c2*(T_old(nx, 1:ny-2, 1) + T_old(nx, 3:ny, 1) + 2*T_old(nx-1, 2:ny-1,1) + 2*T_old(nx, 2:ny-1, 2)) + c3*temp_amb;
    T_new(1, 2:ny-1, 1) = c2*(T_old(1, 1:ny-2, 1) + T_old(1, 3:ny, 1) + 2*T_old(2, 2:ny-1, 1) + 2*T_old(1, 2:ny-1, 2)) + c3*temp_amb;
    
    T_new(nx, ny, 2:nz-1) = c2*(T_old(nx, ny, 1:nz-2) + T_old(nx, ny, 3:nz) + 2*T_old(nx, ny-1,2:nz-1) + 2*T_old(nx-1, ny, 2:nz-1)) + c3*temp_amb;
    T_new(1, ny, 2:nz-1) = c2*(T_old(nx, 1, 1:nz-2) + T_old(nx, 1, 3:nz) + 2*T_old(nx, 2,2:nz-1) + 2*T_old(nx-1, 1, 2:nz-1)) + c3*temp_amb;
    T_new(nx, 1, 2:nz-1) = c2*(T_old(1, ny, 1:nz-2) + T_old(1, ny, 3:nz) + 2*T_old(1, ny-1,2:nz-1) + 2*T_old(2, ny, 2:nz-1)) + c3*temp_amb;
    T_new(1, 1, 2:nz-1) = c2*(T_old(1, 1, 1:nz-2) + T_old(1, 1, 3:nz) + 2*T_old(1, 2 , 2:nz-1) + 2*T_old(2, 1, 2:nz-1)) + c3*temp_amb;
    
    % equations for the corners
    T_new(nx, ny, nz) = c4*(T_old(nx-1, ny, nz) + T_old(nx, ny-1, nz) + T_old(nx, ny, nz-1)) + c5*temp_amb;
    T_new(nx, 1, nz) = c4*(T_old(nx-1, 1, nz) + T_old(nx, 2, nz) + T_old(nx, 1, nz-1)) + c5*temp_amb;
    T_new(nx, ny, 1) = c4*(T_old(nx-1, ny, 1) + T_old(nx, ny-1, 1) + T_old(nx, ny, 2)) + c5*temp_amb;
    T_new(nx, 1, 1) = c4*(T_old(nx-1,1,1) + T_old(nx, 2, 1) + T_old(nx, 1, 2)) + c5*temp_amb;
    T_new(1, ny, nz) = c4*(T_old(2, ny, nz) + T_old(1, ny-1, nz) + T_old(1, ny, nz-1)) + c5*temp_amb;
    T_new(1, 1, nz) = c4*(T_old(2, 1, nz) + T_old(1, 2, nz) + T_old(1, 1, nz-1)) + c5*temp_amb;
    T_new(1, ny, 1) = c4*(T_old(2, ny, 1) + T_old(1, ny-1, 1) + T_old(1, ny, 2)) + c5*temp_amb;
    T_new(1, 1, 1) = c4*(T_old(2, 1, 1) + T_old(1, 2, 1) + T_old(1, 1, 2)) + c5*temp_amb;
    
    % equations for the surfaces
    T_new(2:nx-1, 2:ny-1, nz) = c6*(T_old(2:nx-1, 1:ny-2, nz) + T_old(2:nx-1, 3:ny, nz) + T_old(1:nx-2, 2:ny-1, nz) + T_old(3:nx, 2:ny-1, nz) + 2*T_old(2:nx-1, 2:ny-1, nz-1)) + c7*temp_amb;
    T_new(2:nx-1, 2:ny-1, 1) = c6*(T_old(2:nx-1, 1:ny-2, 1) + T_old(2:nx-1, 3:ny, 1) + T_old(1:nx-2, 2:ny-1, 1) + T_old(3:nx, 2:ny-1, 1) + 2*T_old(2:nx-1, 2:ny-1, 2)) + c7*temp_amb;
    T_new(2:nx-1, ny, 2:nz-1) = c6*(T_old(2:nx-1, ny, 1:nz-2) + T_old(2:nx-1, ny, 3:nz) + T_old(1:nx-2, ny, 2:nz-1) + T_old(3:nx, ny, 2:nz-1) + 2*T_old(2:nx-1, ny-1, 2:nz-1)) + c7*temp_amb;
    T_new(2:nx-1, 1, 2:nz-1) = c6*(T_old(2:nx-1, 1, 1:nz-2) + T_old(2:nx-1, 1, 3:nz) + T_old(1:nx-2, 1, 2:nz-1) + T_old(3:nx, 1, 2:nz-1) + 2*T_old(2:nx-1, 2, 2:nz-1)) + c7*temp_amb;
    T_new(nx, 2:ny-1, 2:nz-1) = c6*(T_old(nx, 2:ny-1, 1:nz-2) + T_old(nx, 2:ny-1, 3:nz) + T_old(nx, 1:ny-2, 2:nz-1) + T_old(nx, 3:ny, 2:nz-1) + 2*T_old(nx-1, 2:ny-1, 2:nz-1)) + c7*temp_amb;
    T_new(1, 2:ny-1, 2:nz-1) = c6*(T_old(1, 2:ny-1, 1:nz-2) + T_old(1, 2:ny-1, 3:nz) + T_old(1, 1:ny-2, 2:nz-1) + T_old(1, 3:ny, 2:nz-1) + 2*T_old(2, 2:ny-1, 2:nz-1)) + c7*temp_amb;

    % enforce fixed temperature at the sink
    T_new(sx,sy,sz)=temp_sink;
    
    % calculate mean change across the matrix
    errc=mean2(abs(T_new-T_old));
    if it==1
        erro=errc;
    end
    
    % swap old and new values for next iteration and increment count
    T_old=T_new;
    it=it+1;
    % continue while loop (unless change ratio drops below tolerance)
end
%Plotting visual data
% provide a visual aid of the data. Slices through the block of data
slice(T_new, [ny/4 2*ny/3 ny], [nx/4 2*nx/3, nx], [1, nz/4 2*nz/3])

% Heat Transfer Calculations
% calculate average temperature changes
h1=mean2(T_new(sx(1)-1, sy, sz)-T_new(sx(1), sy, sz));
h2=mean2(T_new(sx(end)+1, sy, sz)-T_new(sx(end), sy, sz));
h3=mean2(T_new(sx, sy(1)-1, sz)-T_new(sx, sy(1), sz));
h4=mean2(T_new(sx, sy(end)+1, sz)-T_new(sx, sy(end), sz));
h5=mean2(T_new(sx, sy, sz(1)-1)-T_new(sx, sy, sz(1)));
h6=mean2(T_new(sx, sy, sz(end)+1)-T_new(sx, sy, sz(end)));

% total temp change
total_heat = (h1 + h2 + h3 + h4+ h5 + h5)*k;

% Analysis of results
% calculate internal dimensions of sink
inx = (length_x-t)*scalar; % conversion to m
iny = (length_y-t)*scalar;
inz = (length_z-t)*scalar;
internal_volume = inx*iny*inz; % calculating total internal volume

% Calculate energy needed to raise temperature of box
vaccination_volume = 0.185*0.095*0.06*3; % 3 boxes. Dimensions provided externally
latency = 3340; %specific latency for ice
density = 920; % density of ice
ice_volume = internal_volume - vaccination_volume; % Space available for ice
mass = ice_volume * density; % total mass of ice
Latent_energy = mass * latency; % energy required to change state of all ice to water

seconds = Latent_energy / total_heat;
minutes = seconds/60;
hours = minutes/60
print("finished")

toc


