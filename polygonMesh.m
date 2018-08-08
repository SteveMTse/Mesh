%%  No faults ( voronoi mesh)
%5
%%  rotate the perm for a certain degree
%%
%%  output data file format:
%%  The first line is "number of Node" 
%%  N             (# of Node) followed by N lines
%%  1 xCoor, yCoor
%%  2  xCoor, yCoor
%%  ...
%%  N xCoor, yCoor
%%  line "ENDOFSECTION"
%%  "Number of elements"
%%  M ---> number of elements
%% "Number of vertices of elements
%%  L1
%% ...
%% ENDOFSECTION
%%  indices of vertices of elements and permeability
%%  1  i1 i2 i3 ...      K_11 K_22 K_12       ( K is the value of relative permeability on this element)
%%  ...
%%  ENDOFSECTION
%%  Flux type
%%  multipoint flux
%%  N        % number of  elements with multi-point flux.
%%  k1 
%%  ...
%%  k_N
%%  ENDOFSECTION
%%  Boundary Conditions
%%  dirichlet m1
%%  1  eleID  faceID
%%  ...
%%  ENDOFSECTION
%%  neumann  m2    ( if there is Neumann BC)
%%  1 eleID faceID 
%%  ...
%%  ENDOFSECTION

%%  Author: Qian-Yong Chen


close all;
clear all;

NFig =0;

outFileName  = 'polyMesh.cm'



% physical domain
xLeft = -0.5; xRight = 0.5; yBottom = -0.5; yTop = 0.5;

% % % how much should be the verices be moved. 
%  ratio = 0.2;    randomMove = 0 ;
 
 
 Nx = 12; Ny = 12;

% perm: 
homoPerm = [1 ; 1000];
               

%% rotation angle of perm (clockwise direction)
theta = -pi/3;       
perm_rotate = zeros(3, 1);
coorTrans = [ cos(theta)  sin(theta); -sin(theta) cos(theta)];


permOld = [ homoPerm(1) 0; 0  homoPerm(2) ];
permOld_rotate = coorTrans*permOld*coorTrans' 



% initialization of variables
NNode = (Nx+1);   %% number of Nodes
NEle = 2;           %% number of elements

EleToNode = zeros(6, 2);  %% Mapping from element to its nodes

VertX = zeros(NNode,1);
VertY = VertX;              %% coordinates of nodes.


% some elements will have multipoint flux.
MulFlux = zeros(NEle,1);    


% initialize perm
 perm = zeros(3, NEle);

 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % compute the points from left to right, and in each column from bottom to top
 
 % 1st: left most column points
 dy = (yTop-yBottom)/Ny;
 VertX(1:Ny+2) = xLeft;   
 VertY(1) = yBottom;  VertY(2:Ny+1) = (yBottom+dy/2:dy:yTop-dy/2)';
VertY(Ny+2) = yTop;
 
NPts = Ny+2;
 
 
 dx = (xRight-xLeft)/2/Nx;
 % middle points
 for i=1:Nx
     % bottom point
     NPts = NPts + 1;
     VertX(NPts) = xLeft+dx + (i-1)*2*dx;  VertY(NPts) = yBottom;
     
     % interior points
     for j=1:Ny*2
         NPts = NPts + 1;         
         k = fix(j/2);
         if j-k*2 > 0
             VertX(NPts) = xLeft + (i-1)*2*dx + dx*2/3;
             VertY(NPts) = yBottom + (k)*dy + dy/3;         
         else
             VertX(NPts) = xLeft + (i-1)*2*dx + 2*dx*2/3;
             VertY(NPts) = yBottom + (k-1)*dy + 2*dy/3;                      
         end
     end
     
     % top point
     NPts = NPts + 1;
     VertX(NPts) = xLeft+dx + (i-1)*2*dx;  VertY(NPts) = yTop;     
 end
 
 
 % 3rd: right most column points
  VertX(NPts+1:Ny+2+NPts) = xRight;   
 VertY(1+NPts) = yBottom;  VertY(NPts+2:NPts+Ny+1) = (yBottom+dy/2:dy:yTop-dy/2)';
VertY(NPts+Ny+2) = yTop;
 
 NPts = NPts + Ny+2;
 
%  
% NFig = NFig + 1;    figure(NFig);
%  plot(VertX, VertY,'-o');
%  axis([xLeft-0.1 xRight+0.1 yBottom-0.1 yTop + 0.1]);
 
 
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % compute EleToNode from left to right, and in each column from bottom to top
 
 
 % 1st: left most column
 %% bottom one
 NEle = 1;
 EleToNode(1:4,NEle) = [1 Ny+3 Ny+4  2]';
 NEleToNode(NEle) = 4;
 

 % middle one
 for j=1:Ny
     NEle = NEle + 1;
     EleToNode(1:5, NEle) = [1+j  Ny+2+2*(j)  Ny+2+2*(j)+1  Ny+2+2*(j)+2   2+j ]';
        
    NEleToNode(NEle) = 5; 
 end
 
 
 
% 2nd: middle columns
for i=1:Nx-1
    % bottom one
	NEle = NEle + 1;
	NEleToNode(NEle) = 5;
    aID =Ny+3+(i-1)*(2*Ny+2);
    bID = Ny+3+(i)*(2*Ny+2);
	EleToNode(1:5, NEle) = [aID bID bID+1 aID+2 aID+1]';
    
    % middle ones:
    for j=1:Ny-1
        NEle=NEle+1;
        NEleToNode(NEle) = 6;
        aID =Ny+2+ (i-1)*(2*Ny+2) + j*2+1;
        bID = aID+2*Ny+2-1;
		EleToNode(1:6, NEle) = [aID  bID bID+1 bID+2 aID+2 aID+1]';    
    end
    
    % top one
    NEle = NEle + 1;
	NEleToNode(NEle) = 5;
    aID =Ny+2+ i*(2*Ny+2);
    bID = Ny+2+(i+1)*(2*Ny+2);
	EleToNode(1:5, NEle) = [aID aID-1 bID-2 bID-1  bID]';    
end


 % 3rd: right most column
 % middle one
 for j=1:Ny
     NEle = NEle + 1;
     NEleToNode(NEle) = 5;
     aID = Ny+2 + (Nx-1)*(2*Ny+2) + 1 + (j-1)*2;
     bID = Ny+2+Nx*(2*Ny+2) + j;
     EleToNode(1:5, NEle) = [aID  bID bID+1  aID+2 aID+1]';
 end
 
 % top one
 NEle = NEle + 1;
 NEleToNode(NEle) = 4;
 aID = Ny+2 + (Nx)*(2*Ny+2) ;
 bID = NPts;
 EleToNode(1:4, NEle) = [aID aID-1 bID-1  bID]';

 
 
%  %% boundary conditions: Dirichlet & Neumann
% [NNBC,NBC, NDBC,DBC] = BCForMeshInRectDomain(xLeft, xRight, ...
%     yBottom, yTop, EleToNode, NEleToNode, VertX, VertY);
% 
% 
%  %% draw 2d mesh
% if Nx <20
%  NFig = NFig + 1; figure(NFig);
% draw2DMeshBC_allMesh(VertX, VertY, EleToNode, NEleToNode, NNBC,NBC, NDBC,DBC)
% %axis equal;
% end
% 
% 


    
% assign the perm
for k=1:NEle
	perm(1, k) = permOld_rotate(1,1);
	perm(2, k) = permOld_rotate(2,2);
	perm(3, k) = permOld_rotate(1,2);
end


        
   


% # of multipoint flux block
MulFlux = ones(NEle,1);
NMulFlux = sum(MulFlux);



%% remove repeating points.
oldPts = [VertX  VertY];
[numPts, pts, ptsMapID] = removeRepeatPts(oldPts);

VertX = pts(:,1); VertY = pts(:,2);
clear oldPts pts;

% update EleToNode
for k=1:NEle
    locNNode = NEleToNode(k);
    oldNodeID = EleToNode(1:locNNode, k);
    EleToNode(1:locNNode, k) = ptsMapID(oldNodeID)';
end







%% boundary conditions: Dirichlet & Neumann
[NNBC,NBC, NDBC,DBC] = BCForMeshInRectDomain(xLeft, xRight, ...
    yBottom, yTop, EleToNode, NEleToNode, VertX, VertY);

    



%% draw 2d mesh
% if Nx <20
 NFig = NFig + 1; figure(NFig);
draw2DMeshBC_allMesh(VertX, VertY, EleToNode, NEleToNode, NNBC,NBC, NDBC,DBC)
%axis equal;
% end













%% output the mesh
umFid = fopen(outFileName, 'wt');

fprintf(umFid, '%s\n','Number of nodes');
fprintf(umFid, '%7i\n', length(VertX));
for i=1:length(VertX)
    fprintf(umFid, '%7i %14.10f %14.10f \n', i, VertX(i), VertY(i));
end
fprintf(umFid, '%s\n','ENDOFSECTION');    


fprintf(umFid, '%s\n', 'Number of elements');
fprintf(umFid, '%7i\n', NEle);

fprintf(umFid, '%s\n', 'Number of vertices of elements');
for i=1:NEle
    fprintf(umFid, '%i \n' ,  NEleToNode(i));
end
fprintf(umFid, '%s\n','ENDOFSECTION'); 

fprintf(umFid, '%s\n', 'indices of vertices of elements and permeability');
for i=1:NEle            
    fprintf(umFid, ' %i \t ', i);

    for j=1:NEleToNode(i)
        fprintf(umFid, ' %i ', EleToNode(j, i));
    end

    fprintf(umFid, '\t  %f   %f   %f\n' ,  perm(1,i), perm(2,i),perm(3,i));
end
fprintf(umFid, '%s\n','ENDOFSECTION');    


%  Flux type
% if NMulFlux >0
    fprintf(umFid, '%s\n','Flux type');
    fprintf(umFid, '%s \n','multipoint flux');
    fprintf(umFid, '%d \n', NMulFlux);
    counter = 0;
    for i=1:NEle
        if MulFlux(i) > 0
            counter = counter + 1;
            fprintf(umFid, '%d \t %d\n', counter, i);
        end
    end
    fprintf(umFid, '%s\n','ENDOFSECTION');    
% end


%% output the boundary conditions

% Neumann BC
if NNBC > 0
    fprintf(umFid, '%s\n','Boundary Conditions');
    fprintf(umFid, '%s  %d \n','neumann', NNBC);
    for i=1:NNBC
        fprintf(umFid, '%d   %d   %d  \n', i, NBC(i,:));
    end
    fprintf(umFid, '%s\n','ENDOFSECTION');    
end

% Dirichlet BC
if NDBC > 0
    fprintf(umFid, '%s\n','Boundary Conditions');
    fprintf(umFid, '%s  %d \n','dirichlet', NDBC);
    for i=1:NDBC
        fprintf(umFid, '%d   %d   %d  \n', i, DBC(i,:));
    end
    fprintf(umFid, '%s\n','ENDOFSECTION');    

end


fclose(umFid);

