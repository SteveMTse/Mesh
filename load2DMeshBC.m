%%  We waste memory a bit to make coding easier. (array is used).
%%
%%  load 2D  mesh data file with format: (some parts rely on special structure of 2D)
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
%%  multipoint flux,
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
%% ...
%%  ENDOFSECTION
%%  Author: Qian-Yong Chen
%%
%% It loads into mesh data: 
%%                         NEle        = number of elements
%%                         NNode     = number of nodes
%%                         VertX      = X coordinates of nodes
%%                         VertY      = Y coordinates of nodes
%%                         EleEdgeBCType = integer flag for different
%%                                              types of boundary conditions
%%                                              at node, 
%%              0- interior, 1- Dirichelet BC, 2- Neumann BC 
%%                         EleToNode   = [aN x NEle] = list of nodes in each element
%%
%% and computes the mapping from global to local indices:
%%                         EleToEle = [aN x NEle] = element to element list
%%                                     [i,j]-->j-th neighbour of i-th element 
%%                        EleToFace = [aN x NNel] = element to face list
%%                                      [i,j] ---> local index of the face shared by i-th element
%%                                      and [i,j]-the element.
%%                          NodeToNode= [M x NNode]: node to its neighbour node list
%%                          NodeToEle = [M x NNode]: node to element lists
%%                          NodeToFace= [M x NNode]: node to index of faces of elements.
%%
%% It also reorders the node and element order so that the graphs of the element 
%% connectivities and node connectivities are bandwidth minimized (using Reverse Cuthill McKee)
%%  Not done for this staircase like mesh yet.


 clear all;
close all;

% Mesh data file name
%  inFileName ='triMesh_layer.cm';
% %  inFileName ='quaMesh_layer.cm';
%  
% inFileName = 'pinchOutMesh.cm';
% 
% 
  inFileName = 'quaMesh_layer.cm'
% 
%  inFileName = 'quaMesh_loadPerm.cm'

%inFileName = 'polyMesh.cm'




NFig = 0;


%% load data file
umFid = fopen(inFileName, 'rt');

line = fgetl(umFid);
% number of nodes
NNode = fscanf(umFid, '%d\n',1);

% read node coordinates 
VertX = zeros(1,NNode);
VertY = zeros(1,NNode);

for i = 1:NNode
  line = fgetl(umFid);
  tmpx = sscanf(line, '%lf');
  VertX(i) = tmpx(2);
  VertY(i) = tmpx(3);
end

line = fgetl(umFid);
line = fgetl(umFid);

% number of elements
NEle = fscanf(umFid, '%d\n',1);

% number of vertices of elements
NEleToNode = zeros(NEle,1);
line = fgetl(umFid);
for k=1:NEle
    NEleToNode(k) = fscanf(umFid, '%d\n',1);
end

NEleToNode_max = max(NEleToNode);

%%@!$^^, Here comes the waste of memory

% read element to node connectivity  & permeability
EleToNode = zeros(NEleToNode_max, NEle);
ElmPerm = zeros(2,2,NEle);

line = fgetl(umFid);
line = fgetl(umFid);
for k = 1:NEle
  line   = fgetl(umFid);
  tmpcon = sscanf(line, '%lf');

  locNEleNode  = NEleToNode(k);
  EleToNode(1:locNEleNode, k ) = tmpcon(2:1+locNEleNode );
  
  ElePerm(1,1,k) = tmpcon(2+locNEleNode);
  ElePerm(2,2,k) = tmpcon(3+locNEleNode);
  ElePerm(1,2,k) = tmpcon(4+locNEleNode);
  ElePerm(2,1,k) = tmpcon(4+locNEleNode);    
%   ElePerm(k,1:3) = tmpcon(2+locNEleNode :4+locNEleNode )';
end


%read into element type, multipoint flux or two point flux
MulFlux = zeros(NEle,1);
line = fgetl(umFid);
line = fgetl(umFid);
line = fgetl(umFid);
NMulFlux = fscanf(umFid,'%d\n',1);
for k=1:NMulFlux
    line = fgetl(umFid);
              
    tmpID = sscanf(line, '%d');
    
    eleID = tmpID(2);
    MulFlux(eleID) = 1;
end



line = fgetl(umFid);

% read boundary conditions
DIRICHLET = 1;
NEUMANN = 2;

BCEdgeType = zeros( NEleToNode_max, NEle);
NDBC=0; NNBC=0; DBC = zeros(1,2); NBC = zeros(1,2);

line = fgetl(umFid);
line = fgetl(umFid);

while line ~= -1
    
    if ~isempty(strmatch('neumann', line))      %% Neumann BC
        line = fgetl(umFid);
        while isempty(strmatch('ENDOFSECTION', line))
            tmpID = sscanf(line, '%d');
            eleID = tmpID(2);
            faceID= tmpID(3);
            BCEdgeType(faceID, eleID) = NEUMANN;
            NNBC = NNBC + 1;
            NBC(NNBC,1) = eleID; NBC(NNBC,2) = faceID;
            
            line = fgetl(umFid);
        end
        line = fgetl(umFid);
        line = fgetl(umFid);
        
    elseif  ~isempty(strmatch('dirichlet', line))        %% Dirichlet BC
        line = fgetl(umFid);
        while isempty(strmatch('ENDOFSECTION', line))
            tmpID = sscanf(line, '%d');
            eleID = tmpID(2);
            faceID= tmpID(3);
            BCEdgeType(faceID, eleID) = DIRICHLET;
            NDBC = NDBC + 1;
            DBC(NDBC,1) = eleID; DBC(NDBC,2) = faceID;
            
            line = fgetl(umFid);
        end
        line = fgetl(umFid);
        line = fgetl(umFid);
        
    end
    
end

fclose(umFid);


% waste of memeory
EleEdgeBCType = zeros(NEleToNode_max, NEle);



% 
% %draw the mesh
% NFig = NFig + 1;    figure(NFig);
% draw2DMesh_patch(VertX', VertY', EleToNode, NEleToNode);
% axis equal;
% % 
% % 
% % draw flux type & BC
% NFig = NFig + 1; figure(NFig);
% draw2DMeshBC_FluxType_allMesh(VertX, VertY, EleToNode, NEleToNode, NNBC,NBC, NDBC,DBC,MulFlux)
% axis equal;






% %% draw the element & node index
% NFig = NFig + 1;
% figure(NFig);
% draw2DMeshID(VertX, VertY, EleToNode);


% % reorder nodes or element to be banded.  
%  umEleToNode = zeros(NEle*3, 3);
%  counter = 0;
%  for k=1:NEle
%      for j=1:NEleToNode(k)
%          counter = counter + 1;
%          umEleToNode(counter, 1) = k;         
%          umEleToNode(counter, 2) = EleToNode(j, k);
%          umEleToNode(counter, 3) = 1;         
%      end
%  end
%  SparseEleToNode = spconvert(umEleToNode);
%  SparseEleToEle = SparseEleToNode*SparseEleToNode';
%  SparseNodeToNode = SparseEleToNode'*SparseEleToNode;
%  
%  NodeMap = symrcm(SparseNodeToNode);
%  EleMap = symrcm(SparseEleToEle);  
%  
%  
% %  % reorder nodes
% %  VertX = VertX(NodeMap); 
% %  VertY = VertY(NodeMap);
% %  [tmpNodeMap, invNodeMap] = sort(NodeMap);      %% not a smart way.
% %  for k=1:NEle
% %      IDs = EleToNode(k,:);
% %      IDs = invNodeMap(IDs);     
% %      EleToNode(k,:) = IDs;    
% %  end
% % clear tmpNodeMap invNodeMap
% 
% 
% %  %% draw the element & node index
% % NFig = NFig + 1;
% % figure(NFig);
% % draw2DMeshID(VertX, VertY, EleToNode); 
%  
% EleToNode = EleToNode(:,EleMap);    % update EleToNode 
% ElePerm = ElePerm(:,:, EleMap);     % update Perm
% MulFlux = MulFlux(EleMap);          % update flux type 
% NEleToNode = NEleToNode(EleMap);    % update NEleToNode
% 
% %% one way to update the boundary edges after reordering the elements.
% [tmpMap, IDs] = sort(EleMap);
% IDs = IDs';
% if NDBC >0
%     DBC(:,1) = IDs(DBC(:,1));
% end
% if NNBC > 0  
%     NBC(:,1) = IDs(NBC(:,1));    
% end  




for k=1:NDBC
    eleID = DBC(k,1);
    faceID= DBC(k,2);
    EleEdgeBCType(faceID, eleID) = 1;
end
for k=1:NNBC
    eleID = NBC(k,1);
    faceID= NBC(k,2);
    EleEdgeBCType(faceID,eleID) = 2;
end








% %draw the mesh
% NFig = NFig + 1;    figure(NFig);
% draw2DMesh_patch(VertX', VertY', EleToNode, NEleToNode);
% axis equal;


% % % MulFlux(:) = 0;
% % draw flux type & BC
% NFig = NFig + 1; figure(NFig);
% draw2DMeshBC_FluxType_allMesh(VertX, VertY, EleToNode, NEleToNode, NNBC,NBC, NDBC,DBC,MulFlux)
% axis equal;


% 
% %% draw the permeability of the elements
% NFig = NFig + 1;
% figure(NFig);
% draw2DMeshPerm(VertX, VertY, EleToNode, NEleToNode, ElePerm);
% axis equal;



% %  %% draw the element & node index
% NFig = NFig + 1;
% figure(NFig);
% draw2DMeshID(VertX, VertY, EleToNode); 


 % now compute the EleToEle and EleToFace arrays 
 EleToEle = zeros(NEleToNode_max, NEle);
 EleToFace = zeros(NEleToNode_max, NEle);
 
 sumNEleToNode = 0;
 for k=1:NEle
     for j=1:NEleToNode(k)
         umEleToNode(sumNEleToNode + j,1) = k;         
         umEleToNode(sumNEleToNode + j,2) = EleToNode(j, k);
         umEleToNode(sumNEleToNode + j,3) = 1;         
     end
     sumNEleToNode = sumNEleToNode + NEleToNode(k);
 end
  
 SparseEleToNode = spconvert(umEleToNode);
 SparseEleToEle = SparseEleToNode*SparseEleToNode';
%  
%  % To compute node to neighbour elements' node list
%  SparseNodeToNode = SparseEleToNode'*SparseEleToNode;
%  NodeBandWidth = max(sum(SparseNodeToNode>0))-1;
 

 tmpEleToEle = zeros(NEleToNode_max, NEle);

 % being careful to follow the column-first storage in the sparse matrix
 for k=1:NEle
  conn = find(SparseEleToEle(:,k)==2);
  dims = size(conn);
  
  tmpEleToEle(1:dims(1), k) = conn;
 end

 
 
% Reorder tmpEleToEle(i,:) for each i, and save as EleToEle(i,:) such that it is in counterclockwise direction
% and EleToEle(i,1) share the face 1-2 with i-th element.

% find EleToEle and Eletoface connectivities
for k=1:NEle        
    locNEleNode  = NEleToNode(k);
    NodeID = EleToNode(1:locNEleNode, k);         % nodes of k-th element
    NodeID(1+locNEleNode ) = NodeID(1);
    
    neighbourEle = tmpEleToEle(:, k);        % neighbour elements of k-th element    
    neighbourEle_M = max(neighbourEle, 1);  % to make sure the  index is always positive
         
    
    %% find the edge share by k-th element and tmpEleToEle(k,j)
    for j=1:locNEleNode  
        tmpID = neighbourEle_M(j);
        %% check if they share edge i (local index) of k-th element.
        for i=1:locNEleNode             
            oneNode = NodeID(i);       % i-th edge consists nodes i & i+1 (local index)
            twoNode = NodeID(i+1);
            
            EleToEle(i, k) = EleToEle(i, k) + SparseEleToNode(tmpID, oneNode)*...
                SparseEleToNode(tmpID, twoNode) * neighbourEle(j);
        end 
    end
        

    %% generate EleToFace   i-face
    for i=1:locNEleNode 
        nbEleID = EleToEle(i,k);
        if (nbEleID)            %% if it has i-th neighbour element.

            nbLocNEleNode  = NEleToNode(nbEleID);
            nbNodeID = EleToNode(1:nbLocNEleNode, nbEleID);
            nbNodeID(1+nbLocNEleNode ) = nbNodeID(1);
            
            for j=1:nbLocNEleNode        %% check if i-th face of k-th element is the j-th face of 
                                %% the i-th neighbour element, ie, EleToEle(k,i)
                EleToFace(i, k) = EleToFace(i, k) + SparseEleToNode(k, nbNodeID(j))...
                    *SparseEleToNode(k, nbNodeID(j+1))*j;
            end
        end
    end
    
 end

 % free some variables
 clear SparseEleToEle SparseEleToNode umEleToNode;
 clear tmpEleToEle SparseNodeToNode NodeMap EleMap;
 
 %  %% draw EleToEle
% NFig = NFig + 1;
% figure(NFig);
%  drawEleToEle2D(EleToEle,EleToNode,NEleToNode,VertX,VertY);
%  
%  
% %% draw EleToFace
% NFig = NFig + 1;
% figure(NFig);
% drawEleToFace2D(EleToFace,EleToNode,NEleToNode, VertX,VertY);

 
 
 
 
 % Build NodeToNode, NodeToFace, NodeToEle, each row is arranged in the counterclockwise direction.
 % THE FOLLOWING algorithm relies on the fact EleToNode is arranged in counterclockwise direction for each
 % element.
 
 % NodeType:    
 %
 %      =1 --  corner boundary node
 %       2 --  normal boundary node
 %       3 --  interior non-matching node
 %       4 -- normal interior node   

 
% cheat
% physical domain,  should match the domain in mesh generation
xLeft = min(VertX); xRight = max(VertX);
yBottom = min(VertY); yTop = max(VertY);


NodeType = compNodeType(VertX, VertY, xLeft, xRight, yBottom, yTop);


% NFig=NFig+1;
% plot(VertX, VertY,'.');
% for i=1:NNode
%    textNNode(i,:) = num2str(NodeType(i));
% end
% text(VertX', VertY',  textNNode,'FontSize',16);
% axis equal;

 
 
 
 NNodeToNode = zeros(NNode,1);
 for k=1:NEle
     locID = EleToNode(1:NEleToNode(k), k);
     
     NNodeToNode(locID) = 1 + NNodeToNode(locID);
 end 
% misNodeBandWidth = max(NNodeToNode);
 

% % draw NNodeToNode
% NFig = NFig + 1;    figure(NFig);
% draw2DMesh_patch(VertX', VertY', EleToNode, NEleToNode);
% hold;
% for i=1:NNode
%    textNNode(i,:) = num2str(NNodeToNode(i));
% end
% text(VertX', VertY',  textNNode,'FontSize',16);
% axis equal;
% 


%% NodeToNode is NOT a connectivity list of all nodes. Instead,  i-row is the list of its next node in
%% the element which has node i (global index).
%%
%  % a little waste of memory
 NodeToNode = cell(NNode,1);
 NodeToEle = cell(NNode,1);
 NodeToFace = cell(NNode,1);
 
 
 tmpNNodeToNode = zeros(NNode,1);
 for k=1:NEle
     locNEleToNode = NEleToNode(k); 
     locID = EleToNode(1:locNEleToNode, k);
     locID(1+locNEleToNode) = locID(1);
     
     for i=1:locNEleToNode
        NodeToEle{locID(i),1} = [NodeToEle{locID(i),1}  k];
        NodeToFace{locID(i),1} = [NodeToFace{locID(i),1}  i];        
        NodeToNode{locID(i),1} = [NodeToNode{locID(i),1}  locID(i+1)];
     end
 end
 
 
 
 
 
 
 
 
 
 
 
 % arrange NodeToNode, NodeToEle, NodeToFace in counterclockwise direction.
 for i=1:NNode
     locNNodeToNode = NNodeToNode(i);
     locNodeID = NodeToNode{i,1};
  

    centerCoor = [VertX(i);  VertY(i)];    
    nbCoor = [VertX(locNodeID);  VertY(locNodeID)];
    [CCOrder] = sort2DPtsAroundCenter(centerCoor, nbCoor);

    % update NodeToNode
    NodeToNode{i,1} = locNodeID(CCOrder); 

    % update NodeToEle
    locEleID = NodeToEle{i,1}; 
    NodeToEle{i,1} = locEleID(CCOrder);
    
    
    % update NodeToFace
    locFaceID = NodeToFace{i,1};
    NodeToFace{i,1} = locFaceID(CCOrder);

    
    
    % rearrange NodeToNode, NodeToEle, NodeToFace for boundary node & abnormal interior node
    % such that it starts  'from  interor and ends at interior'
    % or  the element with 'node angle' = pi. See notes for details.
    if NodeType(i) <= 3        
        newLocNodeID = locNodeID(CCOrder);
   
        newCCOrder = 1:locNNodeToNode;
        if NodeType(i) <= 2     %% boundary node
            for j=1:locNNodeToNode
                nbNodeID = newLocNodeID(j);
                
                tmpEleID = NodeToEle{i,1}(j);
                tmpFaceID = NodeToFace{i,1}(j);

                if EleEdgeBCType(tmpFaceID, tmpEleID) > 0  %---> starting from boundary node
                    startID = j;
                end
            end
                    
        else                    %% abnormal interior node
            stop
        end
        
        % startID,..., locNNodeToNode, 1,...startID-1
        newCCOrder = [startID:locNNodeToNode   1:(startID-1)];
        
        %update NodeToNode...
        NodeToNode{i,1} = newLocNodeID(newCCOrder); 
    
        newLocEleID = locEleID(CCOrder);
        NodeToEle{i,1} = newLocEleID(newCCOrder);
        
        newLocFaceID = locFaceID(CCOrder);
        NodeToFace{i,1} = newLocFaceID(newCCOrder);
    end
end
 







 
%   % draw the element & node index
% NFig = NFig + 1;
% figure(NFig);
% draw2DMeshID(VertX, VertY, EleToNode, NEleToNode); 

%  
% %% to print out data
% pNodeToNode = zeros(NNode, misNodeBandWidth+1);
% pNodeToNode(:,1) = (1:NNode)'; 
% pNodeToNode(:,2:1+misNodeBandWidth) = NodeToNode;
%  
%  
%  %% to print out data
% pNodeToEle = zeros(NNode, misNodeBandWidth+1);
% pNodeToEle(:,1) = (1:NNode)'; 
% pNodeToEle(:,2:1+misNodeBandWidth) = NodeToEle;
% 
% 
% %% to print out data
% pNodeToFace = zeros(NNode, misNodeBandWidth+1);
% pNodeToFace(:,1) = (1:NNode)'; 
% pNodeToFace(:,2:1+misNodeBandWidth) = NodeToFace;
% % 
%  






% %% draw 2d mesh
if NEle <= 400
NFig = NFig + 1; figure(NFig);
draw2DMeshBC_allMesh(VertX, VertY, EleToNode, NEleToNode, NNBC,NBC, NDBC,DBC)
axis equal;
end