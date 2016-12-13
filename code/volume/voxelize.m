function [gridOUTPUT,varargout] = voxelize(gridX,gridY,gridZ,varargin)
% VOXELISE  Voxelise a 3D triangular-polygon mesh
%==========================================================================
% FILENAME:          VOXELISE.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@physics.cr.man.ac.uk
% DATE:              10th May 2010
% PURPOSE:           Voxelise a 3D triangular-polygon mesh
%
% USAGE:             [gridOUTPUT,varargout] = VOXELISE(gridX,gridY,gridZ,varargin)
%                    takes the mesh defined by <varargin> and voxelises it
%                    according to the grid defined by gridX, gridY and
%                    gridZ.
%
% INPUT PARAMETERS:  MANDATORY
%
%       gridX   - A 1xP array    - List of the grid X coordinates. 
%                 OR an integer  - Number of voxels in the grid in the X direction.
%
%       gridY   - A 1xQ array    - List of the grid Y coordinates.
%                 OR an integer  - Number of voxels in the grid in the Y direction.
%
%       gridZ   - A 1xR array    - List of the grid Z coordinates.
%                 OR an integer  - Number of voxels in the grid in the Z direction.
%
% INPUT PARAMETERS:  OPTIONAL (1 of the following is required)
%
%       STLin   - string        - Filename of the STL file.
%
%       meshX   - 3xN array     - List of the mesh X coordinates for the 3 vertices of each of the N triangular patches
%       meshY   - 3xN array     - List of the mesh Y coordinates for the 3 vertices of each of the N triangular patches
%       meshZ   - 3xN array     - List of the mesh Z coordinates for the 3 vertices of each of the N triangular patches
%
%       meshXYZ - Nx3x3 array   - The vertex positions for each facet, with:
%                                 - 1 row for each facet
%                                 - 3 columns for the x,y,z coordinates
%                                 - 3 pages for the three vertices
%
% OUTPUT PARAMETERS: MANDATORY
%
%       gridOUTPUT - 3D array of size (P,Q,R) - Voxelised data
%                                 1 => Inside the mesh
%                                 0 => Outside the mesh
%
% OUTPUT PARAMETERS: OPTIONAL
%
%       gridCOx - A 1xP array   - List of the grid X coordinates.
%       gridCOy - A 1xQ array   - List of the grid Y coordinates.
%       gridCOz - A 1xR array   - List of the grid Z coordinates.
%
% EXAMPLES:
%
%       To voxelise an STL file, the following command is used:
%       >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,STLin)
% 
%       To voxelise a patch-compatible mesh, the following command is used:
%       >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshX,meshY,meshZ)
%
%       To voxelise a mesh defined by a single Nx3x3 array, the following command is used:
%       >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshXYZ)
%
%       To also output the lists of X,Y,Z coordinates, the following command is used:
%       >>  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,STLin)
%
% REFERENCES:
%
%       This code uses a ray intersection method similar to the method
%       described by:
%     - Patil S and Ravi B.  Voxel-based Representation, Display and
%       Thickness Analysis of Intricate Shapes. Ninth International
%       Conference on Computer Aided Design and Computer Graphics (CAD/CG
%       2005)
%==========================================================================

%==========================================================================
% VERSION:  USER:  CHANGES:
% -------   -----  --------
% 100510    AHA    Original version
% 100514    AHA    Now works with non-STL input.  Changes also provide a
%                  significant speed improvement
% 100520    AHA    Now optionally output the grid x,y,z coordinates.
%                  Robustness also improved
%==========================================================================


%======================================================
% CHECK THE REQUIRED NUMBER OF OUTPUT PARAMETERS
%======================================================

if nargout~=1 && nargout~=4
  error('Incorrect number of output arguments.')
end

%======================================================
% READ INPUT PARAMETERS
%======================================================

% Whatever the input mesh format is, it is converted to an Nx3x3 array
% defining the vertex positions for each facet, with 1 row for each facet,
% 3 cols for the x,y,z coordinates, and 3 pages for the three vertices.

if nargin==4
  
  if ischar(varargin{1})
    STLin = varargin{1};
    % LOAD THE STL FILE
    %Identify whether the STL is ascii or binary:
    STLformat = IDENTIFY_stl_format(STLin);
    %Load the STL file:
    if strcmp(STLformat,'ascii')
      [meshCO] = READ_stlascii(STLin);
    elseif strcmp(STLformat,'binary')
      [meshCO] = READ_stlbinary(STLin);
    end %if
  else
    meshCO = varargin{1};
  end %if

elseif nargin==6

  meshX = varargin{1};
  meshY = varargin{2};
  meshZ = varargin{3};
 
  meshCO = zeros( size(meshX,2) , 3 , size(meshX,1) );

  meshCO(:,1,:) = reshape(meshX',size(meshX,2),1,3);
  meshCO(:,2,:) = reshape(meshY',size(meshY,2),1,3);
  meshCO(:,3,:) = reshape(meshZ',size(meshZ,2),1,3);
  
else
  
  error('Incorrect number of input arguments.')
  
end

%======================================================
% IDENTIFY THE MIN AND MAX X,Y,Z COORDINATES OF THE POLYGON MESH
%======================================================

meshXmin = min(min(meshCO(:,1,:)));
meshXmax = max(max(meshCO(:,1,:)));
meshYmin = min(min(meshCO(:,2,:)));
meshYmax = max(max(meshCO(:,2,:)));
meshZmin = min(min(meshCO(:,3,:)));
meshZmax = max(max(meshCO(:,3,:)));

%======================================================
% CHECK THE DIMENSIONS OF THE 3D OUTPUT GRID
%======================================================
% The output grid will be defined by the coordinates in gridCOx, gridCOy, gridCOz

if numel(gridX)>1
  if size(gridX,1)>size(gridX,2)   %gridX should be a row vector rather than a column vector
    gridCOx = gridX';
  else
    gridCOx = gridX;
  end
elseif numel(gridX)==1 && rem(gridX,1)==0   %If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
  voxwidth  = (meshXmax-meshXmin)/(gridX-1);
  gridCOx   = [meshXmin:voxwidth:meshXmax];
end

if numel(gridY)>1
  if size(gridY,1)>size(gridY,2)   %gridY should be a row vector rather than a column vector
    gridCOy = gridY';
  else
    gridCOy = gridY;
  end
elseif numel(gridY)==1 && rem(gridY,1)==0   %If gridX is a single integer (rather than a vector) then automatically create the list of y coordinates
  voxwidth  = (meshYmax-meshYmin)/(gridY-1);
  gridCOy   = [meshYmin:voxwidth:meshYmax];
end

if numel(gridZ)>1
  if size(gridZ,1)>size(gridZ,2)   %gridZ should be a row vector rather than a column vector
    gridCOz = gridZ';
  else
    gridCOz = gridZ;
  end
elseif numel(gridZ)==1 && rem(gridZ,1)==0   %If gridZ is a single integer (rather than a vector) then automatically create the list of z coordinates
  voxwidth  = (meshZmax-meshZmin)/(gridZ-1);
  gridCOz   = [meshZmin:voxwidth:meshZmax];
end

%Count the number of voxxels in each direction:
voxcountX = numel(gridCOx);
voxcountY = numel(gridCOy);
voxcountZ = numel(gridCOz);

% Prepare array to hold the voxelised data:
gridOUTPUT = zeros(voxcountX,voxcountY,voxcountZ);

%======================================================
% PERFORM VARIOUS OTHER CHECKS AND INITIALISATIONS
%======================================================

%Check that the output grid is large enough to cover the mesh:
if min(gridCOx)>meshXmin || max(gridCOx)<meshXmax
  error('The output grid is too small.  It does not completely enclose the mesh in the X direction.')
elseif min(gridCOy)>meshYmin || max(gridCOy)<meshYmax
  error('The output grid is too small.  It does not completely enclose the mesh in the Y direction.')
elseif min(gridCOz)>meshZmin || max(gridCOz)<meshZmax
  error('The output grid is too small.  It does not completely enclose the mesh in the Z direction.')
end

%Identify the min and max x,y coordinates (pixels) of the mesh:
meshXminp = find(abs(gridCOx-meshXmin)==min(abs(gridCOx-meshXmin)));
meshXmaxp = find(abs(gridCOx-meshXmax)==min(abs(gridCOx-meshXmax)));
meshYminp = find(abs(gridCOy-meshYmin)==min(abs(gridCOy-meshYmin)));
meshYmaxp = find(abs(gridCOy-meshYmax)==min(abs(gridCOy-meshYmax)));

%Make sure min < max for the mesh coordinates:
if meshXminp > meshXmaxp
  [meshXminp,meshXmaxp] = deal(meshXmaxp,meshXminp);
end %if
if meshYminp > meshYmaxp
  [meshYminp,meshYmaxp] = deal(meshYmaxp,meshYminp);
end %if

%Identify the min and max x,y,z coordinates of each facet:
meshCOmin = min(meshCO,[],3);
meshCOmax = max(meshCO,[],3);

%======================================================
% TURN OFF DIVIDE-BY-ZERO WARNINGS
%======================================================
%This prevents the Y1predicted, Y2predicted, Y3predicted and YRpredicted
%calculations creating divide-by-zero warnings.  Suppressing these warnings
%doesn't affect the code, because only the sign of the result is important.
%That is, 'Inf' and '-Inf' results are ok.
%The warning will be returned to its original state at the end of the code.
warningrestorestate = warning('query', 'MATLAB:divideByZero');
warning off MATLAB:divideByZero

%======================================================
% VOXELISE THE MESH
%======================================================

correctionLIST = [];   %Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
                       %it ought to be relatively small should not incur too much of a speed penalty.
                       
% Loop through each x,y pixel.
% The mesh will be voxelised by passing rays in the z-direction through
% each x,y pixel, and finding the locations where the rays cross the mesh.
for loopY = meshYminp:meshYmaxp
  for loopX = meshXminp:meshXmaxp
      
    % - 1 - Find which mesh facets could possibly be crossed by the ray:
    possibleCROSSLIST = find( meshCOmin(:,1) <= gridCOx(loopX) );
    possibleCROSSLIST = possibleCROSSLIST(find( meshCOmax(possibleCROSSLIST,1) >= gridCOx(loopX) ));
    possibleCROSSLIST = possibleCROSSLIST(find( meshCOmin(possibleCROSSLIST,2) <= gridCOy(loopY) ));
    possibleCROSSLIST = possibleCROSSLIST(find( meshCOmax(possibleCROSSLIST,2) >= gridCOy(loopY) ));

    if isempty(possibleCROSSLIST)==0  %Only continue the analysis if some nearby facets were actually identified
          
      % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
          
      % GENERAL METHOD:
      % 1. Take each edge of the facet in turn.
      % 2. Find the position of the opposing vertex to that edge.
      % 3. Find the position of the ray relative to that edge.
      % 4. Check if ray is on the same side of the edge as the opposing vertex.
      % 5. If this is true for all three edges, then the ray definitely passes through the facet.
      %
      % NOTES:
      % 1. If the ray crosses exactly on an edge, this is counted as crossing the facet.
      % 2. If a ray crosses exactly on a vertex, this is also taken into account.
      
      facetCROSSLIST = [];   %Prepare to record all facets which are crossed by the ray.  This array is built on-the-fly, but since
                             %it ought to be relatively small should not incur too much of a speed penalty.
                             
      vertexCOzCROSS = [];   %Prepare to record all vertices that are crossed by the ray.  This array is built on-the-fly, but since
                             %it ought to be relatively small should not incur too much of a speed penalty.
      
      for loopCHECKFACET = possibleCROSSLIST'
        
        % Check if the ray crosses any of the facet vertices (but exclude facets where the ray is travels along an edge):
        % The parameter <checkvertexcrossing> contains one value for each of the x and y coordinates.  If the ray is coincident with the vertex,
        % then the ray's x and y coordinates will match the vertex's x and y coordinates, and <checkvertexcrossing> will therefore have a value of 2.
        checkvertexcrossing = (meshCO(loopCHECKFACET,1,:)==gridCOx(loopX)) + (meshCO(loopCHECKFACET,2,:)==gridCOy(loopY));
        
        if max(checkvertexcrossing)==2  &&  sum(checkvertexcrossing==2)==1 % Ray crosses at one vertex of the facet only

          vertexcrossed  = find( meshCO(loopCHECKFACET,1,:)==gridCOx(loopX) & meshCO(loopCHECKFACET,2,:)==gridCOy(loopY) );
          vertexCOzCROSS = [vertexCOzCROSS,meshCO(loopCHECKFACET,3,vertexcrossed)];

        elseif max(checkvertexcrossing)==2  &&  sum(checkvertexcrossing==2)>1 %Ray crosses at two vertices of the facet, ie. along an edge
          %The ray is not counted as meeting this facet.

        else %Ray does not cross any of the vertices of the facet
        
          %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.  If so, let testVn=1
        
          Y1predicted = meshCO(loopCHECKFACET,2,2) - ((meshCO(loopCHECKFACET,2,2)-meshCO(loopCHECKFACET,2,3)) * (meshCO(loopCHECKFACET,1,2)-meshCO(loopCHECKFACET,1,1))/(meshCO(loopCHECKFACET,1,2)-meshCO(loopCHECKFACET,1,3)));
          YRpredicted = meshCO(loopCHECKFACET,2,2) - ((meshCO(loopCHECKFACET,2,2)-meshCO(loopCHECKFACET,2,3)) * (meshCO(loopCHECKFACET,1,2)-gridCOx(loopX))/(meshCO(loopCHECKFACET,1,2)-meshCO(loopCHECKFACET,1,3)));
          if (Y1predicted > meshCO(loopCHECKFACET,2,1) && YRpredicted > gridCOy(loopY)) || (Y1predicted < meshCO(loopCHECKFACET,2,1) && YRpredicted < gridCOy(loopY))
            %The ray is on the same side of the 2-3 edge as the 1st vertex.
            testV1 = 1;
          elseif YRpredicted == gridCOy(loopY) || isnan(YRpredicted)
            %The ray exactly crosses the 2-3 edge of the facet.
            testV1 = 1;
          else
            %The ray is on the opposite side of the 2-3 edge to the 1st vertex.
            testV1 = 0;
          end %if
          
          Y2predicted = meshCO(loopCHECKFACET,2,3) - ((meshCO(loopCHECKFACET,2,3)-meshCO(loopCHECKFACET,2,1)) * (meshCO(loopCHECKFACET,1,3)-meshCO(loopCHECKFACET,1,2))/(meshCO(loopCHECKFACET,1,3)-meshCO(loopCHECKFACET,1,1)));
          YRpredicted = meshCO(loopCHECKFACET,2,3) - ((meshCO(loopCHECKFACET,2,3)-meshCO(loopCHECKFACET,2,1)) * (meshCO(loopCHECKFACET,1,3)-gridCOx(loopX))/(meshCO(loopCHECKFACET,1,3)-meshCO(loopCHECKFACET,1,1)));
          if (Y2predicted > meshCO(loopCHECKFACET,2,2) && YRpredicted > gridCOy(loopY)) || (Y2predicted < meshCO(loopCHECKFACET,2,2) && YRpredicted < gridCOy(loopY))
            %The ray is on the same side of the 3-1 edge as the 2nd vertex.
            testV2 = 1;
          elseif YRpredicted == gridCOy(loopY) || isnan(YRpredicted)
            %The ray exactly crosses the 3-1 edge of the facet.
            testV2 = 1;
          else
            %The ray is on the opposite side of the 3-1 edge to the 2nd vertex.
            testV2 = 0;
          end %if

          Y3predicted = meshCO(loopCHECKFACET,2,1) - ((meshCO(loopCHECKFACET,2,1)-meshCO(loopCHECKFACET,2,2)) * (meshCO(loopCHECKFACET,1,1)-meshCO(loopCHECKFACET,1,3))/(meshCO(loopCHECKFACET,1,1)-meshCO(loopCHECKFACET,1,2)));
          YRpredicted = meshCO(loopCHECKFACET,2,1) - ((meshCO(loopCHECKFACET,2,1)-meshCO(loopCHECKFACET,2,2)) * (meshCO(loopCHECKFACET,1,1)-gridCOx(loopX))/(meshCO(loopCHECKFACET,1,1)-meshCO(loopCHECKFACET,1,2)));
          if (Y3predicted > meshCO(loopCHECKFACET,2,3) && YRpredicted > gridCOy(loopY)) || (Y3predicted < meshCO(loopCHECKFACET,2,3) && YRpredicted < gridCOy(loopY))
            %The ray is on the same side of the 1-2 edge as the 3rd vertex.
            testV3 = 1;
          elseif YRpredicted == gridCOy(loopY) || isnan(YRpredicted)
            %The ray exactly crosses the 1-2 edge of the facet.
            testV3 = 1;
          else
            %The ray is on the opposite side of the 1-2 edge to the 3rd vertex.
            testV3 = 0;
          end %if
        
          %If the ray is on the correct side of all 3 edges then it passes through the facet
          if testV1==1 && testV2==1 && testV3==1
            facetCROSSLIST = [facetCROSSLIST,loopCHECKFACET];
          end %if

        end %if
      
      end %for loopCHECKFACET = possibleCROSSLIST'
    

      % - 3 - Find the z coordinate of the locations where the ray crosses each facet:

      gridCOzCROSS = zeros(size(facetCROSSLIST));
      for loopFINDZ = facetCROSSLIST

        % METHOD:
        % 1. Define the equation describing the plane of the facet.  For a
        % more detailed outline of the maths, see:
        % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
        %    Ax + By + Cz + D = 0
        %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
        %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
        %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
        %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
        % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.

        planecoA = meshCO(loopFINDZ,2,1)*(meshCO(loopFINDZ,3,2)-meshCO(loopFINDZ,3,3)) + meshCO(loopFINDZ,2,2)*(meshCO(loopFINDZ,3,3)-meshCO(loopFINDZ,3,1)) + meshCO(loopFINDZ,2,3)*(meshCO(loopFINDZ,3,1)-meshCO(loopFINDZ,3,2));
        planecoB = meshCO(loopFINDZ,3,1)*(meshCO(loopFINDZ,1,2)-meshCO(loopFINDZ,1,3)) + meshCO(loopFINDZ,3,2)*(meshCO(loopFINDZ,1,3)-meshCO(loopFINDZ,1,1)) + meshCO(loopFINDZ,3,3)*(meshCO(loopFINDZ,1,1)-meshCO(loopFINDZ,1,2)); 
        planecoC = meshCO(loopFINDZ,1,1)*(meshCO(loopFINDZ,2,2)-meshCO(loopFINDZ,2,3)) + meshCO(loopFINDZ,1,2)*(meshCO(loopFINDZ,2,3)-meshCO(loopFINDZ,2,1)) + meshCO(loopFINDZ,1,3)*(meshCO(loopFINDZ,2,1)-meshCO(loopFINDZ,2,2));
        planecoD = - meshCO(loopFINDZ,1,1)*(meshCO(loopFINDZ,2,2)*meshCO(loopFINDZ,3,3)-meshCO(loopFINDZ,2,3)*meshCO(loopFINDZ,3,2)) - meshCO(loopFINDZ,1,2)*(meshCO(loopFINDZ,2,3)*meshCO(loopFINDZ,3,1)-meshCO(loopFINDZ,2,1)*meshCO(loopFINDZ,3,3)) - meshCO(loopFINDZ,1,3)*(meshCO(loopFINDZ,2,1)*meshCO(loopFINDZ,3,2)-meshCO(loopFINDZ,2,2)*meshCO(loopFINDZ,3,1));

        if abs(planecoC) < 1e-14
          planecoC=0;
        end
      
        gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*gridCOx(loopX) - planecoB*gridCOy(loopY)) / planecoC;
        
      end %for

      %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
      gridCOzCROSS = gridCOzCROSS( gridCOzCROSS>=meshZmin-1e-12 & gridCOzCROSS<=meshZmax+1e-12 );
      
      %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
      gridCOzCROSS = round(gridCOzCROSS*1e12)/1e12;
      gridCOzCROSS = unique(gridCOzCROSS);
      
      %Take the unique list of crossed vertices:
      vertexCOzCROSS = unique(vertexCOzCROSS);
      
      
      % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:

      if (rem(numel(gridCOzCROSS),2)==0 && rem(numel(vertexCOzCROSS),2)==0)  ||  (rem(numel(gridCOzCROSS),2)==1 && rem(numel(vertexCOzCROSS),2)==1)     % (Even number of both facets AND vertices crossed)  OR  (Odd number of both facets AND vertices crossed)

        allCOzCROSS = [gridCOzCROSS,vertexCOzCROSS];
        allCOzCROSS = sort(allCOzCROSS);
        for loopASSIGN = 1:(numel(allCOzCROSS)/2)
          voxelsINSIDE = find( (gridCOz>allCOzCROSS(2*loopASSIGN-1) & gridCOz<allCOzCROSS(2*loopASSIGN)) );
          gridOUTPUT(loopX,loopY,voxelsINSIDE) = 1;
        end %for
        
      elseif numel(vertexCOzCROSS)~=0 || numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
          
        correctionLIST = [ correctionLIST; loopX,loopY ];
        
      end %if

    end %if

  end %for
end %for

%======================================================
% RESTORE DIVIDE-BY-ZERO WARNINGS TO THE ORIGINAL STATE
%======================================================

warning(warningrestorestate)

%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(correctionLIST,1);

if countCORRECTIONLIST>0
    
  %Add a one-pixel border around the x and y edges of the array.  This
  %prevents an error if the code tries to interpolate a ray at the edge of
  %the x,y grid.
  gridOUTPUT     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),gridOUTPUT,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
  correctionLIST = correctionLIST + 1;
  
  for loopC = 1:countCORRECTIONLIST
    gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),:) = mean( [ (gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)-1,:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2),:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)+1,:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)-1,:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)+1,:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)-1,:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2),:)) ,...
                                                                            (gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)+1,:)) ,...
                                                                          ] ,2);
    gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),:) = round( gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),:) );
  end %for

  %Remove the one-pixel border surrounding the array.
  gridOUTPUT = gridOUTPUT(2:end-1,2:end-1,:);
  
  warning('voxelize:correction',...
      'The voxelisation of %d rays (%5.1f%% of all rays) did not give a clear result.\nThese rays were therefore computed by interpolating from the surrounding rays.',...
      countCORRECTIONLIST,100*countCORRECTIONLIST/(voxcountX*voxcountY)) 

end %if


%======================================================
% PREPARE THE OUTPUT PARAMETERS
%======================================================

if nargout==4
  
  varargout(1) = {gridCOx};
  varargout(2) = {gridCOy};
  varargout(3) = {gridCOz};
  
end
