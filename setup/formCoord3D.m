function [etpl,coord] = formCoord3D(nelsx,nelsy,nelsz,lx,ly,lz)

%Three dimensional finite element grid generation
%--------------------------------------------------------------------------
% Description:
% Function to generate a 3D finite element grid of linear quadrilateral
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD3D(nelsx,nelsy,nelsz,lx,ly,lz)
%--------------------------------------------------------------------------
% Input(s):
% nelsx - number of elements in the x direction
% nelsy - number of elements in the y direction
% nelsz - number of elements in the z direction
% lx    - length in the x direction
% ly    - length in the y direction
% lz    - length in the y direction
%--------------------------------------------------------------------------
% Ouput(s);
% etpl  - element topology
% coord - nodal coordinates
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------

nels  = nelsx*nelsy*nelsz;                                                 % number of elements
nodes = (nelsx+1)*(nelsy+1)*(nelsz+1);                                      % number of nodes

%% node generation
coord = zeros(nodes,3);                                                     % zero coordinates
node  = 0;                                                                  % zero node counter
for k=0:nelsz
    z=lz*k/nelsz;
    for j=0:nelsy
        y=ly*j/nelsy;
        for i=0:nelsx
            node=node+1;
            x=lx*i/nelsx;
            coord(node,:)=[x y z];
        end
    end
end
%% element generation
etpl = zeros(nels,8);                                                       % zero element topology
nel  = 0;                                                                   % zero element counter
for nelz=1:nelsz
    for nely=1:nelsy
        for nelx=1:nelsx
            nel=nel+1;
            etpl(nel,1)=(nelz-1)*(nelsy+1)*(nelsx+1)+(nely-1)*(nelsx+1)+nelx;
            etpl(nel,2)=(nelz-1)*(nelsy+1)*(nelsx+1)+nely*(nelsx+1)+nelx;
            etpl(nel,3)=etpl(nel,2)+1;
            etpl(nel,4)=etpl(nel,1)+1;
            etpl(nel,5)=nelz*(nelsy+1)*(nelsx+1)+(nely-1)*(nelsx+1)+nelx;
            etpl(nel,6)=nelz*(nelsy+1)*(nelsx+1)+nely*(nelsx+1)+nelx;
            etpl(nel,7)=etpl(nel,6)+1;
            etpl(nel,8)=etpl(nel,5)+1;
        end
    end
end

end