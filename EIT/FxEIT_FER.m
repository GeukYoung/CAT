function [imdl,Proj_Mat] = FxEIT_FER(bd_n,thick,hmax,E_posi,bd_alpha,Nskip,isKyunghee,mask_factor,E_posi2)
if nargin < 3
    bd_n=10;
    thick=0.05; % for boundary elemination factor(%)
    hmax=0.05;  % mesh size (%)
    option = [];
    ch = 16;
end

if nargin < 4
    E_posi(:,1) = cos(linspace(0,2*pi,17));
    E_posi(:,2) = -sin(linspace(0,2*pi,17));
    E_posi(end,:) = [];
else
    if exist('E_posi2') == 1
        m_E_posi = mean(E_posi(:,1));
        E_posi(:,1) = E_posi(:,1) - m_E_posi;
        E_posi2(:,1) = E_posi2(:,1) - m_E_posi;
        
        m_E_posi = mean(E_posi(:,2));
        E_posi(:,2) = E_posi(:,2) - m_E_posi;
        E_posi2(:,2) = E_posi2(:,2) - m_E_posi;
        
        max_E_psi = max(max(abs(E_posi)));
        E_posi = E_posi./max_E_psi;
        E_posi2 = E_posi2./max_E_psi;
        option=E_posi2;
        
        ch = length(E_posi2);
    else
        E_posi(:,1) = E_posi(:,1) - mean(E_posi(:,1));
        E_posi(:,2) = E_posi(:,2) - mean(E_posi(:,2));
        E_posi = E_posi./max(max(abs(E_posi)));
        option=[];
        ch = length(E_posi);
    end
end

    
if exist('Nskip') ~= 0
    if any(Nskip) == 0
        Nskip=0;
    end
else
    Nskip=0;
end
    
if exist('isKyunghee') ~= 0
    if length(isKyunghee) == 0
        isKyunghee = 1;
    end
else
    isKyunghee=1;
end

if exist('bd_alpha') ~= 0
    if any(bd_alpha) == 0
        bd_alpha=1e-9;
    end
else
    bd_alpha=1e-9;
end
% if nargin < 5
% Nskip=0;
% isKyunghee=1;
%% Mesh Gen
% option=[];
[Node,Element,erod,BndIndex,fnum,bnd]=func_MeshGeneration_Bd(E_posi,bd_n,thick,hmax,option);
%
indx_near_bd=find(fnum==1);
indx_interior=find(fnum==2);
% indx_near_bd=func_get_near_bd(Geom.Element,Geom.Node,Geom.BndIndex,8);

% Display Mesh result
pix=zeros(size(Element,1),1);
pix(indx_near_bd)=1;
pix(indx_interior)=2;
figure(2000);clf;
patch('Faces',Element,'Vertices',Node,'FaceVertexCData',pix,'FaceColor','flat','EdgeColor','None');
axis normal image off

%% Construct Reconstruction matrix
mm=size(Element,1);

% mask = [1	2	16	17	18	19	34	35	36	51	52	53	68	69	70	85	86	87	102	103	104	119	120	121	136	137	138	153	154	155	170	171	172	187	188	189	204	205	206	221	222	223	238	239	240	241	255	256];
if exist('mask_factor') ~= 0
    if any(mask_factor) ~= 0
        mask = FxEIT_mask(ch,mask_factor);
    else
        mask = FxEIT_mask(ch);
    end
else
    mask = FxEIT_mask(ch);
end
% Sensitivity=func_sensitivity_skip(Geom.Element,Geom.Node,Geom.erod,ones(mm,1),Nskip,isKyunghee);
[Sensitivity,Area]=func_sensitivity_skip_fast(Element,Node,erod,ones(mm,1),Nskip,isKyunghee);
Sensitivity(mask,:)=[];

%%
alpha=1*1e-5;
inv_Sense=(Sensitivity'*Sensitivity+alpha*eye(size(Sensitivity,2)))\(Sensitivity');

[R_WeightAvg,S_normal]=func_Recon_WeightedAvg(Sensitivity);
inv_Sense_corr_infty=S_normal';
% Kotre's formula
R_Kotre=diag(1./sum(Sensitivity))*Sensitivity';
R_Kotre_abs=diag(1./sum(abs(Sensitivity)))*Sensitivity';

%% Boundary Artifact Reduction
Proj_Mat=func_bd_artifact_elim_newnew(Sensitivity,indx_interior,bd_alpha);

%% Data Export
imdl.fwd_model.name = 'FER_model';
imdl.fwd_model.nodes = Node;
imdl.fwd_model.elems = Element;

for i = 1:length(erod);
    imdl.fwd_model.electrode(1,i).nodes = erod(i);
    imdl.fwd_model.electrode(1,i).z_contact = 1.0e-3;
end
imdl.fwd_model.type = 'fwd_model';
imdl.fwd_model.normalize_measurements = 0;
imdl.fwd_model.solve = 'eidors_default';
imdl.fwd_model.system_mat = 'eidors_default';
imdl.fwd_model.jacobian = 'eidors_default';
imdl.fwd_model.gnd_node = 1;
imdl.fwd_model.boundary = find_boundary2(imdl.fwd_model);
imdl.fwd_model.indx_near_bd = indx_near_bd;
imdl.fwd_model.indx_interior = indx_interior;
imdl.fwd_model.BndIndex = BndIndex;

imdl.solve = @solve_use_matrix;
imdl.reconst_type = 'difference'
imdl.jacobian_bkgnd.value = 1;
imdl.type = 'inv_model';
imdl.solve_use_matrix.RM = R_WeightAvg;
imdl.Proj_Mat = Proj_Mat;
show_fem(imdl.fwd_model,[0 1]);

function [srf, idx] = find_boundary2(simp);
% [srf, idx] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  array of elements on each boundary simplex
%        boundary simplices are of 1 lower dimention than simp
%idx  =  index of simplex to which each boundary belongs
%simp = The simplices matrix

% $Id: find_boundary.m 4926 2015-05-07 23:10:02Z aadler $

if isstr(simp) && strcmp(simp,'UNIT_TEST'); do_unit_test; return; end
if isstruct(simp) && strcmp(simp.type,'fwd_model'); simp= simp.elems; end

wew = size(simp,2) - 1;

if wew==3 || wew==2
   [srf,idx]= find_2or3d_boundary(simp,wew);
elseif wew == 1
   [srf,idx]= find_1d_boundary(simp);
else
   eidors_msg('find_boundary: WARNING: not 1D, 2D or 3D simplices',1);
   srf=[]; return;
end

% sort surfaces. If there is more than one, its not on the boundary
function [srf,idx]= find_2or3d_boundary(simp,dim);
   if size(simp,1) < 4e9 % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   localface = nchoosek(1:dim+1,dim);
   srf_local= simp(:,localface');
   srf_local= reshape( srf_local', dim, []); % D x 3E
   srf_local= sort(srf_local)'; % Sort each row
   [sort_srl,sort_idx] = sortrows( srf_local );

   % Fine the ones that are the same
   first_ones =  sort_srl(1:end-1,:);
   next_ones  =  sort_srl(2:end,:);
   same_srl = find( all( first_ones == next_ones, 2) );

   % Assume they're all different. then find the same ones
   diff_srl = logical(ones(size(srf_local,1),1));
   diff_srl(same_srl) = 0;
   diff_srl(same_srl+1) = 0;

   srf= sort_srl( diff_srl,: );
   idx= sort_idx( diff_srl);
   idx= ceil(idx/(dim+1));

function [srf,idx]= find_1d_boundary(simp);
   if size(simp,1) < 4e9 % max of uint32
      % convert to integer to make sort faster
      simp = uint32( simp );
   end
   % we expect two nodes as a result
   idx = find(isunique(simp(:)) == 1);
   srf = simp(idx);
   idx = rem(idx-1,size(simp,1))+1;

function x = isunique(a);
   u=unique(a);
   n=histc(a,u);
   x=ismember(a,u(n==1));

function do_unit_test

%2D Test:  
mdl = mk_common_model('c2c',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

unit_test_cmp('2D test', bdy, bdyc);

%3D Test:  
mdl = mk_common_model('n3r2',[16,2]);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

unit_test_cmp( '3D test n3r2', bdy,bdyc);

%3D Test:  
mdl = mk_common_model('a3cr',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

unit_test_cmp('3D test a3c2', bdy, bdyc);

%3D Test:  
mdl = mk_common_model('b3cr',16);
bdy = find_boundary(mdl.fwd_model.elems);
bdy = sort_boundary(bdy);
bdyc= sort_boundary(mdl.fwd_model.boundary);

unit_test_cmp('3D test b3c2', bdy, bdyc);

simp = [  10 190; ...
         182 183; ...
         183 184; ...
         184 185; ...
          11 182; ...
         185 186; ...
         187 186; ...
         187 188; ...
         188 189; ...
         189 190];
[bdy, idx] = find_boundary(simp);
unit_test_cmp('1D bdy', bdy,[10;11]);
unit_test_cmp('1D bdy', idx,[1;5]);

function bdy= sort_boundary(bdy)
   bdy = sort(bdy,2);
   bdy = sortrows(bdy);

