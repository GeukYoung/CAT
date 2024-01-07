function [result] = FxEIT_Motion(imdl,Data)
for i = 1:length(imdl.fwd_model.electrode)
    indx_near_Erod_single{i} = func_get_near_bd_index(imdl.fwd_model.elems, imdl.fwd_model.electrode(i).nodes);
end
% find bnd elem idx
mm = size(imdl.fwd_model.elems,1);
indx_near_bd=func_get_near_bd_index(imdl.fwd_model.elems,imdl.fwd_model.BndIndex);
indx_interior = setdiff(1:mm,indx_near_bd);

% calc S mat
Nskip = 0;
Sensitivity = func_sensitivity_skip(imdl.fwd_model.elems,imdl.fwd_model.nodes,[imdl.fwd_model.electrode.nodes]',ones(mm,1),Nskip,true);
Sensitivity(FxEIT_mask(16),:)=[];
regulMotion = 1e-3;
Proj_MatMotion=func_bd_artifact_elim_new(Sensitivity,indx_interior,regulMotion);
Proj_Mat_ATF = eye(208)-Proj_MatMotion;
bdATFtmp = Proj_Mat_ATF * Data;
for i=1:size(indx_near_Erod_single,2)
   bdATF(i,:) = mean(imdl.solve_use_matrix.RM(indx_near_Erod_single{i},:)*bdATFtmp,1);
end
result.sd_bdATF = std(bdATF');
result.bdATF = bdATF;
end