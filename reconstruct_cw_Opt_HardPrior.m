function [fwd_mesh,pj_error] = reconstruct_cw_Opt_HardPrior...
                                                (fwd_fn,data_fn,output_fn,...
                                                 filter_n)
tic;

% load fine mesh for fwd solve
fwd_mesh = load_mesh(fwd_fn);
frequency = 0;
% read data
anom1 = load_data(data_fn);anom1 = anom1.paa;
anom = log(anom1(:,1));
mesh = fwd_mesh;
%%%%%%%%%%%%%%%%%%%%%%   Building qvec term     %%%%%%%%%%%%%%%%%%%%%%%%%%
% Now calculate source vector
% NOTE last term in mex file 'qvec' is the source FWHM
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);

% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);

qvec = spalloc(nnodes,nsource,nsource*100);
if mesh.dimension == 2
  for i = 1 : nsource
      s_ind = mesh.source.num == source(i);
    if mesh.source.fwhm(s_ind) == 0
        qvec(:,i) = gen_source_point(mesh,mesh.source.coord(s_ind,1:2));
    else
      qvec(:,i) = gen_source(mesh.nodes(:,1:2),...
			   sort(mesh.elements')',...
			   mesh.dimension,...
			   mesh.source.coord(s_ind,1:2),...
			   mesh.source.fwhm(s_ind));
    end
  end

end

clear junk i nnodes nsource;

% Check for distributed source
if mesh.source.distributed == 1
    qvec = sum(qvec,2);
end

% Catch zero frequency (CW) here
if frequency == 0
%   MASS = real(MASS);
  qvec = real(qvec);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate projection error
pj_error = [];

% Initiate log file
fid_log = fopen([output_fn '.log'],'w');
fprintf(fid_log,'Forward Mesh   = %s\n',fwd_fn);

fprintf(fid_log,'Frequency      = %f MHz\n',frequency);
fprintf(fid_log,'Data File      = %s\n',data_fn);
% fprintf(fid_log,'Initial Reg    = %d\n',lambda);
fprintf(fid_log,'Filter         = %d\n',filter_n);
fprintf(fid_log,'Output Files   = %s_mua.sol\n',output_fn);
fprintf(fid_log,'               = %s_mus.sol\n',output_fn);

%if exist('region','var')==0
% disp('Jacob region');
% Intializations
flag_mesh = 0;
regions = unique(fwd_mesh.region);
R = zeros(length(fwd_mesh.mua), length(regions));
MusREC = zeros(length(regions),1);
MuaREC = MusREC;

% store the region indices
for i = 1:length(regions)
    R1 = find(fwd_mesh.region == regions(i));
    RL(i) = length(R1);
    R(1:length(R1), i) = R1;
    MusREC(i) = mean(fwd_mesh.mus(R1));
    MuaREC(i) = mean(fwd_mesh.mua(R1));
    clear R1;
end

%for it = 1 : iteration

mu =  [(MuaREC)] ;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[foo,pj_error,exitflag, output] =...
    NelderMeadSimplexMethod(@(mu) ObjectiveFunction_Opt_HardPrior...
    (mu,fwd_mesh,anom,R,RL,qvec),mu, ...
    optimset('TolX', 1e-12));
 output
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(regions)
    fwd_mesh.mua(R(1:RL(i),i)) = (foo(i));      
end
fwd_mesh.kappa =  (1./(3.*(fwd_mesh.mus+fwd_mesh.mua)));
clear foo Hess Hess_norm tmp data_diff G

% % % % We dont like -ve mua or mus! so if this happens, terminate
if (any(fwd_mesh.mua<0) | any(fwd_mesh.mus<0))
    disp('---------------------------------');
%     disp('-ve mua or mus calculated...not saving solution');
    fprintf(fid_log,'---------------------------------\n');
    fprintf(fid_log,'STOPPING CRITERIA REACHED\n');
% %     break
end

it = 1;
%if it == 1
    fid = fopen([output_fn '_mua.sol'],'w');
% else
%     fid = fopen([output_fn '_mua.sol'],'a');
% end
fprintf(fid,'solution %g ',it);
fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
fprintf(fid,'-components=1 ');
fprintf(fid,'-type=nodal\n');
fprintf(fid,'%f ',fwd_mesh.mua);
fprintf(fid,'\n');
fclose(fid);

%if it == 1
    fid = fopen([output_fn '_mus.sol'],'w');
%else
%    fid = fopen([output_fn '_mus.sol'],'a');
%end
fprintf(fid,'solution %g ',it);
fprintf(fid,'-size=%g ',length(fwd_mesh.nodes));
fprintf(fid,'-components=1 ');
fprintf(fid,'-type=nodal\n');
fprintf(fid,'%f ',fwd_mesh.mus);
fprintf(fid,'\n');
fclose(fid);


% close log file!
time = toc
fprintf(fid_log,'Computation Time = %f\n',time);
fclose(fid_log);
end
