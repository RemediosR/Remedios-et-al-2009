clear
% script to compute the correlation for each of the Natsounds category
% computes correlations across all repeats fo the sounds, by looping all
% files 
global DPATH

SesName = 'J02QL1';
RoiNum = 1;
goto(SesName);

Ses = goto(SesName);
Ind = Ses.grp.natsounds.exps;
for j=1:length(Ind)
  FileNum(j) = Ses.expp(Ind(j)).scanreco(1);
end
SesName2 = [SesName(1:3) '.' SesName(4:6)];

% use reduced set of files with better activation:
FileNum = [23:27,30:40]; % 21 41
Nfile = length(FileNum);
j=1;
FileName = sprintf('%s/%s/%s_%03d.mat',DPATH,SesName2,SesName,FileNum(j));
load(FileName);
%% create models
x = roiTs{1}.mdl{1}';
x = x(5:end);% -< BASELINE ++++++++++++++=
Model.stim = repmat(x,[1,Nfile]);
% Avoc, cage, mvoc, nsnd, rttl
x =[];
for c=1:5
  x = cat(2,x,zeros(1,4));
  x = cat(2,x,ones(1,4)*c);
end
x = cat(2,x,zeros(1,4));
% remove initial baseline 
x = x(5:end); % -< BASELINE ++++++++++++++=

Model_cat = repmat(x,[1,Nfile]);

% model for individual sounds
for s=1:5
  x = zeros(size(Model_cat));
  % get each stimulus plus pre- and post-baseline
  j = find(Model_cat==s);
  j = [j(1) j(find(diff(j)>1)+1)];
  y =[];
  for l=1:length(j)
    y = [y j(l)+[-4,-3,-2,-1]];
    y = [y j(l)+[0:7]];
  end
  y = y(find(y>0));
  Model_each.m{s} = Model_cat(y);
  Model_each.j{s} = y;
end


% Movc > {avoc,nsnd,cage}
Model.mvoc = zeros(size(Model_cat));
Model.mvoc(find(Model_cat==3)) = 1;
Model.mvoc(find(Model_cat==1)) = -1/2;
Model.mvoc(find(Model_cat==4)) = -1/2;
% Model.mvoc(find(Model_cat==2)) = -1/3;

% rttl > {avoc,nsnd,cage}
Model.rttl = zeros(size(Model_cat));
Model.rttl(find(Model_cat==5)) = 1;
Model.rttl(find(Model_cat==1)) = -1/2;
Model.rttl(find(Model_cat==4)) = -1/2;
% Model.rttl(find(Model_cat==2)) = -1/3;

% cage > {avoc,nsnd,rttl}
Model.cage = zeros(size(Model_cat));
Model.cage(find(Model_cat==2)) = 1;
Model.cage(find(Model_cat==1)) = -1/2;
Model.cage(find(Model_cat==4)) = -1/2;
% Model.cage(find(Model_cat==5)) = -1/3;

% concatenate all data
X = [];
for j=1:Nfile
  j
  FileName = sprintf('%s/%s/%s_%03d.mat',DPATH,SesName2,SesName,FileNum(j));
  load(FileName);
  X = cat(1,X,roiTs{RoiNum}.dat(5:end,:));
end


%% compute correlation for each individual stimulus
Contrast.test = 0;
% Avoc, cage, mvoc, nsnd, rttl
for s=[1:5]
  model = Model_each.m{s};
  y = Model_each.j{s};
  rp = zeros(2,size(X,2));
  for N=1:size(X,2)
    [tmpr,tmpp] = corrcoef(model,X(y,N));
    rp(:,N) = [tmpr(1,2),tmpp(1,2)];
  end
  U.p = rp(2,:)';
  U.r = rp(1,:)';
  Contrast = setfield(Contrast,sprintf('stim%d',s),U);
end
%% compute correlation for each difference model
Mod = fieldnames(Model);
for k=1:length(Mod)
  model = getfield(Model,Mod{k})';
  if k>1
    % remove baseline for correlation
    J = find(model);
    model = model(J);
  else
    J = [1:length(model)];
  end
  rp = zeros(2,size(X,2));
  for N=1:size(X,2)
    [tmpr,tmpp] = corrcoef(model,X(J,N));
    rp(:,N) = [tmpr(1,2),tmpp(1,2)];
  end

  RES.p{k} = (rp(2,:)');
  RES.r{k} = (rp(1,:)');
end
k=1;
U.p = RES.p{k};
U.r = RES.r{k};
Contrast = setfield(Contrast,Mod{k},U);
% add constraints for positive total response
for k=2:length(Mod)
  U.p = RES.p{k}.*(RES.r{1}>0.01);
  U.r = RES.r{k}.*(RES.r{1}>0.01);
  Contrast = setfield(Contrast,Mod{k},U);
end

% put back into roiTs like structure for display
roiTs = ck_insert_anatomy(roiTs,'J02QL1');
if isfield(roiTs{1},'dat')
  roiTs{1} = rmfield(roiTs{1},'dat');
end
roiTs{RoiNum}.dat = X;

% to display the results, use the following functions:

% Avoc, cage, mvoc, nsnd, rttl
R = insert_contrast(roiTs,Contrast,'stim',RoiNum);
% cla_dsproits(R,'brain',0.06,[2,6],'r')




%% get all acivation maps in one matrix
THR = 0.06;
CLUS = [2,8];
clear Allmaps
BigAna = mgetcollage(R{1}.anabig,[],[],[]);

% compute activation maps for overlay
for s=1:5
  R = insert_contrast(roiTs,Contrast,sprintf('stim%d',s),RoiNum);
  % construct contrasts
  R = mroitsget(R,[],'brain');
  R = mroitssel(R,THR,[],'r','mean',CLUS);
  % map
  xyz = R{1}.coords;
  cormap = NaN*ones(size(R{1}.ana));
  for N=1:length(R{1}.r{1}),
    cormap(xyz(N,1),xyz(N,2),xyz(N,3)) = R{1}.r{1}(N);
  end;
  fscan = mgetcollage(cormap,[],[],[]);
  Allmaps2(:,:,s) = fscan;
  fscan2 = fscan;
  fscan2(find(isnan(fscan2(:)))) = 0;
  fscan2 = imresize(fscan2,size(BigAna),'nearest');
  fscan2(find(fscan2(:)==0)) = NaN;
  Allmaps(:,:,s) = fscan2;
end
save('Dataset_natsounds.mat','roiTs','Contrast','RoiNum','Model','BigAna','Allmaps','Model_cat','Model_each');


ck_dspfused(BigAna,Allmaps(:,:,5));

% mvoc & rttl
map_mvoc = ((Allmaps(:,:,3)+Allmaps(:,:,5)))/2;
map_mvoc(find(map_mvoc(:)<THR)) = NaN;

% mvoc & rttl > max(nsnd,avoc)
map_nsnd = max(Allmaps(:,:,[1,4]),[],3);
map_nsnd(find(map_nsnd(:)<THR)) = NaN;

map_diff =  map_mvoc;
map_diff(find(map_mvoc<map_nsnd)) = NaN;


