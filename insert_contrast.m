
function R = insert_contrast(roiTs,Contrast,name,RoiNum)

R = roiTs;
X = getfield(Contrast,name);
R{RoiNum}.r{1} = X.r;
R{RoiNum}.p{1} = X.p;
R{RoiNum}.grpname = name;
