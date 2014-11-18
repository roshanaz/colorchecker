clear,
load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\MacbethsVIGerler';
%MacbethsVI = MacbethsVI-(129/4095);
% svd
vecs=zeros(72,568);

for i=1:568
    vecs(:,i)=reshape(MacbethsVI(:,:,i),72,1);
   % vecs(:,i)=vecs(:,i)-.0315; %min(vecs(:,i));
end

[u s v]=svd(vecs*vecs');
pca1=zeros(24,3,3);

for i=1:3
    pca1(:,:,i)=reshape(u(:,i),24,3);
end


%[eRed, eGreen, eBlue]=MacbethEigen(MacbethsVI);
%white_checker_3pca =( [(eRed(19,1:3)); (eGreen(19,1:3)); (eBlue(19,1:3))]);
white_checker_3pca = squeeze(pca1(19,:,:));


% Gw
load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\colorchecker_shi_greyworld';
%load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\colorchecker_shi_exemplarCC';
%load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\colorchecker_shi_shadesofgrey';
%load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\colorchecker_shi_whitepatch';
%load 'C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\colorchecker_shi_whitepatch';
%Gw = mean(estimated_illuminants,1);
Gw = estimated_illuminants(1,:);
%Gw = min(reshape(mean(estimated_illuminants,1),10,3));
%Gw = min(reshape(estimated_illuminants(1,:,:),10,3));
Gw = reshape(Gw,1,3)./norm(Gw);
abc = white_checker_3pca\(Gw).';
 SynthChecker_GW(:,1) =  eRed(:,1:3)* abc;
 SynthChecker_GW(:,2) =  eGreen(:,1:3)* abc;
 SynthChecker_GW(:,3) =  eBlue(:,1:3)* abc;


% white [1 1 1]
w = [1 1 1];
abc_w = white_checker_3pca \(w).';
% SynthChecker_W(:,1) =  eRed(:,1:3)* abc_w;
% SynthChecker_W(:,2) =  eGreen(:,1:3)* abc_w;
% SynthChecker_W(:,3) =  eBlue(:,1:3)* abc_w;
SynthChecker_W = zeros(24,3);
for i=1:3
   SynthChecker_W = SynthChecker_W + squeeze(reshape(pca1(:,:,i),24,3)).*abc_w(i);
end


% minimize error between Gw and white
T =  inv(SynthChecker_GW'*SynthChecker_GW)*SynthChecker_GW'*SynthChecker_W;


% True colorchecker from spectral data
colorchk = xlsread('C:\Users\mcr13hfu\Dropbox\UEA\colour_constancy\Code\AIC\MacbethColorChecker','Sheet1','B4:Y84');
[XYZ_values,lightXYZ] = XYZ_Compute(colorchk); 
%XYZ_values = XYZ_values*inv(diag(sum(XYZ_values)));
% XYZ_values = bsxfun(@rdivide,XYZ_values,sum(XYZ_values,2));
% XYZ_values = XYZ_values./max(XYZ_values(:));
% XYZ to RGB
%colorchecker_C1 = [abs(eRed(:,1)), abs(eGreen(:,1)), abs(eBlue(:,1))];
RGBtoXYz = inv(SynthChecker_W'*SynthChecker_W) * SynthChecker_W'*  XYZ_values;
%im=permute(reshape(XYZtoRGBconver(XYZ_values,inv(RGBtoXYz)),6,4,3),[2 1 3]);
%figure,imshow(imresize(im./max(im(:)),100,'nearest'));


%
SynthChecker_GW_XYZ = RGBtoXYZconver((MacbethsVI(:,:,1)* T),RGBtoXYz);
% the following calc effects SynthChecker_GW_XYZ goes between 0 and 1
T2 = reshape(XYZ_values,size(XYZ_values,1)*size(XYZ_values,2),1)'/reshape(SynthChecker_GW_XYZ,size(XYZ_values,2)*size(XYZ_values,1),1)';
SynthChecker_GW_XYZ = SynthChecker_GW_XYZ.*T2;
SynthChecker_W_XYZ =  RGBtoXYZconver(SynthChecker_W,RGBtoXYz);
im1=imresize(permute(reshape(XYZtoRGBconver(reshape(SynthChecker_GW_XYZ, 24, 3),inv(RGBtoXYz)),6,4,3),[2 1 3]),100,'nearest');
im2=imresize(permute(reshape(XYZtoRGBconver(reshape(SynthChecker_W_XYZ,24,3),inv(RGBtoXYz)),6,4,3),[2 1 3]),100,'nearest');
figure,imshow([im1./max(im1(:));im2./max(im2(:))]);

% RGBtoXYZ is not right
% XYZ to lab
Gw_XYZ = RGBtoXYZconver(Gw, RGBtoXYz)*T2;
%w_XYZ = RGBtoXYZconver(w, RGBtoXYz);
% remove sum to get the least error
lab_Synth_GW = xyz2lab( reshape(SynthChecker_GW_XYZ, 24, 3), ([1 1 1]).*1);
lab_Synth_W = xyz2lab( reshape(SynthChecker_W_XYZ, 24, 3), ([1 1 1]).*1);
% calculate error (deltaE) for each color patch 
% average deltaE
deltaE_GW = getDelab (lab_Synth_GW, lab_Synth_W);
de00 = deltaE2000(lab_Synth_GW, lab_Synth_W);
mean(de00)




% im1=imresize(permute(reshape((SynthChecker_GW * T),6,4,3),[2 1 3]),100,'nearest');
% im2=imresize(permute(reshape((SynthChecker_W),6,4,3),[2 1 3]),100,'nearest');
% figure,imshow([im1./max(im1(:));im2./max(im2(:))].^0.5);