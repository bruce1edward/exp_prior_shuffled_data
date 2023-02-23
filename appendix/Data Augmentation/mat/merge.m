FileList = dir(fullfile('C:\Users\ZhenbangWang\OneDrive - George Mason University - O365 Production\Paper(EM)\AISTAT 2021\DA example\mat', '*.mat'));  % List of all MAT files
iters = 100;iter = 1000;
track_bias = zeros(iters, 2, 2);
track_SE = zeros(iters, 2, 2);
track_CP = zeros(iters, 2, 2);
track_naive = zeros(iters, 2);
track = zeros(iters,iter,2);
for iFile = 1:numel(FileList)               % Loop over found files
  Data   = load(fullfile('C:\Users\ZhenbangWang\OneDrive - George Mason University - O365 Production\Paper(EM)\AISTAT 2021\DA example\mat', FileList(iFile).name));
  Fields = fieldnames(Data);
    aField1 = Fields{1};
    aField2 = Fields{2};
    aField3 = Fields{3};
    track_para = Data.(aField2);
    track_para1 = Data.(aField3);
    track(iFile,:,:) = track_para;
    track(iFile,:,:) = track_para1;
    track_bias(iFile,:,1) = (mean(track_para) - beta_oracle')';
    track_bias(iFile,:,2) = (mean(track_para1) - beta_oracle')';
    track_SE(iFile,:,1) = std(track_para)/sqrt(iter);
    track_SE(iFile,:,2) = std(track_para1)/sqrt(iter);
    lower = mean(track_para) - quantile(track_para,0.025).*std(track_para)/sqrt(iter);
    Upper = mean(track_para) + quantile(track_para,0.975).*std(track_para)/sqrt(iter);
    lower1 = mean(track_para1) - quantile(track_para1,0.025).*std(track_para1)/sqrt(iter);
    Upper1 = mean(track_para1) + quantile(track_para1,0.975).*std(track_para1)/sqrt(iter);
    track_CP(iFile,:,1) = beta_oracle' >= lower & beta_oracle' <= Upper;
    track_CP(iFile,:,2) = beta_oracle' >= lower1 & beta_oracle' <= Upper1;
    track_naive(iFile,:) =  Data.(aField1)';
end

mean(track_bias)
mean(track_SE)
mean(track_CP)
mean(track_naive)' - beta_oracle
std(track_naive)/sqrt(iters)