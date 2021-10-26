FileList = dir(fullfile(pwd, '*.mat'));  % List of all MAT files
error_est = zeros(100,7,4);
for iFile = 1:numel(FileList)               % Loop over found files
  Data = load(fullfile(pwd, FileList(iFile).name));
  index1 = sscanf(FileList(iFile).name,'results_experiment_FTP_%d');
  Fields = fieldnames(Data);
  for iField = 1:numel(Fields)              % Loop over fields of current file
    aField = Fields{iField};
    error_est(ceil(index1/7), index1 - (ceil(index1/7) - 1)*7, :) = (Data.(aField))';  
  end
end
save(fullfile(pwd, 'hamming_FTP.mat'), 'error_est');