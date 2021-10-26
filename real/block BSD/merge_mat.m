FileList = dir(fullfile(pwd, '*.mat'));  % List of all MAT files
allData  = struct();
for iFile = 1:numel(FileList)               % Loop over found files
  Data   = load(fullfile(pwd, FileList(iFile).name));
  Fields = fieldnames(Data);
  for iField = 1:numel(Fields)              % Loop over fields of current file
    aField = Fields{iField};
    if isfield(allData, aField)             % Attach new data:
       allData.(aField) = [allData.(aField) ; Data.(aField)];
    else
       allData.(aField) = Data.(aField);
    end
  end
end
save(fullfile(pwd, 'block_BSD.mat'), '-struct', 'allData');