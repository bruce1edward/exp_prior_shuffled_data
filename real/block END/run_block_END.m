range = 1:700;
jobnames = cell(length(range), 1);
filenames = cell(length(range), 1);
savedir = 'C:\Users\zwang39\OneDrive - George Mason University - O365 Production\Paper(EM)\AISTAT 2021\real2';

pauseinterval=8;

for j=range
    jobnames{j} = ['experiment_block_END' '(' num2str(j) ')']; 
    filenames{j} = ['results_experiment_END_' num2str(j) '.mat'];
end

batchsize = 10;
nrunning = 0;
filenamesrunning = {};
cur = 0;
while cur < max(range)
    
    cur = cur+1; %-nodesktop
    [status,cmdout]= system(['matlab -minimize -r ' jobnames{cur} ' &']);
    nrunning = nrunning + 1;
    filenamesrunning = [filenamesrunning;filenames{cur}];
    
    if nrunning < batchsize
        
    else
        anyfinished=false;
        
        while ~anyfinished
            
            if exist('savedir', 'var')
                dirData = dir(savedir);
            else
                dirData = dir(pwd);      %# Get the data for the current directory
            end
            dirIndex = [dirData.isdir];  %# Find the index for directories
            fileList = {dirData(~dirIndex).name}';
            
            isfinished = ismember(filenamesrunning, fileList);
            
            if any(isfinished)
                filenamesrunning(isfinished) = [];
                nrunning = nrunning - sum(isfinished);
                anyfinished = true;
            else
                pause(pauseinterval)
            end
            
        end
    end
end