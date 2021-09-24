% define array of polar angles examined in the experiment
theta = [0:5:100,120:20:180];

% total number of data files collected
N = length(theta);

% create arrays to store the mean and standard deviation
pmean = zeros(size(theta));
pstd = zeros(size(theta));

% loop through all data files
for (i=1:N)
    % filename of ith file
    FileName=['Pcyl_deg',num2str(theta(i)),'.csv'];

    % read in data from the file
    data = readmatrix(FileName);
    
    % calculate mean and std 
    pmean(i) = mean(data);
    pstd(i) = std(data);
end