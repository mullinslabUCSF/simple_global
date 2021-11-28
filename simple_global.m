function [s,m]=simple_global(data_dir,refsp,extens)

% original version: 4/7/2010
% revised: 11/16/2021
% R. Dyche Mullins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% please cite: Hansen SD, Mullins RD (2010) VASP is a processive actin 
% polymerase that requires monomeric actin for barbed end association. 
% J Cell Biol 191(3):571–584.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs a global, least-squares, fitting of equilibrium 
% analytical ultracentrifuge data. The function assumes that data have been
% collected at three different speeds in a six-channel centerpiece with
% three sample wells. The concentrations of the macromolecules can vary but
% the model assumes that the absorbance profile represents a single 
% sedimenting species. This type of analysis is most useful for 
% macromolecules known to be monomeric or for estimating the size of high
% affinity complexes. In general, this is a good first-pass analysis that 
% provides a sense for whether a macromolecule oligomerizes and how big
% the oligomer might be. Homo-oligomers should be analyzed by more powerful
% fitting programs such as Winnln, Sedphat, or Ultrascan. For heterologous
% associating systems, we use multi_global_kb2(). 

% OUTPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s -- best-fit sigma factor
% m -- vector containing the amplitudes and offsets of the best fit

% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datadir -- string containing the name of data directory 
% refsp -- reference speed for sigma factor (used to calc. speed factors)
% extens -- string with the data file extension (e.g. "RA1" or "txt")

% Subfunctions called by this function (appended below): 
% get_datafile_names()
% get_data() 
% lin_nonlin()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the datafile names from the data directory
[dfiles,fcontent,nds] = get_datafile_names(data_dir,extens);
% nds is the number of datasets to analyze
% dfiles is a cell array holding the names of the files to analyze

% create a cell array to hold the linearized data and parameters
dset=cell(nds,6);
% dset{*,1} is an mx2 matrix containing the x and y values of the dataset
% dset{*,2} is the number of points in the dataset (scalar)
% dset{*,3} is the position of the meniscus (scalar)
% dset{*,4} is the position of the base of the cell (scalar)
% dset{*,5} is the speed factor associated with the data
% dset{*,6} are expected values from the linear fitting

% create a cell array to hold the raw data and final fit
rawdat=cell(nds,1);
% rawdat{1} is an mx4 matrix containing raw x, y, fitted y, and residuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   First we will load all of the data with a given file extension
%   found in a given directory and display the raw data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the datasets
for i=1:nds
    % open a dataset and read out the data and headers
    [x,y,err,speed,temp,sh1,sh2] = getdata(dfiles{i},data_dir);
    rawdat{i}(:,1)=x;
    rawdat{i}(:,2)=y;
    npts = length(x);
    dset{i,2}=npts;
    dset{i,3}=rawdat{i}(1,1);     % meniscus set to first point
    dset{i,4}=rawdat{i}(npts,1);  % base set to last point
    dset{i,5}=speed*speed/(refsp*refsp);    % speed factor for data set
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Next, we linearize the data and use linear least squares methods
%   to solve for the amplitudes and the r.e. molecular weight (sigma)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linearize the data and store in dset array
% checking for negative data
negflag = 0;
for i=1:nds
    for j=1:dset{i,2}
        % calculate r^2/2 - rm^2/2
        dset{i,1}(j,1)=rawdat{i}(j,1)*rawdat{i}(j,1)/2 - dset{i,3}*dset{i,3}/2;
        % check for negative numbers and remove them
        if rawdat{i}(j,2) < 0
            dset{i,1}(j,2) = log(-rawdat{i}(j,2))/dset{i,5};
            negflag=1;
        else
            % take natural log of y data and normalize by the speed factor
            dset{i,1}(j,2) = log(rawdat{i}(j,2))/dset{i,5};
        end
    end
end
% alert the user that some negative numbers have been reflected around 0
if negflag
    negmsg = 'NOTE: some negative numbers removed from data'
end

% build the design matrix
% first count how many data points to catenate
tpts=0;
for i=1:nds
    tpts=tpts+dset{i,2};
end

% then catenate the data and build the catenated design matrix
Y=zeros(tpts,1);        % this holds the catenated measured Y values
A=zeros(tpts,nds+1);    % A helps build the design matrix
k=0;
for i=1:nds
    for j=1:dset{i,2}
        k=k+1;          % keep track of rows. k goes from 1 to tpts
        A(k,i)=1;       % load up the constant coefficients
        A(k,nds+1) = dset{i,1}(j,1);    % now load the linear part
        Y(k)=dset{i,1}(j,2);     % finally, catenate the y values
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here are the best-fit parameters for the log-transformed data %%%%%%%%%%
% solve the linear least squares problem:
M=A'*A;
B=M\(A'*Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the initial values of return values
m = B(1:nds);
s = B(nds+1);

figure(1)
% plot the linearized data and fits
subplot(2,2,1);
hold on
for i=1:nds
    dset{i,6}=B(i)+B(nds+1).*dset{i,1}(:,1);    % this is the set of expected values from the fit
    
    plot(dset{i,1}(:,1),dset{i,1}(:,2));
    plot(dset{i,1}(:,1),dset{i,6});
end 
sttl = 'transformed data and initial linear fits';
    title(sttl);
    xlabel('r^2/2  -  r_m^2/2  (cm^2)');
    ylabel('ln(absorbance)');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Finally, we will test values around the initial sigma to see whether
%   there is a nearby value with a lower MSE in the non-transformed space.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick a set of sigma values that are +/- 20% of the initial value
sigmas = [0:20]';
sigmas = sigmas/50;
sigmas = sigmas-0.2;
sigmas = sigmas*s + s;

% MSEbfx is a vector of mean square errors with fixed baselines
% MSEbfl is a vector of mean square errors with floating baselines
MSEbfx = zeros(21,1);
MSEbfl = zeros(21,1);

% Bmfl is a cell array that contains vectors of all the best fit
% parameters tested for optimal fitting
Bmfl=cell(nds,1);

% run through the sigma values and calculate a linear least squares best 
% fit for the raw (non-linearized) data
for j=1:21
    % send out the number of data sets, the test sigma value and the data
    [MSEbfl(j),MSEbfx(j),Bmfl{j}]=lin_nonlin(nds,sigmas(j),dset,rawdat);
end

% plot the mean square errors of the tested sigma values
subplot(2,2,2);
hold on
%plot(sigmas,MSEbfx,'o');
plot(sigmas,MSEbfl,'k--o');
plot(sigmas(11),MSEbfl(11),'ob','MarkerFaceColor','b');
% select the minimum of the MSEfl vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vdum,mind]=min(MSEbfl);
plot(sigmas(mind),MSEbfl(mind),'or','MarkerFaceColor','r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sttl = 'MSE of adjacent sigma values';
    title(sttl);
    xlabel('sigma');
    ylabel('mean-square error');
hold off

% plot the raw data and the fits
subplot(2,2,3);
hold on
for i=1:nds
    plot(rawdat{i}(:,1),rawdat{i}(:,2),'o');
    % calculate the new best fit from the minimized MSE
    rawdat{i}(:,3)=Bmfl{mind}(i)+Bmfl{mind}(nds+i)*exp(sigmas(mind)*dset{i,5}*dset{i,1}(:,1));
    % plot the raw fit on the same graph
    plot(rawdat{i}(:,1),rawdat{i}(:,3));
end
sttl = 'untransformed data and final fits';
    title(sttl);
    xlabel('radius (cm)');
    ylabel('absorbance (AU)');
hold off

% calculate and plot the residuals
subplot(2,2,4);
hold on
for i=1:nds
    % calculate the residuals of the fit
    rawdat{i}(:,4)=rawdat{i}(:,2)-rawdat{i}(:,3);
    % plot the residuals
    plot(rawdat{i}(:,1),rawdat{i}(:,4),'o');
end
sttl = 'residuals of final fit (untransformed)';
    title(sttl);
    xlabel('radius (cm)');
    ylabel('absorbance (AU)');
hold off

% reset the return values based on secondary grid search
s = sigmas(mind);
m = Bmfl{mind};

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function lin_nonlin calculates a linear least-squares fit of the 
%   non-transformed data and then returns the mean square error of 
%   the fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msql,msqx,Bfl]=lin_nonlin(nds,sval,dset,rawdat)
% msql is the Mean Square error of the computed fit with floating baeline
% msqx is the Mean Square error of the computed fit with fixed baeline
% Bfl is a vector of best fit parameters for floating baseline
% nds is the number of data sets
% sval is the sigma value to test
% dset contains metadata about the datasets
% rawdat contains the raw x and y data

% build the design matrix
% first count the total number of points to catenate
tpts=0;
for i=1:nds
    tpts=tpts+dset{i,2};
end

% then catenate the data and build the catenated design matrix
Y=zeros(tpts,1);        % this holds the catenated measured Y values
Afl=zeros(tpts,2*nds);  % Afx helps build the design matrix for floating BL
Afx=zeros(tpts,nds);  % Afx helps build the design matrix for fixed BL
k=0;
for i=1:nds
    for j=1:dset{i,2}
        k=k+1;          % keep track of rows. k goes from 1 to tpts
        % for the floating baseline calculation %%%%%%%%%%%%%%%%%%%%%%%%%%
        Afl(k,i)=1;     % load up the constant coefficients
        Afl(k,nds+i) = exp(sval*dset{i,5}*dset{i,1}(j,1));    % load exponential part
        %%%% don't forget the SPEED FACTOR!!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for the fixed baseline calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Afx(k,i) = exp(sval*dset{i,5}*dset{i,1}(j,1));    % exponential only
        %%%% don't forget the SPEED FACTOR!!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y(k)=rawdat{i}(j,2);     % finally, catenate y values
    end
end

% solve the linear least squares problem for floating baseline
Mfl=Afl'*Afl;
Bfl=Mfl\(Afl'*Y);
% solve the linear least squares problem for fixed baseline
Mfx=Afx'*Afx;
Bfx=Mfx\(Afx'*Y);

% now compute the mean square errors
msql = 0;
msqx = 0;
for i=1:nds
    for j=1:dset{i,2}
        k=k+1;          % keep track of rows. k goes from 1 to tpts
        % mean square error for the floating baseline calculation
        msql = msql + (Bfl(i)+Bfl(nds+i)*exp(sval*dset{i,5}*dset{i,1}(j,1))-rawdat{i}(j,2))^2;
        % mean square error for the fixed baseline calculation
        msqx = msqx + (Bfx(i)*exp(sval*dset{i,5}*dset{i,1}(j,1))-rawdat{i}(j,2))^2;
    end
end

return
end

function [namelist,fcontent,numfiles] = get_datafile_names(data_dir,extens)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_datafile_names - goes to a directory; identifies all the .rax files
% and loads their names into a cell string array. 
% data_dir - where to look for the data files
% extens - string containing the extension for the files of interest (e.g.
% 'RA1', 'RA2', etc.
% numfiles - is the number of files whose names are retrieved
% fcontent - is an array containing headers and speeds
% RDM - 10/23/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the return directory to the current one
ret_dir = pwd;
% go to the specified data directory
cd(data_dir)
% do a linux 'ls' command and store the results in 'txtfiles'. This
% makes txtfiles a list of strings corresponding to the file names.
extens = ['.' extens];
searstr = ['*' extens];    % first make a search string
raxfiles = ls(searstr);

% how many files?
ddum = size(raxfiles);
numfiles = ddum(1);

% create an indexable list of namess that can be used to read the data
% files or rearranged to construct new names, linked to he originals. 
% This list of name parts gets passed on to other functions
namelist = cell(numfiles,1);

for i=1:numfiles
    namelist{i}=raxfiles(i,:);
end

% % make sure raxfiles is a single row vector
% raxfiles=raxfiles';
% raxfiles=reshape(raxfiles,1,[]);

% get info about the contents of the data files
% column 1 is the file names; column 2 is the list of speeds
fcontent = cell(numfiles,2);
for k=1:numfiles
    % open the data file
    fileID = fopen(namelist{k},'r');
    dum = fgetl(fileID);
    fcontent{k,1} = fgetl(fileID);
    fclose(fileID);
    
    % find the speed information
    % find all the spaces in the string
    SPind = findstr(char(32),fcontent{k,1}); % ascii code for space is 32
    dumspeed = fcontent{k,1}(SPind(3)+1:SPind(4)-1);
    %fcontent{k,2} = str2num(dumspeed);
    fcontent{k,2} = str2double(dumspeed);
    
end

cd(ret_dir)

return
end

function [x,y,err,speed,temp,sh1,sh2] = getdata(datafile,datadir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open an analytical ultracentrifuge data file and copy out all the data
%
% x - x positional data
% y - absorbance at each position in the cell
% err - calculated standard deviation at each point
% speed - centrifuge speed (determined from header information
% sh1 - first line of the header
% sh2 - second line of the file header
%
% datafile - is the filename to open
% datadir - is the directory that holds all the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the return directory to the current one
ret_dir = pwd;
% go to the specified data directory
cd(datadir)

% get info about the contents of the data files
% column 1 is the file names; column 2 is the list of speeds

% open the data file
fileID = fopen(datafile,'r');
sh1 = fgetl(fileID);    % read the first line of the header
sh2 = fgetl(fileID);    % read the second line of the header
    
% find the speed information
% find all the spaces in the string
SPind = findstr(char(32),sh2); % ascii code for space is 32
dumparam = sh2(SPind(3)+1:SPind(4)-1);
speed = str2num(dumparam);
% the code above fails when the speed is below 10000 RPM so let's fix it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(speed)
    dumparam = sh2(SPind(4)+1:SPind(5)-1);
    speed = str2num(dumparam);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dumparam = sh2(SPind(2)+1:SPind(3)-1);
temp = str2num(dumparam);
% CAUTION: the code above will probably fail if the temperature is <10 C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% count the number of data points
fcontent = fgetl(fileID);
k = 0;
while(fcontent ~= -1)
    k=k+1;
    fcontent = fgetl(fileID);
end

% set up data arrays
numpts = k;
x = zeros(numpts,1);
y = zeros(numpts,1);
err = zeros(numpts,1);

fclose(fileID);
fileID = fopen(datafile,'r');
% skip header
dum = fgetl(fileID);
dum = fgetl(fileID);
for k=1:numpts
    fcontent = fgetl(fileID);
    A = sscanf(fcontent,'  %f  %f   %f');
    x(k) = A(1);
    y(k) = A(2);
    err(k) = A(3);
end

cd(ret_dir)

end
