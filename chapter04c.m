%% Analyzing Neural Time Series Data
% Matlab code for Chapter 4 script C
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% 

% The function 'clear' removes all data from the matlab buffer
clear

% you can also clear specific variables.
leave_me_alone   = 10;
remove_me_please = 20;
clear remove_me_please
whos

remove_me_please % gives an error because the variable doesn't exist anymore!

%% basic importing text data

% The most basic way to import data is to copy and paste. This is the best
% option for small amounts of numeric data that you'll need to put into matlab
% only once. (Hint: use square brackets.)

% The simplest way to read in text data is if all data in the next are
% numbers (no text). Open a text editor and make a small matrix (say, 3x4).
% Next, type: 
data = load('chapter04_datafile.txt');

% slightly more advanced:
[file_name,file_path]=uigetfile('*.txt'); % ui = user-interface
data = load([ file_path file_name ]);

% you can also read in data from excel files, but BE CAREFUL because this
% function can act in unexpected ways, e.g., by removing empty columns and
% rows without warning (this can be seen in comparing "numberdata" to "raw_data"). 
% Therefore, it might be best to use the "raw" data output. 
[numberdata,textdata,raw_data] = xlsread('chapter04_excel_data.xls');

%% advanced importing text data

% Here we borrow from C language to flexibly read in mixed data. Let's say
% you have some poorly organized behavioral data files to read in, but at 
% least you know what text strings to look for: 

fid = fopen('chapter04_headache_data.txt','r');
% fid is a pointer to a location on the physical hard disk (similar to how
% we used variables as handles to axes when plotting). The 'r' means read
% (later we'll use 'w' for write).

% In this particular example, we will extract the trial number, subject
% choice, reaction time (RT), and accuracy for each trial. Fields are separated by tabs.

behavioral_data=[]; % initialize... we can't initialize the full matrix, because we don't know how big this will be.

% The following code will remain inside a loop, reading in and processing new
% lines of data, until we reach the end of the file.
datarow=1;

while ~feof(fid) % feof tests whether we're at the end of the file.
    
    dataline = fgetl(fid); % read a line ("file get line")
    
    dataline = regexp(dataline,'\t','split');
    % regexp can be used to cut data according to delimiters. Here we will
    % cut this string of characters into a cell array in which elements of
    % the array are separated by tabs.
    
    % here we use strcmpi to compare strings. The "i" means to ignore case.
    if ~any(strcmpi('trial',dataline))
        continue % continue means to skip to the next iteration of the loop.
    end
    
    trial_column    = find(strcmpi('trial',   dataline));
    choice_column   = find(strcmpi('choice',  dataline));
    rt_column       = find(strcmpi('rt',      dataline));
    accuracy_column = find(strcmpi('accuracy',dataline));
    
    behavioral_data(datarow,1) = str2double(dataline{trial_column+1});      % Note that we didn't initialize the size of the variable "behavioral_data" so matlab gives a warning.
    behavioral_data(datarow,2) = str2double(dataline{choice_column+1});     % If the variable is relatively small, it doesn't matter. 
    behavioral_data(datarow,3) = str2double(dataline{rt_column+1});         % If the variable is large, however, it's best to initialize it to something really big, and then cut it down to size afterwards.
    behavioral_data(datarow,4) = str2double(dataline{accuracy_column+1});   % See chapter 4 in the book for further discussion of matrix initializations.
    
    datarow=datarow+1; % increment row
end

fclose(fid); % don't forget to close the file after you finish it!

%% initializing variables

num_rows = 10;
num_cols = 35;

% initialize with zeros (typical approach)
largematrix = zeros(num_rows,num_cols);

for rowi=1:num_rows
    for coli=1:num_cols
        % processing here...
    end % end row-loop
end % end column-loop

% note that you can increase the size of a matrix without initializing it
largematrix(num_rows+1,1) = 10;

% similarly, you can increase the dimensionality of a matrix
largematrix(1,round(num_cols/2),3) = 100;

% these last two options (adding elements and dimensions to an existing
% matrix) should be avoided whenever possible. They can create confusion and errors.

% you can also decrease matrix sizes/dimensions, ether by re-assignment:
largematrix = largematrix(:,:,1);

% or by setting parts of the matrix to be empty:
size(largematrix)
largematrix(:,end-4:end) = []; % this removes the last 5 columns
size(largematrix)

% Again, changing matrix sizes and dimensions should be avoided when
% possible, and done carefully when necessary.

%% basic saving data

% save as a .mat file (only matlab can read these files):
save('my_matlab_variables.mat','data','amsterdam','x'); % Question: Why does matlab crash on this line?

% The function 'dlmwrite' is useful if you have a matrix of numbers
% and want to write a text file of only numbers:
dlmwrite('data_written_from_matlab.txt',data,'\t');
% the final argument is the delimieter. This can be tab (\t), space ( ), comma (,), the letter X (X), etc. 

%% advanced saving data

fid = fopen('data_output_SPSS_format.txt','w');

% we want the first row to be variable labels, then rows of mixed string-number data

% variable labels
variable_labels = {'Name';'trial';'choice';'rt';'accuracy'};

% let's add subject names
subject_names={'billy';'bob'};

for vari=1:length(variable_labels)
    fprintf(fid,'%s\t',variable_labels{vari});
    % the %s is for string; %g is for number.
end

% insert a new-line character
fprintf(fid,'\n');

for datarowi=1:size(behavioral_data,1)
    
    % print subject name
    fprintf(fid,'%s\t',subject_names{datarowi});
    
    % now loop through columns (variables)
    for columni=1:size(behavioral_data,2)
        fprintf(fid,'%g\t',behavioral_data(datarowi,columni));
    end
    fprintf(fid,'\n'); % end-of-line 
    
    % You could also do this in one line:
    % fprintf(fid,'%s\t%g\t%g\t%g\t%g\n',subject_names{datarowi},behavioral_data(datarowi,1),behavioral_data(datarowi,2),behavioral_data(datarowi,3),behavioral_data(datarowi,4));
    
    fprintf('Finished writing line %g of %g\n',datarowi,size(behavioral_data,1));
end

fclose(fid);

% Now you can easily import these data into SPSS or Excel.

% Note that if you omit the first argument to fprintf, it puts the output
% in the command instead of the text file, as in the final line of this for-loop.

%% end.
