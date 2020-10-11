%% Analyzing Neural Time Series Data
% Matlab code for Chapter 4 script A
% Mike X Cohen
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% Mike X Cohen assumes no responsibility for inappropriate or incorrect use of this code. 

%% variables, part I
% Let's start with variables. A variable is a place-holder for
% information. To create a variable, you simply assign it information. 
% For example,

mike = 10;
bob  = 20;

% If you type these into the command window (or highlight them here and hit
% F9, or highlight them, right-click, and choose "Evaluate Selection"), you 
% will create the variables.
% Because the variables refer to numbers, you can add them, multiple them,
% etc: 
mike + bob; % Note: You can also put comments after code
(mike+bob)/(mike-bob)

% Notice that when the line ends with a semicolon, the output is suppressed, 
% and when there is no semicolon, it outputs the results in the command window.

% variables can also be strings:
mike = 'mike'; 
% now we've re-assigned mike from a number to a character array. Type
% 'whos' into the matlab command. 
whos

% You can also assign matrices to variables:
a_simple_matrix=[ 3 4 5; 1 2 3; 9 8 7 ];
% type this into the command to see how the semicolon was used to delineate
% separate lines. 

% Square brackets concatinate: 
mikes_full_name = [ mike ' cohen' ];
mikes_two_fav_numbers = [ 7 23 ];

% type whos in the command to see the properties of our variables. Note the 
% difference between double and char (character), and note the sizes of different variables.

%% Variables, part II

% Variables can be more sophisticated. Variables can be cells, which
% are like blocks that may contain different kinds of information. 
var1{1} = [ 1 2 3 4 5 6 7 ];
var1{2} = 'hello world';
var1{3} = [ 1 3 6 7 4 3 5 6 7 87 76 43 4 5 6 767 ];

var1
var1{2}

% The most flexible type of variable is a structure. Structures contain fields that 
% are used for different kinds of data. For example:

ANTS.name = 'mike'; % ANTS = Analyzing Neural Time Series
ANTS.position = 'author';
ANTS.favorite_toothpaste_flavor = 'cinnamon';
ANTS.number_of_watches = 18;
ANTS.favorite_color = [ .8 .1 .8 ]; % RGB values

% You can also have an array of structures
ANTS(2).name = 'Your name here';
ANTS(2).position = 'reader';
ANTS(2).favorite_toothpaste_flavor = 'bratworst'; % gross, but true
ANTS(2).number_of_watches = 1;
ANTS(2).favorite_color = [ 1 1 1 ];

% now you can get information about all fields from one specific element of
% the structure:
ANTS(1)

% or information about one field within one member:
ANTS(1).number_of_watches

% or information about one field from all members:
ANTS.favorite_toothpaste_flavor

% note that this last result came out as two separate answers. If you want
% to combine them into a single output (e.g., a cell array), use curly
% brackets:
{ANTS.favorite_toothpaste_flavor}

%% functions

% functions are modular pieces of code stored in a separate file. Most
% functions can be opened and you can see/modify the code. Some functions
% are compiled and not viewable or editable.

% Functions may take inputs:
randperm(4) % randperm is a function that randomly permutes integers. 4 is the input. 

% or a vector:
mean([1 3 2 4 3 5 4 6])

% to see the guts of this function, highlight "mean" and right-click,
% Open File (or type "edit mean" in the command, or highlight and Ctrl-D)

% IMPORTANT! Do not modify matlab functions unless you really know what
% you're doing! A better idea is to copy the function into a different file
% and use a different name. 

% Most functions also give outputs:
permuted_integers = randperm(4); % now the output of the function is stored in a new variable

whos permuted_in* % Note that you can also use the * character for whos

% some functions can have multiple inputs: 
random_number_matrix = rand(4,6); % Here, we asked for a 4 x 6 matrix of random numbers

% some functions have multiple outputs:
[max_value max_value_index] = max([1 2 3 9 8 7 6]);
% This also shows matrix indexing, which we'll get to soon. 

% IMPORTANT: You can use the output of one function as the input to another
% function. This is a powerful way to make your matlab programming fast and
% efficient. On the other hand, if you embed functions to an extreme you
% might unreadable code.
[max_value max_value_index] = max( randperm( round( rand(1)*10 ) ) );
% Note that when you put the cursor on a parenthesis, it underlines the
% corresponding other parenthesis. 

% type 'help <function_name>' in the matlab command to
% read about a function. 
help max % also try: doc max

%% indexing

% Indexing is a powerful tool in matlab to access particular parts of a
% variable. Indexing is very simple: Imagine you have 100 twinkies arranged
% in a 10 x 10 square: 
twinkies = rand(10,10)

% If you wanted to eat the twinkie in the 4th row, 8th column, you write:
the_twinkie_i_will_eat = twinkies(4,8);

%% The colon operator

% By default, the colon operator increments in integer units from the 
% first to the second number. Observe:
1:10
% You can also increment by a certain amount:
1:2:10
count2ten=1:.23956:10;

% the colon operator is also useful when indexing. Let's say you want to
% eat several twinkies:
twinkies_i_will_eat = twinkies(4:8,2:7);

% Question: How big should the variable twinkies_i_will_eat be? 
whos twin* % answer

% To count backwards, you must specify that you want to skip with a negative number:
rookie_mistake     = 10:1;
how_the_pros_do_it = 10:-1:1;

%% determining the sizes of variables

% You can see variable sizes using whos, but there are also ways 
% to output matrix sizes into variables, which will be useful in many situations.

length(random_number_matrix)
% important! "length" returns the length of the longest dimension, regardless of how many dimensions there are!

% you can use size to find the sizes of all dimensions:
size(twinkies)
size(twinkies,1) % or only specific dimensions...

numel(twinkies) % numel stands for 'total number of elements'

% these functions also produce outputs:
twinkie_array_size = size(twinkies);
prod(twinkie_array_size) % prod returns the product of input numbers

% of course, these functions work on non-numeric variables, e.g.,
length(ANTS)

%% for-loops

% A for-loop is a way to iterate repeatedly:

for counting_variable = 1:10
    disp(counting_variable); % disp stands for display, which prints information in the command window
end

% another example:
for counting_variable = 1:2:10
    disp([ 'The ' num2str(counting_variable) 'th iteration value times 2 divided by 3 and added to 7 is ' num2str(counting_variable*2/3+7) '.' ])
end

% You can embed loops within loops
for i = 1:5
    for j = 3:7
        product_matrix(i,j) = i*j; % Matlab produces a warning here because product_matrix is not initialized. See text in Chapter 4.
    end
end

% Two important things to note here: (1) You can use the same numbers as
% indices AND as variables; (2) Unspecified elements in a matrix are
% automatically created and set to zero. 

number_rows    = 5; % having many spaces is allowed and facilitates code aesthetics
number_columns = 7;
% initialize matrix with zeros
product_matrix = zeros(number_rows,number_columns);

for i=1:number_rows
    for j=1:number_columns
        product_matrix(i,j)=i*j;
    end % end j-loop
end % end i-loop

% note the comments following end-statements. When you have multiple long
% loops, this kind of commenting will be helpful. Also note that when you 
% click on one of the "for" or "end" statements, its pair will be underlined. 

%% if-statements

% Exactly how it sounds
if 4>5
    disp('Something has gone awry in the universe')
end

if 4>5
    disp('Something is still very wrong')
else
    disp('Whew! Everything''s normal.') % note the two single-parenthesis marks inside the string
end

% the 'switch/case' statement is similar to 'if/else'
for counting_variable = 1:2:10
    switch counting_variable
        case 1 % compares 'counting_variable' to '1'
            disp([ 'The ' num2str(counting_variable) 'st iteration value times 2 divided by 3 and added to 7 is ' num2str(counting_variable*2/3+7) '.' ])
        case 2
            disp([ 'The ' num2str(counting_variable) 'nd iteration value times 2 divided by 3 and added to 7 is ' num2str(counting_variable*2/3+7) '.' ])
        case 3
            disp([ 'The ' num2str(counting_variable) 'rd iteration value times 2 divided by 3 and added to 7 is ' num2str(counting_variable*2/3+7) '.' ])
        otherwise
            disp([ 'The ' num2str(counting_variable) 'th iteration value times 2 divided by 3 and added to 7 is ' num2str(counting_variable*2/3+7) '.' ])
    end % end switch
end

%% boolean (true/false)

% sometimes, it's useful to know if a statement is TRUE or FALSE. 
% You can use a double-equals sign:

5==5

% The answer is '1', which in this case means 'true'. The oppose of true (1) is false (0):
5==4

% You can assign these answers to variables:
fourIsFive = 4==5; % see below for disambiguating = and ==

whos fourIsFive
% note that 'fourIsFive' is type "logical", or boolean. 

% Other related useful functions:
matrix0 = false(10,3);
matrix1 = true(4,12);

% you can use the tilde "~" to negate a statement:
~(1==1) % false because the negation of 1==1 is false
~(4>3)  % also false because the negation of 4>3 is false
~(4<3)  % true because it is false (~) that 4 is greater than 3. Tricky!

% things can sometimes get tricky:
truthtest = 1 == 2;

% Remember:
%    One equals sign is a statement ("you have this value").
%    Two equals signs means you are asking a question ("are these the same?").
% Mnemonic: "Say with one (=), ask with two (==), then you know just who is who!" 
%           (I didn't say it was a good mnemonic...)

%% repmat

% repmat is short for "replicate matrix." It is a handy function that
% is often used in analyses. The problem that repmat solves is matrix
% operations with other matrices. First, consider that adding a scalar 
% (a single number) to a matrix is no problem. 

mat = rand(4);

mat+23

% Now imagine that mat is EEG data with 4 electrodes and 10 time points:
mat = rand(4,10);

% now you want to subtract the mean over time:
meanOverTime = mean(mat,2); % second input to mean function is the dimension along which to compute the mean

mat = mat - meanOverTime;

% the previous line crashes because the matrices are of unequal size
whos mat meanOverTime

% repmat will help:
meanOverTimeRepmat = repmat(meanOverTime,1,size(mat,2));
whos mat meanOverTime*
% repmat takes 3 inputs: the matrix you want to replicate, the number of
% times to replicate it over rows, and the number of times to replicate it
% over columns (you can also input more dimensions). We want to replicate
% this matrix only over time points (the size of the second dimension of
% mat). Now the subtraction works:

matmean = mat - meanOverTimeRepmat;


% Other examples of repmat:
mat = [1 2 3; 10 20 30];

repmat(mat,1,1)
repmat(mat,1,2)
repmat(mat,2,1)

%% bsxfun

% bsxfun is a useful function for fast and easy array and matrix manipulations.
% It was introduced to Matlab fairly recently, so older versions of Matlab
% do not have this utility. 

% for example, the following function will add 4 to a random matrix:
bsxfun(@plus,randn(10),4)

% this might not seem any better than "randn(10)+4" and for this small
% case, it isn't. bsxfun is more useful because it performs
% singleton-expansion, which means you may be able to avoid using repmat.
% For example, imagine a dataset with 100 channels and 100,000 time points: 
a  = rand(100,100000);

% To subtract the mean of the entire time series: 
am = a - repmat(mean(a,2),1,size(a,2));

% Notice that the repmat is necessary:
am = a - mean(a,2);

% The previous line crashes because the sizes of a and its mean are not the
% same. However, bsxfun expands this automatically
am = bsxfun(@minus,a,mean(a,2));

% let's do a timing test...
tic
for i=1:100
    am = a - repmat(mean(a,2),1,size(a,2));
end
t(1)=toc;

tic
for i=1:100
    am = bsxfun(@minus,a,mean(a,2));
end
t(2)=toc;

figure
bar(t)
set(gca,'xtick',1:2,'xticklabel',{'repmat';'bsxfun'},'xlim',[.5 2.5])
title([ 'bsxfun took ' num2str(100*t(2)/t(1)) '% of the computation time.' ])
% you'll learn more about the above lines in part B of this code


% Thus, bsxfun is a bit faster but also more convenient to use, and more elegant. 
% There are other similar functions to bsxfun, including arrayfun and cellfun. 
% In addition to speed, the *fun functions allow you to avoid using loops. 
% For example, let's say you want to know the length of items in a cell array. 

% create cell array whose elements have variable lengths
c = cell(1,40);
for i=1:length(c)
    c{i} = randn(1,round(rand*100));
end

% now you want to know how many elements are in each cell. 
% Normally you need a loop:
cell_lengths = zeros(size(c));
for i=1:length(c)
    cell_lengths(i) = numel(c{i});
end

% But cellfun is more efficient:
cell_lengths = cellfun(@length,c);

%% end

% continue on the part b of this code...
