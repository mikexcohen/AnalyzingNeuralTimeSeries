% 4.7.1
% 1)
random_matrix=rand(4,8)

% 2)
random_matrix_inspection= cell(size(random_matrix));
random_matrix_inspection(:)={NaN}

for i = 1:size(random_matrix,1)
    for j = 1:size(random_matrix,2)
        if random_matrix(i,j)>0.5
            random_matrix_inspection{i,j}='GREATER';
            disp("The value at row " + num2str(i) + " and column " +num2str(j)+ " is equal to " +num2str(random_matrix(i,j)) +" and thus GREATER than 0.5")
        else
            random_matrix_inspection{i,j}='LOWER';
            disp("The value at row " + num2str(i) + " and column " +num2str(j)+ " is equal to " +num2str(random_matrix(i,j)) +" and thus LOWER than 0.5")
        end
    end
end

disp(random_matrix_inspection)

% 3)

function random_matrix_comparison(u,v)

    random_matrix=rand(u,v)
    
    random_matrix_inspection= cell(size(random_matrix));
    random_matrix_inspection(:)={NaN}
    
    for i = 1:size(random_matrix,1)
        for j = 1:size(random_matrix,2)
            if random_matrix(i,j)>0.5
                random_matrix_inspection{i,j}='GREATER';
                disp("The value at row " + num2str(i) + " and column " +num2str(j)+ " is equal to " +num2str(random_matrix(i,j)) +" and thus GREATER than 0.5")
            else
                random_matrix_inspection{i,j}='LOWER';
                disp("The value at row " + num2str(i) + " and column " +num2str(j)+ " is equal to " +num2str(random_matrix(i,j)) +" and thus LOWER than 0.5")
            end
        end
    end
    disp(random_matrix_inspection)
end
random_matrix_comparison(6,4)

