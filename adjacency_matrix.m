A=full(createRandRegGraph(900,5));
fileID = fopen('deg5_900.txt','w');
fprintf(fileID,'%d\n',A);
fclose(fileID);

%for finding eigen values we change adjacency matrix
%to probability transition matrix
%{
for c = 1:1000
    for r = 1:1000

        if (A(r,c) == 1)
            B(r,c) = 1/5;
       else
            B(r,c) = 0;
        end

    end
end
e=eig(B);
%}