function Write_file( Result )
       [rows, cols] = size (Result);    
       fid=fopen (  'LP_rank.xls', 'w');
       for i = 1:rows
           for j = 1:cols
               fprintf(fid, '%s\t', Result{i,j});
           end
           fprintf(fid, '\n');
       end
       fclose(fid);      


end

