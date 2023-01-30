function write_xyz(filename,data,element_type,element_name)
fid = fopen(filename,'w+');
    if fid==-1
        error(['Error opening ' filename]); 
    end
     L_element_type=length(element_type);
     fprintf(fid,[num2str(L_element_type), '\n']);
     fprintf(fid, '\n');
     for i=1:L_element_type
         ele_type= element_type(i);
         ele_name = element_name{ele_type};
         fprintf(fid, '%s  ', ele_name);
         fprintf(fid, ' %10.6f %10.6f %10.6f \n', data(i,:));
     end
   
  fprintf(fid, '\n');
%   fprintf(fid, '%19.0f %2.0f %f %f %f\n', new_data');
  fclose(fid);