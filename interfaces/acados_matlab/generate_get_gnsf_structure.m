function generate_get_gnsf_structure(model, output_dir)

fileID = fopen(fullfile(output_dir, 'get_gnsf_structure.m'), 'w');

fprintf(fileID, 'function model = get_gnsf_structure(model)\n');
fprintf(fileID, 'model.dim_gnsf_nx1 = %d;\n', model.dim_gnsf_nx1);
%fprintf(fileID, 'model.dim_gnsf_nx2 = %d;\n', model.dim_gnsf_nx2);
fprintf(fileID, 'model.dim_gnsf_nz1 = %d;\n', model.dim_gnsf_nz1);
%fprintf(fileID, 'model.dim_gnsf_nz2 = %d;\n', model.dim_gnsf_nz2);
fprintf(fileID, 'model.dim_gnsf_nuhat = %d;\n', model.dim_gnsf_nuhat);
fprintf(fileID, 'model.dim_gnsf_ny = %d;\n', model.dim_gnsf_ny);
fprintf(fileID, 'model.dim_gnsf_nout = %d;\n', model.dim_gnsf_nout);

fclose(fileID);
