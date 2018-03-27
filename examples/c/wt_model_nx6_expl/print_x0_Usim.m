% load data
load testSim.mat

% open file
myfile = fopen("u_x0.c", "w");

% x0
%fprintf(myfile, "double x0[] = {\n%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e\n};\n", statesFAST(1,1), statesFAST(1,2), statesFAST(1,3), statesFAST(1,4), statesFAST(1,5), statesFAST(1,6));

%fprintf(myfile, "\n");

% u
fprintf(myfile, "int nsim = %d;\n", size(Usim,1));
fprintf(myfile, "\n");
fprintf(myfile, "double u_sim[] = {\n");
for ii=1:size(Usim,1)-1
	fprintf(myfile, "%1.15e, %1.15e, %1.15e,\n", Usim(ii,1), Usim(ii,2), Usim(ii,3));
end
ii = size(Usim,1);
fprintf(myfile, "%1.15e, %1.15e, %1.15e\n", Usim(ii,1), Usim(ii,2), Usim(ii,3));
fprintf(myfile, "};\n");

fprintf(myfile, "\n");

% x_ref
fprintf(myfile, "\n");
fprintf(myfile, "double x_ref[] = {\n");
for ii=1:size(Usim,1)-1
	fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e,\n", statesFAST(ii,1), statesFAST(ii,2), statesFAST(ii,3), statesFAST(ii,4), statesFAST(ii,5), statesFAST(ii,6));
end
ii = size(Usim,1);
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e\n", statesFAST(ii,1), statesFAST(ii,2), statesFAST(ii,3), statesFAST(ii,4), statesFAST(ii,5), statesFAST(ii,6));
fprintf(myfile, "};\n");

% close file
fclose(myfile);
