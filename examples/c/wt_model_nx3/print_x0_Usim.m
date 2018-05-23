% load data
load u_x0.mat

% open file
myfile = fopen("u_x0.c", "w");

% x0
fprintf(myfile, "double x0[] = {\n%1.15e, %1.15e, %1.15e\n};\n", x0(1), x0(2), x0(3));

fprintf(myfile, "\n");

% u
fprintf(myfile, "int nsim = %d;\n", size(u,1));
fprintf(myfile, "\n");
fprintf(myfile, "double u_sim[] = {\n");
for ii=1:size(u,1)-1
	fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e,\n", u(ii,1), u(ii,2), u(ii,3), u(ii,4));
end
ii = size(u,1);
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e\n", u(ii,1), u(ii,2), u(ii,3), u(ii,4));
fprintf(myfile, "};\n");

% close file
fclose(myfile);
