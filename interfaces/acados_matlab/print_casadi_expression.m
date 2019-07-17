function print_casadi_expression(f)
    for ii = 1:length(f)
        disp(f(ii,:));
    end
    disp(' ');
end