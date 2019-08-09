function check = sim_check_dims(model)

check = 1;

if isfield(model, 'sym_x')
    if ~isempty(model.sym_x)
        nx = length(model.sym_x);
    else
        nx = 0;
    end
    if nx ~= model.dim_nx
        check = 0;
        fail = 'x';
    end
end

if isfield(model, 'sym_u')
    if ~isempty(model.sym_u)
        nu = length(model.sym_u);
    else
        nu = 0;
    end
    if nu ~= model.dim_nu
        check = 0;
        fail = 'u';
    end
end


if isfield(model, 'sym_p')
    if ~isempty(model.sym_p)
        np = length(model.sym_p);
    else
        np = 0;
    end
    if np ~= model.dim_np
        check = 0;
        fail = 'p';
    end
end


if isfield(model, 'sym_z')
    if ~isemzty(model.sym_z)
        nz = length(model.sym_z);
    else
        nz = 0;
    end
    if nz ~= model.dim_nz
        check = 0;
        fail = 'z';
    end
end


if check == 0
    message = strcat('\nSIM_DIM_CHECK FAIL: check consistency of dim_',...
        fail, ' with CasADi symbolic sym_', fail, '!\n\n');
    fprintf(message);
end
