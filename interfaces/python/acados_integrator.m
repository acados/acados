classdef acados_integrator_model:
	

	properties
		acados
		py_model
		model_name
		type
		ode_expr
		x
		u
		xdot
		z
		nx
		nu
		nz
	end


	methods
		
		function obj = acados_integrator_model()
%			obj.numpy = py.importlib.import_module('numpy');
			obj.acados = py.importlib.import_module('acados');
			obj.py_model = obj.acados.acados_integrator_model();
			obj.model_name = []; % check with isempty()
			obj.type = [];
			obj.ode_expr = [];
			obj.x = [];
			obj.xdot = [];
			obj.u = [];
			obj.z = [];
			obj.nx = 0;
			obj.nu = 0;
			obj.nz = 0;
		end
		
		function set(obj, field, value)
			if (field=='model_name')
				obj.model_name = value;
				obj.py_model.set(field, value);
			elseif (field=='type')
				obj.type = value;
				obj.py_model.set(field, value);
			elseif (field=='ode_expr')
				obj.ode_expr = value;
			elseif (field=='x')
				obj.x = value;
				obj.nx = value.size()[1];
				obj.py_model.set(field, int32(obj.nx));
			elseif (field=='xdot')
				obj.xdot = value;
			elseif (field=='u')
				obj.u = value;
				obj.nu = value.size()[1];
				obj.py_model.set(field, int32(obj.nu));
			elseif (field=='z')
				obj.z = value;
				obj.nz = value.size()[1];
				obj.py_model.set(field, int32(obj.nz));
			end
		end
	
	end

end





classdef acados_integrator_opts
	
	properties
		acados
		py_opts
		scheme
		sens_forw
		codgen_model
	end


	methods
		
		function obj = acados_integrator_obj()
			obj.acados = py.importlib.import_module('acados');
			obj.py_opts = obj.acados.acados_integrator_opts();
			obj.scheme = [];
			obj.sens_forw = 'false';
			obj.codgen_model = 'true';
			obj.py_opts.set('codgen_model', 'false');
		end

		function set(obj, field, value)
			if (field=='scheme')
				obj.scheme = value;
				obj.py_opts.set(field, value);
			elseif (field=='sens_forw')
				obj.sens_forw = value;
				obj.py_opts.set(field, value);
			elseif (field=='codgen_model')
				obj.codgen_model = value;
			end
		end
	
	end

end





classdef acados_integrator
	


	properties
		numpy
		acados
		py_model
		py_opts
		py_sim
	end



	methods
		

		function obj = acados_intagrator(model, opts)

			obj.numpy = py.importlib.import_module('numpy');
			obj.acados = py.importlib.import_module('acados');

			% checks

			if ~ ((model.type=='explicit') | (model.type=='implicit'))
				fprintf('\nwrong model type value!\n');
				exit();
			end

			if not ((opts.scheme=='erk') | (opts.scheme=='irk')):
				fprintf('\nwrong opts scheme value!\n');
				exit();
			end

			if not ((opts.sens_forw=='true') | (opts.sens_forw=='false')):
				fprintf('\nwrong opts sens_forw value!\n');
				exit();
			end

			if ((model.type=='explicit') & (opts.scheme!='erk')):
				fprintf('\nwrong opts scheme value for explicit model!\n');
				exit();
			end

			if ((model.type=='implicit') & (opts.scheme!='irk')):
				fprintf('\nwrong opts scheme value for implicit model!\n');
				exit();
			end

			obj.py_model = model.py_model;
			obj.py_opts = model.py_opts;

			obj.codgen_model(model, opts);

			obj.py_sim = obj.acados.acados_integrator(obj.py_model, obj.py_opts);

		end


		function codgen_model(model, opts)

			import casadi.*

			% x
			if isempty(model.x)
				x = SX.sym('x', 0, 1);
			elseif
				x = model.x;
			% xdot
			if isempty(model.xdot)
				xdot = SX.sym('xdot', 0, 1);
			elseif
				xdot = model.xdot;
			% u
			if isempty(model.u)
				u = SX.sym('u', 0, 1);
			elseif
				u = model.u;
			% z
			if isempty(model.z)
				z = SX.sym('z', 0, 1);
			elseif
				z = model.z;

			% fun
			fun = model.ode_expr;

			% sizes
			nx = model.nx;
			nu = model.nu;
			nz = model.nz;

			% define functions & generate C code
			casadi_opts = struct('mex', false, 'casadi_int', 'int', 'casadi_real', 'double');
			c_sources = ' ';

			if (opts.scheme=='erk')
				
				if(opts.sens_forw=='false')
					
					fun_name = 'expl_ode_fun';
					casadi_fun = Function(fun_name, {x, u}, {fun});
					casadi_fun.generate(casadi_opts);
					c_sources = [c_sources, fun_name, '.c'];

				else

					fun_name = 'expl_vde_for';
					Sx = SX.sym('Sx', nx, nx);
					Su = SX.sym('Su', nx, nu);
					vde_x = jtimes(fun, x, Sx);
					vde_u = jacobian(fun, u) + jtimes(fun, x, Su);
					casadi_fun = Function(fun_name, {x, Sx, Su, u}, {fun, vde_x, vde_u});
					casadi_fun.generate(casadi_opts);
					c_sources = [c_sources, fun_name, '.c'];

			elseif (opts.scheme=='irk')
				
				fun_name = 'impl_ode_fun';
				casadi_fun = Function(fun_name, {x, xdot, u, z}, {fun});
				casadi_fun.generate(casadi_opts);
				c_sources = [c_sources, fun_name, '.c'];

				fun_name = 'impl_ode_fun_jac_x_xdot_z';
				jac_x = jacobian(fun, x);
				jac_xdot = jacobian(fun, xdot);
				jac_z = jacobian(fun, z);
				casadi_fun = Function(fun_name, {x, xdot, u, z}, {fun, jac_x, jac_xdot, jac_z});
				casadi_fun.generate(casadi_opts);
				c_sources = [c_sources, fun_name, '.c'];

				if(opts.sens_forw=='true')

					fun_name = 'impl_ode_jac_x_xdot_u_z'
					jac_x = jacobian(fun, x);
					jac_xdot = jacobian(fun, xdot);
					jac_u = jacobian(fun, u);
					jac_z = jacobian(fun, z);
					casadi_fun = Function(fun_name, {x, xdot, u, z}, {jac_x, jac_xdot, jac_u, jac_z});
					casadi_fun.generate(casadi_opts);
					c_sources = [c_sources, fun_name, '.c'];

				end

			end
						
			% create model library
			lib_name = model.model_name + '.so';
			system(['gcc -fPIC -shared ', c_sources, ' -o ', lib_name]);

		end


		function set(field, value)
			
			if (field=='x' | field=='xdot')
				% TODO
			elseif (field=='t')
				% TODO
			end

		end


		function flag = solve()
			
			flag = int64(obj.py_sim.solve());

		end


		function value = get(field)
			
			if (field=='xn' | field=='Sxn')
				py_value = obj.py_sim.get(field);
				value = py2m(py_value, obj.numpy);
			end

		end


	end

end
