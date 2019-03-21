% function to check if running in octave or matlab
function r = is_octave()
	persistent x;
	if (isempty(x))
		x = exist( 'OCTAVE_VERSION', 'builtin');
	end
	r = x;
end

