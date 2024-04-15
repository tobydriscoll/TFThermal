function struct_string(p)
    str = ""
	for s in propertynames(p)
		v = getproperty(p, s)
		if v isa Quantity
			str = str * String(s) * " = " * "$(round(typeof(v), v, sigdigits=5))" * ", "
		else
			str = str * String(s) * " = " * "$(round(v, sigdigits=5))" * ", "
		end
	end
    return isempty(str) ? "" : str[1:end-2]
end
