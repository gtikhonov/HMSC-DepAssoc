function setY(m, Y)

if any(size(Y) ~= [m.ny, m.ns])
	error('HMSC: Y dimentions are not compatible with specified');
end
m.Y = Y;

end