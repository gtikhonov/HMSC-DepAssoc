function setLevelNames(m, levelNames)

if length(levelNames) ~= m.nr
	error('HMSC: Length of level names vector must be equal to number of levels');
end;
m.levelNames = levelNames;

end