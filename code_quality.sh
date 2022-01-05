pylama -l radon volmdlr
if [[ $(pylama -l radon volmdlr) ]];
	then echo "Error in code quality check">&2;
fi;
