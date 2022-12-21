function val = getValue(filename, var)
% getValue search for the value for var in the file filename

    switch var
        case 'nnz'
            [s0, r0] = system(['echo `grep -m 1 maxnz ' filename...
        ' | awk ''{print $12}''`']);
        case 'nnx'
            [s0, r0] = system(['echo `grep -m 1 maxnx ' filename...
        ' | awk ''{print $6}'' | sed ''s/,//''`']);
        case 'nny'
            [s0, r0] = system(['echo `grep -m 1 maxny ' filename...
        ' | awk ''{print $9}'' | sed ''s/,//''`']);
        case 'nscl'
            [s0, r0] = system(['echo `grep -m 1 nscl ' filename...
        ' | awk ''{print $6}'' | sed ''s/,//''`']);
        case 'xl'
            [s0, r0] = system(['echo `grep "domain size x" ' filename...
                ' | awk ''{print $3}''`']);
        case 'yl'
            [s0, r0] = system(['echo `grep "domain size y" ' filename...
                ' | awk ''{print $3}''`']);
        case 'zl'
            [s0, r0] = system(['echo `grep "domain size z" ' filename...
                ' | awk ''{print $3}''`']);
        case 'rlat'
            [s0, r0] = system(['echo `grep "rlat " ' filename ...
                ' | awk ''{print $3}''`']);
    end

    if ~s0
        val = str2double(r0);
    end
end
