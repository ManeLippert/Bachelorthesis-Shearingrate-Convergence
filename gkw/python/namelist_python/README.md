# A Python module to parse Fortran namelist files

Read in a namelist file:
```
import namelist_python
namelist = namelist_python.read_namelist_file('config.dat')
namelist.groups['foo']['bar']
```

This creates an instance of `namelist_python.Namelist` whose attribute
`groups` holds the data in a nested ordered case-insensitive dictionary
structure.

Write a `Namelist` object back to a file:
```
with open('new_file.dat', 'w') as f:
	f.write(namelist.dump())
```

`dump` takes an optional argument `array_inline` a boolean which sets whether
arrays should be inline or given in index notation.

If you use ipython or a another interactive REPL prompt you may want to use 
the `data` attribute which allows you to do tab completion on the group and
variable names:

```
In [7]: namelist.data.ATHAM_SETUP.dt
namelist.data.ATHAM_SETUP.dt
namelist.data.ATHAM_SETUP.dtmax
namelist.data.ATHAM_SETUP.dtmin
In [7]: namelist.data.ATHAM_SETUP.dt
Out[7]: 3.0

In [8]: namelist.data.ATHAM_SETUP.dt = 4.0

In [9]: namelist.data.ATHAM_SETUP.dt
Out[9]: 4.0
```

## Features
 - Parses ints, floats, booleans, escaped strings and complex numbers.
 - Can output in namelist format.
 - Tab-completion and variable assignment in interactive console

## Missing features
 - Currently can't handle line continuations
 - Currently can't handle lines with several parameters
 - Comments are not kept, and so won't exist in output.
 - Module does not help to create a Namelist object from scratch

## Contribute
Please send any namelist files that don't parse correctly or fix the code
yourself and send a pull request :)

Thanks
