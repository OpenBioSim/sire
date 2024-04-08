Generating the wrappers for this module is slightly different to the other modules.
You need to use the custom scanheaders.py script here, and run things from this directory.

Assuming you are in this directory, run:

$ python scanheaders.py /path/to/sire/wrapper/Convert/SireOpenMM .
$ python ../../AutoGenerate/create_wrappers.py

This should work. The `module_info` file describes the module. The contents should be similar to

```
Module SireOpenMM
Source SireOpenMM
Root ../../../corelib/src/libs
```

Note that (currently) the file `vector_less__OpenMM_scope_Vec3__greater_.pypp.cpp`
does not correctly generate with the correct `#include <OpenMM.h>` line.

You can either add this manually after regeneration, or simply revert back
to the working version of this file using

```
git checkout vector_less__OpenMM_scope_Vec3__greater_.pypp.cpp
```

after generating wrappers.
