### How to use the library in python:

Import the module with 

```
import awesome_utils.awesome_pyutils as utils

```

Initialize model interface with

```
urdf_path = "/path/to/urdf/file.urdf"
add_floating_joint = False
model = utils.Model(urdf_path, add_floating_joint)

```

Get robot mass with

```
mass = model.get_robot_mass()
```

Get states dimensions

```
nq = model.get_nq()
nv = model.get_nv()
```

Get number of joints (usually Pinocchio adds a "universe" additional joint)
```
n_jnts = model.get_jnt_number()

```

Set a state vector with 
```
import numpy as np

q = np.zeros(nq)

model.set_q(q)

```
