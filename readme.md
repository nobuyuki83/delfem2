[![wercker status](https://app.wercker.com/status/03b6d924ec82270e22a04c3584fbf4de/s/master "wercker status")](https://app.wercker.com/project/byKey/03b6d924ec82270e22a04c3584fbf4de)   [![Documentation Status](https://readthedocs.org/projects/delfem2/badge/?version=master)](https://delfem2.readthedocs.io/en/master/?badge=master)


# DelFEM2

A handy toolset for coding geometry processing and fem simulation

The implementation is based on the [DelFEM](https://github.com/nobuyuki83/DelFEM) library


## How to Build

### Build Examples
```
cd examples
mkdir buildMake
cd buildMake
cmake ..
make
```

### Build Python Binding
```
cd python
mkdir buildMake
cd buildMake
cmake ..
make
```


## Examples

- working directory is in the folder `examples/`
- binary files are put in the folder `examples/bin/`
- input files are put in the folder `test_inputs/`
- The scripts to run the executable is put in `examples/script/`


| Name | Screen Shot |
| ------------- | ------------- |
| triangulation | ![triangulation](docs/screenshot_triangulation.png) |
| transform_handler  | ![handler](docs/screenshot_handler.png) |
| cloth_internal | ![cloth_internal](docs/screenshot_clothinternal.png) |
| subdiv | ![subdiv](docs/screenshot_subdiv.png) |
| read_bvh | ![read_bvh](docs/screenshot_readbvh.png) |
| exponential_map | ![expmap](docs/screenshot_expmap.png) |
| selfcollision_bvh | ![selfcollisionbvh](docs/screenshot_selfcollisionbvh.png) |
| edge_collapse | ![edgecollapse](docs/screenshot_edgecollapse.png) |


