# Matlab integrator template migration
- [ ] make matlab `acados_sim_template` class similar to acados_sim.py
- [ ] create instance in `acados_sim`
- [ ] render templates
- [ ] build using CMake
- [ ] remove: `sim_create.c`, `sim_precompute.c`, `sim_destroy.c`, `sim_set_ext*.c`
- [ ] make tests pass
- [ ] borrow `sim_set`, `sim_get` or better directly call the shared library!

Possible issues:
- when creating integrator and ocp, files might overwrite each other.
-> we might put all files of the integrator into another folder `build_sim`