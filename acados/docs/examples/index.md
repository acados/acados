# Examples

Examples can be found in the [`examples` folder of the `acados` repository](https://github.com/acados/acados/tree/master/examples).

## Matlab and Octave
Examples for `Matlab` and `Octave` can be found in the [folder `examples/acados_matlab_octave` of the `acados` repository](https://github.com/acados/acados/tree/master/examples/acados_matlab_octave).

In order to run and understand these examples, we refer to the documentation of the [`Matlab` and `Octave` interface](../matlab_octave_interface/index.md).

## Python
Examples for `Python` can be found in the [folder `examples/acados_python` of the `acados` repository](https://github.com/acados/acados/tree/master/examples/acados_python).

In order to run and understand these examples, we refer to the documentation of the [`Python` interface](../python_interface/index.md).

## C
Examples for `C` can be found in the [folder `examples/c` of the `acados` repository](https://github.com/acados/acados/tree/master/examples/c).

In order to run and understand these examples, we refer to the documentation of the [`C` interface](../c_interface/index.md).

Note that the `Matlab` and `Python` interfaces can be used to codegenerate `C` examples which are cleaner than some of the examples in [`examples/c`](https://github.com/acados/acados/tree/master/examples/c).
A recommended workflow is thus to prototype an NMPC controller from one of the high-level interfaces and deploy the generated code with minor modifications in a `C`, `C++` or `ROS` framework.
