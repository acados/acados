### Control input rates

This example implements a controller with constraints and cost on the control input rates. In order to do so, the original control input is included as a state and the control input rate becomes the new control input.

The plant model is based on [this paper](https://www.sciencedirect.com/science/article/abs/pii/000510987790070X) and represents an F8 Crusader aircraft.

The original code was posted [here](https://discourse.acados.org/t/penalties-and-constraints-for-control-input-rates-via-augmented-states/1070), while additional information about control input rates can be found [here](https://discourse.acados.org/t/implementing-rate-constraints-and-rate-costs/197).