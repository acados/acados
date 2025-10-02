# Generic reference message
std_msgs/Header header

# Reference vectors
{%- if dims.ny_0 > 0 %}
float64[{{ dims.ny_0 }}] yref_0
{%- endif %}
{%- if dims.ny > 0 %}
float64[{{ dims.ny }}] yref
{%- endif %}
{%- if dims.ny_e > 0 %}
float64[{{ dims.ny_e }}] yref_e
{%- endif %}
