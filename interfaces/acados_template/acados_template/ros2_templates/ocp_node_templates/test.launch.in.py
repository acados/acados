import re
from typing import Union
from unittest import result
import rclpy
import unittest
import launch
import time
import launch_testing
import pytest
import subprocess
from launch_ros.actions import Node

from {{ ros_opts.package_name }}_interface.msg import State, ControlInput, References
{%- if dims.np > 0 -%}
, Parameters
{%- endif %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
{%- set control_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.parameters_topic %}
{%- else %}
{%- set control_topic = "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.parameters_topic %}
{%- endif %}

@pytest.mark.launch_test
def generate_test_description():
    """Generate launch description for node testing."""
    start_{{ ros_opts.node_name }} = Node(
        package='{{ ros_opts.package_name }}',
        executable='{{ ros_opts.node_name }}',
        name='{{ ros_opts.node_name }}'
    )

    return launch.LaunchDescription([
        start_{{ ros_opts.node_name }},
        launch.actions.TimerAction(
                    period=5.0, actions=[launch_testing.actions.ReadyToTest()]),
    ])


class GeneratedNodeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rclpy.init()

    @classmethod
    def tearDownClass(cls):
        rclpy.shutdown()

    def setUp(self):
        self.node = rclpy.create_node('generated_node_test')

    def tearDown(self):
        self.node.destroy_node()

    def test_set_constraints(self, proc_info):
        """
        Test if constraints compile-time declared default parameters.
        """
        {%- for field, param in constraints %}
        {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
        param_name = "constraints.{{ field }}"
        expected_value = [{{- param | join(sep=', ') -}}]
        self.__check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}

    def test_set_cost(self, proc_info):
        """
        Test if cost compile-time declared default parameters.
        """
        # --- Weights ---
        {%- for field, param in cost %}
        {%- if param and (field is starting_with('W')) %}
        param_name = "cost.{{ field }}"
        expected_value = [
        {%- set n_diag = param | length -%}
        {%- for i in range(end=n_diag) -%}
            {{- param[i][i] -}}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
        ]
        self.__check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}
        {%- if has_slack %}

        # --- Slacks ---
        {%- for field, param in cost %}
        {%- set field_l = field | lower %}
        {%- if param and (field_l is starting_with('z')) %}
        param_name = "cost.{{ field }}"
        expected_value = [{{- param | join(sep=', ') -}}]
        self.__check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}
        {%- endif %}

    def test_set_solver_options(self, proc_info):
        """
        Test if solver options compile-time declared default parameters.
        """
        # --- Solver Options ---
        param_name = "ts"
        expected_value = {{ solver_options.Tsim }}
        self.__check_parameter_set(param_name, expected_value)

    def test_subscribing(self, proc_info):
        """Test if the node subscribes to all expected topics."""
        self.wait_for_subscription('{{ state_topic }}')
        self.wait_for_subscription('{{ references_topic }}')
        {%- if dims.np > 0 %}
        self.wait_for_subscription('{{ parameters_topic }}')
        {%- endif %}

    def test_publishing(self, proc_info):
        """Test if the node publishes to all expected topics."""
        self.wait_for_publisher('{{ control_topic }}')

    def wait_for_subscription(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            subs = self.node.get_subscriptions_info_by_topic(topic)
            if subs:
                return True
            time.sleep(0.1)
        self.fail(f"Node has NOT subscribed to '{topic}'.")

    def wait_for_publisher(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            pubs = self.node.get_publishers_info_by_topic(topic)
            if pubs:
                return True
            time.sleep(0.1)
        self.fail(f"Node has NOT published to '{topic}'.")

    def __check_parameter_get(self, param_name: str, expected_value: Union[list[float], float]):
        """Run a subprocess command and return its output."""
        output = get_parameter(param_name)
        numbers = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", output)]
        if isinstance(expected_value, list):
            self.assertListEqual(numbers, expected_value, f"Parameter {param_name} has the wrong value! Got {numbers}")
        else:
            self.assertEqual(numbers[0], expected_value, f"Parameter {param_name} has the wrong value! Got {numbers[0]}")

    def __check_parameter_set(self, param_name: str, new_value: Union[list[float], float]):
        """Run a subprocess command and return its output."""
        try:
            set_parameter(param_name, new_value)
            self.__check_parameter_get(param_name, new_value)
        except subprocess.CalledProcessError as e:
            self.fail(f"Failed to set parameter {param_name}.\n"
                      f"Exit-Code: {e.returncode}\n"
                      f"Stderr: {e.stderr}\n"
                      f"Stdout: {e.stdout}")


def get_parameter(param_name: str):
    """Run a subprocess command and return its output."""
    cmd = ['ros2', 'param', 'get', '{{ ros_opts.node_name }}', param_name]
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    return result.stdout

def set_parameter(param_name: str, value: Union[list[float], float]):
    """Run a subprocess command to set a parameter."""
    if isinstance(value, list):
        value_str = "[" + ",".join(map(str, value)) + "]"
    else:
        value_str = str(value)

    cmd = ['ros2', 'param', 'set', '{{ ros_opts.node_name }}', param_name, value_str]
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    return result.stderr