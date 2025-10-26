import re
import os
import rclpy
import unittest
import launch
import time
import launch_testing
import pytest
import subprocess
import numpy as np

from typing import Union
from unittest import result
from launch_ros.actions import Node


from {{ ros_opts.package_name }}_interface.msg import State, Control, References
{%- if ros_opts.publish_control_sequence -%}
, ControlSequence
{%- endif %}
{%- if dims.np > 0 -%}
, Parameters
{%- endif %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
{%- set control_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.parameters_topic %}
{%- set control_sequence_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic ~ "_sequence" %}
{%- else %}
{%- set control_topic = "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.parameters_topic %}
{%- set control_sequence_topic = "/" ~ ros_opts.control_topic ~ "_sequence" %}
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
        {%- if ros_opts.publish_control_sequence %}
        self.wait_for_publisher('{{ control_sequence_topic }}')
        {%- endif %}
    
    {% if ros_opts.publish_control_sequence %}
    def test_control_sequence_values(self, proc_info):
        """
        Test if the node's published control sequence matches
        the one pre-computed by the Python solver.
        """
        try:
            test_script_dir = os.path.dirname(os.path.realpath(__file__))
            expected_u_file = os.path.abspath(os.path.join(test_script_dir, '..', '..', 'expected_control_sequence.npy'))
            expected_u_sequence = np.load(expected_u_file)
        
        except FileNotFoundError:
            self.skipTest(f"Expected control sequence file not found: {expected_u_file}")

        self.received_control_sequence = None  # Reset before test

        # 1. Subscriber erstellen
        sub = self.node.create_subscription(
            ControlSequence,
            '{{ control_sequence_topic }}',
            self.__control_sequence_callback,
            10
        )

        pub = self.node.create_publisher(
            State,
            '{{ state_topic }}',
            10
        )
        time.sleep(1.0)

        state_msg = State()
        state_msg.x = [0.0, 3.1415926535, 0.0, 0.0] # ocp.constraints.x0 TODO: change accordingly
        state_msg.u = [0.0] * {{ dims.nu }}
        pub.publish(state_msg)

        end_time = time.time() + 10.0
        while time.time() < end_time and self.received_control_sequence is None:
            rclpy.spin_once(self.node, timeout_sec=0.1)

        self.assertIsNotNone(
            self.received_control_sequence,
            "TEST FAILED: Keine Nachricht auf '{{ control_sequence_topic }}' empfangen."
        )
        self.assertEqual(self.received_control_sequence.status, 0, "Solver-Status war nicht erfolgreich (0).")

        expected_length = {{ solver_options.N_horizon }}
        self.assertEqual(len(self.received_control_sequence.control_sequence), expected_length)
        self.assertEqual(expected_u_sequence.shape[0], expected_length)

        for i, control_msg in enumerate(self.received_control_sequence.control_sequence):
            expected_u = expected_u_sequence[i] 
            
            self.assertEqual(
                len(control_msg.u),
                len(expected_u),
                f"Steuervektor bei Schritt {i} hat falsche Dimension."
            )

            for j in range(len(expected_u)):
                self.assertAlmostEqual(
                    control_msg.u[j],
                    expected_u[j],
                    places=2,
                    msg=(f"Wert-Abweichung bei sequence[{i}].u[{j}]. "
                         f"Erhalten: {control_msg.u[j]}, Erwartet: {expected_u[j]}")
                )
    {% endif %}

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
            
    def __control_sequence_callback(self, msg):
        """Callback to store the received control sequence."""
        self.received_control_sequence = msg


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