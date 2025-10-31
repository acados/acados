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
from rclpy.parameter import Parameter, parameter_value_to_python
from rcl_interfaces.srv import GetParameters, SetParameters
from rclpy.publisher import Publisher

from {{ ros_opts.package_name }}_interface.msg import State, Control, References
{%- if ros_opts.publish_control_sequence -%}
, ControlSequence
{%- endif %}
{%- if dims.np > 0 -%}
, Parameters
{%- endif %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
    {%- set topic_prefix = '/' ~ ns ~ '/' %}
{%- else %}
    {%- set topic_prefix = '/' %}
{%- endif %}
{%- set control_topic = topic_prefix ~ ros_opts.control_topic %}
{%- set state_topic = topic_prefix ~ ros_opts.state_topic %}
{%- set references_topic = topic_prefix ~ ros_opts.reference_topic %}
{%- set parameters_topic = topic_prefix ~ ros_opts.parameters_topic %}
{%- set control_sequence_topic = control_topic ~ '_sequence' %}


class AsyncParametersClient:
    """Simple async parameter client."""
    def __init__(self, node, remote_node_name):
        self.node = node
        self.remote_node_name = remote_node_name
        self._get_client = self.node.create_client(
            GetParameters,
            f'/{remote_node_name}/get_parameters')
        self._set_client = self.node.create_client(
            SetParameters,
            f'/{remote_node_name}/set_parameters')

    def wait_for_service(self, timeout_sec=1.0):
        ok1 = self._get_client.wait_for_service(timeout_sec=timeout_sec)
        ok2 = self._set_client.wait_for_service(timeout_sec=timeout_sec)
        return ok1 and ok2

    def get_parameters(self, names):
        request = GetParameters.Request()
        request.names = names
        return self._get_client.call_async(request)
    
    def set_parameters(self, params: list[tuple[str, object]]):
        """Set parameters synchron via service (returns Future)."""
        request = SetParameters.Request()
        for name, value in params:
            if isinstance(value, (list, tuple)):
                value = [float(v) for v in value]
            if isinstance(value, (int, float)):
                value = float(value)
            p = Parameter(name, value=value)
            request.parameters.append(p.to_parameter_msg())
        return self._set_client.call_async(request)


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
        self.param_client = AsyncParametersClient(self.node, '{{ ros_opts.node_name }}')

    def tearDown(self):
        self.node.destroy_node()

    def test_set_constraints(self, proc_info):
        """Test if constraints compile-time declared default parameters."""
        {%- for field, param in constraints %}
        {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
        param_name = "constraints.{{ field }}"
        expected_value = [{{- param | join(sep=', ') -}}]
        self._check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}

    def test_set_cost(self, proc_info):
        """Test if cost compile-time declared default parameters."""
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
        self._check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}
        {%- if has_slack %}

        # --- Slacks ---
        {%- for field, param in cost %}
        {%- set field_l = field | lower %}
        {%- if param and (field_l is starting_with('z')) %}
        param_name = "cost.{{ field }}"
        expected_value = [{{- param | join(sep=', ') -}}]
        self._check_parameter_set(param_name, expected_value)
        {%- endif %}
        {%- endfor %}
        {%- endif %}

    def test_set_solver_options(self, proc_info):
        """Test if solver options compile-time declared default parameters."""
        param_name = "ts"
        expected_value = {{ solver_options.Tsim }}
        self._check_parameter_set(param_name, expected_value)

    def test_subscribing(self, proc_info):
        """Test if the node subscribes to all expected topics."""
        self._wait_for_subscription('{{ state_topic }}')
        self._wait_for_subscription('{{ references_topic }}')
        {%- if dims.np > 0 %}
        self._wait_for_subscription('{{ parameters_topic }}')
        {%- endif %}

    def test_publishing(self, proc_info):
        """Test if the node publishes to all expected topics."""
        self._wait_for_publisher('{{ control_topic }}')
        {%- if ros_opts.publish_control_sequence %}
        self._wait_for_publisher('{{ control_sequence_topic }}')
        {%- endif %}
    
    {% if ros_opts.publish_control_sequence %}
    def test_control_sequence_values(self, proc_info):
        """
        Test if the node's published control sequence matches
        the one pre-computed by the Python solver.
        """
        expected_u_sequence = self._load_expected_sequence()

        self.received_control_sequence = None
        sub = self.node.create_subscription(
            ControlSequence,
            '{{ control_sequence_topic }}',
            self._control_sequence_callback,
            10
        )

        pub = self.node.create_publisher(State, '{{ state_topic }}', 10)
        self._publish_initial_state(pub)

        condition_met = self._wait_for_condition(
            condition_check=lambda: self.received_control_sequence is not None,
            timeout=10.0
        )

        self.assertTrue(
            condition_met,
            "TEST FAILED: '{{ control_sequence_topic }}' message not received (Timeout)."
        )
        self.assertEqual(self.received_control_sequence.status, 0, "Solver status was not successful (0).")

        expected_length = {{ solver_options.N_horizon }}
        self.assertEqual(len(self.received_control_sequence.control_sequence), expected_length)
        self.assertEqual(expected_u_sequence.shape[0], expected_length)

        for i, control_msg in enumerate(self.received_control_sequence.control_sequence):
            expected_u = expected_u_sequence[i] 
            
            self.assertEqual(
                len(control_msg.u),
                len(expected_u),
                f"Control vector at step {i} has wrong dimension."
            )

            for j in range(len(expected_u)):
                self.assertAlmostEqual(
                    control_msg.u[j],
                    expected_u[j],
                    places=2,
                    msg=(f"Value deviation at sequence[{i}].u[{j}]. "
                         f"Received: {control_msg.u[j]}, Expected: {expected_u[j]}")
                )
        
    def _control_sequence_callback(self, msg):
        """Callback to store the received control sequence."""
        self.received_control_sequence = msg

    def _load_expected_sequence(self):
        test_script_dir = os.path.dirname(os.path.realpath(__file__))
        expected_u_file = os.path.abspath(
            os.path.join(test_script_dir, '..', '..', 'expected_control_sequence.npy')
        )
        if not os.path.exists(expected_u_file):
            self.skipTest(f"Expected control sequence file not found: {expected_u_file}")
        return np.load(expected_u_file)
    
    def _publish_initial_state(self, publisher: Publisher, timeout_sec: float = 2.0):
        state_msg = State()
        {%- if constraints.has_x0 %}
        state_msg.x = [float(v) for v in ({{ constraints.lbx_0 | join(sep=', ') }})]
        {%- else %}
        state_msg.x = [0.0] * {{ dims.nx }}
        {%- endif %}

        end_time = time.time() + timeout_sec
        while time.time() < end_time and publisher.get_subscription_count() < 1:
            rclpy.spin_once(self.node, timeout_sec=0.1)

        if publisher.get_subscription_count() < 1:
            self.fail(f"Publisher on Topic '{publisher.topic_name}' "
                      f"could not find a Subscriber (Timeout).")

        publisher.publish(state_msg)
        rclpy.spin_once(self.node, timeout_sec=0.1)
    {% endif %}

    def _wait_for_subscription(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            rclpy.spin_once(self.node, timeout_sec=0.1)
            subs = self.node.get_subscriptions_info_by_topic(topic)
            if subs:
                return True
        self.fail(f"Node has NOT subscribed to '{topic}'.")

    def _wait_for_publisher(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            rclpy.spin_once(self.node, timeout_sec=0.1)
            pubs = self.node.get_publishers_info_by_topic(topic)
            if pubs:
                return True
        self.fail(f"Node has NOT published to '{topic}'.")

    def _wait_for_condition(self, condition_check: callable, timeout: float = 10.0) -> bool:
        """Spin the node until a condition is met or a timeout occurs."""
        end_time = time.time() + timeout
        while time.time() < end_time:
            rclpy.spin_once(self.node, timeout_sec=0.1)
            if condition_check():
                return True
        return False

    def _check_parameter_get(self, param_name: str, expected_value: Union[list[float], float]):
        """Run a subprocess command and return its output."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail(f"Parameter service for '{{ ros_opts.node_name }}' not available.")

        future = self.param_client.get_parameters([param_name])
        rclpy.spin_until_future_complete(self.node, future, timeout_sec=2.0)
        if not future.done():
            self.fail(f"Timeout while getting parameter {param_name}")
        resp = future.result()
        if not resp.values:
            self.fail(f"Parameter {param_name} does not exist.")
        actual_value = parameter_value_to_python(resp.values[0])

        if isinstance(expected_value, list):
            self.assertListEqual(
                list(actual_value), expected_value, 
                f"Parameter {param_name} has the wrong value! Got {actual_value}")
        else:
            self.assertEqual(
                actual_value, expected_value, 
                f"Parameter {param_name} has the wrong value! Got {actual_value}")

    def _check_parameter_set(self, param_name: str, new_value: Union[list[float], float]):
        """Run a subprocess command and return its output."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail(f"Parameter service for '{{ ros_opts.node_name }}' not available.")

        future = self.param_client.set_parameters([(param_name, new_value)])
        rclpy.spin_until_future_complete(self.node, future, timeout_sec=3.0)
        if not future.done() or future.result() is None:
            self.fail(f"Timeout while setting parameter {param_name}")

        results = future.result().results
        if not results or not results[0].successful:
            self.fail(f"Failed to set parameter {param_name}: {results[0].reason if results else 'unknown'}")

        self._check_parameter_get(param_name, new_value)