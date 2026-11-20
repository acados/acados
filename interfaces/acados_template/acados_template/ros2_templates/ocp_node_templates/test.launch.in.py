import math
import os
import rclpy
import unittest
import launch
import time
import launch_testing
import pytest

from typing import Union
from launch_ros.actions import Node
from rclpy.parameter import Parameter, parameter_value_to_python
from rcl_interfaces.srv import GetParameters, SetParameters
{%- set has_slack = dims.ns > 0 or dims.ns_0 > 0 or dims.ns_e > 0 %}
{%- set has_yref = dims.ny_0 > 0 or dims.ny > 0 or dims.ny_e > 0 %}
{%- set has_W_0 = dims.ny_0 > 0 and ((cost.cost_type_0 == "LINEAR_LS") or (cost.cost_type_0 == "NONLINEAR_LS")) %}
{%- set has_W = dims.ny > 0 and ((cost.cost_type == "LINEAR_LS") or (cost.cost_type == "NONLINEAR_LS")) %}
{%- set has_W_e = dims.ny_e > 0 and ((cost.cost_type_e == "LINEAR_LS") or (cost.cost_type_e == "NONLINEAR_LS")) %}
{%- if ros_opts.publish_control_sequence %}
import numpy as np
from {{ ros_opts.package_name }}_interface.msg import ControlSequence{% if constraints.has_x0 %}, State{% endif %}
{%- if constraints.has_x0 %}
from rclpy.publisher import Publisher
{%- endif %}
{%- endif %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
    {%- set topic_prefix = '/' ~ ns ~ '/' %}
    {%- set node_fqn = '/' ~ ns ~ '/' ~ ros_opts.node_name %}
{%- else %}
    {%- set topic_prefix = '/' %}
    {%- set node_fqn = '/' ~ ros_opts.node_name %}
{%- endif %}
{%- set raw_control_topic = ros_opts.control_topic | trim(chars='/') %}
{%- set raw_state_topic = ros_opts.state_topic | trim(chars='/') %}
{%- set raw_reference_topic = ros_opts.reference_topic | trim(chars='/') %}
{%- set raw_parameters_topic = ros_opts.parameters_topic | trim(chars='/') %}
{%- set control_topic = topic_prefix ~ raw_control_topic %}
{%- set state_topic = topic_prefix ~ raw_state_topic %}
{%- set references_topic = topic_prefix ~ raw_reference_topic %}
{%- set parameters_topic = topic_prefix ~ raw_parameters_topic %}
{%- set control_sequence_topic = control_topic ~ '_sequence' %}


class AsyncParametersClient:
    """Simple async parameter client."""
    def __init__(self, node, remote_node_name):
        self.node = node
        self.remote_node_name = remote_node_name
        node_prefix = remote_node_name if remote_node_name.startswith('/') else f'/{remote_node_name}'
        self._get_client = self.node.create_client(
            GetParameters,
            f'{node_prefix}/get_parameters')
        self._set_client = self.node.create_client(
            SetParameters,
            f'{node_prefix}/set_parameters')

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
            elif isinstance(value, bool):
                value = bool(value)
            elif isinstance(value, int):
                value = int(value)
            elif isinstance(value, float):
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
        name='{{ ros_opts.node_name }}'{%- if ns %},
        namespace='{{ ns }}'{%- endif %}
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
        self.param_client = AsyncParametersClient(self.node, '{{ node_fqn }}')

    def tearDown(self):
        self.node.destroy_node()

    def test_set_constraints(self, proc_info):
        """Test declared constraints parameters and boundary/negative updates."""
        {%- for field, param in constraints %}
        {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (('bx_0' not in field) or (constraints.has_x0 == false)) %}
        param_name = "constraints.{{ field }}"
        expected_value = [float(v) for v in [{{- param | join(sep=', ') -}}]]
        self._check_parameter_set(param_name, expected_value)
        self._check_parameter_set_rejected(param_name, expected_value + [0.0], unchanged_expected=expected_value)
        {%- endif %}
        {%- endfor %}

    def test_set_cost(self, proc_info):
        """Test declared cost parameters, diagonal shorthand, and invalid-length rejection."""
        # --- Weights ---
        {%- for field, param in cost %}
        {%- if param and ((field == "W_0" and has_W_0) or (field == "W" and has_W) or (field == "W_e" and has_W_e)) %}
        param_name = "cost.{{ field }}"
        expected_value = [float(v) for v in [
        {%- for row in param -%}
            {%- for val in row -%}
                {{- val -}}
                {%- if not loop.last %}, {% endif -%}
            {%- endfor -%}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
        ]]
        weight_dim = math.isqrt(len(expected_value))
        diag_value = [expected_value[i * weight_dim + i] for i in range(weight_dim)]
        self._check_parameter_set(param_name, diag_value)
        self._check_parameter_set(param_name, expected_value)
        self._check_parameter_set_rejected(param_name, expected_value + [0.0], unchanged_expected=expected_value)
        {%- endif %}
        {%- endfor %}
        {%- if has_slack %}

        # --- Slacks ---
        {%- for field, param in cost %}
        {%- set field_l = field | lower %}
        {%- if param and (field_l is starting_with('z')) %}
        param_name = "cost.{{ field }}"
        expected_value = [float(v) for v in [{{- param | join(sep=', ') -}}]]
        self._check_parameter_set(param_name, expected_value)
        self._check_parameter_set_rejected(param_name, expected_value + [0.0], unchanged_expected=expected_value)
        {%- endif %}
        {%- endfor %}
        {%- endif %}

    def test_set_solver_options(self, proc_info):
        """Test solver/ROS parameter updates and invalid-value rejection."""
        default_ts = float({{ solver_options.Tsim }})
        self._check_parameter_set("ts", default_ts)
        self._check_parameter_set("ts", default_ts * 2.0)
        self._check_parameter_set("ts", default_ts)

        for bad_ts in (0.0, -0.1, float('nan'), float('inf')):
            self._check_parameter_set_rejected("ts", bad_ts, unchanged_expected=default_ts)

        default_print_level = int({{ solver_options.print_level | default(value=0) }})
        self._check_parameter_set("solver_options.print_level", 1 if default_print_level == 0 else 0)
        self._check_parameter_set("solver_options.print_level", default_print_level)

        self._check_parameter_set("verbose", True)
        self._check_parameter_set("verbose", False)

        self._check_parameter_set_rejected("nonexistent_parameter", 1.0)

    def test_subscribing(self, proc_info):
        """Test if the node subscribes to all expected topics."""
        {%- if constraints.has_x0 %}
        self._wait_for_subscription('{{ state_topic }}')
        {%- endif %}
        {%- if has_yref %}
        self._wait_for_subscription('{{ references_topic }}')
        {%- endif %}
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

        {%- if constraints.has_x0 %}
        pub = self.node.create_publisher(State, '{{ state_topic }}', 10)
        self._publish_initial_state(pub)
        {%- endif %}

        # Wait until the solver converges to the expected control sequence (e.g. across multiple ticks for SQP_RTI or SQP).
        self._wait_for_condition(
            condition_check=lambda: self._control_sequence_matches(
                self.received_control_sequence, expected_u_sequence),
            timeout=10.0
        )

        msg = self.received_control_sequence
        self.assertIsNotNone(
            msg,
            "TEST FAILED: '{{ control_sequence_topic }}' message not received (Timeout)."
        )
        self.assertEqual(msg.status, 0, "Solver status was not successful (0).")

        expected_length = {{ solver_options.N_horizon }}
        self.assertEqual(len(msg.control_sequence), expected_length)
        self.assertEqual(expected_u_sequence.shape[0], expected_length)

        for i, control_msg in enumerate(msg.control_sequence):
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

    def _control_sequence_matches(self, msg, expected_u_sequence, places: int = 2) -> bool:
        """Check whether a received control sequence matches the expected values within tolerance."""
        if msg is None or msg.status != 0 or len(msg.control_sequence) != expected_u_sequence.shape[0]:
            return False
        for control_msg, expected_u in zip(msg.control_sequence, expected_u_sequence):
            if len(control_msg.u) != len(expected_u):
                return False
            for u_val, exp_val in zip(control_msg.u, expected_u):
                if round(abs(float(u_val) - float(exp_val)), places) != 0:
                    return False
        return True

    def _load_expected_sequence(self):
        test_script_dir = os.path.dirname(os.path.realpath(__file__))
        expected_u_file = os.path.abspath(
            os.path.join(test_script_dir, '..', '..', 'expected_control_sequence.npy')
        )
        if not os.path.exists(expected_u_file):
            self.skipTest(f"Expected control sequence file not found: {expected_u_file}")
        return np.load(expected_u_file)
    
    {%- if constraints.has_x0 %}

    def _publish_initial_state(self, publisher: Publisher, timeout_sec: float = 2.0):
        state_msg = State()
        state_msg.x = [float(v) for v in ({{ constraints.lbx_0 | join(sep=', ') }})]

        end_time = time.time() + timeout_sec
        while time.time() < end_time and publisher.get_subscription_count() < 1:
            rclpy.spin_once(self.node, timeout_sec=0.1)

        if publisher.get_subscription_count() < 1:
            self.fail(f"Publisher on Topic '{publisher.topic_name}' "
                      f"could not find a Subscriber (Timeout).")

        publisher.publish(state_msg)
        rclpy.spin_once(self.node, timeout_sec=0.1)
    {%- endif %}
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
        """Verify parameter value via ROS2 GetParameters service."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail(f"Parameter service for '{{ node_fqn }}' not available.")

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
        """Set and verify parameter value via ROS2 SetParameters service."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail(f"Parameter service for '{{ node_fqn }}' not available.")

        future = self.param_client.set_parameters([(param_name, new_value)])
        rclpy.spin_until_future_complete(self.node, future, timeout_sec=3.0)
        if not future.done() or future.result() is None:
            self.fail(f"Timeout while setting parameter {param_name}")

        results = future.result().results
        if not results or not results[0].successful:
            self.fail(f"Failed to set parameter {param_name}: {results[0].reason if results else 'unknown'}")

        self._check_parameter_get(param_name, new_value)

    def _check_parameter_set_rejected(self, param_name: str, invalid_value, unchanged_expected=None):
        """Verify that an invalid parameter update is rejected without altering the parameter."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail("Parameter service for '{{ node_fqn }}' not available.")

        future = self.param_client.set_parameters([(param_name, invalid_value)])
        rclpy.spin_until_future_complete(self.node, future, timeout_sec=3.0)
        if not future.done() or future.result() is None:
            self.fail(f"Timeout while testing rejected parameter {param_name}")

        results = future.result().results
        self.assertTrue(results, "SetParameters returned empty results list.")
        self.assertFalse(
            results[0].successful,
            f"Expected parameter '{param_name}' update with {invalid_value!r} to be rejected."
        )
        self.assertTrue(results[0].reason, f"Expected rejection reason for '{param_name}'.")
        if unchanged_expected is not None:
            self._check_parameter_get(param_name, unchanged_expected)