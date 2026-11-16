from typing import Union
import time
import rclpy
import unittest
import launch
import launch_testing
import pytest
from launch_ros.actions import Node
from rclpy.parameter import Parameter, parameter_value_to_python
from rcl_interfaces.srv import GetParameters, SetParameters

from {{ ros_opts.package_name }}_interface.msg import State, ControlInput
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- set raw_control_topic = ros_opts.control_topic | trim(chars='/') %}
{%- set raw_state_topic = ros_opts.state_topic | trim(chars='/') %}
{%- if ns %}
{%- set topic_prefix = "/" ~ ns ~ "/" %}
{%- set node_fqn = "/" ~ ns ~ "/" ~ ros_opts.node_name %}
{%- else %}
{%- set topic_prefix = "/" %}
{%- set node_fqn = ros_opts.node_name %}
{%- endif %}
{%- set control_topic = topic_prefix ~ raw_control_topic %}
{%- set state_topic = topic_prefix ~ raw_state_topic %}


class AsyncParametersClient:
    """Simple async parameter client."""
    def __init__(self, node, remote_node_name):
        self.node = node
        node_prefix = remote_node_name if remote_node_name.startswith('/') else f'/{remote_node_name}'
        self._get_client = self.node.create_client(
            GetParameters,
            f'{node_prefix}/get_parameters')
        self._set_client = self.node.create_client(
            SetParameters,
            f'{node_prefix}/set_parameters')

    def wait_for_service(self, timeout_sec=2.0):
        ok1 = self._get_client.wait_for_service(timeout_sec=timeout_sec)
        ok2 = self._set_client.wait_for_service(timeout_sec=timeout_sec)
        return ok1 and ok2

    def get_parameters(self, names):
        request = GetParameters.Request()
        request.names = names
        return self._get_client.call_async(request)

    def set_parameters(self, params: list[tuple[str, object]]):
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
        {%- if ns %}
        namespace='{{ ns }}',
        {%- endif %}
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
        self.param_client = AsyncParametersClient(self.node, '{{ node_fqn }}')

    def tearDown(self):
        self.node.destroy_node()

    def test_parameters_set(self, proc_info):
        """Test valid, boundary, and rejected parameter updates."""
        default_tsim = float({{ solver_options.Tsim }})
        self.__check_parameter_set("solver_options.Tsim", default_tsim)
        self.__check_parameter_set("solver_options.Tsim", default_tsim * 2.0)
        self.__check_parameter_set("solver_options.Tsim", default_tsim)

        for bad_tsim in (0.0, -0.1, float('nan'), float('inf')):
            self.__check_parameter_set_rejected("solver_options.Tsim", bad_tsim, unchanged_expected=default_tsim)

        self.__check_parameter_set_rejected("nonexistent_parameter", 1.0)

    def test_subscribing(self, proc_info):
        """Test if the node subscribes to all expected topics."""
        self.wait_for_subscription('{{ control_topic }}')

    def test_publishing(self, proc_info):
        """Test if the node publishes to all expected topics."""
        self.wait_for_publisher('{{ state_topic }}')

    def wait_for_subscription(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            rclpy.spin_once(self.node, timeout_sec=0.1)
            subs = self.node.get_subscriptions_info_by_topic(topic)
            if subs:
                return True
        self.fail(f"Node has NOT subscribed to '{topic}'.")

    def wait_for_publisher(self, topic: str, timeout: float = 2.0):
        end_time = time.time() + timeout
        while time.time() < end_time:
            rclpy.spin_once(self.node, timeout_sec=0.1)
            pubs = self.node.get_publishers_info_by_topic(topic)
            if pubs:
                return True
        self.fail(f"Node has NOT published to '{topic}'.")

    def __check_parameter_get(self, param_name: str, expected_value: Union[list[float], float]):
        """Verify parameter value via ROS2 GetParameters service."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail("Parameter service for '{{ node_fqn }}' not available.")
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

    def __check_parameter_set(self, param_name: str, new_value: Union[list[float], float]):
        """Set and verify parameter value via ROS2 SetParameters service."""
        if not self.param_client.wait_for_service(timeout_sec=2.0):
            self.fail("Parameter service for '{{ node_fqn }}' not available.")
        future = self.param_client.set_parameters([(param_name, new_value)])
        rclpy.spin_until_future_complete(self.node, future, timeout_sec=3.0)
        if not future.done() or future.result() is None:
            self.fail(f"Timeout while setting parameter {param_name}")
        results = future.result().results
        if not results or not results[0].successful:
            self.fail(f"Failed to set parameter {param_name}: {results[0].reason if results else 'unknown'}")
        self.__check_parameter_get(param_name, new_value)

    def __check_parameter_set_rejected(self, param_name: str, invalid_value, unchanged_expected=None):
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
            self.__check_parameter_get(param_name, unchanged_expected)
