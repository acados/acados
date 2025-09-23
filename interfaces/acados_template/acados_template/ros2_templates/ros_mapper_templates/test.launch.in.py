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

{%- for m in out_msgs | concat(with=in_msgs) %}
{%- set parts = m.msg_type | split(pat="::") %}
from {{ parts[0] }}.{{ parts[1] }} import {{ parts[2] }}
{%- endfor %}
{%- set ns = namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}

@pytest.mark.launch_test
def generate_test_description():
    """Generate launch description for node testing."""
    start_{{ node_name }} = Node(
        package='{{ package_name }}',
        executable='{{ node_name }}',
        name='{{  node_name }}'
    )

    return launch.LaunchDescription([
        start_{{ node_name }},
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

    def test_subscribing(self, proc_info):
        """Test if the node subscribes to all expected topics."""
        {%- for m in in_msgs %}
        self.wait_for_subscription('{{ m.topic_name }}')
        {%- endfor %}

    def test_publishing(self, proc_info):
        """Test if the node publishes to all expected topics."""
        {%- for m in out_msgs %}
        self.wait_for_publisher('{{ m.topic_name }}')
        {%- endfor %}

    def wait_for_subscription(self, topic: str, timeout: float = 1.0, threshold: float = 0.5):
        end_time = time.time() + timeout + threshold
        while time.time() < end_time:
            subs = self.node.get_subscriptions_info_by_topic(topic)
            if subs:
                return True
            time.sleep(0.05)
        self.fail(f"Node has NOT subscribed to '{topic}'.")

    def wait_for_publisher(self, topic: str, timeout: float = 1.0, threshold: float = 0.5):
        end_time = time.time() + timeout + threshold
        while time.time() < end_time:
            pubs = self.node.get_publishers_info_by_topic(topic)
            if pubs:
                return True
            time.sleep(0.05)
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
    cmd = ['ros2', 'param', 'get', '{{ node_name }}', param_name]
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

    cmd = ['ros2', 'param', 'set', '{{ node_name }}', param_name, value_str]
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    return result.stderr