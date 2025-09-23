#ifndef {{ node_name | upper }}_H
#define {{ node_name | upper }}_H

#include <rclcpp/rclcpp.hpp>
#include <mutex>
#include <array>
#include <vector>
#include <unordered_map>

// ROS2 message includes
{%- for header_incl in header_includes %}
#include "{{ header_incl }}"
{%- endfor %}

// Package includes
#include "{{ package_name }}/utils.hpp"


namespace {{ package_name }}
{

{%- set ClassName = node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
class {{ ClassName }} : public rclcpp::Node {
private:
    // --- ROS Subscriptions ---
    {%- for in_msg in in_msgs %}
    rclcpp::Subscription<{{ in_msg.msg_type }}>::SharedPtr {{ in_msg.topic_name }}_sub_;
    {%- endfor %}

    // --- ROS Publishers ---
    {%- for out_msg in out_msgs %}
    rclcpp::Publisher<{{ out_msg.msg_type }}>::SharedPtr {{ out_msg.topic_name }}_pub_;
    {%- endfor %}


    // --- Data and States ---
    std::mutex data_mutex_;
    {%- for in_msg in in_msgs %}
        {%- set fields_to_store = in_msg.flat_field_tree | filter(attribute="needs_storage", value=true) %}
        {%- if fields_to_store | length > 0 %}
    // {{ in_msg.topic_name }}
            {%- for field in fields_to_store %}
                {%- if field.is_array and field.array_size == 0 %}
    std::vector<{{ field.cpp_type }}> {{ in_msg.topic_name }}_{{ field.name | replace(from=".", to="_") }}_;
                {%- elif field.is_array and field.array_size > 0 %}
    std::array<{{ field.cpp_type }}, {{ field.array_size }}> {{ in_msg.topic_name }}_{{ field.name | replace(from=".", to="_") }}_;
                {%- else %}
    {{ field.cpp_type }} {{ in_msg.topic_name }}_{{ field.name | replace(from=".", to="_") }}_;
                {%- endif %}
            {%- endfor %}
        {%- endif %}
    {% endfor %}

public:
    {{ ClassName }}();
    ~{{ ClassName }}();

private:
    // --- ROS Callbacks ---
    {%- for in_msg in in_msgs %}
    void {{ in_msg.topic_name }}_callback(const {{ in_msg.msg_type }}::SharedPtr msg);
    {%- endfor %}
};

} // namespace {{ package_name }}

#endif // {{ node_name | upper }}_H