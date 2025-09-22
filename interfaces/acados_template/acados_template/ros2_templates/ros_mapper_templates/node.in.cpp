#include "{{ package_name }}/node.h"


namespace {{ package_name }}
{
{% set ClassName = node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{{ ClassName }}::{{ ClassName }}()
    : Node("{{ node_name }}")
{
    RCLCPP_INFO(this->get_logger(), "Initializing {{ node_name | replace(from="_", to=" ") | title }}...");

    // --- Subscriptions ---
    {%- for in_msg in in_msgs %}
    {{ in_msg.topic_name }}_sub_ = this->create_subscription<{{  in_msg.msg_type }}>(
        "{{ in_msg.topic_name }}", 3,
        std::bind(&{{ ClassName }}::{{ in_msg.topic_name }}_callback, this, std::placeholders::_1));
    {%- endfor %}

    // --- Publishers ---
    {%- for out_msg in out_msgs %}
    {{ out_msg.topic_name }}_pub_ = this->create_publisher<{{ out_msg.msg_type }}>(
        "{{ out_msg.topic_name }}", 3);
    {%- endfor %}
}

{{ ClassName }}::~{{ ClassName }}() {
    RCLCPP_INFO(this->get_logger(), "Shutting down {{ node_name | replace(from="_", to=" ") | title }}.");
}


// --- ROS Callbacks ---
{%- for in_msg in in_msgs %}
void {{ ClassName }}::{{ in_msg.topic_name }}_callback(const {{ in_msg.msg_type }}::SharedPtr msg) {
    {%- set fields_to_store = in_msg.flat_field_tree | filter(attribute="needs_storage", value=true) %}
    {%- if fields_to_store | length > 0 %}
    {
        std::scoped_lock lock(data_mutex_);
        {%- for field in fields_to_store %}
        this->{{ in_msg.topic_name }}_{{ field.name | replace(from=".", to="_") }}_ = msg->{{ field.name }};
        {%- endfor %}
    }
    {%- endif %}

    // --- Publish msgs ---
    {%- for out_msg in out_msgs %}
    {%- if out_msg.exec_topic == in_msg.topic_name %}
    {
        auto out_msg_ptr = std::make_unique<{{ out_msg.msg_type }}>();

        {%- if out_msg.needs_publish_lock %}
        std::scoped_lock lock(data_mutex_);
        {%- endif %}

        {%- for map in out_msg.mapping %}
            {%- set src = map.source %}
            {%- set dest = map.dest %}

            {#- Source accessor #}
            {%- if src.topic == in_msg.topic_name %}
                {%- set source_accessor = "msg->" ~ src.field %}
            {%- else %}
                {%- set source_accessor = "this->" ~ src.topic ~ "_" ~ src.field ~ "_" | replace(from=".", to="_") %}
            {%- endif %}
            {%- if src.index or src.index == 0 %}
                {%- set source_accessor = source_accessor ~ "[" ~ src.index ~ "]" %}
            {%- endif %}

            {#- Destination accessor #}
            {%- set dest_accessor = "out_msg_ptr->" ~ dest.field %}
            {%- if dest.index or dest.index == 0 %}
                {%- set dest_accessor = dest_accessor ~ "[" ~ dest.index ~ "]" %}
            {%- endif %}
        {{ dest_accessor }} = {{ source_accessor }};
        {%- endfor %}
        {{ out_msg.topic_name }}_pub_->publish(std::move(out_msg_ptr));
    }
    {%- endif %}
    {%- endfor %}
}
{% endfor %}

} // namespace {{ package_name }}


// --- Main ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<{{ package_name }}::{{ ClassName }}>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
