# `{{ project().name }}`
BGCFlow Report for: `{{ project().name }}`

## Project Description
- {{ project().description }}
- Sample size **{{ project().sample_size }}**

## Available reports:
{% for i in project().rule_used %}
  - [*{{ i }}*](/{{ i }})
{% endfor %}

## References
{% for i in project().references %}
  - *{{ i }}*
{% endfor %}
