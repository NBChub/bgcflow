# `{{ project().name }}`
Summary report for project `{{ project().name }}`. Generated using [**`BGCFlow v{{ project().bgcflow_version}}`**](https://github.com/NBChub/bgcflow){:target="_blank"}

## Project Description
- {{ project().description }}
- Sample size **{{ project().sample_size }}**
{% endraw %}

## Available reports
{{ rule_table }}

{% raw %}
## References
If you find BGCFlow useful, please cite:

<font size="2">

  - *Nuhamunada, M., B.O. Palsson, O. S. Mohite, and T. Weber. 2022. BGCFlow [Computer software]. GITHUB: [https://github.com/NBChub/bgcflow](https://github.com/NBChub/bgcflow)*
  
  - *Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. [F1000Res 10, 33](https://f1000research.com/articles/10-33/v1).*
  
  - *Nathan C Sheffield, Michał Stolarczyk, Vincent P Reuter, André F Rendeiro, Linking big biomedical datasets to modular analysis with Portable Encapsulated Projects, [GigaScience, Volume 10, Issue 12, December 2021, giab077](https://doi.org/10.1093/gigascience/giab077)*

</font>

Please also cite each tools used in the analysis:

<font size="2">
{% for i in project().references %}
  - *{{ i }}*
{% endfor %}
</font>