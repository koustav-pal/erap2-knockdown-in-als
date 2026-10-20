include: "workflow/Snakefile"

rule manuscript:
    default_target: True
    input:
        rules.all.input
