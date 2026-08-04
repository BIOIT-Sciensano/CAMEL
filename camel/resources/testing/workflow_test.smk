rule all:
    input:
        TXT = 'out_rule_b.txt'

rule a:
    output:
        TXT = 'out_rule_a.txt'
    run:
        with open(output.TXT, 'w') as handle:
            handle.write('rule a done')
            handle.write('\n')

rule b:
    input:
        TXT = rules.a.output.TXT
    output:
        TXT = 'out_rule_b.txt'
    run:
        with open(input.TXT, 'r') as handle_in, open(output.TXT, 'w') as handle_out:
            handle_out.write(handle_in.read())
            handle_out.write('\n')
            handle_out.write('rule b done')
            handle_out.write('\n')
