HEXMAERS = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA',
            'AATATA', 'AATACA', 'AATAGA', 'AAAAAG', 'ACTAAA', 'AAGAAA',
            'AATGAA', 'TTTAAA', 'AAAACA', 'GGGGCT']


filter_style_A = 'A'
for nda in range(15, 100, 10):
    for ltc in range(0, 7, 2):
        for nbr in range(0, 7, 2):
            for maxbrtl in range(0, 7, 2):
                for hxm in [HEXMAERS[:2], HEXMAERS]:
                    hxm = ' '.join(hxm)
                    print('python post-kleat-param-sweep.py UHR/C1/tasrkleat-results/kleat/cba.KLEAT {filter_style_A} {nda} {ltc} {nbr} {maxbrtl} {hxm}'.format(**locals()))

filter_style_B = 'B'
for tbr in range(1, 10):
    print('python post-kleat-param-sweep.py UHR/C1/tasrkleat-results/kleat/cba.KLEAT {filter_style_B} {tbr}'.format(**locals()))

