# -*- coding: utf-8 -*-
import io
HERE = r'C:\Users\user\Desktop\Seungwon\solarsail\2026\SolarSail-MFBO\code'
for xc, src in (('0p2', 'seeds/modes_xc0p2_o2.txt'), ('0p8', 'seeds/modes_xc0p8_o2.txt')):
    with io.open(HERE + '\\' + src, encoding='utf-8', errors='replace') as f:
        header = []
        for line in f:
            if line.startswith('MODE '):
                break
            header.append(line)
    with io.open(HERE + '\\modes_t25_o2.txt', encoding='utf-8', errors='replace') as f:
        body = []
        started = False
        for line in f:
            if line.startswith('MODE '):
                started = True
            if started:
                body.append(line)
    out = HERE + '\\modes_t25_o2_x%s.txt' % xc
    with io.open(out, 'w', encoding='utf-8') as f:
        f.write(''.join(header))
        f.write(''.join(body))
    print('wrote', out, len(body))
