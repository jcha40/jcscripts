import numpy as np
import re
import xml.dom.minidom as dom
import argparse

class DNASymbol:
    path = None
    color = None
    max_bits = 2
    DNA_alphabet = ()

    @classmethod
    def get_symbol(cls, i):
        return cls.DNA_alphabet[i]

class DNA_A(DNASymbol):
    path = 'M 0 100 L 33 0 L 66 0 L 100 100 L 75 100 L 66 75 L 33 75 L 25 100 L 0 100 M 41 55 L 58 55 L 50 25 L 41 55'
    color = '#FF0000'

class DNA_C(DNASymbol):
    path = 'M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 10 75 28 Z'
    color = '#0000FF'

class DNA_G(DNASymbol):
    path = 'M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 100 48 L 55 48 L 55 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 5 75 28 Z'
    color = '#FFA500'

class DNA_T(DNASymbol):
    path = 'M 0 0 L 0 20 L 35 20 L 35 100 L 65 100 L 65 20 L 100 20 L 100 0 Z'
    color = '#228B22'

DNASymbol.DNA_alphabet = (DNA_A, DNA_C, DNA_G, DNA_T)

def pwm2logo(pwm, out_fn, symbol=DNASymbol, glyph_width=100, stack_height=200):
    width = len(pwm)

    document = dom.Document()
    svg = document.appendChild(document.createElement('svg'))
    svg.setAttribute('baseProfile', 'full')
    svg.setAttribute('version', '1.1')
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
    svg.setAttribute('viewBox', '0 0 {} {}'.format(width * glyph_width, stack_height))

    defs = svg.appendChild(document.createElement('defs'))
    for base in symbol.DNA_alphabet:
        path = defs.appendChild(document.createElement('path'))
        path.setAttribute('id', base.__name__)
        path.setAttribute('d', base.path)
        path.setAttribute('fill', base.color)

    for i, pwv in enumerate(pwm):
        n = np.sum(pwv)
        if n == 0:
            continue
        vec = pwv / n
        vec[vec > 0] *= np.log2(vec[vec > 0])
        bits = sum(vec, symbol.max_bits)
        heights = pwv / n * bits / symbol.max_bits * stack_height
        idx = np.argsort(pwv)

        stack = svg.appendChild(document.createElement('g'))
        stack.setAttribute('transform', 'translate({} 0)'.format(i * glyph_width))

        y_offset = 0
        for j in idx:
            if heights[j] == 0:
                continue
            base = symbol.get_symbol(j)
            y_offset += heights[j]
            glyph = stack.appendChild(document.createElement('path'))
            glyph.setAttribute('d', base.path)
            glyph.setAttribute('fill', base.color)
            glyph.setAttribute('transform', 'matrix({} 0 0 {} 0 {})'.format(glyph_width / 100., heights[j] / 100., stack_height - y_offset))

    with open(out_fn, 'w') as f:
        svg.writexml(f, addindent='    ', newl='\n')

def alinged_pwms2logo_stack(aligned_pwms, out_fn, symbol=DNASymbol, glyph_width=100, stack_height=200):
    height, width, _ = aligned_pwms.shape

    document = dom.Document()
    svg = document.appendChild(document.createElement('svg'))
    svg.setAttribute('baseProfile', 'full')
    svg.setAttribute('version', '1.1')
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
    svg.setAttribute('viewBox', '0 0 {} {}'.format(width * glyph_width, height * stack_height))

    for y, pwm in enumerate(aligned_pwms):
        row = svg.appendChild(document.createElement('g'))
        row.setAttribute('transform', 'translate(0 {})'.format(y * stack_height))
        for i, pwv in enumerate(pwm):
            n = np.sum(pwv)
            if n == 0:
                continue
            vec = pwv / n
            vec[vec > 0] *= np.log2(vec[vec > 0])
            bits = sum(vec, symbol.max_bits)
            heights = pwv / n * bits / symbol.max_bits * stack_height
            idx = np.argsort(pwv)
    
            stack = row.appendChild(document.createElement('g'))
            stack.setAttribute('transform', 'translate({} 0)'.format(i * glyph_width))
    
            y_offset = 0
            for j in idx:
                if heights[j] == 0:
                    continue
                base = symbol.get_symbol(j)
                y_offset += heights[j]
                glyph = stack.appendChild(document.createElement('path'))
                glyph.setAttribute('d', base.path)
                glyph.setAttribute('fill', base.color)
                glyph.setAttribute('transform', 'matrix({} 0 0 {} 0 {})'.format(glyph_width / 100., heights[j] / 100., stack_height - y_offset))

    with open(out_fn, 'w') as f:
        svg.writexml(f, addindent='    ', newl='\n')

if __name__ == '__main__':
    name_pattern = re.compile('MOTIF\s+(.+)')

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_fn', help='MEME file')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('--revcomp', '-rc', action='store_true', help='output reverse complement')
    args = parser.parse_args()

    with open(args.meme_fn) as f:
        for line in f:
            if line.startswith('MOTIF'):
                motif = name_pattern.match(line).group(1).strip()
                for line in f:
                    if line.startswith('letter-probability matrix'):
                        break
                pwm = []
                for line in f:
                    if not line.strip():
                        break
                    pwm.append(list(map(float, line.split())))
                pwm = np.flip(np.array(pwm)) if args.revcomp else np.array(pwm)
                pwm2logo(pwm, '{}/{}.svg'.format(args.out_dir, motif))
