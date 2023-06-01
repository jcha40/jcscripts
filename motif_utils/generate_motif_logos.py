import numpy as np
import re
import xml.dom.minidom as dom
import argparse

class DNASymbol:
    path = None
    mask = None
    color = None
    max_bits = 2
    DNA_alphabet = ()

    @classmethod
    def get_symbol(cls, i):
        return cls.DNA_alphabet[i]

class DNA_A(DNASymbol):
    path = 'M 0 100 L 33 0 L 66 0 L 100 100 L 75 100 L 66 75 L 33 75 L 25 100 Z'
    mask = 'M 41 55 L 50 25 L 58 55 Z'
    color = '#FF0000'

class DNA_C(DNASymbol):
    path = 'M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 10 75 28 Z'
    mask = None
    color = '#0000FF'

class DNA_G(DNASymbol):
    path = 'M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 100 48 L 55 48 L 55 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 5 75 28 Z'
    mask = None
    color = '#FFA500'

class DNA_T(DNASymbol):
    path = 'M 0 0 L 0 20 L 35 20 L 35 100 L 65 100 L 65 20 L 100 20 L 100 0 Z'
    mask = None
    color = '#228B22'

DNASymbol.DNA_alphabet = (DNA_A, DNA_C, DNA_G, DNA_T)

def pwm2logo(pwm, out_fn, n=1., symbol=DNASymbol, glyph_width=100, stack_height=200):
    width = len(pwm)

    document = dom.Document()
    svg = document.appendChild(document.createElement('svg'))
    svg.setAttribute('baseProfile', 'full')
    svg.setAttribute('version', '1.1')
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
    svg.setAttribute('viewBox', '0 0 {} {}'.format(width * glyph_width, stack_height))

    for i, pwv in enumerate(pwm):
        vec = pwv / n
        vec[vec > 0] *= np.log2(vec[vec > 0])
        bits = sum(vec, symbol.max_bits)
        heights = pwv / n * bits / symbol.max_bits * stack_height
        idx = np.argsort(pwv)

        stack = svg.appendChild(document.createElement('g'))
        stack.setAttribute('transform', 'translate({} 0)'.format(i * glyph_width))

        y_offset = 0
        for j in idx:
            symbol = symbol.get_symbol(j)
            glyph = stack.appendChild(document.createElement('g'))
            y_offset += heights[j]
            glyph.setAttribute('transform', 'matrix({} 0 0 {} 0 {})'.format(glyph_width / 100., heights[j] / 100., stack_height - y_offset))
            if symbol.mask:
                mask = glyph.appendChild(document.createElement('mask'))
                mask.setAttribute('id', 'mask-{}-{}'.format(i, j))
                mask.setAttribute('fill', 'white')
                background = mask.appendChild(document.createElement('rect'))
                background.setAttribute('x', '0')
                background.setAttribute('y', '0')
                background.setAttribute('width', '100')
                background.setAttribute('height', '100')
                mask_path = mask.appendChild(document.createElement('path'))
                mask_path.setAttribute('d', symbol.mask)
                mask_path.setAttribute('fill', 'black')
            path = glyph.appendChild(document.createElement('path'))
            path.setAttribute('d', symbol.path)
            path.setAttribute('fill', symbol.color)
            if symbol.mask:
                path.setAttribute('mask', 'url(#mask-{}-{})'.format(i, j))

    with open(out_fn, 'w') as f:
        svg.writexml(f, addindent='    ', newl='\n')

if __name__ == '__main__':
    name_pattern = re.compile('MOTIF\s+(.+)')

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_fn', help='MEME file')
    parser.add_argument('out_dir', help='output directory')
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
                pwm = np.array(pwm)
                pwm2logo(pwm, '{}/{}.svg'.format(args.out_dir, motif))
