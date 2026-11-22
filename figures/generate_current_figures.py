"""Reproduce the adopted manuscript figures with their stable export names."""
from pathlib import Path
import argparse

import publication_figure_style as style

ROOT = Path(__file__).resolve().parent
STEMS = {1: 'figure1_depth_state_ghg_framework', 2: 'figure2_depth_aware_validation',
         3: 'figure3_depth_state_framework'}
_finish = style.Sheet.finish


def _current_finish(self, unused_stem, scientific_checks):
    stem = STEMS[self.number]
    result = _finish(self, stem, scientific_checks)
    (ROOT/f'{stem}_qa.json').replace(ROOT/f'figure{self.number}_layout_checks.json')
    return result


style.OUT = ROOT
style.Sheet.finish = _current_finish


def run(number=None):
    for item in (1, 2, 3) if number is None else (number,):
        getattr(style, f'figure{item}')()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure', type=int, choices=(1, 2, 3))
    run(parser.parse_args().figure)
