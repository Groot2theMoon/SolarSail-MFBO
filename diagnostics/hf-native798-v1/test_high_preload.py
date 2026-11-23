"""Fabricated byte-level preload edits; no Abaqus or private inputs."""

import unittest

from make_high_preload import derive, verify_delta


CABLE_ROWS = (
    ('CL', 'CL', '171000.'), ('CR', 'CR', '171000.'),
    ('Left', 'BotL', '1.71e+07'), ('Right', 'BotR', '1.71e+07'), ('Top', 'Top', '1.71e+07'),
)


def fixture():
    lines = ['*Heading', '** comment: 171000. remains unchanged', '*Part, name=MEMBRANE',
             '*Node', '1, 0., 0., 1.7e-16', '*End Part']
    for _, part, _ in CABLE_ROWS:
        lines.extend((f'*Part, name=Cable_{part}', '*Elset, elset=All_Elem', '1,', '*End Part'))
    lines.append('*Assembly, name=Assembly')
    for instance, part, _ in CABLE_ROWS:
        lines.extend((f'*Instance, name=Inst_cable_{instance}, part=Cable_{part}', '*End Instance'))
    lines += ['*End Assembly', '*Material, name=m', '*Elastic', '171000., 0.3',
              '*Initial Conditions, type=STRESS', 'MEMBRANE-1.All_Elem, 500., 500., 0., ']
    for instance, _, stress in CABLE_ROWS:
        lines.extend(('*Initial Conditions, type=STRESS', f'Inst_cable_{instance}.All_Elem, {stress}, 0., 0., '))
    lines += ['*Step, name=Step-GlobalTension', '*Static, allsdtol=0', '0.001, 1., 1e-08, 0.005', '*End Step', '']
    return '\r\n'.join(lines).encode('ascii')


class HighPreloadTests(unittest.TestCase):
    def test_two_literal_stress_edits_preserve_unrelated_matches_and_crlf(self):
        a = fixture()
        b = derive(a)
        expected = (a.replace(b'Inst_cable_CL.All_Elem, 171000.', b'Inst_cable_CL.All_Elem, 8550000.')
                     .replace(b'Inst_cable_CR.All_Elem, 171000.', b'Inst_cable_CR.All_Elem, 8550000.'))
        self.assertEqual(b, expected)
        self.assertEqual(set(verify_delta(a, b)), {'Inst_cable_CL', 'Inst_cable_CR'})

    def test_explicit_initial_stresses_below_and_above_original(self):
        a = fixture()
        for val in (123000, 342000, 855000):
            with self.subTest(clamp_pa=val):
                b = derive(a, clamp_pa=val)
                expected = a
                for target in (b'CL', b'CR'):
                    prefix = b'Inst_cable_' + target + b'.All_Elem, '
                    expected = expected.replace(prefix + b'171000.,', prefix + str(val).encode() + b'.,')
                self.assertEqual(b, expected)
                self.assertEqual(set(verify_delta(a, b, clamp_pa=val)), {'Inst_cable_CL', 'Inst_cable_CR'})
                with self.assertRaisesRegex(ValueError, 'unapproved byte change'):
                    verify_delta(a, b, clamp_pa=val+1)

    def test_explicit_stress_rejects_non_values_and_baseline_no_op(self):
        a = fixture()
        for val in (0, -1, 171000, 2.5, True):
            with self.subTest(clamp_pa=val), self.assertRaisesRegex(ValueError, 'positive changed integer'):
                derive(a, clamp_pa=val)

    def test_wrong_target_or_instance_binding_is_rejected(self):
        a = fixture()
        for before, after in ((b'CL.All_Elem', b'CL.WrongSet'), (b'part=Cable_CL', b'part=Cable_Top'),
                              (b'CR.All_Elem, 171000.', b'CR.All_Elem, 8550000.'),
                              (b'Top.All_Elem, 1.71e+07', b'Top.All_Elem, 171000.')):
            with self.subTest(before=before), self.assertRaisesRegex(ValueError, 'target|binding|value'):
                derive(a.replace(before, after))

    def test_missing_duplicate_extra_assignments_are_rejected(self):
        a = fixture()
        row = b'Inst_cable_CL.All_Elem, 171000., 0., 0., \r\n'
        for replacement in (b'', row + row, row + row.replace(b'_CL.', b'_Other.')):
            with self.subTest(row=replacement), self.assertRaisesRegex(ValueError, 'missing|duplicate|target'):
                derive(a.replace(row, replacement))

    def test_independent_check_rejects_missing_edit_or_any_unintended_diff(self):
        a = fixture()
        b = derive(a)
        for before, after in ((b'1.7e-16', b'0.'), (b'171000., 0.3', b'8550000., 0.3'),
                              (b'allsdtol=0', b'allsdtol=0.05'), (b'Top.All_Elem, 1.71e+07', b'Top.All_Elem, 8550000.'),
                              (b'CR.All_Elem, 8550000.', b'CR.All_Elem, 171000.'), (b'CR.All_Elem, 8550000., 0.', b'CR.All_Elem, 8550000., 1.')):
            with self.subTest(before=before), self.assertRaisesRegex(ValueError, 'unapproved|only two'):
                verify_delta(a, b.replace(before, after))
