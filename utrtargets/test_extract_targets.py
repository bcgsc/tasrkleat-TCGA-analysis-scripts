import unittest
import extract_targets as E


def assertEqual(a, b):
    try:
        assert a == b
    except AssertionError as err:
        print('{0} != {1}'.format(a, b))
        # raise


class TestGetCoords(unittest.TestCase):
    def setUp(self):
        adf = E.load_gtf('./test_data/test.gtf')
        self.groups = adf.groupby('transcript_id')

    def test_get_coords_selective(self):
        # one 5' UTR section and one 3' UTR section
        F = E.get_coords_selective

        # one 5' UTR section and one 3' UTR section
        g = self.groups.get_group('ENST00000246163')
        self.assertEqual(g.seqname.unique()[0], 14)
        self.assertEqual(
            F(g, '-'),
            [(65550462, 65551017), (65560426, 65560533)])

        # three 5' UTR sections and one 3' UTR section
        g = self.groups.get_group('ENST00000232014')
        self.assertEqual(g.seqname.unique()[0], 3)
        self.assertEqual(
            F(g, '-'),
            [(187440186, 187440389), (187442729, 187442866)])

        # two 5' UTR sections and two 3' UTR sections
        g = self.groups.get_group('ENST00000256078')
        self.assertEqual(g.seqname.unique()[0], 12)
        self.assertEqual(
            F(g, '-'),
            [(25362365, 25362845), (25368371, 25368494), (25378548, 25378707)])

        # two 5' UTR sections and one 3' UTR section
        g = self.groups.get_group('ENST00000278616')
        self.assertEqual(g.seqname.unique()[0], 11)
        self.assertEqual(
            F(g, '+'),
            [(108235809, 108235945), (108236052, 108239826)])

        # two 5' UTR sections and one 3' UTR section
        g = self.groups.get_group('ENST00000219476')
        self.assertEqual(g.seqname.unique()[0], 16)
        self.assertEqual(
            F(g, '+'),
            [(2138228, 2138326), (2138447, 2138713)])

        # no UTR
        g = self.groups.get_group('ENST00000299252')
        self.assertEqual(g.seqname.unique()[0], 12)
        self.assertEqual(F(g, '+'), [])

        # this transcript contains a bug reported on
        # https://www.biostars.org/p/206362/#206363 incorrect transcript structure,
        # such a coincidence, one of the 0.15% or 113 transcripts with imbalanced
        # L_{exon} = L_{CDS} + L_{UTR} + L_{stop_codon}
        g = self.groups.get_group('ENST00000367180')
        self.assertEqual(g.seqname.unique()[0], 1)
        self.assertEqual(
            F(g, '+'),
            [(204501897, 204501947), (204506046, 204506169)])

        # this one has only 5' UTR
        g = self.groups.get_group('ENST00000578709')
        self.assertEqual(g.seqname.unique()[0], 17)
        self.assertEqual(F(g, '+'), [])

    def test_get_coords(self):
        # one 5' UTR section and one 3' UTR section
        F = E.get_coords
        buffer_size = 300

        # one 5' UTR section and one 3' UTR section
        g = self.groups.get_group('ENST00000246163')
        self.assertEqual(g.seqname.unique()[0], 14)
        self.assertEqual(
            F(g, buffer_size),
            [(65550462 - buffer_size, 65551017),
             (65560426, 65560533),
             (65568264, 65568290),
             (65569022, 65569227 + buffer_size)])

        # three 5' UTR sections and one 3' UTR section
        g = self.groups.get_group('ENST00000232014')
        self.assertEqual(g.seqname.unique()[0], 3)
        self.assertIn((187440186 - buffer_size, 187440389), F(g, buffer_size))

        # two 5' UTR sections and two 3' UTR sections
        g = self.groups.get_group('ENST00000256078')
        self.assertEqual(g.seqname.unique()[0], 12)
        self.assertIn((25362365 - buffer_size, 25362845), F(g, buffer_size))

        # two 5' UTR sections and one 3' UTR section
        g = self.groups.get_group('ENST00000278616')
        self.assertEqual(g.seqname.unique()[0], 11)
        self.assertIn((108236052, 108239826 + buffer_size), F(g, buffer_size))

        # no UTR
        g = self.groups.get_group('ENST00000299252')
        self.assertEqual(g.seqname.unique()[0], 12)
        self.assertIn(((69202992 - buffer_size, 69203072)), F(g, buffer_size))
        self.assertIn(((69233054, 69233629 + buffer_size)), F(g, buffer_size))
