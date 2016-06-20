from parse_sample_sheet import ParseSampleSheet
import parse_vcfs
import unittest


class UnitTests(unittest.TestCase):
    def test_parse_sample_sheet(self):
        self.parse = ParseSampleSheet('ExampleSampleSheet.csv')
        self.run_dict, self.sample_dict = self.parse.parse_sample_sheet()

        self.assertIsNotNone(self.run_dict, msg="Run details are missing.")
        self.assertIsNotNone(self.sample_dict, msg="Sample details are missing.")
        self.assertEqual(len(self.run_dict), 12, msg="Unexpected length for run dictionary: got %s expected 12."
                                                     % len(self.run_dict))
        self.assertEqual(self.run_dict.get("Assay"), "Nextera",
                         msg="Unexpected value for \'Assay\': got %s expected \'Nextera\'."
                             % self.run_dict.get("Assay"))
        self.assertEqual(len(self.sample_dict), 5, msg="Unexpected length for sample dictionary: got %s expected 5."
                                                       % len(self.sample_dict))

    def test_get_delly_output(self):
        self.delly_dict = parse_vcfs.get_delly_output('ExampleDelly.vcf')
        self.assertIsNotNone(self.delly_dict, msg="Delly details are missing.")
        self.assertIsInstance(self.delly_dict, dict,
                              msg="Unexpected type for delly dictionary: got %s expect \'dict\'."
                                  % type(self.delly_dict))
        self.assertEqual(len(self.delly_dict), 2, msg="Unexpected length for delly dictionary: got %s expected 2."
                                                      % len(self.delly_dict))
        self.assertEqual(self.delly_dict['chr2_4561865_N_<DEL>']['SVTYPE'], 'DEL',
                         msg="Unexpected value for \'SVTYPE\': got %s expected \'DEL\'."
                             % self.delly_dict['chr2_4561865_N_<DEL>']['SVTYPE'])

    def test_get_pindel_output(self):
        self.pindel_dict = parse_vcfs.get_pindel_output('ExamplePindel.vcf')
        self.assertIsNotNone(self.pindel_dict, msg="Pindel details are missing.")
        self.assertIsInstance(self.pindel_dict, dict,
                              msg="Unexpected type for pindel dictionary: got %s expect \'dict\'."
                                  % type(self.pindel_dict))
        self.assertEqual(len(self.pindel_dict), 56, msg="Unexpected length for pindel dictionary: got %s expected 56."
                                                        % len(self.pindel_dict))
        self.assertEqual(self.pindel_dict['chr4_153247486_CTTTTTT_C']['END'], 153247492,
                         msg="Unexpected value for \'END\': got %s expected 153247492."
                             % self.pindel_dict['chr4_153247486_CTTTTTT_C']['END'])

    def test_get_vs2_output(self):
        self.vs2_dict = parse_vcfs.get_vs2_output('ExampleVarScan2.vcf')
        self.assertIsNotNone(self.vs2_dict, msg="VarScan2 details are missing.")
        self.assertIsInstance(self.vs2_dict, dict,
                              msg="Unexpected type for VarScan2 dictionary: got %s expect \'dict\'."
                                  % type(self.vs2_dict))
        self.assertEqual(len(self.vs2_dict), 143, msg="Unexpected length for VarScan2 dictionary: got %s expected 143."
                                                      % len(self.vs2_dict))
        self.assertEqual(self.vs2_dict['chrX_100056736_T_TA']['GT'], '0/1',
                         msg="Unexpected value for \'GT\': got %s expected \'0/1\'."
                             % self.vs2_dict['chrX_100056736_T_TA']['GT'])

    def test_print_vcf(self):
        self.variant_in_vcf = parse_vcfs.print_vcf(self.delly_dict, 'chr2_4561865_N_<DEL>')
        self.assertIsInstance(self.variant_in_vcf, str,
                              msg="Unexpected output for \'print_vcf\': got %s expected \'str\'."
                                  % type(self.variant_in_vcf))

