import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

warnings.filterwarnings('ignore', category=ResourceWarning)

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = ROOT.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from QUEEN.queen import QUEEN, pcr, quine
from QUEEN.qexperiment import digestion, ligation, primerdesign


class QuineExecutionIsolationTests(unittest.TestCase):
    def setUp(self):
        self._state = {
            'dna_dict': QUEEN.dna_dict.copy(),
            'queried_feature_dict': QUEEN.queried_feature_dict.copy(),
            'queried_features_dict': QUEEN.queried_features_dict.copy(),
            'queried_features_name_dict': QUEEN.queried_features_name_dict.copy(),
            '_products': QUEEN._products.copy(),
            '_processes': QUEEN._processes.copy(),
            '_namespace': QUEEN._namespace.copy(),
            '_namespaceflag': QUEEN._namespaceflag,
            '_num_history': QUEEN._num_history,
            '_qnum': QUEEN._qnum,
            '_source': QUEEN._source,
            '_project': QUEEN._project,
        }
        QUEEN.dna_dict = {}
        QUEEN.queried_feature_dict = {}
        QUEEN.queried_features_dict = {}
        QUEEN.queried_features_name_dict = {}
        QUEEN._products = {}
        QUEEN._processes = {}
        QUEEN._namespace = {}
        QUEEN._namespaceflag = 0
        QUEEN._source = None
        QUEEN._project = None

    def tearDown(self):
        QUEEN.dna_dict = self._state['dna_dict']
        QUEEN.queried_feature_dict = self._state['queried_feature_dict']
        QUEEN.queried_features_dict = self._state['queried_features_dict']
        QUEEN.queried_features_name_dict = self._state['queried_features_name_dict']
        QUEEN._products = self._state['_products']
        QUEEN._processes = self._state['_processes']
        QUEEN._namespace = self._state['_namespace']
        QUEEN._namespaceflag = self._state['_namespaceflag']
        QUEEN._num_history = self._state['_num_history']
        QUEEN._qnum = self._state['_qnum']
        QUEEN._source = self._state['_source']
        QUEEN._project = self._state['_project']

    def _build_simple_pcr_product(self):
        template = QUEEN(
            seq='ATGCGTACGTTAGCTAGCTAGGATCCGATCGTACGTAGCTAGCTAACGTTAGC',
            project='template',
        )
        fw = QUEEN(seq='ATGCGTACGTTAGCTAGCTA', ssdna=True, project='fw')
        rv = QUEEN(seq='GCTAACGTTAGCTAGCTACG', ssdna=True, project='rv')
        amplicon = pcr(template, fw, rv, product='amp')
        return template, fw, rv, amplicon

    def _run_exported_script(self, script_path, result_key, output_gbk):
        runner = "\n".join(
            [
                'import importlib.util',
                'import sys',
                'script_path, module_name, result_key, output_gbk = sys.argv[1:5]',
                'spec = importlib.util.spec_from_file_location(module_name, script_path)',
                'if spec is None or spec.loader is None:',
                "    raise RuntimeError(f'failed to load {script_path}')",
                'module = importlib.util.module_from_spec(spec)',
                'spec.loader.exec_module(module)',
                'module.QUEEN.dna_dict[result_key].outputgbk(output_gbk)',
            ]
        )
        return subprocess.run(
            [sys.executable, '-c', runner, str(script_path), script_path.stem, result_key, str(output_gbk)],
            cwd=script_path.parent,
            capture_output=True,
            text=True,
        )

    def test_execution_true_isolated_does_not_replace_parent_objects(self):
        template, fw, rv, amplicon = self._build_simple_pcr_product()
        refs = {key: QUEEN.dna_dict[key] for key in ('template', 'fw', 'rv', 'amp')}

        ok, reconstructed = quine(amplicon, execution=True, qexperiment_only=True)

        self.assertTrue(ok)
        self.assertEqual(str(reconstructed.seq), str(amplicon.seq))
        for key, ref in refs.items():
            self.assertIs(QUEEN.dna_dict[key], ref)

    def test_qexperiment_only_export_script_is_self_contained_for_pcr(self):
        _, _, _, amplicon = self._build_simple_pcr_product()
        with tempfile.TemporaryDirectory() as tmpdir:
            script_path = Path(tmpdir) / 'amp_qexp.py'
            gbk_path = Path(tmpdir) / 'amp_qexp.gbk'
            quine(amplicon, output=str(script_path), qexperiment_only=True, execution=False)
            text = script_path.read_text()

            self.assertIn("pcr(", text)
            self.assertNotIn("cropdna(", text)
            self.assertNotIn("searchsequence(", text)
            self.assertNotIn("queried_features_dict", text)

            proc = self._run_exported_script(script_path, 'amp', gbk_path)
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            self.assertTrue(gbk_path.exists())

    def test_re_partner_primerdesign_uses_digested_fragment_context(self):
        backbone_donor = QUEEN(
            record=str(WORKSPACE_ROOT / 'gbks' / 'new_107036.gbk'),
            project='backbone_donor',
        )
        insert_donor = QUEEN(
            record=str(WORKSPACE_ROOT / 'gbks' / 'new_112213.gbk'),
            project='insert_donor',
        )
        restriction_pair = ['Acc65I', 'HindIII']
        backbone_fragment = digestion(
            backbone_donor,
            *restriction_pair,
            selection='!label:mCherry',
            product='mCherry_removed_backbone_fragment',
        )
        primer_pairs = primerdesign(
            insert_donor,
            insert_donor['EGFP'],
            adapter_mode='RE',
            fw_partner=backbone_fragment,
            rv_partner=backbone_fragment,
            target_tm=62.0,
            design_num=1,
            fw_name='Primer_F1',
            rv_name='Primer_R1',
        )
        insert_amplicon = pcr(
            insert_donor,
            primer_pairs[0]['fw'],
            primer_pairs[0]['rv'],
            product='EGFP_restriction_adapter_amplicon',
        )
        insert_fragment = digestion(
            insert_amplicon,
            *restriction_pair,
            selection='max',
            product='EGFP_insert_fragment',
        )
        final_product = ligation(
            backbone_fragment,
            insert_fragment,
            follow_order=False,
            product='mCherry_to_EGFP_replacement_construct',
        )

        self.assertIn('GGTACC', str(primer_pairs[0]['fw'].seq).upper())
        self.assertIn('AAGCTT', str(primer_pairs[0]['rv'].seq).upper())
        self.assertEqual(final_product.project, 'mCherry_to_EGFP_replacement_construct')
        self.assertGreaterEqual(
            len(final_product.searchfeature(key_attribute='qualifier:label', query='^EGFP$')),
            1,
        )


if __name__ == '__main__':
    unittest.main()
