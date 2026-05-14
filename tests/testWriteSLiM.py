import unittest
import os
import tempfile
import shutil
import numpy as np
from utils.writeSLiM import writeSLiM
from utils.writeSLiMHPC import writeSLiMHPC


# ── shared helpers ──────────────────────────────────────────────────────────

def _make_start_params(base_path, slurm_dir=None, nonWF=False, neutral=False,
                       count_subs=False, calculate_selection=False,
                       backup=False, output_gens=False, polymorphisms=False):
    return {
        "filenames": [base_path, base_path, None, slurm_dir, None, None, None, base_path, base_path],
        "nonWF": nonWF,
        "ancestral_sequence": ["14", "20", "30"],
        "neutral_evolution": neutral,
        "fitness_profiles": {
            "A": [1.0], "C": [0.9], "D": [0.8], "E": [0.7], "F": [0.6],
            "G": [0.5], "H": [0.4], "I": [0.3], "K": [0.2], "L": [0.1],
            "M": [0.9], "N": [0.8], "P": [0.7], "Q": [0.6], "R": [0.5],
            "S": [0.4], "T": [0.3], "V": [0.2], "W": [0.1], "Y": [0.0],
            "X": [0.0], "neutral": [1.0],
        },
        "fitness_profile_nums": [0, 1, 0],
        "coding_seqs": np.array([[0, 2]]),
        "genome_length": 3,
        "count_subs": count_subs,
        "calculate_selection": calculate_selection,
        "backup": backup,
        "output_gens": output_gens,
        "polymorphisms": polymorphisms,
        "dn_denom": 1.5,
        "ds_denom": 3.0,
    }


def _make_root_pop(**overrides):
    params = {
        "pop_name": "p1",
        "parent_pop_name": None,
        "population_size": 100,
        "mutation_rate": 2.5e-6,
        "recombination_rate": 2.5e-8,
        "jukes_cantor": True,
        "mutation_matrix": (None, "matrix(c(0,1,2,3,1,0,4,5,2,4,0,6,3,5,6,0), ncol=4, byrow=T)"),
        "dist_from_start": 0,
        "end_dist": 100,
        "terminal_clade": True,
        "last_child_clade": False,
        "sample_size": "all",
        "count_subs": False,
        "calculate_selection": False,
        "backup": False,
        "output_gens": False,
        "polymorphisms": False,
        "clade_name": "population_A",
        "split_ratio": 0.5,
        "scaling_value": 100.0,
        "fitness_profile_nums": [0, 1, 0],
        "time": "24:00:00",
        "partition": "normal",
        "memory": None,
    }
    params.update(overrides)
    return params


def _make_child_pop(**overrides):
    params = {
        "pop_name": "p2",
        "parent_pop_name": "p1",
        "population_size": 50,
        "mutation_rate": 2.5e-6,
        "recombination_rate": 2.5e-8,
        "jukes_cantor": True,
        "mutation_matrix": (None, "matrix(c(0,1,2,3,1,0,4,5,2,4,0,6,3,5,6,0), ncol=4, byrow=T)"),
        "dist_from_start": 100,
        "end_dist": 200,
        "terminal_clade": False,
        "last_child_clade": False,
        "sample_size": "all",
        "count_subs": False,
        "calculate_selection": False,
        "backup": False,
        "output_gens": False,
        "polymorphisms": False,
        "clade_name": "population_B",
        "split_ratio": 0.5,
        "scaling_value": 100.0,
        "fitness_profile_nums": [0, 1, 0],
        "time": "24:00:00",
        "partition": "normal",
        "memory": None,
    }
    params.update(overrides)
    return params


# ── writeSLiM tests ─────────────────────────────────────────────────────────

class testWriteSLiM(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.base_path = os.path.join(self.temp_dir, "output")
        self.sp = _make_start_params(self.base_path)
        self.root = _make_root_pop()
        self.child = _make_child_pop()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        if os.path.exists("fitness_profile_nums.txt"):
            os.remove("fitness_profile_nums.txt")

    def _read(self, pop_name):
        with open(self.base_path + "_" + pop_name + ".slim") as f:
            return f.read()

    def _writer_with_open_file(self, pop_name, start_params=None):
        w = writeSLiM(start_params or self.sp)
        w.output_file = open(self.base_path + "_" + pop_name + ".slim", "w")
        return w

    # ── write_initialize ──

    def test_write_initialize_root_jukes_cantor(self):
        w = self._writer_with_open_file("p1")
        w.write_initialize(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("initialize()", c)
        self.assertIn("codonsToNucleotides(c(14,20,30)", c)
        self.assertIn("mmJukesCantor(", c)
        self.assertIn("initializeRecombinationRate(2.5e-08)", c)
        self.assertIn("initializeGenomicElement(g1, 0, 8)", c)  # 0*3=0, (2+1)*3-1=8

    def test_write_initialize_root_mutation_matrix(self):
        params = _make_root_pop(jukes_cantor=False)
        w = self._writer_with_open_file("p1")
        w.write_initialize(params)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("matrix(c(", c)
        self.assertNotIn("mmJukesCantor", c)

    def test_write_initialize_child_loads_parent_fasta(self):
        w = self._writer_with_open_file("p2")
        w.write_initialize(self.child)
        w.output_file.close()
        self.assertIn('initializeAncestralNucleotides("p1.fasta")', self._read("p2"))

    def test_write_initialize_nonwf_adds_model_type(self):
        sp = _make_start_params(self.base_path, nonWF=True)
        w = self._writer_with_open_file("p1", sp)
        w.write_initialize(self.root)
        w.output_file.close()
        self.assertIn('initializeSLiMModelType("nonWF")', self._read("p1"))

    def test_write_initialize_multiple_coding_regions(self):
        sp = dict(self.sp)
        sp["coding_seqs"] = np.array([[0, 1], [3, 4]])
        sp["genome_length"] = 6
        w = self._writer_with_open_file("p1", sp)
        w.write_initialize(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("initializeGenomicElement(g1, 0, 5)", c)   # (0*3, (1+1)*3-1)
        self.assertIn("initializeGenomicElement(g2,", c)
        self.assertIn("initializeGenomicElement(g1, 9, 14)", c)  # (3*3, (4+1)*3-1)

    # ── write_fitness ──

    def test_write_fitness_fitness_based_has_all_functions(self):
        w = self._writer_with_open_file("p1")
        w.write_fitness(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("function (void) setup_fitness(void)", c)
        self.assertIn("function (void) get_fitness (void)", c)
        self.assertIn("function (float) get_genome_fitness (object nucs)", c)
        self.assertIn("fitnessEffect()", c)

    def test_write_fitness_scaling_value_in_denominator(self):
        w = self._writer_with_open_file("p1")
        w.write_fitness(self.root)
        w.output_file.close()
        # denominator = 2 * scaling_value = 2 * 100.0 = 200.0
        self.assertIn("200.0", self._read("p1"))

    def test_write_fitness_neutral_omits_fitness_functions(self):
        sp = _make_start_params(self.base_path, neutral=True)
        w = self._writer_with_open_file("p1", sp)
        w.write_fitness(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("function (void) setup_fitness(void)", c)
        self.assertNotIn("get_fitness", c)
        self.assertNotIn("fitnessEffect()", c)

    def test_write_fitness_count_subs_adds_tracking_vars(self):
        sp = _make_start_params(self.base_path, count_subs=True)
        w = self._writer_with_open_file("p1", sp)
        w.write_fitness(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("fixations_counted_p1", c)
        self.assertIn("subs_p1", c)

    def test_write_fitness_calculate_selection_adds_dn_ds_vars(self):
        sp = _make_start_params(self.base_path, count_subs=True, calculate_selection=True)
        w = self._writer_with_open_file("p1", sp)
        w.write_fitness(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("dN_p1", c)
        self.assertIn("dS_p1", c)

    # ── write_subpop (WF) ──

    def test_write_subpop_root_contains_all_sections(self):
        w = writeSLiM(self.sp)
        w.write_subpop(self.root)
        w.close_file()
        c = self._read("p1")
        self.assertIn("initialize()", c)
        self.assertIn("function (void) setup_fitness(void)", c)
        self.assertIn('sim.addSubpop("p1", 100)', c)
        self.assertIn("1 early()", c)

    def test_write_subpop_child_uses_addsubpopsplit(self):
        w = writeSLiM(self.sp)
        w.write_subpop(self.root)
        w.write_subpop(self.child)
        w.close_file()
        c = self._read("p1")
        self.assertIn("sim.addSubpopSplit", c)
        self.assertIn('"p2"', c)

    def test_write_subpop_last_child_removes_parent_pop(self):
        w = writeSLiM(self.sp)
        w.write_subpop(self.root)
        w.write_subpop(_make_child_pop(last_child_clade=True))
        w.close_file()
        self.assertIn("p1.setSubpopulationSize(0)", self._read("p1"))

    # ── write_subpop_nonwf ──

    def test_write_subpop_nonwf_root_has_fitness_scaling_and_reproduction(self):
        sp = _make_start_params(self.base_path, nonWF=True)
        w = writeSLiM(sp)
        w.write_subpop_nonwf(self.root)
        w.close_file()
        c = self._read("p1")
        self.assertIn("fitnessScaling", c)
        self.assertIn("reproduction()", c)

    def test_write_subpop_nonwf_last_child_uses_takemigrants(self):
        sp = _make_start_params(self.base_path, nonWF=True)
        w = writeSLiM(sp)
        w.write_subpop_nonwf(self.root)
        w.write_subpop_nonwf(_make_child_pop(last_child_clade=True, terminal_clade=True))
        w.close_file()
        self.assertIn("takeMigrants", self._read("p1"))

    def test_write_subpop_nonwf_non_last_child_samples_migrants(self):
        sp = _make_start_params(self.base_path, nonWF=True)
        w = writeSLiM(sp)
        w.write_subpop_nonwf(self.root)
        w.write_subpop_nonwf(_make_child_pop(last_child_clade=False, terminal_clade=True))
        w.close_file()
        self.assertIn("migrants = sample(", self._read("p1"))

    # ── write_repeated_commands ──

    def test_write_repeated_commands_correct_time_range(self):
        w = self._writer_with_open_file("p1")
        w.write_repeated_commands(self.root)
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("1:100", c)
        self.assertIn("late ()", c)

    def test_write_repeated_commands_output_gens_adds_cycle_check(self):
        w = self._writer_with_open_file("p1")
        w.write_repeated_commands(_make_root_pop(output_gens=True))
        w.output_file.close()
        self.assertIn("sim.cycle%100 == 0", self._read("p1"))

    def test_write_repeated_commands_backup_writes_output_full(self):
        sp = dict(self.sp)
        sp["filenames"] = list(self.sp["filenames"])
        sp["filenames"][2] = self.temp_dir
        w = self._writer_with_open_file("p1", sp)
        w.write_repeated_commands(_make_root_pop(backup=True))
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("outputFull", c)
        self.assertIn("writeFile", c)

    def test_write_repeated_commands_count_subs_tracks_fixations(self):
        w = self._writer_with_open_file("p1")
        w.write_repeated_commands(_make_root_pop(count_subs=True))
        w.output_file.close()
        self.assertIn("fixations_counted_p1", self._read("p1"))

    def test_write_repeated_commands_calculate_selection_tracks_dn_ds(self):
        w = self._writer_with_open_file("p1")
        w.write_repeated_commands(_make_root_pop(count_subs=True, calculate_selection=True))
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("dN_p1", c)
        self.assertIn("dS_p1", c)

    # ── write_end_pop ──

    def test_write_end_pop_terminal_outputs_fasta(self):
        w = self._writer_with_open_file("p1")
        w.write_end_pop(_make_root_pop(terminal_clade=True, sample_size="all"))
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("100 late()", c)
        self.assertIn("_nuc.fasta", c)
        self.assertIn("_aa.fasta", c)

    def test_write_end_pop_nonterminal_no_file_output(self):
        w = self._writer_with_open_file("p2")
        w.write_end_pop(_make_child_pop(terminal_clade=False))
        w.output_file.close()
        c = self._read("p2")
        self.assertIn("200 late()", c)
        self.assertNotIn("writeFile", c)

    def test_write_end_pop_count_subs_writes_subs_file(self):
        sp = dict(self.sp)
        sp["filenames"] = list(self.sp["filenames"])
        sp["filenames"][5] = self.temp_dir
        w = self._writer_with_open_file("p1", sp)
        w.write_end_pop(_make_root_pop(count_subs=True))
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("fixations_counted_p1", c)
        self.assertIn("_num_subs.txt", c)

    def test_write_end_pop_polymorphisms_writes_poly_file(self):
        sp = dict(self.sp)
        sp["filenames"] = list(self.sp["filenames"])
        sp["filenames"][6] = self.temp_dir
        w = self._writer_with_open_file("p1", sp)
        w.write_end_pop(_make_root_pop(polymorphisms=True))
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("_polymorphisms.txt", c)
        self.assertIn("_fixed_sites.txt", c)

    # ── write_terminal_output ──

    def test_write_terminal_output_all_outputs_all_genomes(self):
        w = writeSLiM(self.sp)
        result = w.write_terminal_output(self.root)
        self.assertIn("genomes = p1.genomes", result)
        self.assertIn("_nuc.fasta", result)
        self.assertIn("_aa.fasta", result)

    def test_write_terminal_output_consensus_uses_whichmax(self):
        w = writeSLiM(self.sp)
        result = w.write_terminal_output(_make_root_pop(sample_size="consensus"))
        self.assertIn("consensus", result)
        self.assertIn("whichMax", result)
        self.assertNotIn("genomes = p1.genomes", result)

    def test_write_terminal_output_n_samples_uses_sample(self):
        w = writeSLiM(self.sp)
        result = w.write_terminal_output(_make_root_pop(sample_size="10"))
        self.assertIn("sample(", result)
        self.assertIn("10", result)

    # ── write_reproduction ──

    def test_write_reproduction_contains_addcrossed(self):
        w = self._writer_with_open_file("p1")
        w.write_reproduction()
        w.output_file.close()
        c = self._read("p1")
        self.assertIn("reproduction()", c)
        self.assertIn("addCrossed", c)


# ── writeSLiMHPC tests ──────────────────────────────────────────────────────

class testWriteSLiMHPC(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.base_path = os.path.join(self.temp_dir, "output")
        self.slurm_dir = os.path.join(self.temp_dir, "slurm")
        os.makedirs(self.slurm_dir)
        self.sp = _make_start_params(self.base_path, slurm_dir=self.slurm_dir)
        self.root = _make_root_pop()
        self.child = _make_child_pop()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)
        if os.path.exists("fitness_profile_nums.txt"):
            os.remove("fitness_profile_nums.txt")

    def _read_slim(self, pop_name):
        with open(self.base_path + "_" + pop_name + ".slim") as f:
            return f.read()

    def _read_sh(self, pop_name):
        with open(self.base_path + "_" + pop_name + ".sh") as f:
            return f.read()

    # ── create_scripts ──

    def test_create_scripts_sh_has_sbatch_headers(self):
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.root)
        w.output_file.close()
        sh = self._read_sh("p1")
        self.assertIn("#!/bin/sh", sh)
        self.assertIn("#SBATCH", sh)
        self.assertIn("-p normal", sh)
        self.assertIn("-t 24:00:00", sh)
        self.assertIn("slim", sh)

    def test_create_scripts_sh_references_slim_file(self):
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.root)
        w.output_file.close()
        self.assertIn("_p1.slim", self._read_sh("p1"))

    def test_create_scripts_opens_slim_file_for_writing(self):
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.root)
        self.assertFalse(w.output_file.closed)
        w.output_file.close()
        self.assertTrue(os.path.exists(self.base_path + "_p1.slim"))

    # ── write_start_pop ──

    def test_write_start_pop_root_adds_subpop_and_fixations(self):
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.root)
        w.write_start_pop(self.root)
        w.output_file.close()
        c = self._read_slim("p1")
        self.assertIn('sim.addSubpop("p1", 100)', c)
        self.assertIn("fixations", c)
        self.assertIn("setup_fitness", c)

    def test_write_start_pop_child_reads_parent_population_file(self):
        open(self.base_path + "_p1.slim", "w").close()  # parent must exist
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.child)
        w.write_start_pop(self.child)
        w.output_file.close()
        c = self._read_slim("p2")
        self.assertIn("readFromPopulationFile", c)
        self.assertIn("p1.txt", c)

    def test_write_start_pop_child_appends_sbatch_to_parent(self):
        open(self.base_path + "_p1.slim", "w").close()
        w = writeSLiMHPC(self.sp)
        w.create_scripts(self.child)
        w.write_start_pop(self.child)
        w.output_file.close()
        parent_content = self._read_slim("p1")
        self.assertIn("sbatch", parent_content)
        self.assertIn("_p2.sh", parent_content)

    # ── write_subpop (full pipeline) ──

    def test_write_subpop_produces_complete_script(self):
        w = writeSLiMHPC(self.sp)
        w.write_subpop(self.root)
        c = self._read_slim("p1")
        self.assertIn("initialize()", c)
        self.assertIn("function (void) setup_fitness(void)", c)
        self.assertIn("late()", c)
        self.assertTrue(os.path.exists(self.base_path + "_p1.sh"))

    def test_write_subpop_nonwf_includes_reproduction_and_early(self):
        sp = _make_start_params(self.base_path, slurm_dir=self.slurm_dir, nonWF=True)
        w = writeSLiMHPC(sp)
        w.write_subpop_nonwf(self.root)
        c = self._read_slim("p1")
        self.assertIn("reproduction()", c)
        self.assertIn("fitnessScaling", c)

    # ── write_end_sim ──

    def test_write_end_sim_terminal_outputs_fasta_and_fixed_muts(self):
        w = writeSLiMHPC(self.sp)
        w.output_file = open(self.base_path + "_p1.slim", "w")
        w.write_end_sim(self.root)
        w.output_file.close()
        c = self._read_slim("p1")
        self.assertIn("_nuc.fasta", c)
        self.assertIn("sim.outputFixedMutations()", c)

    def test_write_end_sim_nonterminal_writes_population_file(self):
        params = _make_child_pop(terminal_clade=False, last_child_clade=False)
        w = writeSLiMHPC(self.sp)
        w.output_file = open(self.base_path + "_p2.slim", "w")
        w.write_end_sim(params)
        w.output_file.close()
        c = self._read_slim("p2")
        self.assertIn("outputFull", c)
        self.assertIn("p2.txt", c)
        self.assertIn("p2.fasta", c)

    def test_write_end_sim_last_child_removes_parent_files(self):
        params = _make_child_pop(terminal_clade=False, last_child_clade=True)
        w = writeSLiMHPC(self.sp)
        w.output_file = open(self.base_path + "_p2.slim", "w")
        w.write_end_sim(params)
        w.output_file.close()
        c = self._read_slim("p2")
        self.assertIn("rm", c)
        self.assertIn("p1.txt", c)
        self.assertIn("p1.fasta", c)

    # ── -M / --memory flag ──

    def test_create_scripts_sh_no_mem_line_when_memory_not_set(self):
        # memory=None (default) → no --mem directive, script identical to before
        w = writeSLiMHPC(self.sp)
        w.create_scripts(_make_root_pop(memory=None))
        w.output_file.close()
        sh = self._read_sh("p1")
        self.assertNotIn("--mem", sh)

    def test_create_scripts_sh_includes_mem_when_memory_set(self):
        # memory="16g" → #SBATCH --mem=16g present in script
        w = writeSLiMHPC(self.sp)
        w.create_scripts(_make_root_pop(memory="16g"))
        w.output_file.close()
        sh = self._read_sh("p1")
        self.assertIn("#SBATCH --mem=16g", sh)


if __name__ == "__main__":
    unittest.main()
