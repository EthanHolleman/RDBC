from pathlib import Path
from handler import DEFAULT_BATCH, DEFAULT_XML, DEFAULT_OPTIONS
from handler.utils import join_protein_ligand_pdbs
import os, sys

class DockJob():

    job_counter = 0
    rosetta_exe = ''

    def __init__(self, output_dir, ligand_descriptor, protein,
                xml_template=DEFAULT_XML, options_template=DEFAULT_OPTIONS,
                batch_template=DEFAULT_BATCH, number_iterations=2000):
        self.output_dir = Path(output_dir)
        self.ligand_descriptor = ligand_descriptor
        self.xml_template = xml_template
        self.batch_template = batch_template
        self.protein = protein
        self.options_template = options_template
        self.number_iterations = number_iterations

        DockJob.job_counter += 1

        self._id = DockJob.job_counter
    
    @property
    def run_name(self):
        return self.ligand_descriptor.name
    
    @property
    def name(self):
        return '{}_{}_{}'.format(self._id, self.output_dir, self.ligand_descriptor.name)

    @property
    def output_dir(self):
        return self._output_dir

    @output_dir.setter
    def output_dir(self, new_dir):
        new_dir = Path(new_dir)
        if not Path(new_dir).is_dir():
            Path(new_dir).mkdir()
        self._output_dir = new_dir
    
    @property
    def protein_ligand_file(self):
        plf = self._docking_dir.joinpath('{}_{}.pdb'.format(
            self.protein.stem, self.ligand_descriptor.ligand_pdb.stem
        ))
        if not plf.is_file():
            join_protein_ligand_pdbs(self.protein, self.ligand_descriptor.ligand_pdb,
            plf)
        
        return plf

    @property
    def _docking_dir(self):
        dd = self.output_dir.joinpath('docking')
        if not dd.is_dir():
            dd.mkdir()
        return dd
        
    @property
    def control_dir(self):
        cd = self.output_dir.joinpath('control')
        if not cd.is_dir():
            cd.mkdir()
        return cd
    
    @property
    def results_dir(self):
        rd = self.output_dir.joinpath('results')
        if not rd.is_dir():
            rd.mkdir()
        return rd
    
    @property
    def options_file(self):
        of = self.control_dir.joinpath('options.txt')
        if not of.is_file():
            self._create_options_file(of)
        return of
    
    @property
    def xml_script_file(self):
        xsf = self.control_dir.joinpath('dock.xml')
        if not xsf.is_file():
            with open(str(xsf), 'w') as handle:
                handle.write(open(self.xml_template).read())
        return xsf
    
    @property
    def sbatch_file(self):
        return self.control_dir.joinpath('job.sbatch')
    
    def set_up_for_submit(self):
        batch_template_string = open(self.batch_template).read()
        batch_template_string = batch_template_string.format(
           self.run_name, self.control_dir, self.control_dir, self._rosetta_cmd()
        )
        with open(str(self.sbatch_file), 'w') as handle:
            handle.write(batch_template_string)

    def submit(self):
        '''Submit this job to the job handler (SLURM).

        Returns:
            [type]: [description]
        '''
        cmd = 'sbatch {}'.format(self.sbatch_file)
        os.system(cmd)
    
    def _create_options_file(self, path):
        options_template_string = open(str(self.options_template)).read()
        options_template_string = options_template_string.format(
            self.protein_ligand_file,
            self.ligand_descriptor.ligand_params,
            self.xml_script_file,
            self.results_dir,
            self.number_iterations
        )
        with open(str(path), 'w') as handle:
            handle.write(options_template_string)
        
        return path

    def _rosetta_cmd(self):
        return '{} @{}'.format(DockJob.rosetta_exe, self.options_file)