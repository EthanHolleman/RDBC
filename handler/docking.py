from pathlib import Path
from handler import DEFAULT_BATCH, DEFAULT_XML
from handler.utils import join_protein_ligand_pdbs
import os

class DockJob():

    job_counter = 0

    def __init__(self, output_dir, ligand_descriptor, 
                xml_template=DEFAULT_XML, options_template=batch_template,
                protein=None
                ):
        self.output_dir = Path(output_dir)
        self.ligand_descriptor = ligand_descriptor
        self.xml_template = xml_template
        self.batch_template = batch_template
        self.protein = protein

        DockJob.job_counter += 1

        self._id = DockJob.job_counter
    
    @property
    def protein_ligand_file(self):
        plf = self._docking_dir.joinpath('{}_{}.pdb'.format(
            self.protein.stem, self.ligand_descriptor.ligand_pdb.stem
        ))
        if not plf.is_file():
            join_protein_ligand_pdbs(self.protein, self.ligand_descriptor.ligand_pdb)
        
        return plf

    @property
    def _docking_dir(self):
        dd = self.output_dir.joinpath('docking')
        if not dd.is_dir():
            dd.mkdir()
        return dd
    
    @property
    def results_dir(self):
        rd = self.output_dir.joinpath('results')
        if not rd.is_dir():
            rd.mkdir()
        return rd
    
    @property
    def control_dir(self):
        cd = self.output_dir.joinpath('results')
        if not cd.is_dir():
            cd.mkdir()
        return cd
    
    @property
    def options_file(self):
        of = self.control_dir.joinpath('options')
        if not of.is_file():
            self._create_options_file()
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
        batch_template_string.format(
           self.options_file
        )
        with open(str(self.sbatch_file)) as handle:
            handle.write(batch_template_string)

    def submit(self):
        '''Submit this job to the job handler (SLURM).

        Returns:
            [type]: [description]
        '''
        os.system('sbatch {}'.format(self.sbatch_file))
    
    def _create_options_file(self, path):
        options_template_string = open(str(self.options_file)).read()
        options_template_string.format(
            self.protein_ligand_file,
            self.ligand_descriptor.ligand_params,
            self.xml_script_file
        )
        with open(str(path)) as handle:
            handle.write(options_template_string)
        
        return path