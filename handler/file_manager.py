from pathlib import Path

class DockJob():

    job_counter = 0
    rosetta_exe = '/bin/rosetta_scripts.linuxgccrelease'

    def __init__(self, output_dir, ligand_descriptor, protein, xml_script,
                options_template, batch_template):
        self.output_dir = Path(output_dir)
        self.ligand_descriptor = ligand_descriptor
        self.protein = protein
        self.xml_script = xml_script
        self.options_template = options_template
        self.batch_template = batch_template

        DockJob.job_counter += 1

        self._id = DockJob.job_counter

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
        
    def generate_batch_file(self):
        pass

    def submit(self):
        '''Submit this job to the job handler (SLURM).

        Returns:
            [type]: [description]
        '''


class FileManager():
    '''The FileManager class should handle the organization of directories for
    each individual ligand docking experiment, all of which will be created
    as children of the base_dir (parent directory)
    '''

    def __init__(self, base_dir, ligands):
        self.base_dir = base_dir
        self.ligands = ligands
    
    def set_up_directories(self):
        pass


    