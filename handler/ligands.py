from pathlib import Path

class LigandDescriptor():
    '''Class which contains (or at least points to the files which do) all data
    used to ultimately describe and individual ligand and conformers (rotomers)
    to Rosetta.  
    '''

    def __init__(self, ligand_pdb, ligand_params, ligand_conformers):
        self.ligand_pdb = Path(ligand_pdb)
        self.ligand_params = Path(ligand_params)
        self.ligand_conformers = Path(ligand_conformers)
    
    @property
    def name(self):
        return self.ligand_pdb.name
    
    @classmethod
    def generate_from_directory(cls, dir):
        '''Creates a generator of LigandDescriptor objects from a directory
        containing the files that describe one or more LigandDescriptor objects.

        Note:
            This method assumes a standard formating for the three files that
            compose a LigandDescriptor object. The ligand_pdb file should be
            the name of the ligand followed by the extension .pdb. The
            ligand_params file should be the ligand name followed by the
            extension .params and the ligand_conformers file should be the
            ligand name followed by the suffix and extension _conformers.pdb.

        Args:
            dir (str): Path to input directory

        Returns:
            list: List of LigandDescriptor objects created from the directory.
        
        '''
        temp_dir = Path(dir)
        ligand_files = list(Path(dir).iterdir())
        ligand_names = [lf.stem for lf in ligand_files if lf.suffix == '.pdb' and '_' not in lf.name]
        for name in ligand_names:
            expected_files = [
                temp_dir.joinpath(name + '.pdb'),
                temp_dir.joinpath(name + '.params'),
                temp_dir.joinpath(name + '_conformers.pdb')
            ]
            assert all(ef.is_file() for ef in expected_files)
            yield cls(*expected_files)