from handler.file_manager import FileManager
from handler.ligands import LigandDescriptor
from handler.docking import DockJob
from handler.args import get_args

l = LigandDescriptor.generate_from_directory(
    '/home/ethan/data/lab_notes/Gino/sdf/ligand_params')
l1 = next(l)

d = DockJob('.', l1, 'my.pdb')


def main():
    args = get_args()
    
    ligand_descriptors = LigandDescriptor.generate_from_directory(args.ligands)
    DockJob.rosetta_exe = args.exe
    jobs = (DockJob(
        args.parent.joinpath(ligand.name), 
        ligand, 
        protein=args.protein) 
    for ligand in ligand_descriptors)

    for j in jobs:
        setup_good = j.set_up_for_submit()
        
        j.submit()


if __name__ == "__main__":
    main()
        


