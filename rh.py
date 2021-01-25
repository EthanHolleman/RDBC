from handler.ligands import LigandDescriptor
from handler.docking import DockJob
from handler.args import get_args
from handler import make_logger
from handler import DEFAULT_XML
from pathlib import Path


logger = make_logger(__name__)


def main():
    logger.info('Started run')
    args = get_args()
    pretty_args = ' '.join(
        ['{}: {}'.format(arg, getattr(args, arg)) for arg in vars(args)])
    logger.info('Submitted following arguments: {}'.format(pretty_args))
    ligand_descriptors = LigandDescriptor.generate_from_directory(args.ligands)
    DockJob.rosetta_exe = args.exe
    if args.xml_template:
        xml_template = args.xml_template
    else:
        xml_template = DEFAULT_XML

    protein = Path(args.protein)
    if protein.is_dir():  # if is a directory find all pdb files in this dir
        protein = [pdb for pdb in protein.iterdir() if pdb.suffix == '.pdb']
    jobs = []
    if isinstance(protein, list):
        # multiple proteins to deal with so need to create directories for
        # each of these proteins
        for pdb_path in protein:
            # this assumes unique filename for each pdb
            pdb_dir = args.parent.joinpath(pdb_path.stem)
            pdb_dir.mkdir()
            for ligand in ligand_descriptors:
                jobs.append(
                    DockJob(
                        pdb_dir.joinpath(ligand.name), ligand, protein=pdb_path,
                        number_iterations=args.num_iters,
                        xml_template=xml_template
                    )
                )
    else:  # only need to create jobs for one protein
        jobs = (DockJob(
            args.parent.joinpath(ligand.name),
            ligand,
            protein=args.protein,
            number_iterations=args.num_iters,
            xml_template=xml_template)
            for ligand in ligand_descriptors)

    if args.moist:  # only run the one job (usually for testing)
        logger.info('Running moist')
        first_job = next(jobs)
        first_job.set_up_for_submit()
        first_job.submit()
    else:
        job_counter = 0
        for j in jobs:
            job_counter += 1
            setup_good = j.set_up_for_submit()
            if not args.dry:
                j.submit()
                logger.info('Submitted job: {}'.format(j.name))
        if job_counter == 0:
            logger.critical('NO JOBS WERE CREATED')


if __name__ == "__main__":
    main()
