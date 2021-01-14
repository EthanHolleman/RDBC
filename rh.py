from handler.ligands import LigandDescriptor
from handler.docking import DockJob
from handler.args import get_args
from handler import make_logger


logger = make_logger(__name__)


def main():
    logger.info('Started run')
    args = get_args()
    pretty_args = ' '.join(['{}: {}'.format(arg, getattr(args, arg)) for arg in vars(args)])
    logger.info('Submitted following arguments: {}'.format(pretty_args))
    ligand_descriptors = LigandDescriptor.generate_from_directory(args.ligands)
    DockJob.rosetta_exe = args.exe
    jobs = (DockJob(
        args.parent.joinpath(ligand.name), 
        ligand, 
        protein=args.protein) 
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
        


