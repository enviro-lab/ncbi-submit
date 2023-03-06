#!/usr/bin/env python3
from pathlib import Path
from ncbi_submit.helpers import ensure_outdir_viable

def write_file(filename,outdir):
    """Writes `filename` to `outdir`.

    Args:
        filename (str): name of file to write out
        outdir (Path): output directory
    """

    infile = Path(__file__).parent / filename
    outfile = outdir / filename
    print(f"Writing out '{outfile}'")
    with infile.open() as fh, outfile.open('w') as out:
        for line in fh: out.write(line)

def get_files(outdir,config=False,template=False):
    """Writes out the desired file(s) to the `outdir`.

    Args:
        outdir (Path | str): The directory to write out file.
        config (bool, optional): A flag to write out the config file to `outdir`. Defaults to False.
        template (bool, optional): A flag to write out the template file to `outdir`. Defaults to False.
    """

    # ensure directory can be used
    outdir = ensure_outdir_viable(outdir)
    
    # write out desired file(s)
    if config:
        write_file("config.py",outdir)
    if template:
        write_file("template.sbt",outdir)
