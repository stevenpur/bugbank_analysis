import os

def get_manhattan_file(pathogen = None, source = None, specimen = None, tax_lev = "species", stem = None, mac100 = True):
    """
    Get the manhattan plot for the given pathogen and source.
    """
    if source == "sgss":
        if specimen is None:
            specimen = "all"
    elif source == "hes":
        if specimen is not None:
            raise ValueError("Specimen should not be provided for HES data")
    elif stem is None:
        raise ValueError("Invalid source, must be 'sgss' or 'hes', or stem must be provided")
    if pathogen is not None:
        pathogen = pathogen.replace(" ", "_")
    manhattan_dir = os.path.expanduser("~/saige_pipe_test/")
    if stem is None:
        if source == "hes":
            stem = "05062023_" + source + "_" + tax_lev + "_regenie." + pathogen + "." + tax_lev
        elif source == "sgss":
            stem = "05062023_" + source + "_" + tax_lev + "." + pathogen + "." + tax_lev + "." + specimen
    if mac100:
        manhattan_file = "Manhattan.100MAC." + stem + "_regenie.png"
    else:
        manhattan_file = "Manhattan." + stem + "_regenie.png"
    manhattan_file = os.path.join(manhattan_dir, manhattan_file)
    # read the png file
    #img = Image.open(manhattan_file)
    return manhattan_file


    