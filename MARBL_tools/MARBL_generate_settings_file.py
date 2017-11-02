#!/usr/bin/env python

"""
Generate a MARBL settings file containing values for all MARBL settings, as defined
in a YAML file (and optional additional settings file).

This file can be run as a command line script or imported as part of the MARBL_tools module.
To use from a module, all arguments are required in the call

generate_settings_file(yaml_file, settings_file_in, grid, settings_file_out)

To run from a command line, the following arguments are all optional otherwise default values are used:

  -y YAML_FILE, --yaml_file YAML_FILE
                        Location of YAML-formatted MARBL configuration file
                        (default:
                        /home/mlevy/codes/marbl/src/default_settings.yaml)
  -g GRID, --grid GRID  Some default values are grid-dependent (default:
                        CESM_x1)
  -i SETTINGS_FILE_IN, --settings_file_in SETTINGS_FILE_IN
                        A file that overrides values in YAML (default: None)
  -o SETTINGS_FILE_OUT, --settings_file_out SETTINGS_FILE_OUT
                        Name of file to be written (default: marbl.input)
  -d SETTINGS_CLASS_DIR, --MARBL_settings_class_override_dir SETTINGS_CLASS_DIR
                        Alternate location of MARBL_settings_file_class.py
                        (default: None)

"""

#######################################

def generate_settings_file(yaml_file, settings_file_in, grid, settings_file_out, settings_class_dir=None):
    """ Produce a valid MARBL input file from a YAML parameter file
    """

    if settings_class_dir == None:
        from MARBL_tools import MARBL_settings_file_class
    else:
        import imp
        import os
        import sys
        import logging
        logger = logging.getLogger(__name__)
        logging.info('Importing MARBL_settings_file_class.py from %s' % settings_class_dir)
        settings_class_module = settings_class_dir+'/MARBL_settings_file_class.py'
        if os.path.isfile(settings_class_module):
            MARBL_settings_file_class = imp.load_source('MARBL_settings_file_class', settings_class_module)
        else:
            logger.error('Can not find %s' % settings_class_module)
            sys.exit(1)

    DefaultParms = MARBL_settings_file_class.MARBL_settings_class(yaml_file, grid, settings_file_in)

    fout = open(settings_file_out,"w")
    # Sort variables by subcategory
    for subcat_name in DefaultParms.get_subcategory_names():
        fout.write("! %s\n" % subcat_name.split('. ')[1])
        for varname in DefaultParms.get_parm_dict_variable_names(subcat_name):
            fout.write("%s = %s\n" % (varname, DefaultParms.parm_dict[varname]))
        if subcat_name != DefaultParms.get_subcategory_names()[-1]:
            fout.write("\n")
    fout.close()

#######################################

def _parse_args(marbl_root):
    """ Parse command line arguments
    """

    import argparse

    parser = argparse.ArgumentParser(description="Generate a MARBL settings file from a YAML file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Command line argument to point to YAML file (default is $MARBLROOT/src/default_settings.yaml)
    parser.add_argument('-y', '--yaml_file', action='store', dest='yaml_file',
                        default=os.path.join(marbl_root, 'src', 'default_settings.yaml'),
                        help='Location of YAML-formatted MARBL configuration file')

    # Command line argument to specify resolution (default is CESM_x1)
    parser.add_argument('-g', '--grid', action='store', dest='grid', default='CESM_x1',
                        help='Some default values are grid-dependent')

    # Command line argument to specify an input file which would override the YAML
    parser.add_argument('-i', '--settings_file_in', action='store', dest='settings_file_in', default=None,
                        help='A file that overrides values in YAML')

    # Command line argument to specify an input file which would override the YAML
    parser.add_argument('-o', '--settings_file_out', action='store', dest='settings_file_out', default='marbl.input',
                        help='Name of file to be written')

    # Path to directory containing MARBL_defaults.py
    parser.add_argument('-d', '--MARBL_settings_class_override_dir', action='store', dest='settings_class_dir',
                        default=None, help='Alternate location of MARBL_settings_file_class.py')

    return parser.parse_args()

#######################################

if __name__ == "__main__":
    # We need marbl_root in python path so we can import MARBL_tools from generate_settings_file()
    import sys, os
    marbl_root = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))
    sys.path.append(marbl_root)

    # Parse command-line arguments (marbl_root is used to set default for YAML file location)
    args = _parse_args(marbl_root)

    # Set up logging
    import logging
    logging.basicConfig(format='%(levelname)s (%(funcName)s): %(message)s', level=logging.DEBUG)

    # Write the input file
    generate_settings_file(args.yaml_file, args.settings_file_in, args.grid, args.settings_file_out, args.settings_class_dir)