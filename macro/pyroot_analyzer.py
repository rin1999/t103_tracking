import ROOT
import argparse
import sys
import logging
from os.path import expanduser

def main(input_filepath: str, output_filepath: str, tree_name: str):
    logging.debug(f'Reading inputfile: {input_filepath}')
    file = ROOT.TFile.Open(input_filepath)
    if (not file) or (file.IsZombie()):
        logging.error(f'Input file [{input_filepath}] cannot be found!')
        sys.exit()
    
    logging.debug(f'Reading TTree: {tree_name}')
    tree = file.Get(tree_name)
    if (not tree):
        logging.error(f'There\'s no TTree: {tree_name}')
        file.Close()
        sys.exit()
    
    branch_names = [branch.GetName() for branch in tree.GetListOfBranches()]

    



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file (rootfiles in root_analyzed)')
    parser.add_argument('-o', '--output', type=str, default=' ~/workspace/t103/Analyzer2/image/test.pdf', help='Path to the output file (default=image/test.pdf)')
    parser.add_argument('-l', '--loglevel', default='INFO', choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help='Set log level (default=INFO) choices=[DEBUG, INFO, WARNING, ERROR, CRITICAL]')
    parser.add_argument('-t', '-tree', type=str, default='tree', help='Name of TTree (default=tree)')
    args = parser.parse_args()
    input_filepath = expanduser(args.input)
    output_filepath= expanduser(args.output)
    tree_name = args.tree
    logging.basicConfig(level=getattr(logging, args.log_level))

    logging.debug(f'Input filepath: {input_filepath}')
    logging.debug(f'Output filepath: {output_filepath}')

    main(input_filepath, output_filepath, tree_name)

