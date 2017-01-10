import utils, poised_filter
import sys, gzip, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

### start field name defintions #########################################
#TODO TIM ??? WHAt are these

field_Reactgroup = "poised"

### start main execution #########################################

# Define SMARTS patterns here


def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('--smiles', help='query structure as smiles (incompatible with -molfile arg)')
    parser.add_argument('--molfile',
                        help='query structure as filename in molfile format (incompatible with -smiles arg)')
    parser.add_argument('-f', '--fragment', choices=['hac', 'mw'],
                        help='Find single fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight )')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    utils.add_default_io_args(parser)
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')

    args = parser.parse_args()
    utils.log("Screen Args: ", args)


    if args.smiles and args.molfile:
        raise ValueError('Cannot specify -smiles and -molfile arguments together')
    elif args.smiles:
        query_rdkitmol = Chem.MolFromSmiles(args.smiles)
    elif args.molfile:
        query_rdkitmol = Chem.MolFromMolFile(args.molfile)
    else:
        raise ValueError('No query structure specified')

    input, output, suppl, writer, output_base = utils.default_open_input_output(args.input, args.informat, args.output,
                                                                                'screen', args.outformat)

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    if args.meta:
        t = open(output_base + '_types.txt', 'w')
        t.write(field_Reactgroup+'\n')
        #TODO TIM I don't understand this??t.write(field_Similarity + '=integer\n')
        t.flush()
        t.close()

    i = 0
    count = 0
    for mol in suppl:
        i += 1
        if mol is None: continue
        if args.fragment:
            mol = filter.fragment(mol, args.fragment, quiet=args.quiet)
        if not filter.filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax,
                             quiet=args.quiet):
            continue
        # Return a dict/class here - indicating which filters passed
        filter_pass = poised_filter.pass_filter(mol)
        if filter_pass:
            count += 1
            if not args.quiet:
                utils.log(i, "passed")
            for name in mol.GetPropNames():
                mol.ClearProp(name)
            mol.SetDoubleProp(field_Reactgroup, "passed")
            writer.write(mol)

    utils.log("Found", count, "similar molecules")

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RDKitScreen': count})

    return count


if __name__ == "__main__":
    main()