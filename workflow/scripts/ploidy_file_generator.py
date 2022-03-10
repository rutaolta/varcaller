#!/usr/bin/env python3
from argparse import ArgumentParser, FileType
from statistics import mean


def get_PAR_coordinates(*bed_files) -> dict:
    bed_files = bed_files[0]
    all_start_coordinates = []
    all_stop_coordinates = []
    for infile in bed_files:
        with open(infile, 'r') as f:
            for line in f:
                chr_name, start_coordinate, stop_coordinate = line.strip().split('\t')
                all_start_coordinates.append(int(start_coordinate))
                all_stop_coordinates.append(int(stop_coordinate))
    if len(all_start_coordinates) != len(bed_files):
        print("Error! Different number of lines in input PAR.bed files!")
    if len(set(all_start_coordinates)) != 1:
        print(f"Warning! The start coordinates are different:\n{all_start_coordinates}")
        print("The mean value will be written to the ploidy.file.")
    if len(set(all_stop_coordinates)) != 1:
        print(f"Warning! The stop coordinates are different:\n{all_stop_coordinates}")
        print("The mean value will be written to the ploidy.file.")
    start_coordinate = round(mean(all_start_coordinates)) + 1 # bcftools needs 1-based coords
    stop_coordinate = round(mean(all_stop_coordinates)) + 1
    return (chr_name, start_coordinate, stop_coordinate)

def get_sex_chromosome_length(lenfile, chr_name):
    with open(lenfile, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if chr_name == line[0]:
                sex_chromosome_length = int(line[1])
                break
    return sex_chromosome_length

def main():
    # get start and stop coordinates (mean values will be returned if they differ)
    chr_name, start_coordinate, stop_coordinate = get_PAR_coordinates(args.input)
    # get length of sex chromosome
    sex_chromosome_length = get_sex_chromosome_length(args.lenfile, chr_name)
    # get resulting ploidy.file
    with open(args.output, 'a') as ploidy_outfile:
        ploidy_outfile.write(f"{chr_name}\t1\t{start_coordinate}\tM\t1\n")
        ploidy_outfile.write(f"{chr_name}\t{stop_coordinate}\t{sex_chromosome_length}\tM\t1\n")
        # ploidy_outfile.write("*\t*\t*\tM\t2\n")
        # ploidy_outfile.write("*\t*\t*\tF\t2\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, nargs='+', required=True,
                             help="BED files with PAR coords (space separated list)")
    group_required.add_argument('-l', '--lenfile', type=str, help="assembly.len file")
    group_required.add_argument('-o','--output', type=str, help="output file name")
    args = parser.parse_args()
    main()

