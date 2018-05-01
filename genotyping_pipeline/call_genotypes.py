#!/usr/bin/env python

import argparse
import subprocess
import os

def parse_config(config_file):
    fam = cel_map = threshold = chip = build = outputname = ""
    birdsuite_config = ""
    with open(config_file) as configs:
        for line in configs:
            line = line.strip()
            if line.startswith('famFile'):
                fam = "-f " + line.split("=",1)[1]
            elif line.startswith('celMap'):
                cel_map = "-m " + line.split("=",1)[1]
            elif line.startswith('threshold'):
                threshold = "-t " + line.split("=",1)[1]
            elif line.startswith('chipType'):
                chip = "-c " + line.split("_",1)[1]
                birdsuite_config += " --" + line
            elif line.startswith('genomeBuild'):
                parameter = line.split("=",1)[1]
                if parameter == "hg17":
                    build = "-b 17"
                elif parameter == "hg18":
                    build = "-b 18"
                birdsuite_config += " --" + line
            elif line.startswith('outputName'):
                outputname = "-n " + line.split("=",1)[1]
            else:
                birdsuite_config += " --" + line
    return fam, cel_map, threshold, chip, build, outputname, birdsuite_config

def generate_cel_gender_files(individuals_file):
    plates = {} # key = plate ID, value = [(cel file, gender)]
    
    with open(individuals_file) as individuals:
        for line in individuals:
            line = line.strip().split()
            cel_file = line[0]
            gender = line[1]
            plateID = line[2]
            try:
                plates[plateID].append((cel_file,gender))
            except KeyError:
                plates[plateID] = []
                plates[plateID].append((cel_file,gender))
    
    for plateID in plates.keys():
        cels = open(plateID + '.cels','w')
        print >> cels, "cel_files"
        genders = open(plateID + '.gender','w')
        print >> genders, "gender"
        individuals = plates[plateID]
        for (cel_file,gender) in individuals:
            print >> cels, cel_file
            print >> genders, gender
        cels.close()
        genders.close()

    return plates.keys()

def call_birdsuite_serial(plates, fam, cel_map, threshold, chip, build, outputname, birdsuite_config, output):
    path = os.path.dirname(os.path.realpath(__file__))
    
    if not os.path.exists(output):
        os.makedirs(output)
    
    for plateID in plates:
        if not os.path.exists(output + "/" + plateID):
            os.makedirs(output + "/" + plateID)
        subprocess.call([path + "/bin/birdsuite.sh", "--basename=" + plateID,"--celFiles=" + plateID + ".cels", "--genderFile=" + plateID + ".gender", "--outputDir=" + output + "/" + plateID, birdsuite_config])
        move_genders_file = "mv " + plateID + ".gender " + output + "/" + plateID
        os.system(move_genders_file)

    print "Calling Birdsuite_to_PLINK..."
    birdsuite_to_PLINK_command = ' '.join([path + "/bin/birdsuite_converter", "-p " + args.output, "-x " + path + "/metadata/", "-d " + args.output, fam, cel_map, threshold, chip, build, outputname])
    print(birdsuite_to_PLINK_command)
    subprocess.call(birdsuite_to_PLINK_command,shell=True)

def call_birdsuite_parallel(plates, fam, cel_map, threshold, chip, build, outputname, birdsuite_config, output):
    path = os.path.dirname(os.path.realpath(__file__))
    
    if not os.path.exists(output):
        os.makedirs(output)
    
    for plateID in plates:
        if not os.path.exists(output + "/" + plateID):
            os.makedirs(output + "/" + plateID)
        script = open("run_" + plateID + "_birdsuite.sh",'w')
        birdsuite_command = ' '.join([path + "/bin/birdsuite.sh", "--basename=" + plateID,"--celFiles=" + plateID + ".cels", "--genderFile=" + plateID + ".gender", "--outputDir=" + output + "/" + plateID, birdsuite_config])
        print >> script, birdsuite_command
        move_genders_command = "mv " + plateID + ".gender " + output + "/" + plateID
        print >> script, move_genders_command
        script.close()
        subprocess.call("chmod a+x run_" + plateID + "_birdsuite.sh",shell=True)

    script = open("run_birdsuite_to_plink.sh",'w')
    birdsuite_to_PLINK_command = ' '.join([path + "/bin/birdsuite_converter", "-p " + args.output, "-x " + path + "/metadata/", "-d " + args.output, fam, cel_map, threshold, chip, build, outputname])
    print >> script, birdsuite_to_PLINK_command
    script.close()
    subprocess.call("chmod a+x run_birdsuite_to_plink.sh",shell=True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Arguments to call genotypes. See XWAS manual for full usage details.')

    parser.add_argument('-c', '--config', action="store", dest="config", type=str, required = True, help='Config file for Birdsuite, Birdsuite_to_PLINK options')
    parser.add_argument('-i', '--individuals', action="store", dest="individuals", type=str, required = True, help='Individuals file, tab deliminated: cel file, gender, plate ID')
    parser.add_argument('-o', '--output', action="store", dest="output", type=str, required = True, help='Output directory name')
    parser.add_argument('-p', '--parallel', action='store_true', dest='parallel', help = 'Use this flag to create parallelizable scripts')

    args = parser.parse_args()

    print "Parsing " + args.config + "..."
    fam, cel_map, threshold, chip, build, outputname, birdsuite_config = parse_config(args.config)

    print "Generating cel files and gender files for each plate..."
    plates = generate_cel_gender_files(args.individuals)

    if args.parallel:
        print "Creating Birdsuite scripts for each plate..."
        call_birdsuite_parallel(plates, fam, cel_map, threshold, chip, build, outputname, birdsuite_config, args.output)
    else:
        print "Calling Birdsuite on each plate..."
        call_birdsuite_serial(plates, fam, cel_map, threshold, chip, build, outputname, birdsuite_config, args.output)
