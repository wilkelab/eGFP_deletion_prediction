import os, re
import numpy as np


def extract_scores(score_filename):
    all_scores = []
    names = []
    print score_filename
    in_file = open(score_filename, "r")
    file_lines = in_file.readlines()
    file_lines = file_lines[2:]
    for line in file_lines:
        score = line[10:22]
        score = float(score.strip())
        all_scores.append(score)
        prename = line[-38:].strip()
        prename_parts = re.split("\.", prename)
        #print "prename_parts : ", prename_parts
        name_parts = re.split(" ", prename_parts[-2])  
        names.append(name_parts[1])
      
    mean_score = np.mean(all_scores)
    return names, all_scores

def main():
    #main_dir = "../modeller/rosetta_scoring/" #Directory with score files
    main_dir = "../modeller/relax/" #Directory with score files
    #out = open("../data/egfp_model_scores.csv", "w") #File with all of the scores for each model
    out = open("../data/egfp_relax_model_scores.csv", "w") #File with all of the scores for each model
    out.write("mutant,total_score\n")   
    #out2 = open("../data/egfp_model_summary.csv", "w") #File with mean scores for each model
    #out2.write("mutant,mean_total_score\n")
    score_dirs = [] 
    score_files = []
    for f in os.listdir(main_dir):
        if os.path.isdir(main_dir + f):
            score_dirs.append(main_dir + f)
    for dir in score_dirs:
        files = os.listdir(dir)
        #print files
        for file in files:
            if file.endswith(".fasc"):
                total_filename = dir + "/" + file
                names, scores = extract_scores(total_filename)
                mutant = names[0]
                #mean_score = np.mean(scores)
                #out2.write(mutant + "," + str(mean_score) + "\n")
                score_files.append(file) 
                num_mutants = len(names)
                for n in xrange(0, num_mutants):
                    out.write(names[n] + "," + str(scores[n]) + "\n")
    out.close()
    #out2.close()


if __name__ == "__main__":
    main()