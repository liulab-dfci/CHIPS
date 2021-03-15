import os

def checking_targets(wildcards):
    ls = []
    ls.append(output_path + "/logs/empty_file_list.txt")
    return ls


def emptyCheckingInput(wildcards):
    # # Same as all_targets. Should be modified if all_targets modified.
    # # however DO NOT have checking_targets
    ls = all_targets(wildcards)
    for t in checking_targets(wildcards):
    	#print(t)
        if t in ls:
            ls.remove(t)
    return ls


rule emptyChecking:
    input:
        emptyCheckingInput
    output:
        file=output_path + "/logs/empty_file_list.txt"
    message:
        "EMPTYCHECKING: checking whether any files are empty"
    run:
        empty_list = []
        for i in input:
            size = os.path.getsize(i)
            if size == 0:
                empty_list.append(i)
            else:
                continue
        with open(output.file,"w") as op:
            op.write("\n".join(empty_list))














