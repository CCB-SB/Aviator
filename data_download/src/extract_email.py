import sys

with open(sys.argv[1]) as f, open(sys.argv[2], 'w') as of:
    of.write("PMID\tEmail\n")
    record = ""
    cur_tag = ""
    pmid = ""
    em = ""
    for l in f:
        if l[:2] != "  ":
            cur_tag = l[:2]
        if cur_tag == "PM":
            pmid = l[3:].rstrip()
        elif cur_tag == "EM":
            em = l[3:].rstrip()
        elif l == "\n":
            if pmid != "":
                of.write(f"{pmid}\t{em}\n")
            pmid = ""
            em = ""

