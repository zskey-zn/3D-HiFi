import os
import pysam
import argparse
from collections import defaultdict

def get_pos_lst(inass, chr_num, out_tab):
    all_dict, superscaf_list, discard_list = {}, [], []
    superscaf_index = 1
    loc, last_ctg_name = 0, ""
    out_tab_h = open(out_tab, 'w')
    unanchor = 1
    for line in open(inass):
        line = line.strip()
        llst = line.split()
        if line.startswith(">"):
            ctg_name = line.split()[0].split("::")[0].strip(">")
            if ctg_name != last_ctg_name:
                loc, last_ctg_name = 0, ctg_name

            start, end = loc, loc + int(llst[-1])
            all_dict[llst[1]] = [line.split()[0].strip(">"), start, end]
            loc += int(llst[-1])

        else:
            line_info = []
            for i in llst:
                strand = ''
                if i.startswith("-"):
                    strand = "-"
                    info = all_dict[i.strip("-")] + [strand]
                else:
                    strand = '+'
                    info = all_dict[i] + [strand]
                line_info.append(info)

                if superscaf_index <= int(chr_num):
                    print("Superscaffold{}\t{}\t{}\t{}".format(superscaf_index, strand, info[0], info[2]-info[1]), file=out_tab_h)
                else:
                    print("unanchor{}\t{}\t{}\t{}".format(unanchor, strand, info[0], info[2]-info[1]), file=out_tab_h)
                    unanchor += 1

            if superscaf_index <= int(chr_num):
                superscaf_list.append(line_info)
                superscaf_index += 1
            else:
                discard_list.extend(line_info)
    out_tab_h.close()
    return superscaf_list, discard_list

def recom(seq):
    seq = seq.upper()
    seq = seq.replace("A", "t")
    seq = seq.replace("C", "g")
    seq = seq.replace("G", "c")
    seq = seq.replace("T", 'a')
    seq = seq.replace("N", 'n')
    recom_seq = seq.upper()[::-1]
    return recom_seq

def format_print(scafseq, outfile):
    scaflen = len(scafseq)
    linenum = scaflen // 60
    for i in range(linenum):
        start, end = i*60, (i+1)*60
        print(scafseq[start:end], file=outfile)
    if scaflen % 60 != 0:
        print(scafseq[linenum*60:], file=outfile)

def lst2seq(superscaf, discard, fa_seq, outfa):
    ns = "N" * 500
    with open(outfa, 'w') as output:
        for n, scaf in enumerate(superscaf):
            ctg_lst = []
            for contig in scaf:
                ctg_name = contig[0].split("::")[0]
                ctg_seq = fa_seq[ctg_name][contig[1]:contig[2]]
                if contig[-1] == "-":
                    ctg_seq = recom(ctg_seq)
                ctg_lst.append(ctg_seq)
            scaf_seq = ns.join(ctg_lst)
            print(">Superscaffold{}_{}_{}".format(n+1, len(scaf), len(scaf_seq)), file=output)
            format_print(scaf_seq, output)

        for n, contig in enumerate(discard):
            ctg_name = contig[0].split("::")[0]
            ctg_seq = fa_seq[ctg_name][contig[1]:contig[2]]
            if contig[-1] == "-":
                ctg_seq = recom(ctg_seq)
            print(">{}".format(contig[0]), file=output)
            format_print(ctg_seq, output)

def valid_output(out_fa, ref):
    out_dict, ref_dict = defaultdict(int), defaultdict(int)
    base_lst = ["A", "T", "C", "G", "N"]
    for line in open(out_fa):
        if line.startswith(">"): continue
        line = line.rstrip().upper()
        for x in base_lst:
            out_dict[x] += line.count(x)

    for line in open(ref):
        if line.startswith(">"): continue
        line = line.rstrip().upper()
        for x in base_lst:
            ref_dict[x] += line.count(x)

    print("Stat of {}:".format(out_fa))
    print("Ref_fa [ATCG]: {:,}".format(sum([ref_dict[x] for x in base_lst if x != "N"])))
    print("Ref_fa [N   ]: {:,}".format(ref_dict["N"]))
    print("Out_fa [ATCG]: {:,}".format(sum([out_dict[x] for x in base_lst if x != "N"])))
    print("Out_fa [N   ]: {:,}".format(out_dict["N"]))
    print()

def sort_ass(inass, chr_num, outass):
    chrom_id, frag_len, chrom_dict, switch = 0, {}, {}, True
    with open(outass, 'w') as outbuff, open(inass) as inbuff:
        for line in inbuff:
            line = line.strip()
            llst = line.split()
            if line.startswith(">"):
                print(line, file=outbuff)
                frag_len[llst[1]] = int(llst[2])
            elif chrom_id < chr_num:
                chrom_len = 0
                for i in llst:
                    chrom_len += frag_len[i.strip("-")]
                chrom_dict[line] = chrom_len
                chrom_id += 1
            else:
                if switch:
                    sorted_list = sorted(chrom_dict.items(), key=lambda x:x[1], reverse=True)
                    for cline, _ in sorted_list:
                        print(cline, file=outbuff)
                    switch = False
                    print(line, file=outbuff)
                else:
                    print(line, file=outbuff)
    return outass

def ass2fasta(inass, ref, chr_num, output_prefix="input_prefix", sort=False):
    fa_seq = {fa.name: fa.sequence for fa in pysam.FastxFile(ref)}

    if output_prefix == "input_prefix":
        output_prefix = os.path.splitext(os.path.basename(inass))[0]
    out_fa = "{}.fasta".format(output_prefix)
    out_tab = "{}.order".format(output_prefix)
    out_ass = "{}_sort.assembly".format(output_prefix)

    if sort:
        inass = sort_ass(inass, chr_num, out_ass)
    superscaf, discard = get_pos_lst(inass, chr_num, out_tab)
    lst2seq(superscaf, discard, fa_seq, out_fa)

    valid_output(out_fa, ref)

    if sort:
        print("[output] {}, {}, {}".format(out_ass, out_fa, out_tab))
    else:
        print("[output] {}, {}".format(out_fa, out_tab))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly", required=True,
                        help="input assembly file [required]")
    parser.add_argument("-r", '--ref', required=True,
                        help="input reference fasta file [required]")
    parser.add_argument("-n", '--chr_num', required=True, type=int,
                        help='chromosome number [required]')
    parser.add_argument("-o", '--outpfix', default="input_prefix",
                        help='output prefix. if not assign, then guess from input prefix. [default: <inpfix>]')
    parser.add_argument("-s", '--sort', action="store_true",
                        help='whether sort chromosome by chrom length [default: False]')
    args = parser.parse_args()
    ass2fasta(args.assembly, args.ref, args.chr_num, args.outpfix, args.sort)
