#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import pysam
from collections import defaultdict

def cmp(a, b):
    """Python 3 中不再有 cmp 函数，这里自定义一个"""
    return (a > b) - (a < b)

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
    linenum = scaflen // 60  # 使用整数除法
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
                if contig[-1] == "-":  # 直接比较字符串
                    ctg_seq = recom(ctg_seq)
                ctg_lst.append(ctg_seq)
            scaf_seq = ns.join(ctg_lst)
            print(">Superscaffold{}_{}_{}".format(n+1, len(scaf), len(scaf_seq)), file=output)
            format_print(scaf_seq, output)

        for n, contig in enumerate(discard):
            ctg_name = contig[0].split("::")[0]
            ctg_seq = fa_seq[ctg_name][contig[1]:contig[2]]
            if contig[-1] == "-":  # 直接比较字符串
                ctg_seq = recom(ctg_seq)
            print(">unanchor{}_1_{}".format(n+1, len(ctg_seq)), file=output)
            format_print(ctg_seq, output)

def valid_output(out_fa, ref):
    out_dict, ref_dict = defaultdict(int), defaultdict(int)
    base_lst = ["A", "T", "C", "G", "N"]
    for line in open(out_fa):
        if line.startswith(">"): 
            continue
        line = line.rstrip().upper()
        for x in base_lst:
            out_dict[x] += line.count(x)

    for line in open(ref):
        if line.startswith(">"): 
            continue
        line = line.rstrip().upper()
        for x in base_lst:
            ref_dict[x] += line.count(x)

    print("Stat of {}:".format(out_fa))
    print("Ref_fa [ATCG]: {:,}".format(sum([ref_dict[x] for x in base_lst if x != "N"])))
    print("Ref_fa [N   ]: {:,}".format(ref_dict["N"]))
    print("Out_fa [ATCG]: {:,}".format(sum([out_dict[x] for x in base_lst if x != "N"])))
    print("Out_fa [N   ]: {:,}".format(out_dict["N"]))
    print()

def ass2fasta(inass, ref, chr_num, output_prefix="input_prefix"):
    fa_seq = {fa.name: fa.sequence for fa in pysam.FastxFile(ref)}

    if output_prefix == "input_prefix":
        output_prefix = os.path.splitext(os.path.basename(inass))[0]
    out_fa = "{}.fasta".format(output_prefix)
    out_tab = "{}.order".format(output_prefix)

    superscaf, discard = get_pos_lst(inass, chr_num, out_tab)
    lst2seq(superscaf, discard, fa_seq, out_fa)

    valid_output(out_fa, ref)

    print("[output] {}, {}".format(out_fa, out_tab))

def generate_hic(assembly, outdir):
    # generate cprops and ass file
    # this program fix the bug of hicat-view which generate cprops and ass file
    # without EOF, namely without "\n" at last line
    try: 
        os.stat(outdir)
    except: 
        os.mkdir(outdir)
    
    out_cprops = os.path.splitext(os.path.basename(assembly))[0] + '.cprops'
    out_cprops = os.path.join(outdir, out_cprops)
    out_ass = os.path.splitext(os.path.basename(assembly))[0] + ".asm"
    out_ass = os.path.join(outdir, out_ass)

    with open(out_cprops, 'w') as cprops_h, open(out_ass, 'w') as ass_h:
        for line in open(assembly):
            if line.startswith(">"):
                print(line.strip().strip(">"), file=cprops_h)
            else:
                print(line.strip(), file=ass_h)

def main(args):
    try: 
        os.stat(args["outdir"])
    except: 
        os.mkdir(args["outdir"])

    for x in ["ref", "mnd", "assembly"]:
        args[x] = os.path.abspath(args[x])
        
    cwd_path = os.getcwd()
    script_path = os.path.split(os.path.realpath(__file__))[0]
    os.chdir(args["outdir"])
    
    def trim_unanchor(input_ass, chrom_num, output_ass):
        with open(input_ass) as inass, open(output_ass, 'w') as outass:
            count = 1
            for line in inass:
                line = line.strip()
                if line.startswith(">"):
                    print(line, file=outass)
                else:
                    if count > chrom_num: 
                        break
                    print(line, file=outass)
                    count += 1            
        
    ass2fasta(args["assembly"], args["ref"], args["chrom_num"], args["sample"])
    
    assembly = "{}.assembly".format(args["sample"])
    if os.path.exists(assembly): 
        os.remove(assembly)
    
    trim_unanchor(args["assembly"], args["chrom_num"], assembly)    
    ass2fasta(assembly, args['ref'], args["chrom_num"], "{}_chrom".format(args["sample"]))

    generate_hic(assembly, os.getcwd())

    cprops = "{}.cprops".format(os.path.splitext(assembly)[0])
    asm = "{}.asm".format(os.path.splitext(assembly)[0])

    _3ddna_path = args.get("_3ddna_path", "")
    
    with open('new_mnd_and_visualizer.sh', 'w') as outshell:
        print("#!/usr/bin/bash\n", file=outshell)
        cmd = "{} {}/edit/edit-mnd-according-to-new-cprops.sh {} {} > {}".format(
            "/usr/bin/bash", args['_3ddna_path'], cprops,  
            args["mnd"], "{}.new.mnd".format(args["sample"]))
        print(cmd, file=outshell)
        
        cmd = '{} {}/visualize/run-asm-visualizer.sh {} {} {} {} {} {}'.format(
            "/usr/bin/bash", args['_3ddna_path'], "-c",
            "-r", "1000000,500000,250000",
            cprops, asm, "{}.new.mnd".format(args["sample"]))
        print(cmd, file=outshell)

    with open('hic_viewer.sh', 'w') as plotshell:
        cmd = """#!/usr/bin/bash
python {2}/hic_viewer.py --hicfile {0}.hic --ref {0}_chrom.fasta --rslu {1} --outpfix {0} --norm {3}""".format(
            args['sample'], args["resolution"], script_path,args['norm_method'])
        print(cmd, file=plotshell)

    subprocess.call(["/usr/bin/bash", 'new_mnd_and_visualizer.sh'])
    os.system("bash hic_viewer.sh")
    os.chdir(cwd_path)

def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', "--ref", dest="ref", required=True,
        help="raw contig file [required]")
    parser.add_argument("-a", "--assembly", dest="assembly", required=True,
        help="review assembly file, generated by juicerbox [required]")
    parser.add_argument("-m", "--mnd", dest='mnd', required=True,
        help="raw mnd file used to generate hic file [required]")
    parser.add_argument("-c", "--chrom_num", dest="chrom_num", type=int, required=True,
        help='chrom number, used to parse assembly file [required]')
    parser.add_argument('-e', "--resolution", dest='resolution', type=int, default=500000,
        help='resolution to heatmap visualize [default:500000]')
    parser.add_argument('-n', '--norm_method', dest='norm_method', 
        choices={"KR", "VC", "VC_SQRT", "NONE"}, default="KR", 
        help="normalization method, choices from {'KR', 'VC', 'VC_SQRT', 'NONE'}, default:KR")
    parser.add_argument('-s', "--sample", dest='sample', default="sample", type=str, 
        help='sample name used to generate files [default:sample]')
    parser.add_argument('-o', '--outdir', dest='outdir', default='./',
        help="output dir, [default:$PWD]")
    parser.add_argument("--3ddna_path", dest='_3ddna_path',
        help='3ddna software path')
    args = vars(parser.parse_args())
    return args

if __name__ == '__main__':
    myargs = args_parser()
    main(myargs)
