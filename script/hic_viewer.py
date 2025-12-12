import fire 
import straw
from Bio import SeqIO
import numpy as np
import pandas as pd 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
__version__ = "v2.0"
''' change log:
20190712. add `ref' options, because in sometimes
          the bin number of hic file (which genome size > 2G), 
          is less than bin number of hic file (which genome size < 2G).
          so it is invalid to get axis limit by function:
          `(len(uniq_x) - 1) * int(rslu) / 1000 / 1000`

20200826. add vlines & hlines, which position correspond to chrom length in reference
'''

def hic_viewer(hicfile, ref, rslu=500000, norm="KR", xchrom="assembly", ychrom='assembly', 
	outpfix='samplename', outfmt="pdf"):
	''' ############# generate figure from .hic file #############
	[usage] python hic_viewer.py 
	            --hicfile arabidopsis_thaliana.hic
				--ref reference.fasta
	            --rslu 500000
	            --norm KR
	            --xchrom 1:20000:50000
	            --ychrom 1:20000:50000
	            --outpfix arabidopsis_thaliana	
	[note] `ref' used to get x/y axis limit, sometimes (eg. genome size > 2G) 
	       bin number in hic file is less than which genome size < 2G
	'''
	xchrom, ychrom = str(xchrom), str(ychrom)
	hic_matrix = straw.straw(norm, hicfile, xchrom, ychrom, "BP", int(rslu))
	#uniq_x = np.unique(hic_matrix[0])
	#axis_limit = (len(uniq_x) - 1) * int(rslu) / 1000 / 1000
	infa = [i for i in SeqIO.parse(ref,'fasta')]
	reflens, reflen_total = list(), 0
	for s in infa:
		if s.id.startswith("unanchor"): continue
		reflens.append(len(s.seq))
		reflen_total += len(s.seq)

	axis_limit = reflen_total / 1000 / 1000

	x_bin = np.hstack([hic_matrix[0], hic_matrix[1]])
	y_bin = np.hstack([hic_matrix[1], hic_matrix[0]])
	link_bin = np.log(np.hstack([hic_matrix[2], hic_matrix[2]]))
	bin_matrix = pd.DataFrame({"x": x_bin, "y": y_bin, "counts": link_bin}, 
		columns=["x", 'y', 'counts'])
	bin_matrix = bin_matrix.drop_duplicates().fillna(0)
	bin_matrix = bin_matrix.reset_index(drop=True)
	fig_matrix = bin_matrix.pivot("x", "y", "counts").fillna(0)
	fig_matrix.to_csv("{}.matrix".format(outpfix), index=False, header=False, sep='\t')
		
	plt.figure(figsize=(10, 7))
	plt.imshow(fig_matrix, cmap='OrRd', 
			interpolation='nearest', 
			extent=[0, axis_limit, axis_limit, 0])
	plt.colorbar()
	
	real_x = np.array(reflens).cumsum() * plt.xlim()[1] / reflen_total 
	real_y = np.array(reflens).cumsum() * plt.ylim()[0] / reflen_total
	plt.vlines(x=real_x[:-1], ymin=0, ymax=axis_limit, color="grey", 
			linestyles="dashed", linewidth=0.5, alpha=0.6)
	plt.hlines(y=real_y[:-1], xmin=0, xmax=axis_limit, color="grey", 
			linestyles="dashed", linewidth=0.5, alpha=0.6)
	plt.xlabel("Position (Mb)", fontdict={"family":"serif"})
	plt.ylabel("Position (Mb)", fontdict={"family":"serif"})
	plt.savefig("{}.{}".format(outpfix, outfmt), dpi=600, bbox_inches='tight')

if __name__ == '__main__':
	fire.Fire(hic_viewer)
