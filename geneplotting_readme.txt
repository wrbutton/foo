GENEPLOTTING
-----------
usage: geneplotting.py [-h] [-mp] [-mg] [-dp] [-pp] [-c] [-t] [-ds] [-os]
                       [-sl] [-fl] [--fsize FSIZE FSIZE] [--input INPUT]
                       [--output OUTPUT]
                       dims dims

positional arguments:
  dims                 subplot grid dimensions, enter: rows <space> columns

optional arguments:
  -h, --help           show this help message and exit
  -c, --columnorder    plot in column order, row dflt
  -t, --test           test mode, only generates one image to check
  -ds, --dontsep       don't separate vehicles and poscons, plot one series
  -os, --othersep      dont plot poscons, highlight vehicles and provided test
                       name
  -sl, --speclayout    specify layout of suplot images, start at 1, no blanks
  -fl, --filelist      use data files in /plotfiles/filelist.txt with the text
                       file listing absoltue paths of data files
  --fsize FSIZE FSIZE  size in inches for final fig, eg: 15 8

plot layout:
  -mp, --multiplate    one plot per gene, multiple plate subplots
  -mg, --multigene     one plot per plate, multiple gene subplots

plot format:
  -dp, --dotplot       horizontal dotoplot of gene expression
  -pp, --plateplot     plate based plot of gene expression

locations:
  --input INPUT        input data file dir including genefile.txt, default:
                       Desktop/plotfiles/
  --output OUTPUT      output directory for images, default:
                       Desktop/output_images/
                       
----



Code has been updated 11/9/17 to rely on updated gct and gt modules for dataframe
	extraction and the definition of default folders for input and output
	

To Use:
	- default input is Desktop/plotfiles/
		in this folder, place your gct files of interest as well as text
		file named genelist.txt with one gene per line you want to plot
	
	- then run the code, specifying:
		- the dimensions: '3 1' will plot 3 rows one column 
		- dot plot or plate plot 'dp' or 'pp' depending on representation 
		- multiplate or multigene 'mp' or 'mg' to specify the orientation
			of the plots, one gene across multiple plates, or one plate with
			multiple genes (although only with a max of 16 genes in the genelist)
			
	- example:
		python3 geneplotting.py 3 1 -mp -dp -t
			will test one image only, plotting 3 rows of multiplate dot plot per gene