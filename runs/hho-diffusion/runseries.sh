#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

executable="./hho-diffusion";

# Options:
if [[ $1 == "help" ]]; then
	echo -e "\nExecute tests using parameters in data.sh, creates and compile latex file, and calculate rates.\n";
	exit;
fi;

###
# Directories
root="../../";
buildir=$root"build/Schemes/"
meshdir=$root"typ2_meshes/"   # mesh directory from inside $buildir
outdir=$(pwd)"/outputs/"
errorsfile="data_rates.dat"
latexfile="rates.tex"

if [ ! -d $outdir ]; then
	mkdir $outdir
fi
\rm -r $outdir/*

###
# LOAD DATA
# (LATER: Test that each required data exists (files, parameters...))
. data.sh

echo "degrees: (edge) k=$k, (cell) l=$l"
echo "boundary conditions bc=$bc"
echo -e "test case: solution=$tcsol, diffusion=$tcdiff\n"

###
# EXECUTE FOR EACH MESH
cd $buildir
nbmesh=${#mesh[@]}
for i in `seq 1 $nbmesh`; 
do
	meshfile=$meshdir${mesh[$i]}".typ2"
	echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}.typ2"
	# Execute code
	$executable -m $meshfile -k $k -l $l -b $bc -c $tcsol $tcdiff --solver_type $solver_type
	r=$?
	if [ "$r" != "0" ]; then
		exit
	fi
	# Move outputs
	mv results.txt $outdir/results-$i.txt
	if [ -f T-solution.vtu ]; then
	mv T-solution.vtu $outdir/T-solution-$i.vtu
	fi
	mv solution.vtu $outdir/solution-$i.vtu
	mv exact-solution.vtu $outdir/exact-solution-$i.vtu
done;

# CREATE DATA FILE FOR LATEX
cd $outdir
echo -e "meshsize L2error H1error EnergyError NbEdgeDOFs" > $errorsfile
for i in `seq 1 $nbmesh`; 
do
	meshsize=$(awk '/MeshSize:/ {print $NF}' results-$i.txt)
	L2error=$(awk '/L2error:/ {print $NF}' results-$i.txt)
	H1error=$(awk '/H1error:/ {print $NF}' results-$i.txt)
	EnergyError=$(awk '/EnergyError:/ {print $NF}' results-$i.txt)
	NbEdgeDOFs=$(awk '/NbEdgeDOFs:/ {print $NF}' results-$i.txt)
	echo -e "$meshsize $L2error $H1error $EnergyError $NbEdgeDOFs" >> $errorsfile
done;

## CREATE AND COMPILE LATEX

# Look for minimal meshsize, that we reduce a bit
a=$(cat $errorsfile | awk ' END{print} ' | cut -d ' ' -f 1);
xmin=$(echo ${a/[eE]/*10^}*0.6 | bc -l);

echo -e "\\\documentclass{article}

\\\usepackage{amsfonts,latexsym,graphicx,pgfplots,pgfplotstable}

\\\setlength{\\\textwidth}{16cm}
\\\setlength{\\\textheight}{23cm}
\\\setlength{\\\oddsidemargin}{0cm}
\\\setlength{\\\evensidemargin}{0cm}
\\\setlength{\\\topmargin}{-1cm}
\\\parindent=0pt

\\\begin{document}


\\\begin{figure}
\t \\\centering
\t \\\begin{tikzpicture}[scale=1]
\t\t	 \\\begin{loglogaxis}[
\t\t			 xmin = $xmin, 
\t\t			 legend style = {
\t\t			   legend pos = south east
\t\t			 }
\t\t		 ]
\t\t		% L2 error
\t\t     \\\addplot table[x=meshsize,y={create col/linear regression={y=L2error}}] {$errorsfile}
\t\t     coordinate [pos=0.75] (A)
\t\t     coordinate [pos=1.00] (B);
\t\t     \\\xdef\\\slopea{\\\pgfplotstableregressiona}
\t\t     \\\draw (A) -| (B) node[pos=0.75,anchor=east] {\\\pgfmathprintnumber{\\\slopea}};
\t\t		% H1 error
\t\t     \\\addplot table[x=meshsize,y={create col/linear regression={y=H1error}}] {$errorsfile}
\t\t     coordinate [pos=0.75] (A)
\t\t     coordinate [pos=1.00] (B);
\t\t     \\\xdef\\\slopeb{\\\pgfplotstableregressiona}
\t\t     \\\draw (A) -| (B) node[pos=0.75,anchor=east] {\\\pgfmathprintnumber{\\\slopeb}};
\t\t		% Energy error
\t\t     \\\addplot table[x=meshsize,y={create col/linear regression={y=EnergyError}}] {$errorsfile}
\t\t     coordinate [pos=0.75] (A)
\t\t     coordinate [pos=1.00] (B);
\t\t     \\\xdef\\\slopeb{\\\pgfplotstableregressiona}
\t\t     \\\draw (A) -| (B) node[pos=0.75,anchor=east] {\\\pgfmathprintnumber{\\\slopeb}};
\t\t     \\\legend{L2 error,H1 error,Energy error};
\t\t   \\\end{loglogaxis}
\t\t \\\end{tikzpicture}
\t\t\\\end{figure}

Degrees: (edge) \$k=$k\$, (cell) \$l=$l\$

Boundary conditions: bc=$bc

Test case: tcsol=$tcsol, tcdiff=$tcdiff \n\n" > $latexfile;

for i in `seq 1 $nbmesh`; 
do
	echo -e "mesh[$i]=\\\verb!${mesh[$i]}!\n\n" >> $latexfile;
done


echo -e "\\\end{document}" >> $latexfile;


# Tests data at the end of the file
for i in `seq 1 $nbmesh`; 
do
	echo -e "Test $i:\n" >> $latexfile
	cat results-$i.txt >> $latexfile
	pdflatex $latexfile > /dev/null
done

##
# COMPUTATION OF CONVERGENCE RATES
#
echo -e "\n-----------------------------------\nData:"
echo "degrees: (edge) k=$k, (cell) l=$l"
echo "boundary conditions bc=$bc"
echo -e "test case: solution=$tcsol, diffusion=$tcdiff\n"
cd ..
make compute_rates_run


